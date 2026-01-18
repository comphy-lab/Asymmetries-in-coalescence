#!/bin/bash
# ==============================================================================
# parse_params.sh - Shell Library for Parameter File Parsing
#
# A reusable shell library for parsing parameter files and generating parameter
# sweep combinations. Designed for use with simulation workflows that need
# configurable parameters.
#
# Usage:
#   source src-local/parse_params.sh
#
# After sourcing, the following functions are available:
#   - parse_param_file <file>        Parse a parameter file
#   - get_param <key> [default]      Get a parameter value
#   - generate_sweep_cases <file>    Generate sweep case files
#   - validate_required_params ...   Check required parameters exist
#   - print_params                   Print all loaded parameters
#
# Parameter File Format:
#   key=value
#   # Comments start with #
#   key2=value2  # Inline comments are supported
#
# Example parameter file (default.params):
#   # Simulation parameters
#   OhOut=1e-2
#   RhoIn=1e-3
#   MAXlevel=12
#   tmax=40.0
#
# Author: Vatsal Sanjay
# Last updated: Jan 2026
# ==============================================================================

# ==============================================================================
# Function: parse_param_file
#
# Description:
#   Parse a parameter file and export all key=value pairs as environment
#   variables with a PARAM_ prefix. Comments and empty lines are skipped.
#
# Parameters:
#   $1 - Path to the parameter file
#
# Returns:
#   0 on success, 1 if file not found
#
# Side Effects:
#   Exports environment variables named PARAM_<key> for each parameter
#
# Example:
#   parse_param_file "default.params"
#   echo $PARAM_OhOut  # Prints: 1e-2
# ==============================================================================
parse_param_file() {
    local param_file=$1

    if [ ! -f "$param_file" ]; then
        echo "ERROR: Parameter file $param_file not found" >&2
        return 1
    fi

    # Read parameters (skip comments and empty lines)
    while IFS='=' read -r key value || [ -n "$key" ]; do
        # Skip comments and empty lines
        [[ "$key" =~ ^[[:space:]]*# ]] && continue
        [[ -z "$key" ]] && continue

        # Remove inline comments and whitespace using sed and xargs
        value=$(echo "$value" | sed 's/#.*//' | xargs)
        key=$(echo "$key" | xargs)

        # Skip if key or value is empty
        [ -z "$key" ] && continue
        [ -z "$value" ] && continue

        # Export as environment variable with PARAM_ prefix
        export "PARAM_${key}=${value}"
    done < "$param_file"

    return 0
}

# ==============================================================================
# Function: get_param
#
# Description:
#   Retrieve a parameter value by key. Returns default value if not set.
#
# Parameters:
#   $1 - Parameter key (without PARAM_ prefix)
#   $2 - Default value (optional, defaults to empty string)
#
# Returns:
#   Prints the parameter value to stdout
#
# Example:
#   OhOut=$(get_param "OhOut" "1e-2")
#   MAXlevel=$(get_param "MAXlevel" "10")
# ==============================================================================
get_param() {
    local key=$1
    local default=${2:-}
    local var_name="PARAM_${key}"
    echo "${!var_name:-$default}"
}

# ==============================================================================
# Function: generate_sweep_cases
#
# Description:
#   Generate parameter files for all combinations in a parameter sweep.
#   Creates a Cartesian product of all SWEEP_* variables defined in the
#   sweep file.
#
# Parameters:
#   $1 - Path to the sweep configuration file
#
# Returns:
#   Prints path to temporary directory containing generated case files
#   Returns 1 on error
#
# Sweep File Format:
#   BASE_CONFIG=path/to/base.params    # Required: base parameter file
#   OUTPUT_TEMPLATE=output/{Rr}/{Oh}   # Output directory template
#   SWEEP_Rr=0.5,0.7,1.0               # Sweep variable (comma-separated)
#   SWEEP_Oh=1e-2,1e-3                 # Another sweep variable
#
# Generated Files:
#   Creates case_0000.params, case_0001.params, etc. in a temporary directory.
#   Each file contains the base configuration with sweep values overridden.
#
# Example:
#   cases_dir=$(generate_sweep_cases "sweep.params")
#   for case_file in "$cases_dir"/*.params; do
#       parse_param_file "$case_file"
#       # Run simulation with these parameters...
#   done
# ==============================================================================
generate_sweep_cases() {
    local sweep_file=$1

    if [ ! -f "$sweep_file" ]; then
        echo "ERROR: Sweep file $sweep_file not found" >&2
        return 1
    fi

    # Source the sweep file to get variables
    source "$sweep_file"

    # Check if BASE_CONFIG is defined
    if [ -z "$BASE_CONFIG" ]; then
        echo "ERROR: BASE_CONFIG not defined in sweep file" >&2
        return 1
    fi

    # Parse base configuration
    parse_param_file "$BASE_CONFIG"

    # Create temporary directory for generated cases
    local temp_dir=$(mktemp -d "${TMPDIR:-/tmp}/sweep.XXXXXX")

    # Extract sweep variables from the sweep file
    # Match lines like: SWEEP_varname=value1,value2,value3
    local sweep_vars=()
    local sweep_values=()

    while IFS='=' read -r line; do
        # Remove comments
        line=$(echo "$line" | sed 's/#.*//')
        [ -z "$line" ] && continue

        # Match SWEEP_* variables using regex
        if [[ "$line" =~ ^[[:space:]]*SWEEP_([^=]+)=(.+)$ ]]; then
            var_name="${BASH_REMATCH[1]}"
            var_values="${BASH_REMATCH[2]}"
            sweep_vars+=("$var_name")
            sweep_values+=("$var_values")
        fi
    done < "$sweep_file"

    if [ ${#sweep_vars[@]} -eq 0 ]; then
        echo "ERROR: No SWEEP_* variables found in $sweep_file" >&2
        rm -rf "$temp_dir"
        return 1
    fi

    # Generate cartesian product of sweep values
    local case_num=0

    # --------------------------------------------------------------------------
    # Recursive function: generate_recursive
    #
    # Recursively generates all combinations of sweep variable values.
    # At each level of recursion, iterates through values for one variable.
    # When all variables are assigned, creates a parameter file.
    #
    # Parameters:
    #   $1      - Current recursion depth (0 to num_sweep_vars-1)
    #   $2...   - Accumulated values for variables 0 to depth-1
    # --------------------------------------------------------------------------
    generate_recursive() {
        local depth=$1
        shift
        local current_values=("$@")

        if [ $depth -eq ${#sweep_vars[@]} ]; then
            # Base case: all variables assigned, create parameter file
            local case_file="${temp_dir}/case_$(printf "%04d" $case_num).params"
            local output_dir="$OUTPUT_TEMPLATE"

            # Copy base config
            cp "$BASE_CONFIG" "$case_file"

            # Override with sweep values and build output directory path
            for i in "${!sweep_vars[@]}"; do
                local var="${sweep_vars[$i]}"
                local val="${current_values[$i]}"

                # Replace in parameter file (or append if not present)
                if grep -q "^${var}=" "$case_file"; then
                    sed -i'.bak' "s|^${var}=.*|${var}=${val}|" "$case_file"
                else
                    echo "${var}=${val}" >> "$case_file"
                fi
                rm -f "${case_file}.bak"

                # Replace placeholder in output template: {varname} -> value
                output_dir="${output_dir//\{${var}\}/${val}}"
            done

            # Set output directory in the case file
            if grep -q "^output_dir=" "$case_file"; then
                sed -i'.bak' "s|^output_dir=.*|output_dir=${output_dir}|" "$case_file"
            else
                echo "output_dir=${output_dir}" >> "$case_file"
            fi
            rm -f "${case_file}.bak"

            ((case_num++))
            return
        fi

        # Recursive case: iterate through values for current variable
        local values="${sweep_values[$depth]}"
        IFS=',' read -ra value_array <<< "$values"

        for val in "${value_array[@]}"; do
            val=$(echo "$val" | xargs)  # Trim whitespace
            generate_recursive $((depth + 1)) "${current_values[@]}" "$val"
        done
    }

    # Start recursion with depth 0 and no accumulated values
    generate_recursive 0

    echo "$temp_dir"
    return 0
}

# ==============================================================================
# Function: validate_required_params
#
# Description:
#   Check that all required parameters have been loaded. Prints error messages
#   for any missing parameters.
#
# Parameters:
#   $@ - List of required parameter names (without PARAM_ prefix)
#
# Returns:
#   0 if all parameters are present, 1 if any are missing
#
# Example:
#   if ! validate_required_params OhOut RhoIn MAXlevel tmax; then
#       echo "Missing required parameters, exiting"
#       exit 1
#   fi
# ==============================================================================
validate_required_params() {
    local missing=0

    for param in "$@"; do
        local var_name="PARAM_${param}"
        if [ -z "${!var_name}" ]; then
            echo "ERROR: Required parameter '$param' not found" >&2
            missing=1
        fi
    done

    return $missing
}

# ==============================================================================
# Function: print_params
#
# Description:
#   Print all currently loaded parameters for debugging purposes.
#   Lists all environment variables with the PARAM_ prefix.
#
# Parameters:
#   None
#
# Returns:
#   Prints formatted parameter list to stdout
#
# Example:
#   parse_param_file "default.params"
#   print_params
#   # Output:
#   # Loaded parameters:
#   #   MAXlevel = 12
#   #   OhOut = 1e-2
#   #   ...
# ==============================================================================
print_params() {
    echo "Loaded parameters:"
    env | grep "^PARAM_" | sort | while IFS='=' read -r key value; do
        key="${key#PARAM_}"
        echo "  $key = $value"
    done
}

# ==============================================================================
# Function: validate_restart_file
#
# Description:
#   Validate a restart file exists, is non-empty, and is readable.
#   Provides clear error messages for each failure case.
#
# Parameters:
#   $1 - Path to restart file (optional, defaults to "restart")
#
# Returns:
#   0 if valid, 1 if invalid
#
# Example:
#   if ! validate_restart_file "restart"; then
#       echo "Cannot proceed without valid restart file"
#       exit 1
#   fi
# ==============================================================================
validate_restart_file() {
    local restart_file="${1:-restart}"

    if [ ! -f "$restart_file" ]; then
        echo "ERROR: Restart file not found: $restart_file" >&2
        return 1
    fi

    if [ ! -s "$restart_file" ]; then
        echo "ERROR: Restart file is empty: $restart_file" >&2
        return 1
    fi

    if [ ! -r "$restart_file" ]; then
        echo "ERROR: Restart file not readable: $restart_file" >&2
        return 1
    fi

    return 0
}
