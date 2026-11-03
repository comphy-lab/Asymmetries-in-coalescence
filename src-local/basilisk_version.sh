#!/bin/bash
# basilisk_version.sh - Legacy contour/sweep Basilisk version configuration
#
# Description:
#   Defines the Basilisk version and provides a helper function for installation.
#   Legacy contour/sweep launchers retain their historical January ref.
#
# Usage:
#   source src-local/basilisk_version.sh
#   install_basilisk              # Install with default settings
#   install_basilisk --hard       # Force reinstall
#
# Author: Vatsal Sanjay
# Organization: CoMPhy Lab, Durham University

# ============================================================
# Basilisk Version Configuration
# ============================================================
# Legacy contour/sweep workflows use their own committed historical pin.
# Do not replace a differently pinned installation or silently migrate it.
_basilisk_version_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
if ! source "${_basilisk_version_dir}/src-local/basilisk_pin.sh"; then
    return 1 2>/dev/null || exit 1
fi
if ! BASILISK_REF=$(basilisk_read_pin "${_basilisk_version_dir}/basilisk-legacy.ref"); then
    return 1 2>/dev/null || exit 1
fi
if [[ $BASILISK_REF != v2026-01-13 ]]; then
    echo "ERROR: legacy workflows require v2026-01-13; refusing an unreviewed version change." >&2
    return 1 2>/dev/null || exit 1
fi

# Installation script URL (from comphy-lab/basilisk-C)
BASILISK_INSTALL_URL="https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh"

# ============================================================
# Installation Helper
# ============================================================
# Install Basilisk with the pinned version
# Usage: install_basilisk [options]
# Options are passed directly to the install script (e.g., --hard)
install_basilisk() {
    local pin
    if [[ $(pwd -P) != "$_basilisk_version_dir" ]]; then
        echo "ERROR: legacy installer must run from its own component root; installation refused." >&2
        return 1
    fi
    pin=$(basilisk_read_pin "${_basilisk_version_dir}/basilisk-legacy.ref") || return 1
    if [[ $pin != v2026-01-13 || $pin != "$BASILISK_REF" ]]; then
        echo "ERROR: legacy Basilisk pin changed or is not v2026-01-13; installation refused." >&2
        return 1
    fi
    # The installer targets the current directory. Never let --hard replace
    # an existing installation whose lock is missing or names another ref.
    if [[ -e "$PWD/basilisk" || -L "$PWD/basilisk" ]]; then
        if [[ -L "$PWD/basilisk" || ! -d "$PWD/basilisk" ]]; then
            echo "ERROR: existing Basilisk installation is not a local directory; installation refused." >&2
            return 1
        fi
        basilisk_require_lock "$PWD/basilisk/.comphy-lock" "$BASILISK_REF" || return 1
    fi
    if ! command -v curl >/dev/null 2>&1; then
        echo "ERROR: curl not found (required to install Basilisk)" >&2
        return 1
    fi

    echo "Installing Basilisk (ref: ${BASILISK_REF})..."
    curl -sL "$BASILISK_INSTALL_URL" | bash -s -- --ref="$BASILISK_REF" "$@"
}
