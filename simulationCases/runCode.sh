#!/bin/bash

# Check if .project_config exists
if [ ! -f ../.project_config ]; then
    echo "Error: .project_config not found!"
    echo "Please run the installation script first:"
    echo "  cd .. && ./reset_install_requirements.sh"
    exit 1
fi

source ../.project_config
echo "BASILISK: $BASILISK"

# Check if qcc is available
if ! command -v qcc &> /dev/null; then
    echo "Error: qcc command not found!"
    echo "Please run the installation script first:"
    echo "  cd .. && ./reset_install_requirements.sh"
    exit 1
fi

# Default parameters for droplet under shear simulation
MAXlevel="8"
L="10"
Re="1.0"
Bi="0.0"
Ca="0.3"
lambda="1.0"
tend="2.0"

FILE_NAME=${1:-droplet_shear} # pass from command line, default to droplet_shear

mkdir -p $FILE_NAME

qcc -O2 -Wall -disable-dimensions -I$PWD/src-local -I$PWD/../src-local $FILE_NAME.c -o $FILE_NAME/$FILE_NAME -lm
cd $FILE_NAME
./$FILE_NAME $MAXlevel $L $Re $Bi $Ca $lambda $tend
cd ..
