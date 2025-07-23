#!/bin/bash

source ../.project_config
echo "BASILISK: $BASILISK"

OhOut="1e-2"
RhoIn="1e-3"
Rr="1e0"
LEVEL="12"
tmax="40.0"
zWall="0.01"

FILE_NAME=$1 # pass from command line

mkdir -p $FILE_NAME

qcc -O2 -Wall -disable-dimensions -I$PWD/src-local -I$PWD/../src-local $FILE_NAME.c -o $FILE_NAME/$FILE_NAME -lm
cp -r DataFiles $FILE_NAME/
cd $FILE_NAME
./$FILE_NAME $OhOut $RhoIn $Rr $LEVEL $tmax $zWall # pass from command line
cd ..
