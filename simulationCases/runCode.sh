#!/bin/bash

source ../.project_config
echo "BASILISK: $BASILISK"

OhOut="1e-2"
RhoIn="1e-3"
Rr="1e0"
LEVEL="12"
tmax="40.0"
zWall="0.01"

qcc -O2 -Wall -disable-dimensions -I$PWD/src-local -I$PWD/../src-local coalescenceBubble.c -o coalescenceBubble -lm
./coalescenceBubble $OhOut $RhoIn $Rr $LEVEL $tmax $zWall
