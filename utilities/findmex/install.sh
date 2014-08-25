#!/bin/bash

# set environmental variables
source ../../bin/setsgenvars.sh

echo Compiling findmex
echo $SGENFC findmex.f90 -o findmex.x $SGENFLAG $SGENLIB -g -traceback -check bounds -check all
$SGENFC findmex.f90 -o findmex.x $SGENFLAG $SGENLIB -g -traceback -check bounds -check all
echo Copying executable to surfgen directory
cp findmex.x $SGENDIR 
echo Done 
