#!/bin/bash

# set environmental variables
source ../../bin/setsgenvars.sh

echo Compiling testpoints
echo $SGENFC testpoints.f90 -o testpoints.x $SGENFLAG $SGENLIB
$SGENFC testpoints.f90 -o testpoints.x -g -fbacktrace $SGENFLAG $SGENLIB
echo Copying executable to surfgen directory
cp testpoints.x $SGENDIR
echo Done 
