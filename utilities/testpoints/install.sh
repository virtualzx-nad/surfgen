#!/bin/bash

# set environmental variables
source ../../bin/setsgenvars.sh

echo Compiling testpoints
echo $SGENFC testpoints.f90 -o testpoints.x $SGENFLAG $SGENLIB
$SGENFC testpoints.f90 -o testpoints.x $SGENFLAG $SGENLIB
echo Done 
