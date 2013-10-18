#!/bin/bash

# set environmental variables
source ../../bin/setsgenvars.sh

progname=analyse
echo Compiling analyse
echo $SGENFC ${progname}.f90 -o ${progname}.x $SGENFLAG $SGENLIB 
$SGENFC ${progname}.f90 -o ${progname}.x $SGENFLAG $SGENLIB
echo Copying executable to surfgen directory
cp ${progname}.x $SGENDIR 
echo Done 
