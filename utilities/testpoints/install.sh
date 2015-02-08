#!/bin/bash

# set environmental variables
source ../../bin/setsgenvars.sh

echo Compiling testpoints
echo $SGENFC testpoints.f90 -o testpoints.x -I$SURFJAY/source/progdata.mod -I$SURFJAY/source/hddata.mod $SURFJAY/source/localcoord.o $SGENFLAG $SJAYLIB $BLASLIB
$SGENFC testpoints.f90 -o testpoints.x -I$SURFJAY/source/ $SURFJAY/source/localcoord.o -g -traceback -check uninit -check bound -check pointers $SGENFLAG $SJAYLIB $BLASLIB
echo Copying executable to surfgen directory
cp testpoints.x $SGENDIR
echo Done 
