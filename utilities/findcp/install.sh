#!/bin/bash

# set environmental variables
source ../../bin/setsgenvars.sh

echo Compiling findcp
echo $SGENFC findcp.f90 -o findcp.x $SJAYFLAG $SJAYLIB $BLASLIB
$SGENFC findcp.f90 -o findcp.x $SGENFLAG $SJAYLIB $BLASLIB
echo Cleaning up
rm opttools.mod
echo Copying executable to surfgen directory
cp findcp.x $SGENDIR 
echo Done 
