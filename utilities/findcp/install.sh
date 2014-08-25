#!/bin/bash

# set environmental variables
source ../../bin/setsgenvars.sh

echo Compiling findcp
echo $SGENFC findcp.f90 -o findcp.x $SJAYFLAG $SJAYLIB
$SGENFC findcp.f90 -o findcp.x $SURFJAY/source/progdata.mod $SURFJAY/source/hddata.mod $SGENFLAG $SJAYLIB
echo Cleaning up
rm opttools.mod
echo Copying executable to surfgen directory
cp findcp.x $SGENDIR 
echo Done 
