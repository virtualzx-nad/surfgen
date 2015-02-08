#!/bin/bash

# set environmental variables
source ../../bin/setsgenvars.sh

echo Compiling findcp
echo $SJAYFC findcp.f90 -o findcp.x $SJAYFLAG $SJAYLIB $BLAS_LIB
$SJAYFC findcp.f90 -o findcp.x $SJAYFLAG $SJAYLIB $BLAS_LIB
echo Cleaning up
rm opttools.mod
echo Copying executable to surfgen directory
cp findcp.x $SJAYDIR 
echo Done 
