#!/bin/bash

# set environmental variables
source ../../bin/setsgenvars.sh

echo Compiling findmex
echo $SJAYFC findmex.f90 -o findmex.x $SJAYFLAG $SJAYLIB -g -traceback -check bounds -check all
$SJAYFC findmex.f90 -o findmex.x $SJAYFLAG $SJAYLIB -g -traceback -check bounds -check all
#echo Copying executable to surfgen directory
#cp findmex.x $SGENDIR 
echo Done 
