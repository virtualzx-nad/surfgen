#!/bin/bash
ulimit -s unlimited
#Set surfgen compilation flags and link line options
export SJAYFLAG='-auto -lpthread -assume byterecl -i8 -no-opt-matmul -assume nounderscore'
export SJAYLIB='/global/u2/m/malbon/Workspace/Surfgen/surfJay/utilities/polyhes/Build/Lib/Surfjay/lib/libsurfgen.a  -Wl,--start-group  /usr/common/usg/intel/13.0.028/composer_xe_2013.1.117/mkl/lib/intel64/libmkl_intel_ilp64.a /usr/common/usg/intel/13.0.028/composer_xe_2013.1.117/mkl/lib/intel64/libmkl_intel_thread.a /usr/common/usg/intel/13.0.028/composer_xe_2013.1.117/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -auto -lpthread -parallel '
export SJAYFC='ifort'
export SJAYDIR='/global/u2/m/malbon/Workspace/Surfgen/surfJay/utilities/polyhes/Build/Lib/Surfjay/bin'
if [[ ":$PATH:" =~ ":$SJAYDIR:" ]]; then
  echo surfgen found in PATH variable
else
  export PATH=$PATH:$SJAYDIR
fi
