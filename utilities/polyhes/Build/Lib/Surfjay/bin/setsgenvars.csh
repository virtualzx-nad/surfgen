#!/bin/tcsh
set stacksize unlimited
#Set surfgen compilation flags and link line options
setenv SJAYFLAG '-auto -lpthread -assume byterecl -i8 -no-opt-matmul -assume nounderscore'
setenv SJAYLIB '/global/u2/m/malbon/Workspace/Surfgen/surfJay/utilities/polyhes/Build/Lib/Surfjay/lib/libsurfgen.a  -Wl,--start-group  /usr/common/usg/intel/13.0.028/composer_xe_2013.1.117/mkl/lib/intel64/libmkl_intel_ilp64.a /usr/common/usg/intel/13.0.028/composer_xe_2013.1.117/mkl/lib/intel64/libmkl_intel_thread.a /usr/common/usg/intel/13.0.028/composer_xe_2013.1.117/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -auto -lpthread -parallel '
setenv SJAYFC 'ifort'
setenv SJAYDIR '/global/u2/m/malbon/Workspace/Surfgen/surfJay/utilities/polyhes/Build/Lib/Surfjay/bin'
set found=`echo "$PATH" | tr ":" "\n" | grep -x "$SJAYDIR"`
if ( ${?found} == 0 ) then
   setenv PATH "${PATH}:${SJAYDIR}"
else 
   echo surfgen found in PATH 
endif
