#!/bin/bash
#
# Install libs
make clean
make libs
cp lib/lib*64.a ../../
echo "Installed"