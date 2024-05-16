#!/bin/bash
mkdir -p build
cd build
f2py -c ../lib/_dlsodes.pyf --opt='-std=legacy -O3' -m ../lib/blkdta000.f ../lib/opkda1.f ../lib/opkda2.f ../lib/dlsodes.f
