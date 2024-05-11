#!/bin/bash

f2py -c odepack/_dlsodes.pyf --opt='-std=legacy' -m odepack/blkdta000.f odepack/opkda1.f odepack/opkda2.f odepack/dlsodes.f