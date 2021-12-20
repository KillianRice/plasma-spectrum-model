#!/bin/bash
# need to give execution priviledges with chmod u+x execute.sh

# compile command
/usr/bin/g++-10 -I source/header/ -std=c++20 -fopenmp main.cpp MDQT.cpp source/md/* source/qt/* source/utils.cpp source/PlasmaSettings.cpp source/wigner-symbols/* -o main -O3 -lm -larmadillo -lstdc++fs

# run command: i should iterate over valid task array values
for ((i=0;i<=0;i++))
do
    ./main -p test -a $i -s ucnp.settings
done