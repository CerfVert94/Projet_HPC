#!/bin/bash
make --file=Makefile
echo $?
if [ $? ]
then
    exe/project
    if [ $? ]
    then
        ./tools/plot_benchmark.py output/benchmark_dilation.dat
    fi
fi