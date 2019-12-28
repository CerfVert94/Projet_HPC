#!/bin/bash
make --file=Makefile
if [ $? -eq 0 ]
then
    exe/project
    if [ $? -eq 0 ]
    then
        ./tools/plot_benchmark.py output/benchmark_dilation.dat
    fi
fi