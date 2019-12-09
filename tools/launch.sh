#!/bin/bash
make --file=Makefile
exe/project
./tools/plot_benchmark.py output/benchmark_dilation.dat
