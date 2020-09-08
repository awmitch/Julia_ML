#!/bin/bash

gcc julia_gen.c -c -fopenmp -lm
gcc julia_gen.o main.c -o main -fopenmp -lm
