#!/bin/sh
set -x

clang -O3 -Wall -march=native mc.c -o mc -lm -lpthread -static
