#!/bin/sh
set -x

srcs="mc.c"

clang -O3 -Wall -march=native $srcs -o mc -lm
