#!/bin/sh
set -x

srcs="mcml.c"

clang -O3 -Wall -march=native $srcs -o mcml
