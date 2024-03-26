#!/bin/sh

srcs="main.cpp Photon.cpp utils.cpp"

c++ -O3 -Wall -march=native $srcs -o mcml
