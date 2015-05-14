
CXX = g++
CC = gcc
LAPACK =/Users/mortaza/lib
OPTS = -O3 -ftree-vectorize
mc: code3.c stringlib.c stringlib.h
$(CXX) -c code3.c stringlib.c $(OPTS)
$(CXX) code3.o stringlib.o $(OPTS) -o code3.x
