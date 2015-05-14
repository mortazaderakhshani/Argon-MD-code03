
CXX = g++
CC = gcc
LAPACK =/Users/mortaza/lib 
OPTS = -O3 -ftree-vectorize

mc: code3.c stringlib.c stringlib.h  
	$(CC) -c  code3.c stringlib.c $(OPTS) 
	$(CC) code3.o stringlib.o $(OPTS) -o code3.x
