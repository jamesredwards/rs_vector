P=rs_vector
OBJECTS=
CFLAGS= -g -Wall -O3
LDLIBS=-lgsl -lgslcblas
CC=gcc

$(P) : $(OBJECTS)