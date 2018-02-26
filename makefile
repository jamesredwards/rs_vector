P=rs_vector
OBJECTS=$(patsubst %.c, %.o, $(SOURCES))
SOURCES=$(wildcard *.c)
CFLAGS=-g -Wall -Wextra -O3
LDLIBS=-lm -lgsl -lgslcblas
CC=gcc

$(P) : $(OBJECTS)

clean:
	rm -f *.o
	rm -f rs_vector