TARGET=rs_vector
OBJECTS=$(patsubst %.c, %.o, $(SOURCES))
SOURCES=$(wildcard *.c)
CFLAGS=-g -Wall -Wextra -Ofast 
LDLIBS=#-lm -lgsl -lgslcblas
CC=gcc

$(TARGET) : $(OBJECTS)

clean:
	rm -f *.o
	rm -f rs_vector