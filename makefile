target = rs_vector
src = $(wildcard src/*.c)
obj = $(src:.c=.o)
LDFLAGS = -lm -lgsl -lgslcblas
CFLAGS = -Wall -Wextra -Wpedantic -Ofast -std=c99
CC = gcc

$(target): $(obj)
	$(CC) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(obj) target