CC = gcc
CFLAGS = -std=c23 -Wall -Wextra -O2

.PHONY: all clean

all: optics

optics: optics.c
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f optics