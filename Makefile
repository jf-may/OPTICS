CC      = gcc
CFLAGS  = -std=c11 -Wall -Wextra -Wpedantic \
          -g -O1 \
          -fsanitize=address,undefined

SRC = helpers.c heap.c rtree.c optics.c csv_parser.c main.c
OBJ = $(SRC:.c=.o)
BIN = optics

.PHONY: all clean

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ -lm

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ) $(BIN)
