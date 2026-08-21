CC      = gcc
MPICC 	= mpicc
CFLAGS  = -std=c11 -Wall -Wextra -Wpedantic -D_POSIX_C_SOURCE=200809L \
          -g -O1 -fsanitize=address,undefined
LDFLAGS = -lm

SHARED_OBJS = csv_parser.o helpers.o rtree.o heap.o optics.o

all: optics_seq optics_mpi

optics_seq: main.o $(SHARED_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

optics_mpi: main_mpi.o doptics.o $(SHARED_OBJS)
	$(MPICC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

main_mpi.o: main_mpi.c
	$(MPICC) $(CFLAGS) -c $< -o $@

doptics.o: doptics.c
	$(MPICC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o optics_seq optics_mpi
