# Compiler and flags
CC = mpicc
CFLAGS = -fopenmp -O3 -Wall
LDFLAGS = -lm

# Source files
SOURCES = source/graph_utils.c source/bitset.c
HEADERS = source/graph_utils.h source/bitset.h

# Executables
TARGETS = csrMaker.out BFS_basic.out BFS_hybrid.out

# Default target
all: $(TARGETS)

# Build CSR maker
csrMaker.out: source/csrMaker.c $(SOURCES) $(HEADERS)
	${CC} $(CFLAGS) source/csrMaker.c $(SOURCES) -o $@ $(LDFLAGS)

# Build basic distributed BFS
BFS_basic.out: source/BFS_basic.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) source/BFS_hybrid.c $(SOURCES) -o $@ $(LDFLAGS)


# Build hybrid BFS
BFS_hybrid.out: source/BFS_hybrid.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) source/BFS_hybrid.c $(SOURCES) -o $@ $(LDFLAGS)


# Clean build artifacts
clean:
	rm -f $(TARGETS) *.o

# Rebuild everything
rebuild: clean all

.PHONY: all clean rebuild
