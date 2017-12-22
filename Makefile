# Warnings
WFLAGS	:= -Wall -Wextra -Wsign-conversion -Wsign-compare

# Optimization and architecture
OPT		:= -O3
ARCH   	:= -march=native

# Language standard
CCSTD	:= -std=c99
CXXSTD	:= -std=c++11

# Linker options
LDOPT 	:= $(OPT)
LDFLAGS := 
BIN = "/usr/local/gcc/6.4.0/bin/gcc"
.DEFAULT_GOAL := all

.PHONY: debug
debug : OPT  := -O0 -g -fno-omit-frame-pointer -fsanitize=address
debug : LDFLAGS := -fsanitize=address
debug : ARCH :=
debug : $(EXEC)

all : problem1A problem3 problem2 

problem1A: problem1A.cpp
	g++ -fopenmp -fopt-info-vec -Wall -Wextra -Wsign-conversion -Wsign-compare -O3 -o problem1A problem1A.cpp --std=c++14

problem3: problem3.cu
	module load cuda;nvcc -o problem3 $(OPT) problem3.cu -ccbin $(BIN)
problem2: problem2.cu
	module load cuda;nvcc -o problem2 $(OPT) problem2.cu -ccbin $(BIN)
# TODO: add targets for building executables


.PHONY: clean
clean:
	rm -f *.o *.exe *.dat *.out
