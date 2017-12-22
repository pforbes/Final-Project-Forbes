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

all : myODE1 myODE_seq myODE_vec myODE_openMP myODE_cuda myODE_cuda_prof 

myODE1: myODE1.cpp
	g++ -O3 -o myODE1 myODE1.cpp --std=c++11
myODE_seq: myODE_seq.cpp
	g++ -O3 -o myODE_seq myODE_seq.cpp --std=c++11
myODE_vec: myODE_seq.cpp
	g++ -O3 -mavx -fopt-info-vec -o myODE_vec myODE_seq.cpp --std=c++11
myODE_openMP: myODE_openMP.cpp
	g++ -O3 -fopenmp -o myODE_openMP myODE_openMP.cpp --std=c++11
myODE_cuda: myODE_cuda.cu
	module load cuda;nvcc -o myODE_cuda $(OPT) myODE_cuda.cu -ccbin $(BIN)
myODE_cuda_prof: myODE_cuda_prof.cu
	module load cuda;nvcc -o myODE_cuda_prof $(OPT) myODE_cuda_prof.cu -ccbin $(BIN)
# TODO: add targets for building executables


.PHONY: clean
clean:
	rm -f *.o *.exe *.dat *.out
