CC = g++
CCFLAGS_PAR = -fopenmp -O3  -funroll-loops -ftree-vectorize -march=native -fopt-info-vec
CCFLAGS_SEQ = -O3  -funroll-loops -ftree-vectorize -march=native -fopt-info-vec

all: jacobi jacobi_seq

jacobi:
	$(CC) $(CCFLAGS_PAR) main_jacobi.cc -o jacobi

jacobi_seq:
	$(CC) $(CCFLAGS_SEQ) main_jacobi.cc  -o jacobi_seq

# need to change flag to compare sequential/parallel
benchmark:
	$(CC) $(CCFLAGS_PAR) bm_jacobi.cc -o bm_jacobi

plot:
	python3 run.py 

test:
	python3 test.py

clean:
	rm jacobi jacobi_seq bm_jacobi