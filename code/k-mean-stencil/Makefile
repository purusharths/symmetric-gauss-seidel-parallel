CC = g++
# check without loop unrolling as well
CCFLAGS_PAR = -fopenmp -O3  -funroll-loops -ftree-vectorize -march=native -fopt-info-vec
CCFLAGS_SEQ = -O3  -funroll-loops -ftree-vectorize -march=native -fopt-info-vec

all: benchmark 

benchmark:
	$(CC) $(CCFLAGS_SEQ) bm_stencil_mean.cc -o bm_stencil

blocked:
	$(CC) $(CCFLAGS_SEQ) main_sgs.cc -o sgs_seq

simd4d:
	$(CC) $(CCFLAGS_SEQ) vanilla/main_sgs_vanilla.cc -o sgs_vanilla

simd_omp:
	$(CC) $(CCFLAGS_PAR) vanilla/main_sgs_vanilla.cc -o sgs_vanilla

# need to change flag to compare sequential/parallel
simd8d:
	$(CC) $(CCFLAGS_PAR) bm_sgs.cc -o bm_sgs

plot:
	python3 run.py 

test:
	python3 test.py

clean:
	rm vanilla blocked simd4d simd_omp simd8d
