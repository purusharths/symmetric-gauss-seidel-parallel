CC = g++
CCFLAGS_PAR = -fopenmp -O3  -funroll-loops -ftree-vectorize -march=native -fopt-info-vec
CCFLAGS_SEQ = -O3  -funroll-loops -ftree-vectorize -march=native -fopt-info-vec
CCFLAGS_TEST =  -fopenmp -O3 -funroll-loops -march=native
all: sgs sgs_seq

sgs:
	$(CC) $(CCFLAGS_PAR) main_sgs.cc -o sgs

sgs_seq:
	$(CC) $(CCFLAGS_SEQ) main_sgs.cc -o sgs_seq

sgs_vanilla:
	$(CC) $(CCFLAGS_SEQ) vanilla/main_sgs_vanilla.cc -o sgs_vanilla

# need to change flag to compare sequential/parallel
benchmark:
	$(CC) $(CCFLAGS_PAR) bm_sgs.cc -o bm_sgs

plot:
	python3 run.py

test:
	$(CC) $(CCFLAGS_TEST) test.cc -o unit_test

clean:
	rm sgs sgs_seq sgs_vanilla bm_sgs
