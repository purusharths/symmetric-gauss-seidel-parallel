// TODO: int -> size_t
#include <iostream>
#include <ostream>

#include "fw_bw.hh"
#include "sgs-omp-bw.hh"
#include "sgs-omp-fw.hh"

#include "time_experiment.hh"
#include "utilities.hh"

#define BENCHMARK

int main(int argc, char **argv) {
  int n; // 0 -> n = n+1
  int k;
  int blocksize;
  int iterations;
  if (argc < 2 || argc > 4) {
    std::cout << "usage: ./sgs <grid size> <stencil size (k)> <iterations>"
              << std::endl;
    exit(1);
  }
  sscanf(argv[1], "%d", &n);
  sscanf(argv[2], "%d", &k);
  // sscanf(argv[3], "%d", &blocksize);
  if (argc > 3) {
    sscanf(argv[3], "%d", &iterations);
  } else {
    iterations = 100;
  }

  if (k > n) {
    std::cout << "Stencil Size must be smaller than Matrix Size" << std::endl;
    exit(1);
  }
  
  // std::cout << "Iterations: " << iterations;

  int stencil_size = 2 * k + 1;
  double *grid = new (std::align_val_t(64)) double[n * n];
  // double *local_mean = new (std::align_val_t(64)) double[n * n];
  double *alpha = new double[stencil_size * stencil_size];

  // fill boundary and initial values
  auto g = [&](int i0, int i1) {
    return (i0 > 0 && i0 < n - 1 && i1 > 0 && i1 < n - 1)
               //  ? i0+i1 : ((double)(i0+i1))/n; };
               ? i0 + i1
               : i0 * i1 + 5 * i1 + 5 * i0;
  }; //

  // fill alpha (coefficient stencil)
  coeff_stencil(k, alpha);

  //  fill in grid / local mean warmup
  for (int i1 = 0; i1 < n; i1++) {
    for (int i0 = 0; i0 < n; i0++) {
      grid[i1 * n + i0] = g(i0, i1);
    }
  }

  // print the grid
  // print_grid(grid, n);
  // print_grid(local_mean, n);
  for(int i = 0; i < iterations; i++){
    forward_gauss_sidel_omp(n, k, grid, alpha);
    backward_gauss_sidel_omp(n, k, grid, alpha);
  }
// gauss_sidel_omp(n, k, iterations, grid, local_mean, alpha);
  // if (local_mean != grid) {
  //   std::swap(grid, local_mean);
  // }

  print_grid(grid, n);
  // std::cout << "Final Matrix:";
  // print_grid(local_mean, n);
 
  delete[] grid;
  // delete[] local_mean;
  delete[] alpha;
}
