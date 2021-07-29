// TODO: int -> size_t
#include <iostream>
#include <ostream>

#include "jacobi.hh"
#include "time_experiment.hh"
#include "utilities.hh"

#define BENCHMARK

int main(int argc, char **argv) {
  int n; 
  int k;
  int blocksize;
  int iterations = 1;
  if (argc < 2 || argc > 4) {
    std::cout
        << "usage: ./jacobi <grid size (n)> <stencil size (k)> <iterations>"
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
  // std::cout << "iterations: " << iterations << std::endl;

  if (k > n) {
    std::cout << "Stencil Size must be smaller than Matrix Size" << std::endl;
    exit(1);
  }

  int stencil_size = 2 * k + 1;
  double *grid = new (std::align_val_t(64)) double[n * n];
  double *local_mean = new (std::align_val_t(64)) double[n * n];
  double *alpha = new double[stencil_size * stencil_size];

  // fill alpha (coefficient stencil)
  coeff_stencil(k, alpha);
  // fill the grid
  fill_array(n, grid);
  // fill the local_mean
  fill_array(n, local_mean, true);

  // print the grid
  // print_grid(grid, n);

  // jacobi_vanilla(n, k, iterations, grid, local_mean, alpha);
  jacobi_omp(n, k, iterations, grid, local_mean, alpha);
  // print_grid(local_mean, n);
  
  delete[] grid;
  delete[] local_mean;
  delete[] alpha;
}
