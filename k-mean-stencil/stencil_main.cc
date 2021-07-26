// TODO: int -> size_t
#include <iostream>
#include <ostream>

// #include "stencil_benchmark.hh"
#include "stencil_large_simd.hh"
#include "stencil_vanilla.hh"
#include "stencil_blocked.hh"
#include "stencil_omp_simd.hh"
#include "stencil_simd.hh"
#include "time_experiment.hh"
#include "utilities.hh"

// #define BENCHMARK  

int main(int argc, char **argv) {
  int n; // 0 -> n = n+1
  int k;
  int blocksize;
  if (argc < 2 || argc > 4) {
    std::cout
        << "usage: ./stencil_vanilla <grid size> <stencil size (k)>\n"
        << "\t./stencil_vanilla <grid size> <stencil size (k)> <block size>"
        << std::endl;
    exit(1);
  }
  sscanf(argv[1], "%d", &n);
  sscanf(argv[2], "%d", &k);
  sscanf(argv[3], "%d", &blocksize);

  if ((n) % blocksize != 0) {
    std::cout << "Blocksize must be a multiple of n" << std::endl;
    exit(1);
  }

  if (k > n) {
    std::cout << "Stencil Size must be smaller than n" << std::endl;
    exit(1);
  }

  double *grid = new (std::align_val_t(64)) double[n * n];
  double *local_mean = new (std::align_val_t(64)) double[n * n];

  // fill boundary and initial values
  auto g = [&](int i0, int i1) {
    return (i0 > 0 && i0 < n - 1 && i1 > 0 && i1 < n - 1) ? 5 : 1;
    // change after questionmark: ((double)(i0 + i1))
  };

  // warmup
  for (int i1 = 0; i1 < n; i1++)
    for (int i0 = 0; i0 < n; i0++)
      grid[i1 * n + i0] = g(i0, i1);

#if defined(BENCHMARK)
  benchmark(k, blocksize);
#endif

  // print the grid
  print_grid(grid, n);

  // stencil_mean_vanilla(n, k, grid, local_mean);
  // stencil_mean_blocked(n, k, grid, local_mean, blocksize);
  // stencil_simd(n, k, grid, local_mean, blocksize);
  stencil_simd_large_vec(n, k, grid, local_mean, blocksize);
  // stencil_omp_simd(n, k, grid, local_mean, blocksize);
  print_grid(local_mean, n);
}
