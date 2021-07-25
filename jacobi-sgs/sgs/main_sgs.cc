// TODO: int -> size_t
#include <iostream>
#include <ostream>

// #include "sgs-bw.hh"
#include "sgs-omp-bw.hh"
// #include "debug-sws-fw.hh"
// #include "sgs-omp-fw.hh"
// #include "fw-sgs-thread.hh"
// #include "sws-omp-vec.hh"
#include "time_experiment.hh"
#include "utilities.hh"

#define BENCHMARK

int main(int argc, char **argv) {
  int n; // 0 -> n = n+1
  int k;
  int blocksize;
  int iterations = 1;
  if (argc < 2 || argc > 5) {
    std::cout << "usage: ./jacobi <grid size> <stencil size (k)>\n"
              << "\t./jacobi <grid size> <stencil size (k)> <>"
              << "\t./jacobi <grid size> <stencil size (k)> <> <iterations>"
              << std::endl;
    exit(1);
  }
  sscanf(argv[1], "%d", &n);
  sscanf(argv[2], "%d", &k);
  // sscanf(argv[3], "%d", &blocksize);
  sscanf(argv[3], "%d", &iterations);

  // if ((n) % blocksize != 0) {
  //   std::cout << "Blocksize must be a multiple of n" << std::endl;
  //   exit(1);
  // }

  if (k > n) {
    std::cout << "Stencil Size must be smaller than Matrix Size" << std::endl;
    exit(1);
  }

  int stencil_size = 2 * k + 1;
  double *grid = new (std::align_val_t(64)) double[n * n];
  double *local_mean = new (std::align_val_t(64)) double[n * n];
  double *alpha = new double[stencil_size * stencil_size];

  // fill boundary and initial values
  auto g = [&](int i0, int i1) {
    return (i0 > 0 && i0 < n - 1 && i1 > 0 && i1 < n - 1)
              //  ? i0+i1 : ((double)(i0+i1))/n; };
               ? i0 + i1
               : i0 * i1 + 5 * i1 + 5 * i0;
  }; // 

  // fill alpha (coefficient stencil)
  for (int i = 0; i < 2 * k + 1; i++) {
    for (int j = 0; j < 2 * k + 1; j++) {
      if (i == k && j == k) {
        alpha[i * stencil_size + j] = 50;
      } else {
        alpha[i * stencil_size + j] =
            std::pow(std::max(std::abs(k - i), std::abs(k - j)), -2);
      }
    }
  }

  //  fill in grid / local mean warmup
  for (int i1 = 0; i1 < n; i1++) {
    for (int i0 = 0; i0 < n; i0++) {
      grid[i1 * n + i0] = g(i0, i1);
      local_mean[i1 * n + i0] = g(i0, i1);
    }
  }

  for (int i1 = 0; i1 < n; i1++) {
    for (int i0 = 0; i0 < n; i0++) {
      local_mean[i1 * n + i0] = 0.0;
    }
  }

  // alpha (coefficient stencil)
  // for (int i = 0; i < 2 * k + 1; i++) {
  //   for (int j = 0; j < 2 * k + 1; j++) {
  //     std::cout << alpha[i * stencil_size + j] << " ";
  //   }
  //   std::cout << std::endl;
  // }


  // print the grid
  print_grid(grid, n);
  // print_grid(local_mean, n);

  // forward_gauss_sidel(n, k, grid, local_mean);
  // print_grid(local_mean, n);
  // backward_gauss_sidel(n, k, local_mean, local_mean);
  backward_gauss_sidel_omp(n, k, grid, local_mean);
  // forward_gauss_sidel_omp(n, k, grid, local_mean);
  print_grid(local_mean, n);
  delete [] grid;
  delete [] local_mean;
  delete [] alpha;
}