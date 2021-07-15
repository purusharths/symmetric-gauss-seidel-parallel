// TODO: int -> size_t
#include <iostream>
#include <ostream>

#include "time_experiment.hh"
#include "stencil_mean_kernel.hh"

void print_grid(double *array, int size) {
  for (int i1 = 0; i1 < size; i1++) {
    for (int i0 = 0; i0 < size; i0++) {
      std::cout << array[i0 * size + i1] << " ";
    }
    std::cout << std::endl;
  }
    std::cout << "\n\n --------------- \n\n " << std::endl;

}

int main(int argc, char **argv) {
  int n;
  int k;
  if (argc < 2 || argc > 4) {
    std::cout << "usage: ./stencil_vanilla <grid size> <stencil size (k)>"
              << std::endl;
    exit(1);
  }
  sscanf(argv[1], "%d", &n);
  sscanf(argv[2], "%d", &k);

  double *grid = new (std::align_val_t(64)) double[n * n];
  double *local_mean = new (std::align_val_t(64)) double[n * n];

  // fill boundary and initial values
  auto g = [&](int i0, int i1) {
    return (i0 > 0 && i0 < n - 1 && i1 > 0 && i1 < n - 1) ? ((double)(i0 + i1))
                                                          : 1;
  };

  // warmup
  for (int i1 = 0; i1 < n; i1++)
    for (int i0 = 0; i0 < n; i0++)
      grid[i1 * n + i0] = g(i0, i1);

  // print the grid
  print_grid(grid, n);
  
  stencil_mean_vanilla(n, k, grid, local_mean);

  print_grid(local_mean, n);
}
