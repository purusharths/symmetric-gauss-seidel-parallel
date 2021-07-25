#include "utilities.hh"
#include <iostream>
#include <omp.h>
#include <unistd.h>

int main(int argc, char **argv) {
  int n, k, blocked;
  if (argc < 2) {
    std::cout << "enter n and k" << std::endl;
    exit(1);
  }

  sscanf(argv[1], "%d", &n);
  sscanf(argv[2], "%d", &k);
  //   std::vector<std::tuple<int, int>> loop_domains;
  //   std::vector<std::vector<std::tuple<int, int>>> vals;

  double *grid = new (std::align_val_t(64)) double[n * n];

  fill_array(n, grid);
  print_grid(grid, n);

  for (int i = 0; i < n; i++) {
    for (int c = 0; c < k + 1; c++) {

      int p = i;
      int q = c;

      int loop_count = std::min(p, (n - 1 - q) / (k + 1)) + 1;
      // #pragma omp parallel for shared(p, q)
      for (int l = 1; l <= loop_count; l++) {
        usleep(100); // usleep stands in for real work
#pragma omp critical
        {
          std::cout << grid[p * n + q] << " "; //<< "  (" << tid << ")\n";
          p = p - 1;
          q = q + (k + 1);
        }
      }
      // #pragma omp barrier
      std::cout << "---" << std::endl;
    }
  }

  std::cout << "xxxx" << std::endl;

  for (int c = k + 1; c < n; c++) {
    int p = n - 1, q = c;
    int loop_count = 1 + std::floor((n - 1 - q) / (k + 1));
#pragma omp parallel for shared(p, q)
    for (int l = 1; l <= loop_count; l++) {
#pragma omp critical
      {
        usleep(100);                         // usleep stands in for real work
        std::cout << grid[p * n + q] << " "; //<< "  (" << tid << ")\n";
        p = p - 1;
        q = q + (k + 1);
      }
    }
    std::cout << "---" << std::endl;
  }

  print_grid(grid, n);
}
