#pragma once
#include <cmath>
#include <iostream>

void print_grid(double *array, int size) {
  std::cout << "\n\n --------------- \n\n " << std::endl;
  for (int i1 = 0; i1 < size; i1++) {
    for (int i0 = 0; i0 < size; i0++) {
      std::cout << array[i0 * size + i1] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "\n\n --------------- \n\n " << std::endl;
}

void fill_array(int n, double *grid, bool zero_grid = false) {

  if (zero_grid) {
    auto g = [&](int i0, int i1) {
      return (i0 > 0 && i0 < n - 1 && i1 > 0 && i1 < n - 1) ? 0 : 0;
    };
    for (int i1 = 0; i1 < n; i1++) {
      for (int i0 = 0; i0 < n; i0++) {
        grid[i1 * n + i0] = g(i0, i1);
      }
    }
  } else {
    auto g = [&](int i0, int i1) {
      return (i0 > 0 && i0 < n - 1 && i1 > 0 && i1 < n - 1)
                 ? i0 + i1
                 : i0 * i1 + 5 * i1 + 5 * i0;
    };
    for (int i1 = 0; i1 < n; i1++) {
      for (int i0 = 0; i0 < n; i0++) {
        grid[i1 * n + i0] = g(i0, i1);
      }
    }
  }
}

void coeff_stencil(int k, double *alpha) {
  int stencil_size = 2 * k + 1;
  double m = 50 + ((stencil_size * stencil_size) - 1);

  for (int i = 0; i < stencil_size; i++) {
    for (int j = 0; j < stencil_size; j++) {
      if (i == k && j == k) {
        alpha[i * stencil_size + j] = 50 / m;
      } else {
        alpha[i * stencil_size + j] =
            std::pow(std::max(std::abs(k - i), std::abs(k - j)), -2) / m;
      }
    }
  }
}