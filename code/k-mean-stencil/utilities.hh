#pragma once
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

void fill_array(int n, double *grid){
  auto g = [&](int i0, int i1) {
    return (i0 > 0 && i0 < n - 1 && i1 > 0 && i1 < n - 1)
               ? i0+i1 : ((double)(i0+i1))/n; };

  for (int i1 = 0; i1 < n; i1++) {
    for (int i0 = 0; i0 < n; i0++) {
      grid[i1 * n + i0] = g(i0, i1);
    }
  }               
}