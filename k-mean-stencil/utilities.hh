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