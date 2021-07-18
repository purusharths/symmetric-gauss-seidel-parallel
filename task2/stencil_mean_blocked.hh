#pragma once
#include <algorithm>
#include <iostream>

void stencil_mean_blocked(int n, int k, double *grid, double *local_mean, int blocksize = 4) {
  //   double max_points = (2 * k +goo 1) * (2 * k + 1);

  for (int i = 0; i < n; i = i + (blocksize)) {
    for (int j = 0; j < n; j = j + (blocksize)) {

      // std::cout << " \nBlock: " << points << "\n";
      for (int blockRow = i; blockRow < i + blocksize; blockRow++) {
        for (int blockCol = j; blockCol < j + blocksize; blockCol++) {
          double sum = 0, points = 0;
        //   std::cout << "(" << blockRow << ", " << blockCol << "): ";
          for (int mm = std::max(blockCol - k, 0);
               mm <= std::min(blockCol + k, n - 1); mm++) {
            for (int ll = std::max(blockRow - k, 0);
                 ll <= std::min(blockRow + k, n - 1); ll++) {
              sum += grid[ll * n + mm];
              // std::cout << ll << ", " << mm << " | ";
              points++;
            }
          }
          //   std::cout << sum << std::endl;
          local_mean[blockCol * n + blockRow] =
              static_cast<double>(sum / 1 );//points); // /
          sum = 0;
        }
      }
    }
    // std::cout << std::endl;
  }
}
//
// std::cout << "Block " << points++ << ": (" << std::max(i - k, 0) << ", "
//     << std::min(i + k, n - 1) << ") -> "
//     << std::min(i + blocksize + k, n - 1) << ", "
//     << std::min(i + blocksize + k, n - 1) << "\n";