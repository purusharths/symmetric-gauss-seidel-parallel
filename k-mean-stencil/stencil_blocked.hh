#pragma once
#include <algorithm>
#include <iostream>

void stencil_mean_blocked(int n, int k, double *grid, double *local_mean, int blocksize = 4) {
  //   double max_points = (2 * k +goo 1) * (2 * k + 1);
  // std::cout << " Blocked" << std::endl;
  for (int i = 0; i < n; i = i + (blocksize)) {
    for (int j = 0; j < n; j = j + (blocksize)) {

      // std::cout << " \nBlock: " << points << "\n";
      for (int blockRow = i; blockRow < i + blocksize; blockRow++) {
        for (int blockCol = j; blockCol < j + blocksize; blockCol++) {
          int col_max = std::min(blockCol + k, n - 1);
          int col_min = std::max(blockCol - k, 0);
          int row_max = std::min(blockRow + k, n - 1);
          int row_min = std::max(blockRow - k, 0);
          int m = col_max - col_min;
          int o = row_max - row_min;
          int mn = (m + 1) * (o + 1);

          double sum = 0;
        //   std::cout << "(" << blockRow << ", " << blockCol << "): ";
          for (int mm = col_min; mm <= col_max; mm++) {
            for (int ll = row_min; ll <= row_max; ll++) {
              sum += grid[ll * n + mm];
              // std::cout << ll << ", " << mm << " | ";
              // points++;
            }
          }
          //   std::cout << sum << std::endl;
          local_mean[blockCol * n + blockRow] =
              static_cast<double>(sum / mn );//points); // /
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