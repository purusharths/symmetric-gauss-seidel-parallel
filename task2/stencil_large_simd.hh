#pragma once
#include "vcl/vectorclass.h"

#include <algorithm>
#include <iostream>
#include <vcl/vectorf128.h>

void stencil_simd_large_vec(int n, int k, double *grid, double *local_mean,
                            int blocksize = 4) {
                              // std::cout << "SIMD 8d" << std::endl;
  if (k < 4) {
    std::cout << "WARNING: The function call uses 8 vector units and may give "
                 "unambigious results for smaller stencil size \n ";
  }
  for (int i = 0; i < n; i = i + (blocksize)) {
    for (int j = 0; j < n; j = j + (blocksize)) {

      for (int blockRow = i; blockRow < i + blocksize; blockRow++) {
        for (int blockCol = j; blockCol < j + blocksize; blockCol++) {

          double sum = 0, points = 0;
          int col_max = std::min(blockCol + k, n - 1);
          int col_min = std::max(blockCol - k, 0);
          int row_max = std::min(blockRow + k, n - 1);
          int row_min = std::max(blockRow - k, 0);
          // int m = (std::min(blockCol + k, n - 1) - std::max(blockCol - k,
          // 0));
          int m = col_max - col_min;
          // int o = (std::min(blockRow + k, n - 1) - std::max(blockRow - k,
          // 0));
          int o = row_max - row_min;
          int mn = (m + 1) * (o + 1);

          int datasize = o + 1;
          int vectorsize = 8;
          int regularpart = datasize & (-vectorsize);

          int ll = std::max(blockRow - k, 0);
          // for (int mm = std::max(blockCol - k, 0); mm <= std::min(blockCol +
          // k, n - 1); mm++) {
          for (int mm = col_min; mm <= col_max; mm++) {
            int a;
            Vec8d sum1(0);

            for (a = 0; a < regularpart; a += vectorsize) {
              Vec8d temp(grid[((row_min + a) * n) + mm],
                         grid[((row_min + (a + 1)) * n) + mm],
                         grid[((row_min + (a + 2)) * n) + mm],
                         grid[((row_min + (a + 3)) * n) + mm],
                         grid[((row_min + (a + 4)) * n) + mm],
                         grid[((row_min + (a + 5)) * n) + mm],
                         grid[((row_min + (a + 6)) * n) + mm],
                         grid[((row_min + (a + 7)) * n) + mm]);
              sum1 += temp;
            }

            if (datasize - a >= 4) {
              // get four more numbers
              Vec4d sum2(grid[((row_min + (a + 8)) * n) + mm],
                         grid[((row_min + (a + 9)) * n) + mm],
                         grid[((row_min + (a + 10)) * n) + mm],
                         grid[((row_min + (a + 11)) * n) + mm]);
              // sum2.load(mydata + i);
              a += 4;
              sum += horizontal_add(sum2);
            }

            for (; a < datasize; a++) {
              sum += grid[(row_min + a) * n + mm];
            }
            sum += horizontal_add(sum1);
          }
          local_mean[blockCol * n + blockRow] =
              static_cast<double>(sum / mn); // /mn
        }
      }
    }
  }
}
