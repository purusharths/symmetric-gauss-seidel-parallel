#pragma once
#include <cmath>
#include <iostream>
#include <omp.h>
#include <thread>
#include <vector>

void division_among_threads(int n, int total_threads,
                            std::vector<std::tuple<int, int>> &loop_domains) {
  const int SIZE = n;
  int p = total_threads;
  int p_size = SIZE > p ? p : SIZE;

  std::vector<int> iter(p_size);

  for (int i = 0; i < SIZE; i++) {
    iter[i % p_size] += 1;
  }
  //   std::reverse(iter.begin(), iter.end());
  //   std::vector<std::tuple<int, int>> loop_domains;

  int val_start = 0;
  int val_end = 0;

  for (unsigned int j = 0; j < iter.size(); ++j) {
    val_end += iter[j];
    // printf("[%d] = %d | %d, %d\n", j, (int)iter[j], val_start, val_end);
    loop_domains.emplace_back(val_start, val_end);
    val_start += iter[j];
  }

  //   for (const auto &i : loop_domains) {
  //     std::cout << std::get<0>(i) << ", " << std::get<1>(i) << ", " <<
  //     std::endl;
  //   }
}

void calculate(int n, int k, int p, int q, double *grid, double *local_mean) {
  double sum = 0;
  std::cout << std::this_thread::get_id() << std::endl;

  int row_min = std::max(p - k, 0);
  int row_max = std::min(p + k, n - 1);
  int col_min = std::max(q - k, 0);
  int col_max = std::min(q + k, n - 1);
  // int tid = omp_get_thread_num();

  // loop for k - box
  for (int mm = row_min; mm <= row_max; mm++) {
    for (int ll = col_min; ll <= col_max; ll++) {
      if (mm < p || (mm == p && ll < q)) { // mm >= k && ll > i) {
                                           // B+

        sum += local_mean[ll * n + mm];
      } else {
// A+
#pragma omp critical
        sum += grid[ll * n + mm];
      }
    }
  }

  local_mean[q * n + p] = static_cast<double>(sum); // / change to multiply
}

void threading_kernel(){

}

void forward_gauss_sidel_omp(int n, int k, double *grid, double *local_mean) {
  int j = 0;
  unsigned int numThreads = std::thread::hardware_concurrency();
  for (int i = 0; i < n; i++) {       // i loop
    for (int c = 0; c < k + 1; c++) { // c loop

      double points = 1; // max_points;
        
      int p = i;
      int q = c, l;
      int loop_count = std::min(p, (n - 1 - q) / (k + 1)) + 1;
      std::vector<std::thread> threads;
      std::vector<std::tuple<int, int>> &loop_domains;
      division_among_threads(loop_count, numThreads, loop_domains)

          for (int l = 1; l <= loop_count; l++) {

        // std::cout << "(" << p << ", " << q << ") ";

        threads.push_back(std::thread{calculate, n, k, p, q, grid, local_mean});
        p = p - 1;
        q = q + (k + 1);
      }
      for (auto &t : threads) {
        t.join();
      }
      threads.clear();

      std::cout << "\n----" << std::endl;
    }
  }
}