#include <cmath>
#include <iostream>
#include <vector>

#include "sgs-omp-bw.hh"
#include "sgs-omp-fw.hh"
#include "vanilla/sgs-bw.hh"
#include "vanilla/sgs-fw.hh"

#include "time_experiment.hh"
#include "utilities.hh"

using NUMBER = double;
const int N = 32 * 1024 * 1024; // problem size

class SGS_Sequential { //: public Experiment {
  int n, k, iterations;
  double *grid = new (std::align_val_t(64)) double[n * n];
  // double *local_mean = new (std::align_val_t(64)) double[n * n];
  double *alpha = new double[(2 * k + 1) * (2 * k + 1)];

public:
  // construct an experiment
  SGS_Sequential(int n_, int k_, int iterations)
      : n(n_), k(k_), iterations(iterations) {
    coeff_stencil(k, alpha);
    fill_array(n, grid);
  }
  // run an experiment; can be called several times
  void run() const {
    for (int i = 0; i < iterations; i++) {
      forward_gauss_sidel_omp(n, k, grid, alpha);
      backward_gauss_sidel_omp(n, k, grid, alpha);
    }
  }
  // report number of operations
  double operations() const {
    // return (n * n) * ((std::pow(2 * k + 1, 2) - 1) * (2 * k + 1));
    return iterations * std::pow((n * n) * std::pow(2 * k + 1, 2) + 1, 2);
  } // no operations
  ~SGS_Sequential() {
    delete[] grid;
    delete[] alpha;
  }
};
class SGS_Vanilla { //: public Experiment {
  int n, k, iterations;
  double *grid = new (std::align_val_t(64)) double[n * n];
  // double *local_mean = new (std::align_val_t(64)) double[n * n];
  double *alpha = new double[(2 * k + 1) * (2 * k + 1)];

public:
  // construct an experiment
  SGS_Vanilla(int n_, int k_, int iterations)
      : n(n_), k(k_), iterations(iterations) {
    coeff_stencil(k, alpha);
    fill_array(n, grid);
  }
  // run an experiment; can be called several times
  void run() const {
    for (int i = 0; i < iterations; i++) {
      forward_gauss_sidel(n, k, grid, alpha);
      backward_gauss_sidel(n, k, grid, alpha);
    }
  }
  // report number of operations
  double operations() const {
    return iterations * std::pow((n * n) * std::pow(2 * k + 1, 2) + 1, 2);
  } // no operations
  ~SGS_Vanilla() {
    delete[] grid;
    delete[] alpha;
  }
};

// ---------------------------------------------------------------------------------

class SGS_OMP { //: public Experiment {
  int n, k, iterations;
  double *grid = new (std::align_val_t(64)) double[n * n];
  // double *local_mean = new (std::align_val_t(64)) double[n * n];
  double *alpha = new double[(2 * k + 1) * (2 * k + 1)];

public:
  // construct an experiment
  SGS_OMP(int n_, int k_, int iterations)
      : n(n_), k(k_), iterations(iterations) {
    coeff_stencil(k, alpha);
    fill_array(n, grid);
    // fill_array(n, local_mean);
  }
  // run an experiment; can be called several times
  void run() const {
    for (int i = 0; i < iterations; i++) {
      forward_gauss_sidel_omp(n, k, grid, alpha);
      backward_gauss_sidel_omp(n, k, grid, alpha);
    }
  }
  // report number of operations
  double operations() const {
    return iterations * std::pow((n * n) * std::pow(2 * k + 1, 2) + 1, 2);
  } // no operations
  ~SGS_OMP() {
    delete[] grid;
    // delete[] local_mean;
    delete[] alpha;
  }
};
// ---------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
void benchmark(int k, int iterations = 100) {
  // std::cout << N * sizeof(NUMBER) / 1024 / 1024 << " MByte per vector"
  //           << std::endl;
  double time_factor = 1e6;
  std::vector<int> sizes = {10,  20,  40,   60,   100,  200, 300,
                            500, 800, 1000, 1500, 2000, 2500};

  std::cout << "Stencil Size: " << k << ", " << std::endl;

  std::cout
      << "experiment, n, time (us), iterations, repetitions, Gflops/s, GByte/s"
      << std::endl;

  // ---------------------------------------------------------------------------------
  // std::cout << "\n**** SEQ ****\n";
  // for (auto n : sizes) {

  //   SGS_Sequential e(n, k, iterations);
  //   auto d = time_experiment(e);
  //   double flops = d.first * e.operations() / d.second * 1e6 / 1e9;

  //   std::cout << "sequential, " << n << ", " << d.second / time_factor << ",
  //   "
  //             << iterations << ", " << d.first << ", " << flops << ", "
  //             << flops * sizeof(NUMBER) << std::endl;
  // }
  // std::cout << "\n\n----Experiment Over---- \n\n";

  // ---------------------------------------------------------------------------------

  // std::cout << "**** OMP ****\n";

  // for (auto n : sizes) {

  //   SGS_OMP e(n, k, iterations);
  //   auto d = time_experiment(e);
  //   double flops = d.first * e.operations() / d.second * 1e6 / 1e9;

  //   std::cout << "omp parllel, " << n << ", " << d.second / time_factor << ",
  //   "
  //             << iterations << ", " << d.first << ", " << flops << ", "
  //             << flops * sizeof(NUMBER) << std::endl;
  // }
  // std::cout << "\n\n----Experiment Over---- \n\n";

  // ---------------------------------------------------------------------------------

  // ---------------------------------------------------------------------------------

  // std::cout << "**** Vanilla ****\n";

  for (auto n : sizes) {

    SGS_Vanilla e(n, k, iterations);
    auto d = time_experiment(e);
    double flops = d.first * e.operations() / d.second * 1e6 / 1e9;

    std::cout << "Vanilla, " << n << ", " << d.second / time_factor << ", "
              << iterations << ", " << d.first << ", " << flops << ", "
              << flops * sizeof(NUMBER) << std::endl;
  }
  std::cout << "\n\n----Experiment Over---- \n\n";

  // ---------------------------------------------------------------------------------
}

int main(int argc, char **argv) {
  int k, iterations;
  sscanf(argv[1], "%d", &k);
  sscanf(argv[2], "%d", &iterations);
  benchmark(k, iterations);
}