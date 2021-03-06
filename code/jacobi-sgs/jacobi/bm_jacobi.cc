#include <cmath>
#include <iostream>
#include <vector>

#include "jacobi.hh"
#include "time_experiment.hh"
#include "utilities.hh"

using NUMBER = double;
const int N = 32 * 1024 * 1024; // problem size

//  stencil_mean_vanilla(n, k, grid, local_mean);
// stencil_mean_blocked(n, k, grid, local_mean, blocksize);
// stencil_simd(n, k, grid, local_mean, blocksize);

class JacobiVanilla { //: public Experiment {
  int n, k, iterations;
  double *grid = new (std::align_val_t(64)) double[n * n];
  double *local_mean = new (std::align_val_t(64)) double[n * n];
  double *alpha = new double[(2 * k + 1) * (2 * k + 1)];

public:
  // construct an experiment
  JacobiVanilla(int n_, int k_, int iterations)
      : n(n_), k(k_), iterations(iterations) {
    coeff_stencil(k, alpha);
    fill_array(n, grid);
    fill_array(n, local_mean);
  }
  // run an experiment; can be called several times
  void run() const {
    jacobi_vanilla(n, k, iterations, grid, local_mean, alpha);
  }
  // report number of operations
  double operations() const {
    return (std::pow(2 * k + 1, 2) - 1) * 2 * (n * n) * iterations;
  } // no operations
  ~JacobiVanilla() {
    delete[] grid;
    delete[] local_mean;
    delete[] alpha;
  }
};
// ---------------------------------------------------------------------------------

class Jacobi_OMP { //: public Experiment {
  int n, k, iterations;
  double *grid = new (std::align_val_t(64)) double[n * n];
  double *local_mean = new (std::align_val_t(64)) double[n * n];
  double *alpha = new double[(2 * k + 1) * (2 * k + 1)];

public:
  // construct an experiment
  Jacobi_OMP(int n_, int k_, int iterations)
      : n(n_), k(k_), iterations(iterations) {
    coeff_stencil(k, alpha);
    fill_array(n, grid);
    fill_array(n, local_mean);
  }
  // run an experiment; can be called several times
  void run() const { jacobi_omp(n, k, iterations, grid, local_mean, alpha); }
  // report number of operations
  double operations() const {
    return (std::pow(2 * k + 1, 2) - 1) * 2 * (n * n) * iterations;
  } // no operations
  ~Jacobi_OMP() {
    delete[] grid;
    delete[] local_mean;
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
  std::vector<int> sizes = {100, 150, 200, 250, 300, 350, 400, 500};

  std::cout << "# Stencil Size: " << k << ", ";

  // std::cout << "\n**** Vanilla ****\n";
  std::cout << "experiment,n,time,repetitions,Gflops/s,GByte/s" << std::endl;
  for (auto n : sizes) {

    JacobiVanilla e(n, k, iterations);
    auto d = time_experiment(e);
    double flops = d.first * e.operations() / d.second * 1e6 / 1e9;

    std::cout << "vanilla, " << n << ", " << d.second / time_factor << ", "
              << d.first << ", " << flops << ", " << flops * sizeof(NUMBER)
              << std::endl;
  }
  // std::cout << "\n\n----Experiment Over---- \n\n";

  // ---------------------------------------------------------------------------------

  // std::cout << "**** OMP ****\n";

  for (auto n : sizes) {

    Jacobi_OMP e(n, k, iterations);
    auto d = time_experiment(e);
    double flops = d.first * e.operations() / d.second * 1e6 / 1e9;

    std::cout << "omp parllel, " << n << ", " << d.second / time_factor << ", "
              << d.first << ", " << flops << ", " << flops * sizeof(NUMBER)
              << std::endl;
  }
  // std::cout << "\n\n----Experiment Over---- \n\n";

  // ---------------------------------------------------------------------------------
}

int main(int argc, char **argv) {
  int k;
  if (argc < 2 || argc > 3) {
    std::cout << "usage: ./bm_jacobi <stencil size (k)>\n" << std::endl;
    exit(1);
  }
  sscanf(argv[1], "%d", &k);
  benchmark(k);
}
