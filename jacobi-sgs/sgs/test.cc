#define BOOST_TEST_MODULE tolerance_02
#include <boost/test/included/unit_test.hpp>
#include <iostream>

#include "sgs-omp-bw.hh"
#include "sgs-omp-fw.hh"

#include "vanilla/sgs-bw.hh"
#include "vanilla/sgs-fw.hh"

#include "utilities.hh"

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

BOOST_AUTO_TEST_CASE(test, *utf::tolerance(0.02)) {
  // framework::master_test_suite().argv[0] 
  int n = 200; // 0 -> n = n+1
  int k = 10;
  int blocksize;
  int iterations = 200;

  int stencil_size = 2 * k + 1;
  double *grid_omp = new (std::align_val_t(64)) double[n * n];
  double *grid_van = new (std::align_val_t(64)) double[n * n];
  double *alpha = new double[stencil_size * stencil_size];
  coeff_stencil(k, alpha);
  fill_array(n, grid_omp);
  fill_array(n, grid_van);

  for (int i = 0; i < iterations; i++) {
    forward_gauss_sidel(n, k, grid_van, alpha);
    backward_gauss_sidel(n, k, grid_van, alpha);
  }

  for (int i = 0; i < iterations; i++) {
    forward_gauss_sidel_omp(n, k, grid_omp, alpha);
    backward_gauss_sidel_omp(n, k, grid_omp, alpha);
  }

  for (int i = 0; i < n*n; i++) {
    BOOST_TEST(grid_omp[i] == grid_van[i], tt::tolerance(0.2));
  }
  //   BOOST_TEST(x == y); // irrelevant by default

  //   BOOST_TEST(x == z); // relevant by default
  //   BOOST_TEST(x == z, tt::tolerance(0.1));
}
