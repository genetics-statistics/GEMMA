#include <catch.hpp>
#include <iostream>
#include "gsl/gsl_matrix.h"
#include "mathfunc.h"
#include <algorithm>
#include <limits>
#include <numeric>

using namespace std;

TEST_CASE( "Math functions", "[math]" ) {
  double data[] = { 2,-1, 0,
                   -1, 2,-1,
                    0,-1, 2};
  gsl_matrix *m = gsl_matrix_alloc(3,3);
  copy(data, data+9, m->data);
  REQUIRE( isMatrixPositiveDefinite(m) );
  REQUIRE( isMatrixSymmetric(m) );
  // REQUIRE( checkMatrixEigen(m,0.001) );

  double data1[] = {1.0,0,0,
                    0,3.0,0,
                    0,0,2.0};
  copy(data1, data1+9, m->data);
  REQUIRE( isMatrixPositiveDefinite(m) );
  // REQUIRE( checkMatrixEigen(m) );

  double data2[] = {1,1,1,
                    1,1,1,
                    1,1,0.5};
  copy(data2, data2+9, m->data);
  REQUIRE( !isMatrixPositiveDefinite(m));
  // REQUIRE( !checkMatrixEigen(m) );

  double data3[] = {1.0,  0,  0,
                    3.0,3.0,  0,
                      0,  0,2.0};
  copy(data3, data3+9, m->data);
  REQUIRE( !isMatrixPositiveDefinite(m) );
  REQUIRE( !isMatrixSymmetric(m) );
  // REQUIRE( checkMatrixEigen(m) );

  // ---- NaN checks
  vector<double> v = {1.0, 2.0};
  REQUIRE (!std::isnan(std::accumulate(v.begin(), v.end(), 0)));
  vector<double> v2 = {1.0, 2.0, std::numeric_limits<double>::quiet_NaN()};
  REQUIRE (std::isnan(v2[2]));
  REQUIRE(has_nan(v2));
  // test minus nan
  vector<double> v3 = {1.0, 2.0, -std::numeric_limits<double>::quiet_NaN()};
  REQUIRE (std::isnan(v3[2]));
  REQUIRE(has_nan(v3));
}
