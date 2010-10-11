#include "vZ.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

typedef vZ::EquationSystem<2> Y;

// y'' = y (y == C*exp(t))
Y
f(double t, Y y)
{
  Y r;
  r[0] = y[1];
  r[1] = y[0];
  return r;
}

int
main()
{
  Y y;
  y[0] = 1.0;
  y[1] = 1.0;
  vZ::GenericDP45Integrator<Y> integrator(f);
  integrator.tol(1e-6).y(y).x(0.0).h(0.06);

  integrator.integrate(2.0);

  double actual   = integrator.y()[0];
  double expected = std::exp(2.0);

  std::cout << std::setprecision(10)
            << "Numerical:  " << actual << std::endl
            << "Expected:   " << expected  << std::endl
            << "h:          " << integrator.h() << std::endl
            << "Iterations: " << integrator.iterations() << std::endl
            << "Rejections: " << integrator.rejections() << std::endl;

  double error = std::fabs(expected - actual)/expected;
  if (error > 6.0e-7) {
    std::cerr << "Error:      " << 100.0*error << "%" << std::endl;
    return EXIT_FAILURE;
  } else {
    std::cout << "Error:      " << 100.0*error << "%" << std::endl;
    return EXIT_SUCCESS;
  }
}
