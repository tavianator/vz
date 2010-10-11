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
  vZ::GenericEulerIntegrator<Y> integrator(f);
  integrator.y(y).x(0.0).h(0.01);

  integrator.integrate(2.0);

  double actual   = integrator.y()[0];
  double expected = std::exp(2.0);

  std::cout << std::setprecision(10)
            << "Numerical:  " << actual << std::endl
            << "Expected:   " << expected  << std::endl
            << "Iterations: " << integrator.iterations() << std::endl;

  double error = std::fabs(expected - actual)/expected;
  if (error > 0.01) {
    std::cerr << "Error:      " << 100.0*error << "%" << std::endl;
    return EXIT_FAILURE;
  } else {
    std::cout << "Error:      " << 100.0*error << "%" << std::endl;
    return EXIT_SUCCESS;
  }
}
