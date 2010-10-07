#include "vZ.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

// y' = y (y == C*exp(t))
double
f(double t, double y)
{
  return y;
}

int
main()
{
  vZ::HeunIntegrator integrator(f);
  integrator.y(1.0).x(0.0).h(0.02);

  integrator.integrate(2.0);

  double actual   = integrator.y();
  double expected = std::exp(2.0);

  std::cout << std::setprecision(10)
            << "Numerical:  " << actual << std::endl
            << "Expected:   " << expected  << std::endl
            << "iterations: " << integrator.iterations() << std::endl;

  double error = std::fabs(expected - actual)/expected;
  if (error > 1.4e-4) {
    std::cerr << "Error:      " << 100.0*error << "%" << std::endl;
    return EXIT_FAILURE;
  } else {
    std::cout << "Error:      " << 100.0*error << "%" << std::endl;
    return EXIT_SUCCESS;
  }
}
