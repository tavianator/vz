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
  double tol = 1e-6;
  vZ::HE12Integrator integrator(f);
  integrator.atol(tol).rtol(tol).y(1.0).x(0.0).h(0.02);

  integrator.integrate(2.0);

  double actual   = integrator.y();
  double expected = std::exp(2.0);

  std::cout << std::setprecision(10)
            << "Numerical: " << actual << std::endl
            << "Expected:  " << expected  << std::endl
            << "h:         " << integrator.h() << std::endl;

  double error = std::fabs(expected - actual)/expected;
  if (error > 8.7e-7) {
    std::cerr << "Error:     " << 100.0*error << "%" << std::endl;
    return EXIT_FAILURE;
  } else {
    std::cout << "Error:     " << 100.0*error << "%" << std::endl;
    return EXIT_SUCCESS;
  }
}
