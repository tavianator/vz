#include "vZ.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

// y' = -y (y == C*exp(-t))
vZ::Vector<3>
f(double t, vZ::Vector<3> y)
{
  return -y;
}

int
main()
{
  vZ::GenericDP45Integrator<vZ::Vector<3> > integrator(f);
  integrator.tol(1e-6).y(vZ::Vector<3>(1.0, 1.0, 1.0)).x(0.0).h(0.06);

  integrator.integrate(2.0);

  double actual   = integrator.y().x();
  double expected = std::exp(-2.0);

  std::cout << std::setprecision(10)
            << "Numerical:  " << actual << std::endl
            << "Expected:   " << expected  << std::endl
            << "h:          " << integrator.h() << std::endl
            << "Iterations: " << integrator.iterations() << std::endl
            << "Rejections: " << integrator.rejections() << std::endl;

  double error = std::fabs(expected - actual)/expected;
  if (error > 1.5e-6) {
    std::cerr << "Error:      " << 100.0*error << "%" << std::endl;
    return EXIT_FAILURE;
  } else {
    std::cout << "Error:      " << 100.0*error << "%" << std::endl;
    return EXIT_SUCCESS;
  }
}
