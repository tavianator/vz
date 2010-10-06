#include "vZ.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>

// y' = y (y == C*exp(t))
double
f(double t, double y)
{
  return y;
}

int
main()
{
  vZ::RK4Integrator integrator(f);
  integrator.y(1.0).x(0.0).h(0.04);

  integrator.integrate(2.0);

  std::cout << integrator.y() << std::endl
            << std::exp(2.0)  << std::endl;

  return EXIT_SUCCESS;
}
