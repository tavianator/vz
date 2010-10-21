/*************************************************************************
 * Copyright (C) 2010 Tavian Barnes <tavianator@gmail.com>               *
 *                                                                       *
 * This file is part of The vZ Test Suite.                               *
 *                                                                       *
 * The vZ Test Suite is free software; you can redistribute it and/or    *
 * modify it under the terms of the GNU General Public License as        *
 * published by the Free Software Foundation; either version 3 of the    *
 * License, or (at your option) any later version.                       *
 *                                                                       *
 * The vZ Test Suite is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 * General Public License for more details.                              *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 *************************************************************************/

#include "vZ.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

// y' = y (y == C*exp(x))
double
f(double x, double y)
{
  return y;
}

int
main()
{
  vZ::RKF45Integrator integrator(f);
  integrator.tol(1e-6)
            .y(1.0)
            .x(0.0)
            .h(0.06);

  integrator.integrate(2.0);

  double actual   = integrator.y();
  double expected = std::exp(2.0);

  std::cout << std::setprecision(10)
            << "Numerical:  " << actual << std::endl
            << "Expected:   " << expected  << std::endl
            << "h:          " << integrator.h() << std::endl
            << "Iterations: " << integrator.iterations() << std::endl
            << "Rejections: " << integrator.rejections() << std::endl;

  double error = std::abs(expected - actual)/expected;
  if (error > 1.7e-6) {
    std::cerr << "Error:      " << 100.0*error << "%" << std::endl;
    return EXIT_FAILURE;
  } else {
    std::cout << "Error:      " << 100.0*error << "%" << std::endl;
    return EXIT_SUCCESS;
  }
}
