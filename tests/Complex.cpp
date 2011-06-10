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
#include <complex>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

// y' = x*y (y == C*exp(x^2/2))
std::complex<double>
f(double x, std::complex<double> y)
{
  return x*y;
}

int
main()
{
  vZ::GenericDP45Integrator<std::complex<double> > integrator(f);
  integrator.tol(1e-6)
            .y(std::complex<double>(1.0, 0.0))
            .x(0.0)
            .h(0.06);

  integrator.integrate(2.0);

  std::complex<double> actual = integrator.y();
  std::complex<double> expected(std::exp(2.0), 0.0);

  std::cout << std::setprecision(10)
            << "Numerical:  " << actual << std::endl
            << "Expected:   " << expected  << std::endl
            << "h:          " << integrator.h() << std::endl
            << "Iterations: " << integrator.iterations() << std::endl
            << "Rejections: " << integrator.rejections() << std::endl;

  double error = std::abs(expected - actual)/std::abs(expected);
  if (error > 5.9e-7 || !std::isfinite(error)) {
    std::cerr << "Error:      " << 100.0*error << "%" << std::endl;
    return EXIT_FAILURE;
  } else {
    std::cout << "Error:      " << 100.0*error << "%" << std::endl;
    return EXIT_SUCCESS;
  }
}
