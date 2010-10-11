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

typedef vZ::EquationSystem<2, vZ::Vector<2> > Y;

// y'' = -(|y'|^2/|y|)*(y/|y|) (circular motion)
Y
f(double t, Y y)
{
  Y r;
  r[0] = y[1];
  r[1] = -(dot(y[1], y[1])*y[0])/(dot(y[0], y[0]));
  return r;
}

int
main()
{
  Y y;
  y[0] = vZ::Vector<2>(0.0, 1.0);
  y[1] = vZ::Vector<2>(8*std::atan(1.0), 0.0);
  vZ::GenericDP45Integrator<Y> integrator(f);
  integrator.tol(1e-6)
            .y(y)
            .x(0.0)
            .h(0.06);

  integrator.integrate(2.0);

  vZ::Vector<2> actual = integrator.y()[0];
  vZ::Vector<2> expected(0.0, 1.0);

  std::cout << std::setprecision(10)
            << "Numerical:  " << actual << std::endl
            << "Expected:   " << expected  << std::endl
            << "h:          " << integrator.h() << std::endl
            << "Iterations: " << integrator.iterations() << std::endl
            << "Rejections: " << integrator.rejections() << std::endl;

  double error = norm(expected - actual)/norm(expected);
  if (error > 4.0e-6) {
    std::cerr << "Error:      " << 100.0*error << "%" << std::endl;
    return EXIT_FAILURE;
  } else {
    std::cout << "Error:      " << 100.0*error << "%" << std::endl;
    return EXIT_SUCCESS;
  }
}
