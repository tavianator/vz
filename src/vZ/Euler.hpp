/*************************************************************************
 * Copyright (C) 2009-2010 Tavian Barnes <tavianator@gmail.com>          *
 *                                                                       *
 * This file is part of The vZ Library.                                  *
 *                                                                       *
 * The vZ Library is free software; you can redistribute it and/or       *
 * modify it under the terms of the GNU Lesser General Public License as *
 * published by the Free Software Foundation; either version 3 of the    *
 * License, or (at your option) any later version.                       *
 *                                                                       *
 * The vZ Library is distributed in the hope that it will be useful, but *
 * WITHOUT ANY WARRANTY; without even the implied warranty of            *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 * Lesser General Public License for more details.                       *
 *                                                                       *
 * You should have received a copy of the GNU Lesser General Public      *
 * License along with this program.  If not, see                         *
 * <http://www.gnu.org/licenses/>.                                       *
 *************************************************************************/

#ifndef VZ_EULER_HPP
#define VZ_EULER_HPP

namespace vZ
{
  // Euler method
  //
  // Simplest Runge-Kutta method
  // First order
  // Its tableau is:
  //
  //   0 |
  //   --+--
  //     | 1
  //
  // y[n + 1] = y[n] + dt*f(y[n])
  template <typename Y>
  class GenericEulerIntegrator : public GenericSimpleIntegrator<Y>
  {
  public:
    typedef typename GenericSimpleIntegrator<Y>::Scalar   Scalar;
    typedef typename GenericSimpleIntegrator<Y>::Function Function;

    GenericEulerIntegrator(Function f)
      : GenericSimpleIntegrator<Y>(f, s_a, s_b) { }
    ~GenericEulerIntegrator() { }

  private:
    typedef typename GenericSimpleIntegrator<Y>::ACoefficients ACoefficients;
    typedef typename GenericSimpleIntegrator<Y>::BCoefficients BCoefficients;

    static ACoefficients s_a;
    static BCoefficients s_b;

    static Scalar s_bArr[1];
  };

  // Type alias
  typedef GenericEulerIntegrator<double> EulerIntegrator;

  // Implementation

  template <typename Y>
  typename GenericEulerIntegrator<Y>::Scalar
  GenericEulerIntegrator<Y>::s_bArr[1] = {
    Scalar(1)
  };

  template <typename Y>
  typename GenericEulerIntegrator<Y>::ACoefficients
  GenericEulerIntegrator<Y>::s_a;

  template <typename Y>
  typename GenericEulerIntegrator<Y>::BCoefficients
  GenericEulerIntegrator<Y>::s_b(s_bArr, s_bArr + 1);
}

#endif // VZ_EULER_HPP
