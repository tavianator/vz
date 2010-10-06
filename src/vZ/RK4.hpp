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

#ifndef VZ_RK4_HPP
#define VZ_RK4_HPP

namespace vZ
{
  // RK4 method
  //
  // Fourth order Runge-Kutta method
  // Four function evaluations per step
  // Its tableau is:
  //
  //   0  |
  //   1/2|1/2
  //   1/2|0   1/2
  //   1  |0   0   1
  //   ---+---------------
  //      |1/6 1/3 1/3 1/6
  //
  //    k1    = dt*f(y[n])
  //    k2    = dt*f(y[n] + (dt/2)*k1)
  //    k3    = dt*f(y[n] + (dt/2)*k2)
  //    k4    = dt*f(y[n] + dt*k3)
  // y[n + 1] = y[n] + (k1 + 2*k2 + 2*k3 + k4)/6
  template <typename Y>
  class GenericRK4Integrator : public GenericSimpleIntegrator<Y>
  {
  public:
    typedef typename GenericSimpleIntegrator<Y>::Scalar   Scalar;
    typedef typename GenericSimpleIntegrator<Y>::Function Function;

    GenericRK4Integrator(Function f)
      : GenericSimpleIntegrator<Y>(f, s_a, s_b) { }
    ~GenericRK4Integrator() { }

  private:
    typedef typename GenericSimpleIntegrator<Y>::ACoefficients ACoefficients;
    typedef typename GenericSimpleIntegrator<Y>::BCoefficients BCoefficients;

    static ACoefficients s_a;
    static BCoefficients s_b;

    static Scalar s_a2Arr[1];
    static Scalar s_a3Arr[2];
    static Scalar s_a4Arr[3];
    static std::vector<Scalar> s_aArr[3];
    static Scalar s_bArr[4];
  };

  // Type alias
  typedef GenericRK4Integrator<double> RK4Integrator;

  // Implementation

  template <typename Y>
  typename GenericRK4Integrator<Y>::Scalar
  GenericRK4Integrator<Y>::s_a2Arr[1] = {
    Scalar(1)/Scalar(2)
  };
  template <typename Y>
  typename GenericRK4Integrator<Y>::Scalar
  GenericRK4Integrator<Y>::s_a3Arr[2] = {
    Scalar(0), Scalar(1)/Scalar(2)
  };
  template <typename Y>
  typename GenericRK4Integrator<Y>::Scalar
  GenericRK4Integrator<Y>::s_a4Arr[3] = {
    Scalar(0), Scalar(0), Scalar(1)
  };

  template <typename Y>
  std::vector<typename GenericRK4Integrator<Y>::Scalar>
  GenericRK4Integrator<Y>::s_aArr[3] = {
    std::vector<Scalar>(s_a2Arr, s_a2Arr + 1),
    std::vector<Scalar>(s_a3Arr, s_a3Arr + 2),
    std::vector<Scalar>(s_a4Arr, s_a4Arr + 3)
  };

  template <typename Y>
  typename GenericRK4Integrator<Y>::ACoefficients
  GenericRK4Integrator<Y>::s_a(s_aArr, s_aArr + 3);

  template <typename Y>
  typename GenericRK4Integrator<Y>::Scalar
  GenericRK4Integrator<Y>::s_bArr[4] = {
    Scalar(1)/Scalar(6),
    Scalar(1)/Scalar(3),
    Scalar(1)/Scalar(3),
    Scalar(1)/Scalar(6)
  };

  template <typename Y>
  typename GenericRK4Integrator<Y>::BCoefficients
  GenericRK4Integrator<Y>::s_b(s_bArr, s_bArr + 4);
}

#endif // VZ_RK4_HPP
