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

#ifndef VZ_BS23_HPP
#define VZ_BS23_HPP

namespace vZ
{
  // Heun-Euler method
  //
  // Second-order with embedded first-order
  // Its tableau is:
  //
  //   0 |
  //   1 | 1
  //   --+---------
  //     | 1/2  1/2
  //     |  1    0
  //
  //    k1    = dt*f(y[n])
  //    k2    = dt*f(y[n] + dt*k1)
  // y[n + 1] = y[n] + 1/2*(k1 + k2)
  template <typename Y>
  class GenericBS23Integrator : public GenericAdaptiveIntegrator<Y>
  {
  public:
    typedef typename GenericAdaptiveIntegrator<Y>::Scalar   Scalar;
    typedef typename GenericAdaptiveIntegrator<Y>::Function Function;

    GenericBS23Integrator(Function f)
      : GenericAdaptiveIntegrator<Y>(f, 3, s_a, s_b, s_bStar) { }
    ~GenericBS23Integrator() { }

  private:
    typedef typename GenericAdaptiveIntegrator<Y>::ACoefficients ACoefficients;
    typedef typename GenericAdaptiveIntegrator<Y>::BCoefficients BCoefficients;

    static ACoefficients s_a;
    static BCoefficients s_b;
    static BCoefficients s_bStar;

    static Scalar s_a2Arr[1];
    static Scalar s_a3Arr[2];
    static Scalar s_a4Arr[3];
    static std::vector<Scalar> s_aArr[3];
    static Scalar s_bArr[4];
    static Scalar s_bStarArr[4];
  };

  // Type alias
  typedef GenericBS23Integrator<double> BS23Integrator;

  // Implementation

  template <typename Y>
  typename GenericBS23Integrator<Y>::Scalar
  GenericBS23Integrator<Y>::s_a2Arr[1] = {
    Scalar(1)/Scalar(2)
  };
  template <typename Y>
  typename GenericBS23Integrator<Y>::Scalar
  GenericBS23Integrator<Y>::s_a3Arr[2] = {
    Scalar(0), Scalar(3)/Scalar(4)
  };
  template <typename Y>
  typename GenericBS23Integrator<Y>::Scalar
  GenericBS23Integrator<Y>::s_a4Arr[3] = {
    Scalar(2)/Scalar(9), Scalar(1)/Scalar(3), Scalar(4)/Scalar(9)
  };

  template <typename Y>
  std::vector<typename GenericBS23Integrator<Y>::Scalar>
  GenericBS23Integrator<Y>::s_aArr[3] = {
    std::vector<Scalar>(s_a2Arr, s_a2Arr + 1),
    std::vector<Scalar>(s_a3Arr, s_a3Arr + 2),
    std::vector<Scalar>(s_a4Arr, s_a4Arr + 3)
  };

  template <typename Y>
  typename GenericBS23Integrator<Y>::ACoefficients
  GenericBS23Integrator<Y>::s_a(s_aArr, s_aArr + 3);

  template <typename Y>
  typename GenericBS23Integrator<Y>::Scalar
  GenericBS23Integrator<Y>::s_bArr[4] = {
    Scalar(2)/Scalar(9),
    Scalar(1)/Scalar(3),
    Scalar(4)/Scalar(9),
    Scalar(0)
  };

  template <typename Y>
  typename GenericBS23Integrator<Y>::BCoefficients
  GenericBS23Integrator<Y>::s_b(s_bArr, s_bArr + 4);

  template <typename Y>
  typename GenericBS23Integrator<Y>::Scalar
  GenericBS23Integrator<Y>::s_bStarArr[4] = {
    Scalar(7)/Scalar(24),
    Scalar(1)/Scalar(4),
    Scalar(1)/Scalar(3),
    Scalar(1)/Scalar(8)
  };

  template <typename Y>
  typename GenericBS23Integrator<Y>::BCoefficients
  GenericBS23Integrator<Y>::s_bStar(s_bStarArr, s_bStarArr + 4);
}

#endif // VZ_BS23_HPP
