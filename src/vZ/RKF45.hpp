/*************************************************************************
 * Copyright (C) 2010 Tavian Barnes <tavianator@gmail.com>               *
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

#ifndef VZ_RKF45_HPP
#define VZ_RKF45_HPP

namespace vZ
{
  // Runge-Kutta-Fehlberg method
  //
  // Fifth-order with embedded fourth-order
  // Its tableau is:
  //
  //   0     |
  //   1/4   |  1/4
  //   3/8   |  3/32       9/34
  //   12/13 |  1932/2197 -7200/2197  7296/2197
  //   1     |  439/216   -8          3680/513   -845/4104
  //   1/2   | -8/27       2         -3544/2565   1859/4104   -11/40
  //   ------+-----------------------------------------------------------
  //   b     | 25/216      0          1408/2565   2197/4104   -1/5   0
  //   b*    | 16/135      0          6656/12825  28561/56430 -9/50  2/55
  template <typename Y>
  class GenericRKF45Integrator : public GenericAdaptiveIntegrator<Y>
  {
  public:
    typedef typename GenericAdaptiveIntegrator<Y>::Scalar   Scalar;
    typedef typename GenericAdaptiveIntegrator<Y>::Function Function;

    GenericRKF45Integrator(Function f)
      : GenericAdaptiveIntegrator<Y>(f, 5, s_a, s_b, s_bStar) { }
    ~GenericRKF45Integrator() { }

  private:
    typedef typename GenericAdaptiveIntegrator<Y>::ACoefficients ACoefficients;
    typedef typename GenericAdaptiveIntegrator<Y>::BCoefficients BCoefficients;

    static ACoefficients s_a;
    static BCoefficients s_b;
    static BCoefficients s_bStar;

    static Scalar s_a2Arr[1];
    static Scalar s_a3Arr[2];
    static Scalar s_a4Arr[3];
    static Scalar s_a5Arr[4];
    static Scalar s_a6Arr[5];
    static std::vector<Scalar> s_aArr[5];
    static Scalar s_bArr[6];
    static Scalar s_bStarArr[6];
  };

  // Type alias
  typedef GenericRKF45Integrator<double> RKF45Integrator;

  // Implementation

  template <typename Y>
  typename GenericRKF45Integrator<Y>::Scalar
  GenericRKF45Integrator<Y>::s_a2Arr[1] = {
    Scalar(1)/Scalar(4)
  };
  template <typename Y>
  typename GenericRKF45Integrator<Y>::Scalar
  GenericRKF45Integrator<Y>::s_a3Arr[2] = {
    Scalar(3)/Scalar(32), Scalar(9)/Scalar(32)
  };
  template <typename Y>
  typename GenericRKF45Integrator<Y>::Scalar
  GenericRKF45Integrator<Y>::s_a4Arr[3] = {
     Scalar(1932)/Scalar(2197),
    -Scalar(7200)/Scalar(2197),
     Scalar(7296)/Scalar(2197)
  };
  template <typename Y>
  typename GenericRKF45Integrator<Y>::Scalar
  GenericRKF45Integrator<Y>::s_a5Arr[4] = {
     Scalar(439)/Scalar(216),
    -Scalar(8),
     Scalar(3680)/Scalar(513),
    -Scalar(845)/Scalar(4104),
  };
  template <typename Y>
  typename GenericRKF45Integrator<Y>::Scalar
  GenericRKF45Integrator<Y>::s_a6Arr[5] = {
    -Scalar(8)/Scalar(27),
     Scalar(2),
    -Scalar(3544)/Scalar(2565),
     Scalar(1859)/Scalar(4104),
    -Scalar(11)/Scalar(40)
  };

  template <typename Y>
  std::vector<typename GenericRKF45Integrator<Y>::Scalar>
  GenericRKF45Integrator<Y>::s_aArr[5] = {
    std::vector<Scalar>(s_a2Arr, s_a2Arr + 1),
    std::vector<Scalar>(s_a3Arr, s_a3Arr + 2),
    std::vector<Scalar>(s_a4Arr, s_a4Arr + 3),
    std::vector<Scalar>(s_a5Arr, s_a5Arr + 4),
    std::vector<Scalar>(s_a6Arr, s_a6Arr + 5)
  };

  template <typename Y>
  typename GenericRKF45Integrator<Y>::ACoefficients
  GenericRKF45Integrator<Y>::s_a(s_aArr, s_aArr + 5);

  template <typename Y>
  typename GenericRKF45Integrator<Y>::Scalar
  GenericRKF45Integrator<Y>::s_bArr[6] = {
     Scalar(16)/Scalar(135),
     Scalar(0),
     Scalar(6656)/Scalar(12825),
     Scalar(28561)/Scalar(56430),
    -Scalar(9)/Scalar(50),
     Scalar(2)/Scalar(55),
  };

  template <typename Y>
  typename GenericRKF45Integrator<Y>::BCoefficients
  GenericRKF45Integrator<Y>::s_b(s_bArr, s_bArr + 6);

  template <typename Y>
  typename GenericRKF45Integrator<Y>::Scalar
  GenericRKF45Integrator<Y>::s_bStarArr[6] = {
     Scalar(25)/Scalar(216),
     Scalar(0),
     Scalar(1408)/Scalar(2565),
     Scalar(2197)/Scalar(4104),
    -Scalar(1)/Scalar(5),
     Scalar(0)
  };

  template <typename Y>
  typename GenericRKF45Integrator<Y>::BCoefficients
  GenericRKF45Integrator<Y>::s_bStar(s_bStarArr, s_bStarArr + 6);
}

#endif // VZ_RKF45_HPP
