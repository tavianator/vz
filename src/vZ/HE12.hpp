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

#ifndef VZ_HE12_HPP
#define VZ_HE12_HPP

namespace vZ
{
  // Heun-Euler method
  //
  // Second-order with embedded first-order
  // Its tableau is:
  //
  //   0  |
  //   1  | 1
  //   ---+---------
  //   b  | 1/2  1/2
  //   b* | 1    0
  template <typename Y>
  class GenericHE12Integrator : public GenericAdaptiveIntegrator<Y>
  {
  public:
    typedef typename GenericAdaptiveIntegrator<Y>::Scalar   Scalar;
    typedef typename GenericAdaptiveIntegrator<Y>::Function Function;

    GenericHE12Integrator(Function f)
      : GenericAdaptiveIntegrator<Y>(f, 2, s_a, s_b, s_bStar) { }
    ~GenericHE12Integrator() { }

  private:
    typedef typename GenericAdaptiveIntegrator<Y>::ACoefficients ACoefficients;
    typedef typename GenericAdaptiveIntegrator<Y>::BCoefficients BCoefficients;

    static ACoefficients s_a;
    static BCoefficients s_b;
    static BCoefficients s_bStar;

    static Scalar s_a2Arr[1];
    static std::vector<Scalar> s_aArr[1];
    static Scalar s_bArr[2];
    static Scalar s_bStarArr[2];
  };

  // Type alias
  typedef GenericHE12Integrator<double> HE12Integrator;

  // Implementation

  template <typename Y>
  typename GenericHE12Integrator<Y>::Scalar
  GenericHE12Integrator<Y>::s_a2Arr[1] = {
    Scalar(1)
  };

  template <typename Y>
  std::vector<typename GenericHE12Integrator<Y>::Scalar>
  GenericHE12Integrator<Y>::s_aArr[1] = {
    std::vector<Scalar>(s_a2Arr, s_a2Arr + 1)
  };

  template <typename Y>
  typename GenericHE12Integrator<Y>::ACoefficients
  GenericHE12Integrator<Y>::s_a(s_aArr, s_aArr + 1);

  template <typename Y>
  typename GenericHE12Integrator<Y>::Scalar
  GenericHE12Integrator<Y>::s_bArr[2] = {
    Scalar(1)/Scalar(2), Scalar(1)/Scalar(2)
  };

  template <typename Y>
  typename GenericHE12Integrator<Y>::BCoefficients
  GenericHE12Integrator<Y>::s_b(s_bArr, s_bArr + 2);

  template <typename Y>
  typename GenericHE12Integrator<Y>::Scalar
  GenericHE12Integrator<Y>::s_bStarArr[2] = {
    Scalar(1), Scalar(0)
  };

  template <typename Y>
  typename GenericHE12Integrator<Y>::BCoefficients
  GenericHE12Integrator<Y>::s_bStar(s_bStarArr, s_bStarArr + 2);
}

#endif // VZ_HE12_HPP
