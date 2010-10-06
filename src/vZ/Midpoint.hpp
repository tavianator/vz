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

#ifndef VZ_MIDPOINT_HPP
#define VZ_MIDPOINT_HPP

namespace vZ
{
  // Midpoint method
  //
  // Second order Runge-Kutta method
  // Two function evaluations per step
  // Its tableau is:
  //
  //   0  |
  //   1/2|1/2
  //   ---+-----
  //      |0   1
  //
  //    k1    = dt*f(y[n])
  //    k2    = dt*f(y[n] + (dt/2)*k1)
  // y[n + 1] = y[n] + k2
  template <typename Y>
  class GenericMidpointIntegrator : public GenericSimpleIntegrator<Y>
  {
  public:
    typedef typename GenericSimpleIntegrator<Y>::Scalar   Scalar;
    typedef typename GenericSimpleIntegrator<Y>::Function Function;

    GenericMidpointIntegrator(Function f)
      : GenericSimpleIntegrator<Y>(f, s_a, s_b) { }
    ~GenericMidpointIntegrator() { }

  private:
    typedef typename GenericSimpleIntegrator<Y>::ACoefficients ACoefficients;
    typedef typename GenericSimpleIntegrator<Y>::BCoefficients BCoefficients;

    static ACoefficients s_a;
    static BCoefficients s_b;

    static Scalar s_a2Arr[1];
    static std::vector<Scalar> s_aArr[1];
    static Scalar s_bArr[2];
  };

  // Type alias
  typedef GenericMidpointIntegrator<double> MidpointIntegrator;

  // Implementation

  template <typename Y>
  typename GenericMidpointIntegrator<Y>::Scalar
  GenericMidpointIntegrator<Y>::s_a2Arr[1] = {
    Scalar(1)/2
  };

  template <typename Y>
  std::vector<typename GenericMidpointIntegrator<Y>::Scalar>
  GenericMidpointIntegrator<Y>::s_aArr[1] = {
    std::vector<Scalar>(s_a2Arr, s_a2Arr + 1)
  };

  template <typename Y>
  typename GenericMidpointIntegrator<Y>::ACoefficients
  GenericMidpointIntegrator<Y>::s_a(s_aArr, s_aArr + 1);

  template <typename Y>
  typename GenericMidpointIntegrator<Y>::Scalar
  GenericMidpointIntegrator<Y>::s_bArr[2] = {
    Scalar(0), Scalar(1)
  };

  template <typename Y>
  typename GenericMidpointIntegrator<Y>::BCoefficients
  GenericMidpointIntegrator<Y>::s_b(s_bArr, s_bArr + 2);
}

#endif // VZ_MIDPOINT_HPP
