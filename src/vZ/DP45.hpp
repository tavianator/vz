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

#ifndef VZ_DP45_HPP
#define VZ_DP45_HPP

namespace vZ
{
  // Dormand-Prince method
  //
  // Fifth-order with embedded fourth-order
  // Its tableau is:
  //
  //   0    |
  //   1/5  | 1/5
  //   3/10 | 3/40        9/40
  //   4/5  | 44/45      -56/15      32/9
  //   8/9  | 19372/6561 -25360/2186 64448/6561 -212/729
  //   1    | 9017/3168  -355/33     46732/5247  49/176   -5103/18656
  //   1    | 35/384      0          500/1113    125/192  -2187/6784    11/84
  //   -----+------------------------------------------------------------------------
  //   b    | 35/384      0          500/1113    125/192  -2187/6784    11/84    0
  //   b*   | 5179/57600  0          7571/16695  393/640  -92097/339200 172/2100 1/40
  template <typename Y>
  class GenericDP45Integrator : public GenericAdaptiveIntegrator<Y>
  {
  public:
    typedef typename GenericAdaptiveIntegrator<Y>::Scalar   Scalar;
    typedef typename GenericAdaptiveIntegrator<Y>::Function Function;

    GenericDP45Integrator(Function f)
      : GenericAdaptiveIntegrator<Y>(f, 5, s_a, s_b, s_bStar) { }
    ~GenericDP45Integrator() { }

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
    static Scalar s_a7Arr[6];
    static std::vector<Scalar> s_aArr[6];
    static Scalar s_bArr[7];
    static Scalar s_bStarArr[7];
  };

  // Type alias
  typedef GenericDP45Integrator<double> DP45Integrator;

  // Implementation

  template <typename Y>
  typename GenericDP45Integrator<Y>::Scalar
  GenericDP45Integrator<Y>::s_a2Arr[1] = {
    Scalar(1)/Scalar(5)
  };
  template <typename Y>
  typename GenericDP45Integrator<Y>::Scalar
  GenericDP45Integrator<Y>::s_a3Arr[2] = {
    Scalar(3)/Scalar(40), Scalar(9)/Scalar(40)
  };
  template <typename Y>
  typename GenericDP45Integrator<Y>::Scalar
  GenericDP45Integrator<Y>::s_a4Arr[3] = {
     Scalar(44)/Scalar(45), -Scalar(56)/Scalar(15), Scalar(32)/Scalar(9)
  };
  template <typename Y>
  typename GenericDP45Integrator<Y>::Scalar
  GenericDP45Integrator<Y>::s_a5Arr[4] = {
     Scalar(19372)/Scalar(6561),
    -Scalar(25360)/Scalar(2187),
     Scalar(64448)/Scalar(6561),
    -Scalar(212)/Scalar(729),
  };
  template <typename Y>
  typename GenericDP45Integrator<Y>::Scalar
  GenericDP45Integrator<Y>::s_a6Arr[5] = {
     Scalar(9017)/Scalar(3168),
    -Scalar(355)/Scalar(33),
     Scalar(46732)/Scalar(5247),
     Scalar(49)/Scalar(176),
    -Scalar(5103)/Scalar(18656)
  };
  template <typename Y>
  typename GenericDP45Integrator<Y>::Scalar
  GenericDP45Integrator<Y>::s_a7Arr[6] = {
     Scalar(35)/Scalar(384),
     Scalar(0),
     Scalar(500)/Scalar(1113),
     Scalar(125)/Scalar(192),
    -Scalar(2187)/Scalar(6784),
     Scalar(11)/Scalar(84)
  };

  template <typename Y>
  std::vector<typename GenericDP45Integrator<Y>::Scalar>
  GenericDP45Integrator<Y>::s_aArr[6] = {
    std::vector<Scalar>(s_a2Arr, s_a2Arr + 1),
    std::vector<Scalar>(s_a3Arr, s_a3Arr + 2),
    std::vector<Scalar>(s_a4Arr, s_a4Arr + 3),
    std::vector<Scalar>(s_a5Arr, s_a5Arr + 4),
    std::vector<Scalar>(s_a6Arr, s_a6Arr + 5),
    std::vector<Scalar>(s_a7Arr, s_a7Arr + 6)
  };

  template <typename Y>
  typename GenericDP45Integrator<Y>::ACoefficients
  GenericDP45Integrator<Y>::s_a(s_aArr, s_aArr + 6);

  template <typename Y>
  typename GenericDP45Integrator<Y>::Scalar
  GenericDP45Integrator<Y>::s_bArr[7] = {
     Scalar(35)/Scalar(384),
     Scalar(0),
     Scalar(500)/Scalar(1113),
     Scalar(125)/Scalar(192),
    -Scalar(2187)/Scalar(6784),
     Scalar(11)/Scalar(84),
     Scalar(0)
  };

  template <typename Y>
  typename GenericDP45Integrator<Y>::BCoefficients
  GenericDP45Integrator<Y>::s_b(s_bArr, s_bArr + 7);

  template <typename Y>
  typename GenericDP45Integrator<Y>::Scalar
  GenericDP45Integrator<Y>::s_bStarArr[7] = {
     Scalar(5179)/Scalar(57600),
     Scalar(0),
     Scalar(7571)/Scalar(16695),
     Scalar(393)/Scalar(640),
    -Scalar(92097)/Scalar(339200),
     Scalar(187)/Scalar(2100),
     Scalar(1)/Scalar(40)
  };

  template <typename Y>
  typename GenericDP45Integrator<Y>::BCoefficients
  GenericDP45Integrator<Y>::s_bStar(s_bStarArr, s_bStarArr + 7);
}

#endif // VZ_DP45_HPP
