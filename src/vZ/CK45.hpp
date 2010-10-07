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

#ifndef VZ_CK45_HPP
#define VZ_CK45_HPP

namespace vZ
{
  // Cash-Karp method
  //
  // Fifth-order with embedded fourth-order
  // Its tableau is:
  //
  //   0     |
  //   1/5   |  1/5
  //   3/10  |  3/40       9/40
  //   3/5   |  3/10      -9/10     6/5
  //   1     | -11/54      5/2     -70/27       35/27
  //   7/8   | 1631/55296  175/512  575/13824   44275/110592 253/4096
  //   ------+-----------------------------------------------------------------
  //   b     | 37/378      0        250/621     125/594      0         512/1771
  //   b*    | 2825/27648  0        18575/48384 13525/55296  277/14336 1/4
  template <typename Y>
  class GenericCK45Integrator : public GenericAdaptiveIntegrator<Y>
  {
  public:
    typedef typename GenericAdaptiveIntegrator<Y>::Scalar   Scalar;
    typedef typename GenericAdaptiveIntegrator<Y>::Function Function;

    GenericCK45Integrator(Function f)
      : GenericAdaptiveIntegrator<Y>(f, 5, s_a, s_b, s_bStar) { }
    ~GenericCK45Integrator() { }

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
  typedef GenericCK45Integrator<double> CK45Integrator;

  // Implementation

  template <typename Y>
  typename GenericCK45Integrator<Y>::Scalar
  GenericCK45Integrator<Y>::s_a2Arr[1] = {
    Scalar(1)/Scalar(5)
  };
  template <typename Y>
  typename GenericCK45Integrator<Y>::Scalar
  GenericCK45Integrator<Y>::s_a3Arr[2] = {
    Scalar(3)/Scalar(40), Scalar(9)/Scalar(40)
  };
  template <typename Y>
  typename GenericCK45Integrator<Y>::Scalar
  GenericCK45Integrator<Y>::s_a4Arr[3] = {
     Scalar(3)/Scalar(10),
    -Scalar(9)/Scalar(10),
     Scalar(6)/Scalar(5)
  };
  template <typename Y>
  typename GenericCK45Integrator<Y>::Scalar
  GenericCK45Integrator<Y>::s_a5Arr[4] = {
    -Scalar(11)/Scalar(54),
     Scalar(5)/Scalar(2),
    -Scalar(70)/Scalar(27),
     Scalar(35)/Scalar(27),
  };
  template <typename Y>
  typename GenericCK45Integrator<Y>::Scalar
  GenericCK45Integrator<Y>::s_a6Arr[5] = {
    Scalar(1631)/Scalar(55296),
    Scalar(175)/Scalar(512),
    Scalar(575)/Scalar(13824),
    Scalar(44275)/Scalar(110592),
    Scalar(253)/Scalar(4096)
  };

  template <typename Y>
  std::vector<typename GenericCK45Integrator<Y>::Scalar>
  GenericCK45Integrator<Y>::s_aArr[5] = {
    std::vector<Scalar>(s_a2Arr, s_a2Arr + 1),
    std::vector<Scalar>(s_a3Arr, s_a3Arr + 2),
    std::vector<Scalar>(s_a4Arr, s_a4Arr + 3),
    std::vector<Scalar>(s_a5Arr, s_a5Arr + 4),
    std::vector<Scalar>(s_a6Arr, s_a6Arr + 5)
  };

  template <typename Y>
  typename GenericCK45Integrator<Y>::ACoefficients
  GenericCK45Integrator<Y>::s_a(s_aArr, s_aArr + 5);

  template <typename Y>
  typename GenericCK45Integrator<Y>::Scalar
  GenericCK45Integrator<Y>::s_bArr[6] = {
    Scalar(37)/Scalar(378),
    Scalar(0),
    Scalar(250)/Scalar(621),
    Scalar(125)/Scalar(594),
    Scalar(0),
    Scalar(512)/Scalar(1771),
  };

  template <typename Y>
  typename GenericCK45Integrator<Y>::BCoefficients
  GenericCK45Integrator<Y>::s_b(s_bArr, s_bArr + 6);

  template <typename Y>
  typename GenericCK45Integrator<Y>::Scalar
  GenericCK45Integrator<Y>::s_bStarArr[6] = {
    Scalar(2825)/Scalar(27648),
    Scalar(0),
    Scalar(18575)/Scalar(48384),
    Scalar(13525)/Scalar(55296),
    Scalar(277)/Scalar(14336),
    Scalar(1)/Scalar(4)
  };

  template <typename Y>
  typename GenericCK45Integrator<Y>::BCoefficients
  GenericCK45Integrator<Y>::s_bStar(s_bStarArr, s_bStarArr + 6);
}

#endif // VZ_CK45_HPP
