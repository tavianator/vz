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

#ifndef VZ_ADAPTIVE_HPP
#define VZ_ADAPTIVE_HPP

#include <cmath>

namespace vZ
{
  // Base class for adaptive RK-style algorithms
  template <typename Y>
  class GenericAdaptiveIntegrator : public GenericRKIntegrator<Y>
  {
  public:
    typedef typename GenericRKIntegrator<Y>::Scalar   Scalar;
    typedef typename GenericRKIntegrator<Y>::Function Function;

    GenericAdaptiveIntegrator& atol(Scalar tol) { m_atol = tol; return *this; }
    GenericAdaptiveIntegrator& rtol(Scalar tol) { m_rtol = tol; return *this; }

    Scalar atol() const { return m_atol; }
    Scalar rtol() const { return m_rtol; }

  protected:
    typedef typename GenericRKIntegrator<Y>::ACoefficients ACoefficients;
    typedef typename GenericRKIntegrator<Y>::BCoefficients BCoefficients;
    typedef typename GenericRKIntegrator<Y>::KVector       KVector;

    GenericAdaptiveIntegrator(Function f, unsigned int order,
                              ACoefficients a, BCoefficients b,
                              BCoefficients bStar)
      : GenericRKIntegrator<Y>(f), m_order(order),
        m_a(a), m_b(b), m_bStar(bStar)
    { }
    virtual ~GenericAdaptiveIntegrator() { }

    void step();

  private:
    Scalar m_atol, m_rtol;
    unsigned int m_order;
    ACoefficients m_a;
    BCoefficients m_b, m_bStar;
  };

  // Type alias
  typedef GenericAdaptiveIntegrator<double> AdaptiveIntegrator;

  // Implementations

  template <typename Y>
  void
  GenericAdaptiveIntegrator<Y>::step()
  {
    static const Scalar S = Scalar(19)/Scalar(20); // Arbitrary saftey factor
    Scalar newH;
    Y y;

    // Attempt the integration step in a loop
    bool rejected = true;
    while (rejected) {
      KVector k = calculateK(m_a);
      y         = calculateY(k, m_b);
      Y yStar   = calculateY(k, m_bStar);

      // Get an error estimate
      using std::abs;
      using std::pow;
      Scalar delta = abs(y - yStar);
      Scalar scale = m_atol + std::max(abs(y), abs(this->y()))*m_rtol;

      newH = S*this->h()*pow(scale/delta, Scalar(1)/m_order);

      if (delta > scale) {
        // Reject the step
        this->h(newH);
      } else {
        rejected = false;
      }
    }

    // Update x and y
    this->y(y);
    this->x(this->x() + this->h());

    // Adjust the stepsize for the next iteration
    this->h(newH);
  }
}

#endif // VZ_ADAPTIVE_HPP
