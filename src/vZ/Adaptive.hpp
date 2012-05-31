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

    GenericAdaptiveIntegrator& tol(Scalar tol)
      { m_atol = tol; m_rtol = tol; return *this; }
    GenericAdaptiveIntegrator& atol(Scalar tol) { m_atol = tol; return *this; }
    GenericAdaptiveIntegrator& rtol(Scalar tol) { m_rtol = tol; return *this; }

    Scalar atol() const { return m_atol; }
    Scalar rtol() const { return m_rtol; }

    unsigned int rejections() const { return m_rejections; }

  protected:
    typedef typename GenericRKIntegrator<Y>::ACoefficients ACoefficients;
    typedef typename GenericRKIntegrator<Y>::BCoefficients BCoefficients;
    typedef typename GenericRKIntegrator<Y>::KVector       KVector;

    GenericAdaptiveIntegrator(Function f, unsigned int order,
                              ACoefficients a, BCoefficients b,
                              BCoefficients bStar);
    virtual ~GenericAdaptiveIntegrator() { }

    void step();

  private:
    Scalar m_atol, m_rtol;
    unsigned int m_order;
    unsigned int m_rejections;
    ACoefficients m_a;
    BCoefficients m_b, m_bStar;

    bool m_fsal, m_k1Set;
    Y m_k1;
  };

  // Type alias
  typedef GenericAdaptiveIntegrator<double> AdaptiveIntegrator;

  // Implementations

  template <typename Y>
  GenericAdaptiveIntegrator<Y>::GenericAdaptiveIntegrator(
    Function f, unsigned int order,
    ACoefficients a, BCoefficients b, BCoefficients bStar
  )
    : GenericRKIntegrator<Y>(f), m_order(order), m_rejections(0),
      m_a(a), m_b(b), m_bStar(bStar), m_fsal(true), m_k1Set(false)
  {
    std::size_t i;
    for (i = 0; i < m_a.back().size(); ++i) {
      if (m_a.back()[i] != m_b.at(i)) {
        m_fsal = false;
        return;
      }
    }

    if (m_b.at(i) != Scalar(0))
      m_fsal = false;
  }

  template <typename Y>
  void
  GenericAdaptiveIntegrator<Y>::step()
  {
    static const Scalar S = Scalar(19)/Scalar(20); // Arbitrary saftey factor
    Scalar newH = this->h();
    Y y;
    KVector k;

    // Attempt the integration step in a loop
    while (true) {
      if (m_k1Set) {
        k = this->calculateK(m_k1, y, m_a);
      } else if (m_fsal) {
        k = this->calculateK(y, m_a);
      } else {
        k = this->calculateK(m_a);
        y = this->calculateY(k, m_b);
      }
      Y yStar = this->calculateY(k, m_bStar);

      // Get an error estimate

      using std::abs;
      using std::pow;
      Scalar delta = abs(y - yStar);

      if (delta == Scalar(0)) {
        break;
      }

      Scalar scale = m_atol + std::max(abs(y), abs(this->y()))*m_rtol;
      newH         = S*this->h()*pow(scale/delta, Scalar(1)/m_order);

      if (delta > scale) {
        // Reject the step
        this->h(newH);
        ++m_rejections;
      } else {
        break;
      }
    }

    // Update x and y
    this->y(y);
    this->x(this->x() + this->h());

    // Handle FSAL optimization
    if (m_fsal) {
      m_k1Set = true;
      m_k1    = k.back();
    }

    // Adjust the stepsize for the next iteration
    this->h(newH);
  }
}

#endif // VZ_ADAPTIVE_HPP
