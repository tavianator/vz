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

#ifndef VZ_SIMPLE_HPP
#define VZ_SIMPLE_HPP

namespace vZ
{
  // Base class for non-adaptive RK-style algorithms
  template <typename Y>
  class GenericSimpleIntegrator : public GenericRKIntegrator<Y>
  {
  public:
    typedef typename GenericRKIntegrator<Y>::Scalar   Scalar;
    typedef typename GenericRKIntegrator<Y>::Function Function;

  protected:
    typedef typename GenericRKIntegrator<Y>::ACoefficients ACoefficients;
    typedef typename GenericRKIntegrator<Y>::BCoefficients BCoefficients;
    typedef typename GenericRKIntegrator<Y>::KVector       KVector;

    GenericSimpleIntegrator(Function f, ACoefficients a, BCoefficients b)
      : GenericRKIntegrator<Y>(f), m_a(a), m_b(b) { }
    virtual ~GenericSimpleIntegrator() { }

    void step();

  private:
    ACoefficients m_a;
    BCoefficients m_b;
  };

  // Type alias
  typedef GenericSimpleIntegrator<double> SimpleIntegrator;

  // Implementations

  template <typename Y>
  void
  GenericSimpleIntegrator<Y>::step()
  {
    this->y(this->calculateY(this->calculateK(m_a), m_b));
    this->x(this->x() + this->h());
  }
}

#endif // VZ_SIMPLE_HPP
