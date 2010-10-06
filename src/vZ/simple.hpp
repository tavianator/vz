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

#ifndef VZ_SIMPLE_HPP
#define VZ_SIMPLE_HPP

#include <vector>

namespace vZ
{
  // Base class for non-adaptive RK-style algorithms
  template <typename T>
  class GenericSimpleIntegrator : public GenericIntegrator<T>
  {
  protected:
    // Coefficients in the tableau representation of the RK algorithm
    typedef std::vector<std::vector<T> > ACoefficients;
    typedef std::vector<T>               BCoefficients;

  public:
    typedef typename GenericIntegrator<T>::Function Function;

    GenericSimpleIntegrator(Function f, T dt,
                            ACoefficients a, BCoefficients b)
      : GenericIntegrator<T>(f, dt), m_a(a), m_b(b) { }
    virtual ~GenericSimpleIntegrator() { }

  protected:
    virtual void step(T& t, T& dt);

  private:
    ACoefficients m_a;
    BCoefficients m_b;
  };

  // Type alias
  typedef GenericSimpleIntegrator<double> SimpleIntegrator;
}

#endif // VZ_SIMPLE_HPP
