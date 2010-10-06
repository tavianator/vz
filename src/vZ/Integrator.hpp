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

#ifndef VZ_INTEGRATOR_HPP
#define VZ_INTEGRATOR_HPP

#include <tr1/functional>
#include <algorithm>

namespace vZ
{
  // Base Integrator class
  //
  // All integration methods derive from this class
  // If the initial value problem is specified as
  //   y' = f(x, y); y(x0) = y0
  // then an Integrator could be constructed as Integrator(f, dt).y(y0).x(x0)
  template <typename Y>
  class GenericIntegrator
  {
  public:
    typedef typename Traits<Y>::Scalar Scalar;
    typedef std::tr1::function<Y (Scalar, Y)> Function;

    // By default, y and t start at zero, h starts UNDEFINED
    GenericIntegrator(Function f)
      : m_f(f), m_y(0), m_x(0), m_h() { }
    virtual ~GenericIntegrator() { }

    GenericIntegrator& y(Y y)      { m_y = y; return *this; }
    GenericIntegrator& x(Scalar x) { m_x = x; return *this; }
    GenericIntegrator& h(Scalar h) { m_h = h; return *this; }

    Y      y() const { return m_y; }
    Scalar x() const { return m_x; }
    Scalar h() const { return m_h; }

    // Integrate until x == x_final
    void integrate(Scalar x_final);

  protected:
    virtual void step() = 0;

    const Function& f() const { return m_f; }

  private:
    Function m_f;
    Y m_y;
    Scalar m_x, m_h;
  };

  // Type alias
  typedef GenericIntegrator<double> Integrator;

  // Implementations

  template <typename Y>
  void
  GenericIntegrator<Y>::integrate(Scalar x_final)
  {
    while (m_x < x_final) {
      m_h = std::min(m_h, x_final - m_x);
      step();
    }
  }
}


#endif // VZ_INTEGRATOR_HPP
