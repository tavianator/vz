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
  // All integration methods derrive from this class
  // If the initial value problem is specified as
  //   y' = f(t, y); y(t0) = y0
  // then an Integrator could be constructed as Integrator(f, dt).y(y0).t(t0)
  template <typename T>
  class GenericIntegrator
  {
  public:
    typedef std::tr1::function<T (T, T)> Function;

    // By default, y and t start at zero
    GenericIntegrator(Function f, T dt)
      : m_f(f), m_y(0), m_t(0), m_dt(dt) { }
    virtual ~GenericIntegrator() { }

    GenericIntegrator& y(T y)   { m_y  = y;  return *this; }
    GenericIntegrator& t(T t)   { m_t  = t;  return *this; }
    GenericIntegrator& dt(T dt) { m_dt = dt; return *this; }

    T y()  const { return m_y; }
    T t()  const { return m_t; }
    T dt() const { return m_dt; }

    // Integrate until time t
    void integrate(T t_final);

  protected:
    virtual void step(T& t, T& dt) = 0;
    Function m_f;

  private:
    T m_y;
    T m_t, m_dt;
  };

  // Type alias
  typedef GenericIntegrator<double> Integrator;

  // Implementations

  template <typename T>
  void
  GenericIntegrator<T>::integrate(T t_final)
  {
    while (m_t < t_final) {
      m_dt = std::min(m_dt, t_final - m_t);
      step(m_t, m_dt);
    }
  }
}


#endif // VZ_INTEGRATOR_HPP
