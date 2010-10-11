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

#ifndef VZ_EQUATIONSYSTEM_HPP
#define VZ_EQUATIONSYSTEM_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace vZ
{
  // A class to easily represent a system of ODEs
  template <std::size_t N, typename T = double>
  class EquationSystem
  {
  public:
    typedef typename Traits<T>::Scalar Scalar;

    // EquationSystem();
    // ~EquationSystem();

    T&       operator[](std::size_t i)       { return m_values[i]; }
    const T& operator[](std::size_t i) const { return m_values[i]; }

    EquationSystem& operator+=(const EquationSystem& rhs);
    EquationSystem& operator-=(const EquationSystem& rhs);
    EquationSystem& operator*=(Scalar rhs);
    EquationSystem& operator/=(Scalar rhs);

  private:
    T m_values[N];
  };

  // Disallow 0-sized EquationSystems
  template <typename T>
  class EquationSystem<0, T>;

  // Traits specialization
  template <std::size_t N, typename T>
  class Traits<EquationSystem<N, T> >
  {
  public:
    typedef typename Traits<T>::Scalar Scalar;

  private:
    Traits();
  };

  // Binary operators

  template <std::size_t N, typename T>
  inline EquationSystem<N, T>
  operator+(const EquationSystem<N, T>& lhs, const EquationSystem<N, T>& rhs)
  {
    EquationSystem<N, T> res = lhs;
    res += rhs;
    return res;
  }

  template <std::size_t N, typename T>
  inline EquationSystem<N, T>
  operator-(const EquationSystem<N, T>& lhs, const EquationSystem<N, T>& rhs)
  {
    EquationSystem<N, T> res = lhs;
    res -= rhs;
    return res;
  }

  template <std::size_t N, typename T>
  inline EquationSystem<N, T>
  operator*(typename EquationSystem<N, T>::Scalar lhs,
            const EquationSystem<N, T>& rhs)
  {
    EquationSystem<N, T> res = rhs;
    res *= lhs;
    return res;
  }

  template <std::size_t N, typename T>
  inline EquationSystem<N, T>
  operator*(const EquationSystem<N, T>& lhs,
            typename EquationSystem<N, T>::Scalar rhs)
  {
    EquationSystem<N, T> res = lhs;
    res *= rhs;
    return res;
  }

  template <std::size_t N, typename T>
  inline EquationSystem<N, T>
  operator/(const EquationSystem<N, T>& lhs,
            typename EquationSystem<N, T>::Scalar rhs)
  {
    EquationSystem<N, T> res = lhs;
    res /= rhs;
    return res;
  }

  // Implementation

  template <std::size_t N, typename T>
  EquationSystem<N, T>&
  EquationSystem<N, T>::operator+=(const EquationSystem<N, T>& rhs)
  {
    for (std::size_t i = 0; i < N; ++i) {
      m_values[i] += rhs.m_values[i];
    }
    return *this;
  }

  template <std::size_t N, typename T>
  EquationSystem<N, T>&
  EquationSystem<N, T>::operator-=(const EquationSystem<N, T>& rhs)
  {
    for (std::size_t i = 0; i < N; ++i) {
      m_values[i] -= rhs.m_values[i];
    }
    return *this;
  }

  template <std::size_t N, typename T>
  EquationSystem<N, T>&
  EquationSystem<N, T>::operator*=(typename EquationSystem<N, T>::Scalar rhs)
  {
    for (std::size_t i = 0; i < N; ++i) {
      m_values[i] *= rhs;
    }
    return *this;
  }

  template <std::size_t N, typename T>
  EquationSystem<N, T>&
  EquationSystem<N, T>::operator/=(typename EquationSystem<N, T>::Scalar rhs)
  {
    for (std::size_t i = 0; i < N; ++i) {
      m_values[i] /= rhs;
    }
    return *this;
  }

  template <std::size_t N, typename T>
  typename EquationSystem<N, T>::Scalar
  abs(const EquationSystem<N, T>& es)
  {
    typename EquationSystem<N, T>::Scalar ret(0);
    for (std::size_t i = 0; i < N; ++i) {
      using std::abs;
      ret = std::max(ret, abs(es[i]));
    }
    return ret;
  }
}

#endif // VZ_EQUATIONSYSTEM_HPP
