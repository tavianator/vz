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

#ifndef VZ_VECTOR_HPP
#define VZ_VECTOR_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace vZ
{
  // An N-dimensional vector
  template <std::size_t N, typename T = double>
  class Vector
  {
  public:
    typedef typename Traits<T>::Scalar Scalar;

    Vector()              { }
    explicit Vector(T x)  { m_values[0] = x; }
    Vector(T x, T y)      { m_values[0] = x; m_values[1] = y; }
    Vector(T x, T y, T z) { m_values[0] = x; m_values[1] = y; m_values[2] = z; }
    // Vector(const Vector& v);
    // ~Vector();

    // Vector& operator=(const Vector& v);

    // Component access

    T&       operator[](std::size_t i)       { return m_values[i]; }
    const T& operator[](std::size_t i) const { return m_values[i]; }

    T x() const { return m_values[0]; }
    T y() const { return m_values[1]; }
    T z() const { return m_values[2]; }

    // Operators
    inline Vector& operator+=(const Vector& rhs);
    inline Vector& operator-=(const Vector& rhs);
    inline Vector& operator*=(Scalar rhs);
    inline Vector& operator/=(Scalar rhs);

  private:
    T m_values[N];
  };

  // Disallow 0-sized Vectors
  template <typename T>
  class Vector<0, T>;

  // Traits specialization
  template <std::size_t N, typename T>
  class Traits<Vector<N, T> >
  {
  public:
    typedef typename Traits<T>::Scalar Scalar;

  private:
    Traits();
  };

  // Unary operators

  template <std::size_t N, typename T>
  inline Vector<N, T>
  operator+(const Vector<N, T>& rhs)
  {
    return rhs;
  }

  template <std::size_t N, typename T>
  inline Vector<N, T>
  operator-(const Vector<N, T>& rhs)
  {
    Vector<N, T> res;
    for (std::size_t i = 0; i < N; ++i) {
      res[i] = -rhs[i];
    }
    return res;
  }

  template <std::size_t N, typename T>
  inline typename Vector<N, T>::Scalar
  norm(const Vector<N, T>& v)
  {
    using std::sqrt;
    return sqrt(dot(v, v));
  }

  template <std::size_t N, typename T>
  inline typename Vector<N, T>::Scalar
  abs(const Vector<N, T>& v)
  {
    return norm(v);
  }

  // Binary operators

  template <std::size_t N, typename T>
  inline Vector<N, T>
  operator+(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
  {
    Vector<N, T> res = lhs;
    res += rhs;
    return res;
  }

  template <std::size_t N, typename T>
  inline Vector<N, T>
  operator-(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
  {
    Vector<N, T> res = lhs;
    res -= rhs;
    return res;
  }

  template <std::size_t N, typename T>
  inline Vector<N, T>
  operator*(typename Vector<N, T>::Scalar lhs,
            const Vector<N, T>& rhs)
  {
    Vector<N, T> res = rhs;
    res *= lhs;
    return res;
  }

  template <std::size_t N, typename T>
  inline Vector<N, T>
  operator*(const Vector<N, T>& lhs,
            typename Vector<N, T>::Scalar rhs)
  {
    Vector<N, T> res = lhs;
    res *= rhs;
    return res;
  }

  template <std::size_t N, typename T>
  inline Vector<N, T>
  operator/(const Vector<N, T>& lhs,
            typename Vector<N, T>::Scalar rhs)
  {
    Vector<N, T> res = lhs;
    res /= rhs;
    return res;
  }

  template <std::size_t N, typename T>
  inline typename Vector<N, T>::Scalar
  dot(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
  {
    typename Vector<N, T>::Scalar res(0);
    for (std::size_t i = 0; i < N; ++i) {
      res += lhs[i]*rhs[i];
    }
    return res;
  }

  template <typename T>
  inline Vector<3, T>
  cross(const Vector<3, T>& lhs, const Vector<3, T>& rhs)
  {
    return Vector<3, T>(lhs.y()*rhs.z() - lhs.z()*rhs.y(),
                        lhs.z()*rhs.x() - lhs.x()*rhs.z(),
                        lhs.x()*rhs.y() - lhs.y()*rhs.x());
  }

  // Implementation

  template <std::size_t N, typename T>
  inline Vector<N, T>&
  Vector<N, T>::operator+=(const Vector<N, T>& rhs)
  {
    for (std::size_t i = 0; i < N; ++i) {
      m_values[i] += rhs.m_values[i];
    }
    return *this;
  }

  template <std::size_t N, typename T>
  inline Vector<N, T>&
  Vector<N, T>::operator-=(const Vector<N, T>& rhs)
  {
    for (std::size_t i = 0; i < N; ++i) {
      m_values[i] -= rhs.m_values[i];
    }
    return *this;
  }

  template <std::size_t N, typename T>
  inline Vector<N, T>&
  Vector<N, T>::operator*=(typename Vector<N, T>::Scalar rhs)
  {
    for (std::size_t i = 0; i < N; ++i) {
      m_values[i] *= rhs;
    }
    return *this;
  }

  template <std::size_t N, typename T>
  inline Vector<N, T>&
  Vector<N, T>::operator/=(typename Vector<N, T>::Scalar rhs)
  {
    for (std::size_t i = 0; i < N; ++i) {
      m_values[i] /= rhs;
    }
    return *this;
  }
}

#endif // VZ_VECTOR_HPP
