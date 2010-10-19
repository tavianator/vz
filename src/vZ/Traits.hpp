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

#ifndef VZ_TRAITS_HPP
#define VZ_TRAITS_HPP

#include <complex>

namespace vZ
{
  // Traits class
  //
  // Specialize this class for non-scalar types which are used as template
  // arguments to the Generic* classes
  template <typename T>
  class Traits
  {
  public:
    typedef T Scalar;

  private:
    Traits();
  };

  // Specialization for std::complex<T>
  template <typename T>
  class Traits<std::complex<T> >
  {
  public:
    typedef typename Traits<T>::Scalar Scalar;

  private:
    Traits();
  };
}

#endif // VZ_TRAITS_HPP
