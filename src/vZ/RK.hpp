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

#ifndef VZ_RK_HPP
#define VZ_RK_HPP

#include <vector>

namespace vZ
{
  // Base class for Runge-Kutta type algorithms
  template <typename Y>
  class GenericRKIntegrator : public GenericIntegrator<Y>
  {
  public:
    typedef typename GenericIntegrator<Y>::Scalar   Scalar;
    typedef typename GenericIntegrator<Y>::Function Function;

  protected:
    // Coefficients in the tableau representation of the RK algorithm
    typedef std::vector<std::vector<Scalar> > ACoefficients;
    typedef std::vector<Scalar>               BCoefficients;

    // Result vectors
    typedef std::vector<Y> KVector;

    GenericRKIntegrator(Function f) : GenericIntegrator<Y>(f) { }
    virtual ~GenericRKIntegrator() { }

    // Perform the stages of an RK integration
    KVector calculateK(const ACoefficients& a) const;
    Y calculateY(const KVector& k, const BCoefficients& b) const;
  };

  // Type alias
  typedef GenericRKIntegrator<double> RKIntegrator;

  // Implementation

  template <typename Y>
  typename GenericRKIntegrator<Y>::KVector
  GenericRKIntegrator<Y>::calculateK(const ACoefficients& a) const
  {
    KVector k;
    k.reserve(a.size() + 1);

    // k1
    k.push_back(this->f()(this->x(), this->y()));

    // k2..n
    for (typename ACoefficients::const_iterator i = a.begin();
         i != a.end();
         ++i)
    {
      Scalar c(0);
      Y y = this->y();
      for (typename std::vector<Scalar>::size_type j = 0; j < i->size(); ++j) {
        Scalar aij = i->at(j);
        c += aij;

        y += aij*this->h()*k.at(j);
      }

      k.push_back(this->f()(this->x() + c, y));
    }

    return k;
  }

  template <typename Y>
  Y
  GenericRKIntegrator<Y>::calculateY(const KVector& k, const BCoefficients& b)
    const
  {
    Y y = this->y();

    for (typename std::vector<Scalar>::size_type i = 0; i < k.size(); ++i) {
      y += this->h()*b.at(i)*k.at(i);
    }

    return y;
  }
}

#endif // VZ_RK_HPP
