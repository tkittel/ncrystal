#ifndef NCrystal_Rational_hh
#define NCrystal_Rational_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/core/NCDefs.hh"

namespace NCRYSTAL_NAMESPACE {

  //Exact rational number n/d, for exact arithmetic on values such as atom
  //coordinates (e.g. so the symmetry image 1-9/10 is exactly 1/10). Always
  //normalised (d>0, gcd(|n|,d)=1), and the magnitudes of n and d are limited
  //to 2^53 (kMaxValue), so both are exactly representable as doubles and
  //toDouble() is correctly rounded (a single IEEE division). Operations whose
  //exact result exceeds this limit throw CalcError, so results are never
  //inexact. Not meant for hot loops.

  class Rational final {
  public:
    using int_t = std::int64_t;
    static constexpr int_t kMaxValue = int_t(1) << 53;

    Rational() = default;//zero
    explicit Rational( int_t n );
    Rational( int_t numerator, int_t denominator );//throws if d=0 or too big

    int_t numerator() const noexcept { return m_n; }
    int_t denominator() const noexcept { return m_d; }
    bool isInteger() const noexcept { return m_d == 1; }
    double toDouble() const noexcept;//correctly rounded

    //Largest integer <= value, and value-floor() in [0,1):
    int_t floor() const noexcept;
    Rational wrappedToUnitInterval() const;

    Rational operator-() const noexcept;
    Rational operator+( const Rational& ) const;
    Rational operator-( const Rational& ) const;
    Rational operator*( const Rational& ) const;
    Rational operator/( const Rational& ) const;//throws for division by zero
    Rational& operator+=( const Rational& o ) { return *this = *this + o; }
    Rational& operator-=( const Rational& o ) { return *this = *this - o; }

    bool operator==( const Rational& o ) const noexcept;
    bool operator!=( const Rational& o ) const noexcept;
    bool operator<( const Rational& o ) const noexcept;
    bool operator>( const Rational& o ) const noexcept;
    bool operator<=( const Rational& o ) const noexcept;
    bool operator>=( const Rational& o ) const noexcept;

  private:
    int_t m_n = 0;
    int_t m_d = 1;
    struct NoCheck {};
    Rational( NoCheck, int_t n, int_t d ) noexcept : m_n(n), m_d(d) {}
  };

  std::ostream& operator<<( std::ostream&, const Rational& );//"n/d" or "n"

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  inline double Rational::toDouble() const noexcept
  {
    //Both values are exactly representable, so this is correctly rounded:
    return static_cast<double>( m_n ) / static_cast<double>( m_d );
  }

  inline Rational Rational::operator-() const noexcept
  {
    return Rational( NoCheck(), -m_n, m_d );
  }

  inline bool Rational::operator==( const Rational& o ) const noexcept
  {
    return m_n == o.m_n && m_d == o.m_d;//normalised
  }

  inline bool Rational::operator!=( const Rational& o ) const noexcept
  {
    return !( *this == o );
  }

  inline bool Rational::operator>( const Rational& o ) const noexcept
  {
    return o < *this;
  }

  inline bool Rational::operator<=( const Rational& o ) const noexcept
  {
    return !( o < *this );
  }

  inline bool Rational::operator>=( const Rational& o ) const noexcept
  {
    return !( *this < o );
  }

}

#endif
