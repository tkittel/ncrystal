#ifndef NCrystal_PolyBoundHelper_hh
#define NCrystal_PolyBoundHelper_hh

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

//Private helper of NCIsoMacroXSTable.cc: exact checks whether polynomials of
//low degree are non-negative on an interval.

#include "NCrystal/internal/utils/NCMath.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace XSTable {

    //Polynomial of degree at most 9, with its coefficients in the power basis
    //(lowest order first). This is enough for the checks of the table
    //interpolation (see NCIsoMacroXSTable.cc), and a product of degree above 9
    //throws a LogicError:
    class Poly9 final {
    public:
      static constexpr int maxDegree = 9;

      Poly9() = default;
      explicit Poly9( double c0 ) { m_c[0] = c0; }
      Poly9( double c0, double c1 ) : m_deg( 1 ) { m_c[0] = c0; m_c[1] = c1; }

      int degree() const { return m_deg; }
      double coeff( int i ) const { return i <= m_deg ? m_c[i] : 0.0; }
      double eval( double z ) const;

      Poly9& operator+=( const Poly9& );
      Poly9& operator-=( const Poly9& );
      Poly9& operator*=( double );
      friend Poly9 operator*( const Poly9&, const Poly9& );

    private:
      double m_c[maxDegree+1] = {};
      int m_deg = 0;
    };

    Poly9 operator+( Poly9, const Poly9& );
    Poly9 operator-( Poly9, const Poly9& );
    Poly9 operator*( Poly9, double );
    Poly9 operator*( const Poly9&, const Poly9& );

    //The cubic c[0]+c[1]*u+c[2]*u^2+c[3]*u^3, with u a polynomial:
    Poly9 composeCubic( const double * c, const Poly9& u );

    //The Bernstein coefficients b[0..9] on [0,1], of degree 9 whatever the
    //degree of the polynomial (raised to 9 where it is lower):
    void bernsteinCoefficients( const Poly9&, double * b );

    //Whether the polynomial with the Bernstein coefficients b[0..9] on [0,1] is
    //non-negative there. This is true if all coefficients are non-negative,
    //and false if b[0] or b[9] (the values at the end points) is
    //negative. Otherwise the interval is subdivided (de Casteljau). If this
    //does not decide it within 30 subdivisions, the result is false:
    bool bernsteinNonNegative( const double * b );

  }

}

#endif
