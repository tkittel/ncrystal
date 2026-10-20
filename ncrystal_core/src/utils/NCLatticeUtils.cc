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

#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCString.hh"
namespace NC = NCrystal;

#define NCRYSTAL_EXACT_LATTICEROTS_FOR_SPECIAL_CASES

NC::RotMatrix NC::getLatticeRot( double a, double b, double c,
                                 double alpha, double beta, double gamma )
{
  //sanity check angles in (0,pi) - also confirms that angles are not in degrees by mistake.
  nc_assert_always(alpha<kPi&&alpha>0);
  nc_assert_always(beta<kPi&&beta>0);
  nc_assert_always(gamma<kPi&&gamma>0);
  nc_assert_always(a>0);
  nc_assert_always(b>0);
  nc_assert_always(c>0);
  double cg(2), sg(2), ca(2), cb(2), sb(2);

#ifdef NCRYSTAL_EXACT_LATTICEROTS_FOR_SPECIAL_CASES
  //Explicitly look for angles 90 and 120 that are very common in crystals, so
  //we can avoid some numerical imprecisions in those cases and instead use
  //exact values for the sines and cosines. This is needed because
  //e.g. std::cos(kPiHalf)=6.1e-17 rather than 0.0 exactly, resulting in
  //incorrect tiny non-zero matrix elements where there should be zero. This
  //propagates further, for instance resulting in funny crystal normals like {
  //0, 1, -6.12323e-17 } instead of the correct { 0, 1, 0 }. This is mostly a
  //cosmetic issue of course.
  const double tol = 1e-14;
  bool alpha90 = ncabs(alpha-kPiHalf)<tol;
  bool beta90 = ncabs(beta-kPiHalf)<tol;
  bool gamma90 = ncabs(gamma-kPiHalf)<tol;
  bool gamma120 = ncabs(gamma-kPi*(2/3.))<tol;

  if (gamma90) { cg = 0.0; sg = 1.0; }
  else if (gamma120) { cg = -0.5; sg = 0.86602540378443864676372317075293618347140262690519; }
  if (alpha90) { ca = 0.0; }
  if (beta90) { cb = 0.0; sb = 1.0; }
  if (!ca&&!cb&&!cg) {
    double m[9] = { a, 0., 0.,
                    0., b, 0.,
                    0., 0., c };
    return RotMatrix(m);
  }
#endif

  if (cg==2) cg = std::cos(gamma);
  if (sg==2) sg = std::sin(gamma);
  if (ca==2) ca = std::cos(alpha);
  if (cb==2) cb = std::cos(beta);
  if (sb==2) sb = std::sin(beta);

  nc_assert(sg>0);
  const double m57 = c*(ca-cb*cg)/sg;
  double m[9] = { a, 0., 0.,
                  b*cg, b*sg, 0.,
                  c*cb, m57, 0. };//last entry changed below
  if ( !cb && !m57 ) {
    //Avoid the sqrt and potential introduction of numerical imprecision:
    m[8] = c;
  } else if ( cb && !m57 ) {
    //Avoid the sqrt and potential introduction of numerical imprecision:
    m[8] = c*sb;
  } else {
    //most general case
    double tmp = sb*sb - ncsquare( (ca-cb*cg)/sg );
    nc_assert( tmp >=0.0 );
    m[8]  = c * std::sqrt(tmp);
  }
  return RotMatrix(m);
}


NC::RotMatrix NC::getReciprocalLatticeRot( double a, double b, double c,
                                           double alpha, double beta, double gamma )
{
  //sanity check angles in (0,pi) - also confirms that angles are not in degrees by mistake.
  nc_assert_always(alpha<kPi&&alpha>0);
  nc_assert_always(beta<kPi&&beta>0);
  nc_assert_always(gamma<kPi&&gamma>0);
  nc_assert_always(a>0);
  nc_assert_always(b>0);
  nc_assert_always(c>0);
#ifdef NCRYSTAL_EXACT_LATTICEROTS_FOR_SPECIAL_CASES
  //Like in getLatticeRot we explicitly treat angles 90 and 120.
  const double tol = 1e-14;
  bool alpha90 = ncabs(alpha-kPiHalf)<tol;
  bool beta90 = ncabs(beta-kPiHalf)<tol;
  bool gamma90 = ncabs(gamma-kPiHalf)<tol;
  bool gamma120 = ncabs(gamma-kPi*(2/3.))<tol;
  if ( alpha90 && beta90 && gamma90 ) {
    //Usually this is cubic, tetragonal, orthorombic
    std::array<double,9> m{ k2Pi/a, 0., 0.,
                            0., k2Pi/b, 0.,
                            0., 0., k2Pi/c };
    return RotMatrix(m);
  } else if (alpha90 && beta90 && gamma120 ) {
    //Usually this is trigonal/hexagonal
    const double twopidivsqrt3 = 3.6275987284684357011881565152843114645681324961855;
    const double fourpidivsqrt3 = 7.255197456936871402376313030568622929136264992371;
    std::array<double,9> m{ k2Pi/a, 0., 0.,
                            twopidivsqrt3/a, fourpidivsqrt3/b, 0.,
                            0., 0., k2Pi/c };
    return RotMatrix(m);
  } else if ( alpha90 && gamma90 ) {
    //Usually this is monoclinic (with the b-axis unique setting)

    //Todo: alpha=gamma=90 for monoclinic corresponds to the b-axis unique
    //      setting, which is the most widespread setting. We might consider
    //      also supporting the a-axis unique setting (beta=gamma=90) and the
    //      c-axis unique setting (alpha=beta=90).

    const double sb = std::sin(beta);
    const double tanb = std::tan(beta);
    nc_assert(tanb!=0.0);
    const double cotb = 1.0 / tanb;//slightly more accurate than cos(b)/sin(b)
    std::array<double,9> m{ k2Pi/a, 0., 0.,
                             0., k2Pi/b, 0.,
                            -cotb*k2Pi/a, 0., k2Pi/(c*sb) };
    return RotMatrix(m);
  }
#endif
  //general case (usually triclinic and rhombohedral ends here)
  RotMatrix m = getLatticeRot(a,b,c,alpha,beta,gamma);
  m.inv();
  m *= k2Pi;
  return m;
}

NC::MaxHKL NC::estimateHKLRange( double dcutoff, double a, double b, double c )
{
  //For any lattice, h = a.G/2pi (a being the lattice vector and G the wave
  //vector of the plane), so |h| <= a*|G|/2pi = a/d, with equality when G is
  //parallel to a. Thus |h| <= a/dcutoff is an exact (and tight) bound, and
  //similarly for k and l. The bound is increased slightly as a safeguard
  //against numerical issues (e.g. planes with d exactly equal to dcutoff):
  nc_assert_always( dcutoff > 0.0 && a > 0.0 && b > 0.0 && c > 0.0 );
  auto f = [dcutoff]( double x ) -> int
  {
    const double v = std::floor( std::nextafter( 1.001 * x / dcutoff,
                                                 kInfinity ) );
    if ( !( v >= 1.0 ) )
      return 1;
    if ( v >= static_cast<double>( std::numeric_limits<int>::max() ) )
      return std::numeric_limits<int>::max();
    return static_cast<int>( v );
  };
  return { f( a ), f( b ), f( c ) };
}

bool NC::isRhombohedralSpaceGroup( int sg )
{
  return ( sg==146 || sg==148 || sg==155 || sg==160
           || sg==161 || sg==166 || sg==167 );
}

bool NC::usesRhombohedralAxes( int sg, double alpha )
{
  //(an unspecified angle, i.e. alpha=0, means the default hexagonal axes)
  return isRhombohedralSpaceGroup( sg ) && alpha > 0.0 && alpha != 90.0;
}

NC::KSqCalc::KSqCalc( double a, double b, double c,
                      double alpha, double beta, double gamma )
{
  nc_assert_always(alpha<kPi&&alpha>0);
  nc_assert_always(beta<kPi&&beta>0);
  nc_assert_always(gamma<kPi&&gamma>0);
  nc_assert_always(a>0);
  nc_assert_always(b>0);
  nc_assert_always(c>0);

  //Step 1: the six factors of |k|^2 = F1 h^2 + F2 k^2 + F3 l^2 + F4 hk + F5 hl
  //+ F6 kl from the usual formulas for 1/d^2 (times 4pi^2). Factors which
  //are equal by symmetry must be computed by the same expression, so that
  //they are the same double (step 2 then merges them):
  double F[6];
  const double tol = 1e-14;//(exact angle detection, as in getLatticeRot)
  const bool alpha90 = ncabs(alpha-kPiHalf)<tol;
  const bool beta90 = ncabs(beta-kPiHalf)<tol;
  const bool gamma90 = ncabs(gamma-kPiHalf)<tol;
  const bool gamma120 = ncabs(gamma-kPi*(2/3.))<tol;
  //Fixme: Consider using our custom versions of std::sin and std::cos below,
  //       for added reproducibility.
  if ( alpha90 && beta90 && gamma90 ) {
    //Orthogonal axes (cubic, tetragonal, orthorhombic): 1/d^2 = h^2/a^2 +
    //k^2/b^2 + l^2/c^2 (equal lengths give equal factors):
    F[0] = ncsquare( k2Pi / a );
    F[1] = ncsquare( k2Pi / b );
    F[2] = ncsquare( k2Pi / c );
    F[3] = F[4] = F[5] = 0.0;
  } else if ( alpha90 && beta90 && gamma120 && a == b ) {
    //Hexagonal axes: 1/d^2 = (4/3)(h^2+hk+k^2)/a^2 + l^2/c^2:
    F[0] = F[1] = F[3] = (4.0/3.0) * ncsquare( k2Pi / a );
    F[2] = ncsquare( k2Pi / c );
    F[4] = F[5] = 0.0;
  } else if ( a == b && b == c && alpha == beta && beta == gamma ) {
    //Rhombohedral axes: 1/d^2 = [ (h^2+k^2+l^2) sin^2(alpha)
    //+ 2(hk+hl+kl)(cos^2(alpha)-cos(alpha)) ] / [ a^2 (1-3cos^2(alpha)+2cos^3(alpha)) ]:
    const double ca = std::cos(alpha), sa = std::sin(alpha);
    const double denom = a*a*( 1.0 - 3.0*ca*ca + 2.0*ca*ca*ca );
    nc_assert( denom > 0.0 );
    F[0] = F[1] = F[2] = ncsquare(k2Pi) * sa*sa / denom;
    F[3] = F[4] = F[5] = ncsquare(k2Pi) * 2.0*( ca*ca - ca ) / denom;
  } else {
    //General case (monoclinic, triclinic): the reciprocal metric G from the
    //cell, with exact sines/cosines for 90 and 120 degree angles:
    const double ca = alpha90 ? 0.0 : std::cos(alpha);
    const double cb = beta90 ? 0.0 : std::cos(beta);
    const double cg = gamma90 ? 0.0 : ( gamma120 ? -0.5 : std::cos(gamma) );
    const double sa = alpha90 ? 1.0 : std::sin(alpha);
    const double sb = beta90 ? 1.0 : std::sin(beta);
    constexpr double sqrt3half
      = 0.86602540378443864676372317075293618347140262690519;
    const double sg = gamma90 ? 1.0 : ( gamma120 ? sqrt3half : std::sin(gamma) );
    const double vfact = 1.0 - ca*ca - cb*cb - cg*cg + 2.0*ca*cb*cg;//(V/abc)^2
    nc_assert( vfact > 0.0 );
    const double q = ncsquare(k2Pi) / vfact;
    F[0] = q * ncsquare( sa / a );//G11
    F[1] = q * ncsquare( sb / b );//G22
    F[2] = q * ncsquare( sg / c );//G33
    F[3] = 2.0 * q * ( ca*cb - cg ) / ( a*b );//2*G12
    F[4] = 2.0 * q * ( cg*ca - cb ) / ( a*c );//2*G13
    F[5] = 2.0 * q * ( cb*cg - ca ) / ( b*c );//2*G23
  }

  //Step 2: merge terms with identical factors (zero factors are dropped):
  for ( unsigned t = 0; t < 6; ++t ) {
    if ( F[t] == 0.0 ) {
      m_group_of_term[t] = 6;//(dummy group, never multiplied)
      continue;
    }
    unsigned g = 0;
    while ( g < m_ngroups && m_factor[g] != F[t] )
      ++g;
    if ( g == m_ngroups )
      m_factor[m_ngroups++] = F[t];
    m_group_of_term[t] = g;
  }
  nc_assert( m_ngroups >= 1 );
}

//|k|^2 = F1 h^2 + F2 k^2 + F3 l^2 + F4 hk + F5 hl + F6 kl, with the factors F
//from the reciprocal metric of the cell. For symmetry-equal planes to get
//identical results, the integer parts are summed exactly (64 bit) per distinct
//factor before multiplying, so e.g. cubic 333 and 511 both compute F1*27.
double NC::KSqCalc::ksq( int h, int k, int l ) const
{
  //Only the fma chain rounds. NB: Must stay identical to detail_fillhkl_ksq
  //in NCFillHKL_FMA.hh:
  const std::int64_t H = h, K = k, L = l;//(64 bit: exact for |h| < 1e9)
  const std::int64_t intpart[6] = { H*H, K*K, L*L, H*K, H*L, K*L };
  std::int64_t intsum[7] = { 0, 0, 0, 0, 0, 0, 0 };//[6] collects dropped terms
  for ( unsigned t = 0; t < 6; ++t )
    intsum[m_group_of_term[t]] += intpart[t];
  double acc = 0.0;
  for ( unsigned g = 0; g < m_ngroups; ++g )
    acc = std::fma( m_factor[g], static_cast<double>( intsum[g] ), acc );
  return acc;
}

double NC::KSqCalc::dspacing( int h, int k, int l ) const
{
  if (h==0&&k==0&&l==0)
    NCRYSTAL_THROW(BadInput,"Can not calculate d-spacing for hkl=000");
  const double k2 = ksq( h, k, l );
  if (!(k2>0.0))
    NCRYSTAL_THROW(CalcError,"Created invalid k-vector in d-spacing calculations (bad lattice parameters provided?)");
  return k2Pi / std::sqrt(k2);
}

  // Validate orientation parameters (check against null-vectors, parallel
  // vectors, consistency of angles in the two frames given dirtol,
  // etc.). Throws BadInput exception in case of issues.
namespace NCRYSTAL_NAMESPACE {
  namespace  {
    constexpr const char * vectName( const CrystalAxis& ) { return "CrystalAxis"; }
    constexpr const char * vectName( const HKLPoint& ) { return "HKLPoint"; }
    constexpr const char * vectName( const LabAxis& ) { return "LabAxis"; }
    static constexpr double orient_paralleltol = 1e-6;
    template<class TVect>
    void checkNotNull( const TVect& v, const char * designator = nullptr)
    {
      if (!v.template as<Vector>().mag2())
        NCRYSTAL_THROW2(BadInput,"Specified "
                        <<(designator?designator:"")
                        <<(designator?" ":"")
                        <<vectName(v)<<" is a null-vector.");
    }
    void baseCheckOD( const OrientDir& odir, const char * designator = nullptr )
    {
      checkNotNull(odir.lab,designator);
      const auto& cc = odir.crystal;
      if ( cc.has_value<HKLPoint>() )
        checkNotNull(cc.get<HKLPoint>(),designator);
      else if ( cc.has_value<CrystalAxis>() )
        checkNotNull(cc.get<CrystalAxis>(),designator);
      else
        NCRYSTAL_THROW2(BadInput,"Invalid crystal direction object (empty crystal direction)");
    }
  }

}

void NC::precheckLatticeOrientTol( double dirtol )
{
  if ( !( dirtol > 0.0 ) || dirtol > kPi )
    NCRYSTAL_THROW(BadInput, "Directional tolerance must be in interval (0.0,pi]");
}

void NC::precheckLatticeOrientDir( const OrientDir& odir )
{
  baseCheckOD(odir);
}

void NC::verifyLatticeOrientDef( const LabAxis& l1, const CrystalAxis& c1,
                                 const LabAxis& l2, const CrystalAxis& c2,
                                 double dirtol )
{
  checkNotNull(l1,"primary");
  checkNotNull(c1,"primary");
  checkNotNull(l2,"secondary");
  checkNotNull(c2,"secondary");
  precheckLatticeOrientTol(dirtol);

  if ( l1.as<Vector>().isParallel( l2.as<Vector>(), orient_paralleltol ) )
    NCRYSTAL_THROW(BadInput, "Specified primary and secondary lab directions are parallel");

  if ( c1.as<Vector>().isParallel( c2.as<Vector>(), orient_paralleltol ) )
    NCRYSTAL_THROW(BadInput, "Specified primary and secondary crystal directions are parallel");

  const double anglec = c1.as<Vector>().angle(c2.as<Vector>());
  const double anglel = l1.as<Vector>().angle(l2.as<Vector>());
  if ( ncabs( anglec - anglel ) > dirtol ) {
    NCRYSTAL_THROW2(BadInput,"Chosen orientation defining directions in the lab frame are "<<dbl2shortstr(anglel*kToDeg)
                    <<" deg apart, while the chosen directions in the crystal frame"
                    " are "<<dbl2shortstr(anglec*kToDeg)<<" deg apart. This is not within the specified"
                    " tolerance of "<<dbl2shortstr(dirtol)<<" rad. = "<<dbl2shortstr(dirtol*kToDeg)<<" deg.");
  }

}

void NC::precheckLatticeOrientDef( const OrientDir& dir1, const OrientDir& dir2, double dirtol )
{
  if ( dir1.crystal.has_value<CrystalAxis>() && dir2.crystal.has_value<CrystalAxis>() ) {
    //We can immediately do a full check:
    verifyLatticeOrientDef( dir1.lab, dir1.crystal.get<CrystalAxis>(),
                            dir2.lab, dir2.crystal.get<CrystalAxis>(),
                            dirtol );
    return;
  }
  //Partially check as much as we can:
  baseCheckOD(dir1,"primary");
  baseCheckOD(dir2,"secondary");
  precheckLatticeOrientTol(dirtol);
  if ( dir1.lab.as<Vector>().isParallel( dir2.lab.as<Vector>(), orient_paralleltol ) )
    NCRYSTAL_THROW(BadInput, "Specified primary and secondary lab directions are parallel");
  if ( dir1.crystal.has_value<HKLPoint>() && dir2.crystal.has_value<HKLPoint>() ) {
    if ( dir1.crystal.get<HKLPoint>().as<Vector>().isParallel( dir2.crystal.get<HKLPoint>().as<Vector>(), orient_paralleltol ) )
      NCRYSTAL_THROW(BadInput, "Specified primary and secondary crystal directions (hkl points) are parallel");
  }
}

NC::RotMatrix NC::verifyLatticeOrientDefAndConstructCrystalRotation( const OrientDir& dir1,
                                                                     const OrientDir& dir2,
                                                                     double dirtol,
                                                                     const RotMatrix& reci_lattice )
{
  precheckLatticeOrientDef( dir1, dir2, dirtol );
  nc_assert(!dir1.crystal.empty());
  nc_assert(!dir2.crystal.empty());

  Vector c1 = ( dir1.crystal.has_value<CrystalAxis>()
                ? dir1.crystal.get<CrystalAxis>().as<Vector>()
                : ( reci_lattice * dir1.crystal.get<HKLPoint>().as<Vector>().unit() ) ).unit();
  Vector c2 = ( dir2.crystal.has_value<CrystalAxis>()
                ? dir2.crystal.get<CrystalAxis>().as<Vector>()
                : ( reci_lattice * dir2.crystal.get<HKLPoint>().as<Vector>().unit() ) ).unit();
  Vector l1 = dir1.lab.as<Vector>().unit();
  Vector l2 = dir2.lab.as<Vector>().unit();

  verifyLatticeOrientDef( l1.as<LabAxis>(), c1.as<CrystalAxis>(),
                          l2.as<LabAxis>(), c2.as<CrystalAxis>(),
                          dirtol );

  //We are within the tolerance and all is OK (nothing null, parallel, etc.), so
  //we now ensure exact angle(c1,c2)==angle(l1,l2) by removing components of the
  //secondary direction parallel to the primary direction:

  c2 -= c1 * c2.dot(c1);
  l2 -= l1 * l2.dot(l1);
  c2.normalise();
  l2.normalise();
  return RotMatrix(l1,c1,l2,c2);
}

std::tuple<int,int,int> NC::normalAndDSpacingToHKLIndex( const RotMatrix& lattice_rot,
                                                         double dspacing,
                                                         const Vector& normal )
{
  nc_assert( dspacing > 0.0 );
  auto hkl = lattice_rot * normal;
  hkl /= dspacing;
  if ( hkl < -hkl )
    hkl = -hkl;
  Vector hkl_rounded( std::round(hkl[0]), std::round(hkl[1]), std::round(hkl[2]) );
  if ( ( hkl - hkl_rounded ).mag2() > 1e-10 )
    NCRYSTAL_THROW(CalcError,"HKL point estimated from dspacing+normal is not approximately integral.");
  return { int(hkl_rounded[0]), int(hkl_rounded[1]), int(hkl_rounded[2]) };
}
