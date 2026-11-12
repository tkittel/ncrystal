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

#include "NCrystal/internal/sab/NCSABCellSample.hh"
#include "NCrystal/internal/sab/NCSABUtils.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC = NCrystal;
namespace NCS = NCrystal::SAB;

NCS::LogLinDistSampler::LogLinDistSampler( double a, double fa, double logfa,
                                           double b, double fb, double logfb,
                                           bool force_linlin )
  : m_islinlin(false), m_a(a), m_b(b)
{
  nc_assert( a >= 0.0 );
  nc_assert( b >= a );
  nc_assert( fa >= 0 );
  nc_assert( fb >= 0 );
  nc_assert( std::isfinite(a) );
  nc_assert( std::isfinite(fa) );
  nc_assert( std::isfinite(b) );
  nc_assert( std::isfinite(fb) );
  nc_assert( std::isfinite(logfa) );
  nc_assert( std::isfinite(logfb) );

  const double bma = b-a;

  if ( !( ncmin(fa,fb,bma) > 0.0 ) )
    force_linlin = true;

  if ( !force_linlin ) {
    const double minus_dlnf = logfa-logfb;
    const double df = fb-fa;

    nc_assert(bma>0.0);
    if ( minus_dlnf == 0 || df == 0 ) {
      //fall back to uniform (via linlin code path):
      fa = fb = 1.0;
      force_linlin = true;
    } else {
      m_cache1 = m_b-m_a;
      m_cache2 = fb/fa;
      m_cache3 = logfb-logfa;
    }
  }

  if ( force_linlin ) {
    m_islinlin = true;
    if ( fa!=fb && bma!=0.0 ) {
      //uniform base + triangle top
      m_cache1 = ncmin(fa,fb);//av. height of uniform base
      m_cache2 = 0.5*(fa+fb);//av. height of base + triangle
      m_cache3 = ( fa>fb ? b : a );
      m_cache4 = ( fa>fb ? -bma : bma );
    } else {
      //uniform
      m_cache1 = 1.0;
      m_cache2 = 1.0;
      m_cache3 = a;
      m_cache4 = bma;
    }
  }
}

NCS::FullCellSampler::FullCellSampler( const CellData* cellptr )
  : m_cellptr(cellptr)
{
  //fixme: alternative constructor which simply takes prob1 or (W1,W2), if we
  //already know it?
  nc_assert(m_cellptr!=nullptr);
  const CellData& c = *m_cellptr;
  const double W1 = integrateLogLinOrLinWCL(c.a1,c.S[0],c.a2 , c.S[1],
                                            c.logS[0], c.logS[1]);
  const double W2 = integrateLogLinOrLinWCL(c.a1,c.S[2], c.a2 , c.S[3],
                                            c.logS[2], c.logS[3]);
  nc_assert_always(W1>0.0||W2>0.0);
  m_prob1 = W1/(W1+W2);
}

NC::PairDD NCS::FullCellSampler::sampleAlphaBeta( RNG& rng )
{
  nc_assert(m_cellptr!=nullptr);
  const CellData& c = *m_cellptr;
  //pick if we sample from the edge @ beta=c.b1 or the edge @ beta=c.b2:
  const int offset = ( rng.generate() <= m_prob1 ? 0 : 2 );
  auto& opt_sampler = m_samplers[offset/2];
  if ( !opt_sampler.has_value() ) {
    opt_sampler.emplace( c.a1, c.S[offset], c.logS[offset],
                         c.a2, c.S[offset+1], c.logS[offset+1] );
  }
  const double alpha = opt_sampler.value().sample(rng);

  //interpolation along beta is purely linear, so simply sample a triangle which
  //is zero at the far side:

  const double t = std::sqrt(rng.generate());
  const double beta = ( offset
                        ? nclerp( c.b1, c.b2, t )
                        : nclerp( c.b2, c.b1, t ) );
  return { alpha, ncclamp(beta, c.b1, c.b2) };
}

namespace NCRYSTAL_NAMESPACE {
  namespace SAB {
    namespace {
      Optional<PairDD> findActiveAlphaRange( const CellData& c, EDivKT e )
      {
        double a1 = c.a1;
        double a2 = c.a2;
        const double b1 = c.b1;
        const double b2 = c.b2;
        if (c.b2 <= -e.dbl())
          return NullOpt;//no overlap
        //fixme: reuse sqrt(e) in the next calls?:
        const double ap2 = getAlphaPlus(e,b2);//alpha^+(b2)
        if ( a1 >= ap2 )
          return NullOpt;//no overlap
        const double am2 = getAlphaMinus(e,b2);//alpha^-(b2)
        double am1(-1.0);
        if ( b1 >= -e.dbl() ) {
          am1 = getAlphaMinus(e,b1);//alpha^-(b1)
        }
        //snap along alpha:
        a2 = ncmin(a2,ap2);
        if ( b2 < 0 ) {
          a1 = ncmax(a1,am2);
        } else if ( b1 > 0 ) {
          a1 = ncmax(a1,am1);
        }
        if (!(a2 > a1))
          return NullOpt;//no overlap
        return PairDD{a1,a2};
      }


      class RCSImpl final : private NoCopyMove {
        CellData m_cell;
        EDivKT m_e;
        double m_alow = -1.0, m_aup=-1.0;//active alpha range.
        struct AlphaPtGeom {
          double blow, bup, bmid, bwidth;
        };
        VectD m_abinedges;
        VectD m_bincontrib_commul;
        VectD m_bincontrib;

        AlphaPtGeom alphaPtGem( double a ) const
        {
          nc_assert(m_aup>m_alow);
          nc_assert(valueInInterval(m_alow,m_aup,a));
          AlphaPtGeom res;
          //fixme: avoid repeated sqrt(e):
          const double bl = getBetaMinus(m_e,a);
          const double bu = getBetaPlus(m_e,a);
          res.blow = ncmax(m_cell.b1,bl);
          res.bup = ncmin(m_cell.b2,bu);
          if ( bl >= m_cell.b1 && bu <= m_cell.b2 ) {
            res.bmid = a;
            //Exactly 4*sqrt(e*a) (bu-bl would cancel for narrow bands):
            res.bwidth = 4.0 * std::sqrt( m_e.dbl() * a );
          } else {
            res.bmid = 0.5*(res.blow+res.bup);
            res.bwidth = res.bup - res.blow;
          }
          return res;
        }

        PairDD findSValsAtA( double a ) const
        {
          auto& S = m_cell.S;
          auto& logS = m_cell.logS;
          double s1 = interpLogLinOrLinWCL( m_cell.a1, S[0], m_cell.a2, S[1],
                                            a, logS[0], logS[1] ).first;
          double s2 = interpLogLinOrLinWCL( m_cell.a1, S[2], m_cell.a2, S[3],
                                            a, logS[2], logS[3] ).first;
          return {s1,s2};
        }

        double contribOfA( double a ) const
        {
          //Contribution is integral over beta at S(alpha=a,beta) multiplied
          //with (b2-b1), since it is a common factor everywhere and we are only
          //looking at sampling anyway.
          const AlphaPtGeom pt = alphaPtGem( a );
          const PairDD svals = findSValsAtA( a );
          StableSum ss, ss2;
          ss.add(m_cell.b2);
          ss.add(-pt.bmid);
          ss.mult(svals.first);
          ss2.add(pt.bmid);
          ss2.add(-m_cell.b1);
          ss2.mult(svals.second);
          ss.add(ss2);
          ss.mult(pt.bwidth);
          return ss.sum();
        }

      public:
        RCSImpl( const CellData& cell, EDivKT e )
          : m_cell(cell), m_e(e)
        {
          {
            auto arange = findActiveAlphaRange( cell, e );
            if (!arange.has_value())
              NCRYSTAL_THROW(BadInput,"RefCellSampler can not be used with"
                             " cell with no phasespace overlap");
            m_alow = arange.value().first;
            m_aup = arange.value().second;
          }
          nc_assert_always(m_alow < m_aup );
          nc_assert_always(m_alow >= cell.a1);
          nc_assert_always(m_aup <= cell.a2);

          //The contribution(alpha) might have local maxima at a few special
          //values, so for added safety, we must also evaluate on these if they
          //fall in a given alpha-bin. Otherwise the overlay value determined purely from the
          //alpha-bin edges could be underestimated.
          TinyVector<double,7> special_avals;
          {
            //Candidates. It is no harm to add spurious candidates here if that
            //candidate did not actually exist in the given setup.
            //
            //fixme: check if we can avoid repeated sqrt with the same args:
            auto addIfValid = [&]( double a ) {
              if ( valueInInterval(m_alow,m_aup,a) )
                special_avals.push_back(a);
            };
            addIfValid( m_e.dbl() );
            addIfValid( m_cell.b1/3 );
            addIfValid( m_cell.b2/3 );
            if ( m_cell.b1 + m_e.dbl() >= 0.0 ) {
              addIfValid( getAlphaMinus(m_e,m_cell.b1) );
              addIfValid( getAlphaPlus(m_e,m_cell.b1) );
            }
            if ( m_cell.b2 + m_e.dbl() >= 0.0 ) {
              addIfValid( getAlphaMinus(m_e,m_cell.b2) );
              addIfValid( getAlphaPlus(m_e,m_cell.b2) );
            }
          }

          //Now set up a number of alpha bins in which to sample:
          {
            constexpr std::size_t nbins = 1024;
            m_abinedges = linspace( m_alow, m_aup, nbins+1 );
            m_bincontrib_commul.reserve(nbins);
            StableSumKahan contrib_commul_ssum;
            double bin_alow = m_abinedges.front();
            double contrib_low_edge = contribOfA(bin_alow);
            for ( auto ibin : ncrange(nbins) ) {
              double bin_aup = vectAt(m_abinedges,ibin+1);
              double contrib_high_edge = contribOfA(bin_aup);
              double contrib = ncmax( contrib_low_edge, contrib_high_edge );
              for ( auto aspecial : special_avals ) {
                if ( valueInInterval( bin_alow, bin_aup, aspecial ) )
                  contrib = ncmax( contrib, contribOfA(aspecial));
              }
              contrib_commul_ssum.add( contrib );
              m_bincontrib_commul.push_back(contrib_commul_ssum.sum());
              m_bincontrib.push_back(contrib);
              contrib_low_edge = contrib_high_edge;
              bin_alow = bin_aup;
            }
          }
        }

        RefCellSampler::Result sample( RNG& rng )
        {
          RefCellSampler::Result res;
          res.ntries = 0;
          //Sample alpha value first:
          double a(-1.0);
          while (true) {
            ++res.ntries;
            nc_assert_always(res.ntries<1000);
            auto ibin
              = RandIdxByWeightPicker::oneShotPick( rng, m_bincontrib_commul );
            const double overlay = vectAt( m_bincontrib, ibin );
            a = randInterval( rng,
                              vectAt( m_abinedges, ibin ),
                              vectAt( m_abinedges, ibin+1 ) );
            const double contrib = contribOfA( a );
            constexpr double safety = 1.05;
            nc_assert_always( contrib <= overlay*1.045 );
            if ( overlay*safety*rng.generate() <= contrib )
              break;
          }
          nc_assert_always( valueInInterval(m_alow,m_aup,a) );
          //Sample beta value. This is done linearly from pt.blow to pt.bup
          //according to s-values:
          const AlphaPtGeom pt = alphaPtGem( a );
          const PairDD svals = findSValsAtA( a );
          const double bl = pt.blow;
          const double bu = pt.bup;
          const double b1 = m_cell.b1;
          const double b2 = m_cell.b2;
          auto contrib_of_b = [&svals,b1,b2] ( double b )
          {
            //Returns S(b)*(b2_b1) = (b2-b)*s1+(b-b1)*s2;
            StableSum ss;
            ss.add(b2*svals.first);
            ss.add(-b*svals.first);
            ss.add(b*svals.second);
            ss.add(-b1*svals.second);
            return ss.sum();
          };
          const double f_bl = contrib_of_b( bl );
          const double f_bu = contrib_of_b( bu );
          //Sample from [bl,bu] according to linear fct between (bl,f_bl) to
          //(bu,f_bu). To maximise numerical stability, this is done by
          //splitting into a base rectangle or the triangle top.
          const double c_base = ncmin(f_bl,f_bu);//height of rectangle base
          const double c_tot = 0.5*(f_bl+f_bu);//av. height of base+triangle
          double x;
          if ( c_tot*rng.generate()<=c_base ) {
            x = rng.generate();
          } else {
            //just a triangle, slanted the correct way:
            x = std::sqrt(rng.generate());
            if ( f_bl > f_bu )
              x = 1.0 - x;
          }
          res.alpha = a;
          res.beta = intervalPos( bl, bu, x );
          return res;
        }
      };
    }
  }
}

NCS::RefCellSampler::RefCellSampler( const CellData& cell, EDivKT e )
  : m_impl(static_cast<void*>(new RCSImpl(cell,e)))
{
}

NCS::RefCellSampler::~RefCellSampler()
{
  if ( m_impl ) {
    delete static_cast<RCSImpl*>(m_impl);
  }
}

NCS::RefCellSampler::RefCellSampler( RefCellSampler&& o )
{
  std::swap( m_impl, o.m_impl );
}

NCS::RefCellSampler&
NCS::RefCellSampler::operator=( RefCellSampler&& o )
{
  if ( m_impl ) {
    delete static_cast<RCSImpl*>(m_impl);
    m_impl = nullptr;
  }
  std::swap( m_impl, o.m_impl );
  return *this;
}

NCS::RefCellSampler::Result NCS::RefCellSampler::sampleAlphaBeta( RNG& rng )
{
  nc_assert( m_impl != nullptr );
  return static_cast<RCSImpl*>(m_impl)->sample(rng);
}
