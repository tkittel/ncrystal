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

#include "NCrystal/internal/sab/NCSABRefSampler.hh"
#include "NCrystal/internal/sab/NCSABCellInteg.hh"
#include "NCrystal/internal/sab/NCSABCellSample.hh"
#include "NCrystal/internal/sab/NCSABProcessor.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {

  namespace {

    //TODO: These randShuffle utilities could be moved to somewhere common, if
    //needed.
    template<class TSwapIJ>
    void randShuffleFisherYates(RNG& rng, std::size_t n, TSwapIJ&& swapIJ) {
      //A simple Fisher-Yates (Knuth) shuffle, but written to accomodate our
      //particular RNG interface which typically provides uint32_t values more
      //cheaply than uint64_t values, and with a flexible swap function which
      //for instance allows multiple containers to be shuffled at once in
      //randShuffleMultiple.
      if (n < 2)
        return;

      constexpr std::uint64_t lim_u32 =
        static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max());

      std::uint64_t i64 = static_cast<std::uint64_t>(n);

      // 64-bit part: i64 > lim_u32, so we need the more expensive
      // rng.generateInt64:
      for (; i64 > lim_u32; --i64) {
        std::uint64_t j64 = rng.generateInt64(i64); // [0, i64)
        const std::size_t a = static_cast<std::size_t>(j64);
        const std::size_t b = static_cast<std::size_t>(i64 - 1);
        swapIJ(a, b);
      }

      // 32-bit part: the index now fits in uint32_t so we can use the cheaper
      // rng.generateInt32:
      nc_assert(i64 <= lim_u32);
      std::uint32_t i32 = static_cast<std::uint32_t>(i64);
      for (; i32 > 1; --i32) {
        std::uint32_t j32 = rng.generateInt(i32); // [0, i32)
        swapIJ(static_cast<std::size_t>(j32),
               static_cast<std::size_t>(i32 - 1));
      }
    }

    namespace details {
      template<class V>
      void checkAllSizesEq(std::size_t n, const V& v) {
        if (v.size() != n)
          NCRYSTAL_THROW(BadInput,"randShuffleMultiple: size mismatch");
      }
      template<class V0, class... Vs>
      void checkAllSizesEq(std::size_t n, const V0& v0, const Vs&... vs) {
        if (v0.size() != n)
          NCRYSTAL_THROW(BadInput,"randShuffleMultiple: size mismatch");
        checkAllSizesEq(n, vs...);
      }
      template<class... TVectors>
      void randShuffleMultipleImpl(RNG& rng, std::size_t n, TVectors&... vecs) {
        checkAllSizesEq(n, vecs...);
        randShuffleFisherYates(rng, n, [&](std::size_t a, std::size_t b) {
          using std::swap;
          int dummy[] = {0, (swap(vecs[a], vecs[b]), 0)...};
          (void)dummy;
        });
      }
    }
    // Shuffle one or more containers in lockstep (must have same size):
    template<class TVFirst, class... TVRest>
    void randShuffle(RNG& rng, TVFirst& first, TVRest&... rest) {
      details::randShuffleMultipleImpl(rng, first.size(), first, rest...);
    }

    //If CellData fields a1,a2,b1,b2 have already been set, this can be used to
    //finish up, based on lower corner alpha/beta indices. The S values of the
    //cell are scaled by 2^k (k is returned), such that the largest is in
    //[0.5,1). This keeps them out of the subnormal range, except for values
    //more than ~1e308 below the largest, which are set to 0 (log/exp of
    //subnormal values are unreliable, e.g. with MSVC or FTZ/DAZ modes, exp can
    //return 0). Sampling within the cell does not depend on the scale:
    using CellData = SABUtils::CellData;
    int initCellDataSvals( CellData& c, const VectD& sab,
                           std::size_t na, std::size_t ia, std::size_t ib )
    {
      nc_assert( c.a2 > c.a1 );
      nc_assert( c.a1 >= 0.0 );
      nc_assert( c.b2 > c.b1 );
      std::size_t i = na*ib+ia;
      c.S[0] = vectAt(sab,i++);
      c.S[1] = vectAt(sab,i);
      i += (na-1);
      c.S[2] = vectAt(sab,i++);
      c.S[3] = vectAt(sab,i);
      int k = 0;
      const double smax = ncmax( c.S[0], c.S[1], c.S[2], c.S[3] );
      if ( smax > 0.0 ) {
        std::frexp( smax, &k );
        k = -k;
      }
      for ( int j = 0; j < 4; ++j ) {
        c.S[j] = std::ldexp( c.S[j], k );
        if ( c.S[j] < std::numeric_limits<double>::min() )
          c.S[j] = 0.0;
        c.logS[j] = NC::SABUtils::safeLogOrElse( c.S[j], 0.0 );
      }
      return k;
    }

    //This can initialise a CellData object from scratch based on lower corner
    //alpha/beta indices.
    CellData createCellData( std::size_t ia, std::size_t ib, const SABData& sd )
    {
      std::size_t na = sd.alphaGrid().size();
      CellData c;
      c.a1 = vectAt(sd.alphaGrid(),ia);
      c.a2 = vectAt(sd.alphaGrid(),ia+1);
      c.b1 = vectAt(sd.betaGrid(),ib);
      c.b2 = vectAt(sd.betaGrid(),ib+1);
      initCellDataSvals( c, sd.sab(), na, ia, ib );
      return c;
    }
    static constexpr std::int64_t packfactor = 1048576;
    static_assert( NC::ncconstexpr_ispow2( packfactor ), "" );
  }
}

namespace NCRYSTAL_NAMESPACE {

  namespace SABRef {
    namespace {

      struct EPtContribInfo {
        using CCells = std::vector<std::pair<double,std::int64_t>>;
        CCells contrib_cells;//{ contrib, cellidx }
        VectD cumul;
        //The contributions above are relative to the largest one, and must be
        //multiplied by 2^contribExponent to get the actual values:
        int contribExponent = 0;

        struct Decoded { bool fully_inside; CellData cellData; };
        Decoded decodeCell( const SABData& sab, std::size_t i ) const
        {
          std::int64_t cellidx = vectAt(contrib_cells,i).second;
          Decoded res;
          if ( cellidx < 0 ) {
            res.fully_inside = true;
            cellidx = -cellidx;
          } else {
            res.fully_inside = false;
          }
          nc_assert(cellidx>=1);
          --cellidx;//remove +1, left with packfactor*ib+ia
          nc_assert(cellidx>=0);
          std::size_t ia = static_cast<std::size_t>(cellidx % packfactor);
          std::size_t ib = static_cast<std::size_t>(cellidx / packfactor);
          res.cellData = createCellData( ia, ib, sab );
          return res;
        }
      };

      EPtContribInfo analyseEPt( const SABData& sab, EDivKT e )
      {
        EPtContribInfo res;

        const auto& alphas = sab.alphaGrid();
        const auto& betas = sab.betaGrid();
        const auto na = alphas.size();
        const auto nb = betas.size();
        nc_assert_always( na < packfactor && nb < packfactor );
        nc_assert_always( na>=2 && nb >= 2 && e.dbl() > 0.0
                          && std::isfinite(e.dbl()) );

        //First figure out which cells have which contributions, and if they are
        //completely covered or not (indicated by negative cellidx). The
        //contributions are kept as m*2^e, since the cells are integrated with
        //their S values scaled by 2^kscale (see initCellDataSvals):
        struct Contrib { double m; int e; std::int64_t cellidx; };
        std::vector<Contrib> contribs;
        contribs.reserve(8192);
        auto addContrib = [&contribs]( double c, int kscale,
                                       std::int64_t cellidx )
        {
          Contrib contrib;
          contrib.m = std::frexp( c, &contrib.e );
          contrib.e -= kscale;
          contrib.cellidx = cellidx;
          contribs.push_back( contrib );
        };

        const double foure = 4*e.dbl();
        const double minusE = -e.dbl();
        auto ptInsideKB = [foure]( double a, double b )
        {
          //Strict: a cell with a corner exactly on the boundary is left to
          //the (equally exact) partial-coverage integration:
          return ncsquare( a-b ) < foure * a;
        };

        for ( auto ia : ncrange(na-1) ) {
          const double a1 = vectAt(alphas,ia);
          const double a2 = vectAt(alphas,ia+1);
          for ( auto ib : ncrange(nb-1) ) {
            const double b2 = vectAt(betas,ib+1);
            if ( b2 <= minusE )
              continue;//no overlap (cheap check)
            const double b1 = vectAt(betas,ib);
            //Check if cell is fully above beta+(alpha):
            const double wb2 = std::sqrt( foure * a2 );
            if ( b1 >= a2 + wb2 )
              continue;
            //beta-(alpha) check written as squared comparisons: subtracting a
            //computed sqrt from a raw alpha cancels catastrophically near the
            //boundary. Since dbpm>=0, "b2<=a-dbpm" <=> "a-b2>=0 &&
            //(a-b2)^2>=foure*a":
            if ( a1 >= e.dbl() && a1-b2 >= 0.0 && ncsquare(a1-b2) >= foure*a1 )
              continue;
            if ( a2 <= e.dbl() && a2-b2 >= 0.0 && ncsquare(a2-b2) >= foure*a2 )
              continue;

            //Ok, we have an overlap. Let us check that S is not zero in entire
            //cell:
            CellData cell;
            cell.a1 = a1;
            cell.a2 = a2;
            cell.b1 = b1;
            cell.b2 = b2;
            const int kscale = initCellDataSvals( cell, sab.sab(), na, ia, ib );
            if ( cell.S[0]+cell.S[1]+cell.S[2]+cell.S[3]==0.0 )
              continue;//no actual contribution since S=0 in entire cell

            //We next need to determine if the cell is fully inside the
            //kinematic boundary or not. We simply check if the corners are
            //outside (and we only need to check 3 corners, since (a2,b2) will
            //be inside if the other 3 corners are):

            const bool inside11 = ptInsideKB( a1, b1 );
            const bool inside12 = ptInsideKB( a1, b2 );
            const bool inside21 = ptInsideKB( a2, b1 );
            //const bool inside22 = ptInsideKB( a2, b2 );

            //+1 so the sign of cellidx can encode info:
            const std::int64_t cellidx = ( packfactor*static_cast<std::int64_t>(ib)
                                           + static_cast<std::int64_t>(ia) + 1 );
            nc_assert( cellidx >= 1 );
            if ( inside11 & inside12 & inside21/* & inside22*/ ) {
              //cell is fully inside boundary
              nc_assert(ptInsideKB( a2, b2 ));//fourth corner is also inside
              StableSumKahan sum;
              SABUtils::StdLogLinCellIntegrator::integrateFullCell( cell, sum );
              addContrib( sum.sum(), kscale, -cellidx );//-cellidx
            } else {
              //cell crosses boundary
              StableSumKahan sum;
              SABUtils::StdLogLinCellIntegrator
                ::integrateWithinKB( cell, e,
                                     SABCfg::IntegrationScheme::MaxPrec, sum );
              addContrib( sum.sum(), kscale, cellidx );//+cellidx
            }
          }
        }
        //Contributions relative to the largest one (exactly, as the mantissas
        //are only scaled by powers of two):
        int emax = std::numeric_limits<int>::min();
        for ( auto& contrib : contribs )
          if ( contrib.m > 0.0 )
            emax = std::max( emax, contrib.e );
        if ( emax == std::numeric_limits<int>::min() )
          emax = 0;
        res.contribExponent = emax;
        res.contrib_cells.reserve( contribs.size() );
        for ( auto& contrib : contribs )
          res.contrib_cells.emplace_back( std::ldexp( contrib.m,
                                                      contrib.e - emax ),
                                          contrib.cellidx );
        //Now sort by contribution, small to large. That will make our
        //cumulative sum more robust (otherwise tiny contributions after large
        //ones might get imprecise).
        std::sort( res.contrib_cells.begin(), res.contrib_cells.end() );

        //Create cumul vect:
        res.cumul.reserve(res.contrib_cells.size());
        StableSum contribsum;
        for ( auto& contrib : res.contrib_cells ) {
          nc_assert( contrib.first >= 0.0 );
          contribsum.add( contrib.first );
          res.cumul.push_back( contribsum.sum() );
        }
        return res;
      }


    }
  }
}


namespace NCRYSTAL_NAMESPACE {

  namespace SABRef {
    namespace {
      Optional<RefSampleExtension>
      prepareExtension( shared_obj<const SABData> sab,
                        EDivKT e,
                        Optional<RefSampleExtension> opt_ext )
      {
        if ( !opt_ext.has_value() )
          return NullOpt;//has to be explicitly requested
        RefSampleExtension res = opt_ext.value();
        if ( res.emax.dbl() == 0.0 )
          res.emax = autodetectRefSampleEmax( sab );
        NeutronEnergy ekin( e.dbl() * sab->temperature().kT() );
        if ( ekin <= res.emax )
          return NullOpt;//extension not relevant
        if ( res.extender == nullptr )
          res.extender = makeSO<SAB::SABFGExtender>( sab->temperature(),
                                                     sab->elementMassAMU(),
                                                     SigmaBound{1.0} );
        return res;
      }
    }
  }
}

NC::NeutronEnergy
NC::SABRef::autodetectRefSampleEmax( const shared_obj<const SABData>& sab )
{
  if ( sab->suggestedEmax() > 0.0 )
    return NeutronEnergy( sab->suggestedEmax() );
  //The expensive way:
  using SP = SABUtils::SABProcessor;
  auto cfg =  SABCfg::createConfig( 5 );//6 is too slow, and the focus
                                        //here is sampling, not xs.
  SP sp( cfg, sab, nullptr, SP::SampleSupport::NO );
  return sp.getEMaxInfo().ekin;
}

std::pair<NC::VectD,NC::VectD>
NC::SABRef::refSampleAlphaBeta( RNG& rng,
                                shared_obj<const SABData> sab,
                                EDivKT e,
                                std::uint64_t nsample,
                                Optional<RefSampleExtension> opt_extension )
{
  //Check if we need to use an extension or impose an EMax on the cell-based
  //processing:
  opt_extension = prepareExtension(sab,e,opt_extension);

  const double kT = sab->temperature().kT();
  EDivKT E_div_kT_cell = e;
  if ( opt_extension.has_value() ) {
    E_div_kT_cell = EDivKT{ opt_extension.value().emax.dbl() / kT };
    //(or prepareExtension would not have returned a value):
    nc_assert( E_div_kT_cell < e );
    nc_assert( E_div_kT_cell.dbl() > 0.0 );
  }

  //Analyse requested energy point (but not beyond emax:
  auto analysedEpt = analyseEPt( sab, E_div_kT_cell );
  const std::size_t ncells = analysedEpt.contrib_cells.size();

  //Find the probability for events to be covered by the extender (if any):

  double probabilityExtension(0.0);
  if ( opt_extension.has_value() ) {
    auto& ext = opt_extension.value();
    nc_assert( ext.extender != nullptr );
    const double ext_xs_Ecell
      = ext.extender->crossSection(ext.emax).dbl();
    const double ext_xs_E
      = ext.extender->crossSection(NeutronEnergy(e.dbl()*kT)).dbl();
    double extS = 4.0 * ncfma( ext_xs_E, e.dbl(),
                                  -(ext_xs_Ecell*E_div_kT_cell.dbl()) );
    const double cellS = std::ldexp( analysedEpt.cumul.back(),
                                     analysedEpt.contribExponent );
    probabilityExtension = extS / ( extS + cellS );
  }
  nc_assert( probabilityExtension >= 0.0 && probabilityExtension <= 1.0 );

  //Now count how many are sampled in each cell, so we can afterwards do the
  //actual sampling one cell at a time:
  std::vector<std::size_t> cellcount( ncells, 0 );
  for ( std::uint64_t i = 0; i < nsample; ++i ) {
    if ( probabilityExtension <= 0.0 || rng() >= probabilityExtension )
      ++vectAt( cellcount,
                RandIdxByWeightPicker::oneShotPick( rng, analysedEpt.cumul ) );
  }

  //Now do the actual sampling - first in cells:
  std::pair<NC::VectD,NC::VectD> res;
  auto& a = res.first;
  auto& b = res.second;
  a.reserve( countToSizeT( nsample ) );
  b.reserve( countToSizeT( nsample ) );
  for ( auto i : ncrange( ncells ) ) {
    auto n = vectAt(cellcount,i);
    if (!n)
      continue;//nothing in this cell

    auto decodedCell = analysedEpt.decodeCell( sab, i );
    bool fully_inside = decodedCell.fully_inside;
    auto& cell = decodedCell.cellData;

    //Now sample:
    if ( fully_inside ) {
      //full cell sampling => just use simple alg even for this reference
      //(fixme: we could try to only use RecFellSampler as a check)
      SABUtils::FullCellSampler fc(&cell);
      while ( n ) {
        auto ab = fc.sampleAlphaBeta( rng );
        a.push_back(ab.first);
        b.push_back(ab.second);
        --n;
      }
    } else {
      //crossing => use RefCellSampler
      SABUtils::RefCellSampler rcs( cell, E_div_kT_cell );
      while ( n ) {
        auto ab = rcs.sampleAlphaBeta( rng );
        a.push_back(ab.alpha);
        b.push_back(ab.beta);
        --n;
      }
    }
  }

  //Now do the sampling in the extended region (if any):
  if ( a.size() < nsample ) {
    const NeutronEnergy ekin{ e.dbl() * kT };
    nc_assert(opt_extension.has_value());
    const auto& ext = opt_extension.value();
    nc_assert( ext.extender!=nullptr );
    nc_assert( ekin > ext.emax );
    const double fourecell = 4.0 * E_div_kT_cell.dbl();
    while ( a.size() < nsample ) {
      auto ab = ext.extender->sampleAlphaBeta(rng, ekin );
      if ( ncsquare( ab.first - ab.second ) >= fourecell * ab.first ) {
        a.push_back(ab.first);
        b.push_back(ab.second);
      }
    }
  }

  //Finally, shuffle results to not get them ordered by cell (or extension):
  randShuffle( rng, a, b);
  return res;
}
