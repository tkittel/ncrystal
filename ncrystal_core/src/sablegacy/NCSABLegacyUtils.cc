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

#include "NCrystal/internal/sablegacy/NCSABLegacyUtils.hh"
#include "NCrystal/internal/utils/NCIter.hh"

namespace NC = NCrystal;

void NC::SABLegacy::activeGridRanges( const SABData& data,
                                     EDivKT e,
                                     std::vector<std::pair<std::uint16_t,std::uint16_t>>& out_alpharanges,
                                     std::size_t& ibeta_low  )
{
  const auto& alphaGrid = data.alphaGrid();
  nc_assert(alphaGrid.size()>1);
  nc_assert(nc_is_grid(alphaGrid));
  nc_assert_always(alphaGrid.size()<std::numeric_limits<std::uint16_t>::max());

  ibeta_low = 0;
  out_alpharanges.clear();
  const double agrid_front = alphaGrid.front();
  const double agrid_back = alphaGrid.back();
  //At each beta grid point, we search for grid points marking the front and
  //back of the active grid range. For (large) gains of efficiency, we start off
  //the search for these alpha grid points at a given beta grid point, based on
  //their values at the previous beta grid point. This works well because the
  //kinematic boundaries are smooth curves.
  const auto itBegin = alphaGrid.begin();
  const auto itLast = std::prev(alphaGrid.end());
  auto itLow = alphaGrid.begin();
  auto itUpp = std::prev(alphaGrid.end());

  for (auto&& beta : enumerate(data.betaGrid())) {
    double alow(-1.0),aupp(-2.0);
    if ( beta.val > -e.dbl() ) {
      auto alims = getAlphaLimits( e, beta.val );
      alow = alims.first;
      aupp = alims.second;
    }
    if ( agrid_back <= alow || agrid_front >= aupp || aupp < alow ) {
      //No kinematically accessible alpha grid ranges at this beta point (or
      //energy is so ultra low that numerical imprecision led to aupp=alow) .
      if ( out_alpharanges.empty() ) {
        //Still didn't encounter any beta point with accessible alpha ranges, so
        //simply increment ibeta_low.
        ibeta_low = beta.idx + 1;
      } else {
        //Already encountered at least one beta point with accessible alpha
        //range, so must insert an empty alpha range at this point:
        nc_assert( alphaGrid.size()
                   < static_cast<std::size_t>(std::numeric_limits<std::uint16_t>::max()));
        auto na = static_cast<std::uint16_t>(alphaGrid.size());
        out_alpharanges.emplace_back(na,na);
      }
      continue;
    }

    //Ok, there is a accessible alpha grid range at this beta point, time to
    //find it.

    //Move itLow down or up as needed:
    while ( *itLow > alow && itLow > itBegin )
      --itLow;
    while ( itLow < itLast && *std::next(itLow) <= alow )
      ++itLow;
    //Move itUpp down or up as needed:
    if ( itUpp < itLow )
      itUpp = itLow;
    while ( *itUpp < aupp && itUpp < itLast )
      ++itUpp;
    while ( itUpp > itBegin && *std::prev(itUpp) >= aupp )
      --itUpp;

    //Register result. Notice that in case aupp==alow (for ultra ultra cold
    //neutrons, due to numerical issues) and that value is exactly on an alpha
    //grid point, we could have an empty range (itLow=itUpp). That is OK, we
    //simply register an empty range:
    nc_assert( itUpp < alphaGrid.end() );
    nc_assert( (aupp==alow) ? (itLow <= itUpp) : (itLow < itUpp) );
    out_alpharanges.emplace_back( static_cast<std::uint16_t>(std::distance(itBegin,itLow)),
                                  static_cast<std::uint16_t>(std::distance(itBegin,itUpp)) );

  }

}

NC::SABLegacy::TailedBreakdown NC::SABLegacy::createTailedBreakdown( const Span<const double>& alphaGrid,
                                                                   const Span<const double>& sab,
                                                                   const Span<const double>& logsab,
                                                                   const Span<const double>& alphaIntegrals_cumul,
                                                                   double alpha_low, double alpha_upp,
                                                                   const unsigned aidx_low, const unsigned aidx_upp )
{
  nc_assert( alpha_low <= alpha_upp );
  nc_assert( aidx_low <= aidx_upp );
  nc_assert( aidx_upp < alphaGrid.size() );

  //Constrain ranges to grid (xs outside is modelled as 0):

  alpha_low = ncclamp(alpha_low,alphaGrid.front(),alphaGrid.back());
  alpha_upp = ncclamp(alpha_upp,alphaGrid.front(),alphaGrid.back());

  TailedBreakdown tb;
  if ( aidx_low == aidx_upp || alpha_low == alpha_upp ) {
    //vanishing!
    return tb;
  }

  nc_assert( aidx_upp+1==alphaGrid.size() || alpha_upp <= alphaGrid[aidx_upp] );
  nc_assert( aidx_low==0 || alpha_low >= alphaGrid[aidx_low] );
  nc_assert( aidx_low+1==alphaGrid.size() || alpha_low < alphaGrid[aidx_low+1] );
  nc_assert( aidx_upp==0 || alpha_upp > alphaGrid[aidx_upp-1] );

  auto interpSVal = [&alphaGrid,&sab,&logsab](const std::size_t alphaidx_lowedge, const double alpha)
                    {
                      nc_assert( alphaidx_lowedge + 1 < (unsigned)alphaGrid.size() );
                      const double alpha0(alphaGrid[alphaidx_lowedge]);
                      const double alpha1(alphaGrid[alphaidx_lowedge+1]);
                      nc_assert( valueInInterval(alpha0,alpha1,alpha) );
                      auto i0(alphaidx_lowedge), i1(alphaidx_lowedge+1);
                      return SAB::interpLogLinOrLinWCL( alpha0, sab[i0], alpha1, sab[i1], alpha,
                                                   logsab[i0], logsab[i1] ).first;
                    };
  auto setTailPoint = [&interpSVal](TailedBreakdown::TailPoint& tp,unsigned aidx,double alpha)
                      {
                        tp.alpha = alpha;
                        tp.sval = interpSVal(aidx,alpha);
                        nc_assert(std::isfinite(tp.sval));
                        //We used to have here:
                        // tp.logsval = ( tp.sval > 0.0
                        //                ? std::log(tp.sval)
                        //                : -kInfinity );
                        //However in a rare case (1 test in release builds with
                        //intel one api compiler) that lead to an FPE in the log
                        //call. To move forward we instead limit svalues to the
                        //smallest positive non-subnormal number (2.2e-308
                        //before passing to the log):
                        constexpr double z = std::numeric_limits<double>::min();
                        tp.logsval = std::log( ncmax( z, tp.sval ) );
                      };

  //Enough setting up, time to analyse the tails and how they fall wrt the grid:
  if ( aidx_low + 1 == aidx_upp ) {
    //Special "narrow" case. Only a single alpha bin is touched by the
    //non-expanded region.  Integrate entire alpha-range in one go:
    tb.narrow = true;
    setTailPoint(tb.front,aidx_low,alpha_low);
    setTailPoint(tb.back,aidx_low,alpha_upp);
    tb.xs_front = SAB::integrateLogLinOrLinWCL( tb.front.alpha, tb.front.sval,
                                           tb.back.alpha, tb.back.sval,
                                           tb.front.logsval, tb.back.logsval );
    return tb;
  }

  //Non-narrow case, treat "front" and "back" tails separately.
  tb.imiddle_low = aidx_low;
  tb.imiddle_upp = aidx_upp;

  //Front (not there if alpha_low is outside the grid range):
  if ( alpha_low >= alphaGrid[aidx_low] ) {
    nc_assert( alpha_low <= alphaGrid[aidx_low + 1] );
    setTailPoint(tb.front,aidx_low,alpha_low);
    tb.xs_front = SAB::integrateLogLinOrLinWCL( tb.front.alpha, tb.front.sval,
                                           alphaGrid[aidx_low+1], sab[aidx_low+1],
                                           tb.front.logsval, logsab[aidx_low+1] );
    ++tb.imiddle_low;
  }
  //Back (not there if alpha_upp is outside the grid range):
  if ( alpha_upp <= alphaGrid[aidx_upp] ) {
    nc_assert( aidx_upp != 0 );
    nc_assert( alpha_upp >= alphaGrid[aidx_upp-1] );
    setTailPoint(tb.back,aidx_upp-1,alpha_upp);
    tb.xs_back = SAB::integrateLogLinOrLinWCL( alphaGrid[aidx_upp-1], sab[aidx_upp-1],
                                          tb.back.alpha, tb.back.sval,
                                          logsab[aidx_upp-1], tb.back.logsval );
    --tb.imiddle_upp;
  }
  tb.xs_middle = ( tb.imiddle_upp > tb.imiddle_low ?
                   alphaIntegrals_cumul[tb.imiddle_upp] - alphaIntegrals_cumul[tb.imiddle_low]
                   : 0.0 );
  return tb;
}
