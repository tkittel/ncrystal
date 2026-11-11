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

#include "NCrystal/internal/sab/NCSABUtils.hh"
#include "NCrystal/internal/utils/NCIter.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCMath.hh"
namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {
    namespace {
      std::size_t detail_trimZeroEdgesFromKernel(ScatKnlData& input)
      {
        //Some input kernels do not fill the entire grid on which they are
        //specified. Detect and remove such parts of the grid (no validation).
        const auto nalpha = input.alphaGrid.size();
        const auto nbeta = input.betaGrid.size();

        auto betaRowIsZero = [&input,nalpha] ( std::size_t ibeta )
        {
          auto idx0 = ibeta * nalpha;
          auto idxE = idx0 + nalpha;
          for ( auto idx = idx0; idx != idxE; ++idx ) {
            nc_assert( std::isfinite(vectAt( input.sab, idx ))
                       && vectAt( input.sab, idx )>=0.0 );
            if ( vectAt( input.sab, idx ) )
              return false;
          }
          return true;
        };
        auto alphaColIsZero = [&input,nalpha] ( std::size_t ialpha )
        {
          auto idx0 = ialpha;
          auto idxE = input.sab.size();
          for ( auto idx = idx0; idx < idxE; idx += nalpha ) {
            nc_assert( std::isfinite(vectAt( input.sab, idx ))
                       && vectAt( input.sab, idx )>=0.0 );
            if ( vectAt( input.sab, idx ) )
              return false;
          }
          return true;
        };

        std::size_t nTrimBetaUpper(0);
        std::size_t nTrimBetaLower(0);
        std::size_t nTrimAlphaUpper(0);

        for ( auto i : ncrange( nalpha ) ) {
          const auto ialpha = nalpha - i - 1;
          if ( vectAt(input.alphaGrid,ialpha) > 0.0 && alphaColIsZero( ialpha ) ) {
            ++nTrimAlphaUpper;
          } else {
            break;
          }
        }

        for ( auto i : ncrange( nbeta ) ) {
          const auto ibeta = nbeta - i - 1;
          if ( vectAt(input.betaGrid,ibeta) > 0.0 && betaRowIsZero( ibeta ) ) {
            ++nTrimBetaUpper;
          } else {
            break;
          }
        }

        for ( auto ibeta : ncrange( nbeta ) ) {
          if ( vectAt(input.betaGrid,ibeta) < 0.0 && betaRowIsZero( ibeta ) ) {
            ++nTrimBetaLower;
          } else {
            break;
          }
        }

        //Abort if too agressive (e.g. S=0 everywhere)
        if (  nTrimAlphaUpper >= nalpha ) {
          nTrimAlphaUpper = 0;
          nTrimBetaUpper = 0;
          nTrimBetaLower = 0;
        }

        const auto ntrimtot = nTrimBetaUpper + nTrimBetaLower + nTrimAlphaUpper;

        if ( ntrimtot == 0 )
          return 0;//nothing to trim

        //Trim:
        VectD new_sab;
        new_sab.reserve( ( nalpha - nTrimAlphaUpper ) * ( nbeta - (nTrimBetaUpper + nTrimBetaLower) ) );
        const auto nbetalim = nbeta - nTrimBetaUpper;
        for ( auto ibeta : ncrange( nbeta ) ) {
          if ( ! ( ibeta >= nTrimBetaLower && ibeta < nbetalim ) )
            continue;
          const auto offset = nalpha * ibeta;
          for ( auto ialpha : ncrange( nalpha-nTrimAlphaUpper ) )
            new_sab.push_back( vectAt( input.sab, ialpha + offset ) );
        }
        nc_assert_always( new_sab.size() == ( nalpha - nTrimAlphaUpper ) * ( nbeta - (nTrimBetaUpper + nTrimBetaLower) ) );
        std::swap( input.sab, new_sab );
        if ( nTrimAlphaUpper ) {
          VectD new_alpha( input.alphaGrid.begin(), std::next(input.alphaGrid.begin(),nalpha-nTrimAlphaUpper) );
          std::swap(input.alphaGrid,new_alpha);
        }
        if ( nTrimBetaUpper || nTrimBetaLower ) {
          VectD new_beta( std::next(input.betaGrid.begin(),nTrimBetaLower), std::next(input.betaGrid.begin(),nbeta-nTrimBetaUpper) );
          std::swap(input.betaGrid,new_beta);
        }
        nc_assert_always ( input.sab.size() == input.alphaGrid.size() * input.betaGrid.size() );
        return ntrimtot;
      }
    }
  }
}


double NC::SABUtils::evalSABOneShot( const SABData& sd,
                                     double alpha, double beta )
{
#ifndef NDEBUG
  nc_assert_always( !ncisnan( alpha ) && !ncisnan( beta ) );
#endif
  const auto& A = sd.alphaGrid();
  const auto& B = sd.betaGrid();
  if ( !( alpha >= A.front() && alpha <= A.back()
          && beta >= B.front() && beta <= B.back() ) )
    return 0.0;
  //Index of the cell containing x (the last cell at the upper edge):
  auto cellIdx = []( const VectD& g, double x )
  {
    const auto i = static_cast<std::size_t>
      ( std::upper_bound( g.begin(), g.end(), x ) - g.begin() );
    return ncmin( i, g.size() - 1 ) - 1;
  };
  const std::size_t ia = cellIdx( A, alpha );
  const std::size_t ib = cellIdx( B, beta );
  const std::size_t na = A.size();
  const auto& S = sd.sab();
  double s[4] = { vectAt( S, calcSABIdx( na, ia, ib ) ),
                  vectAt( S, calcSABIdx( na, ia + 1, ib ) ),
                  vectAt( S, calcSABIdx( na, ia, ib + 1 ) ),
                  vectAt( S, calcSABIdx( na, ia + 1, ib + 1 ) ) };
  const double smax = ncmax( s[0], s[1], s[2], s[3] );
  if ( !( smax > 0.0 ) )
    return 0.0;
  //Scale by 2^k so the largest corner is in [0.5,1), keeping the others out
  //of the subnormal range where possible (those still subnormal are set to
  //0, as when integrating kernels):
  int k;
  std::frexp( smax, &k );
  for ( auto& v : s ) {
    v = std::ldexp( v, -k );
    if ( v < std::numeric_limits<double>::min() )
      v = 0.0;
  }
  const double a0 = vectAt( A, ia );
  const double a1 = vectAt( A, ia + 1 );
  const double s0 = interpolate_loglin_fallbacklinlin( a0, s[0], a1, s[1],
                                                       alpha );
  const double s1 = interpolate_loglin_fallbacklinlin( a0, s[2], a1, s[3],
                                                       alpha );
  const double b0 = vectAt( B, ib );
  const double b1 = vectAt( B, ib + 1 );
  //Linear in beta, from the nearest edge (exact at both):
  const double res = ( beta - b0 <= b1 - beta
                       ? nclerp( s0, s1, ( beta - b0 ) / ( b1 - b0 ) )
                       : nclerp( s1, s0, ( b1 - beta ) / ( b1 - b0 ) ) );
  return std::ldexp( res, k );
}

double NC::SABUtils::evalSABCell( const SABData& sd,
                                  std::size_t ia, std::size_t ib,
                                  double alpha, double beta )
{
  const auto& ag = sd.alphaGrid();
  const auto& bg = sd.betaGrid();
  const auto& S = sd.sab();
  const std::size_t na = ag.size();
#ifndef NDEBUG
  nc_assert_always( ia + 1 < na && ib + 1 < bg.size() );
  nc_assert_always( ag[ia] <= alpha && alpha <= ag[ia+1] );
  nc_assert_always( bg[ib] <= beta && beta <= bg[ib+1] );
#endif
  const std::size_t i0 = calcSABIdx( na, ia, ib );
  const std::size_t i1 = calcSABIdx( na, ia, ib + 1 );
  const double s[4] = { vectAt( S, i0 ), vectAt( S, i0 + 1 ),
                        vectAt( S, i1 ), vectAt( S, i1 + 1 ) };
  double ls[4];
  for ( auto i : ncrange( 4 ) )
    ls[i] = safeLogOrElse( s[i], 0.0 );
  return evalSABCellCorners( vectAt( ag, ia ), vectAt( ag, ia + 1 ),
                             vectAt( bg, ib ), vectAt( bg, ib + 1 ),
                             s, ls, alpha, beta );
}

NC::shared_obj<const NC::VectD>
NC::SABUtils::createLogSTable( const SABData& sd )
{
  auto res = makeSO<VectD>();
  res->reserve( sd.sab().size() );
  for ( auto v : sd.sab() )
    res->push_back( safeLogOrElse( v, 0.0 ) );
  return res;
}

NC::SABData NC::SABUtils::transformKernelToStdFormat( ScatKnlData&& input_orig )
{
  validateScatKnlData(input_orig);

  //////////////////////
  // Trim null edges: //
  //////////////////////

  ScatKnlData input{ std::move(input_orig) };
  auto ntrimmed = detail_trimZeroEdgesFromKernel(input);
  if ( ntrimmed )
    NCRYSTAL_WARN("Discarding "<<ntrimmed<<" edges of provided kernel"
                  " data due to missing S values.");

  ////////////////////////////////////////////////////
  // First convert ScatKnlData object to type "SAB" //
  ////////////////////////////////////////////////////

  if ( input.knltype == ScatKnlData::KnlType::SCALED_SYM_SAB ) {
    //Remove symmetry and actually specify full table :
    VectD complete_betagrid, complete_sab;
    expandBetaAndSABToAllBetas( input.betaGrid, input.alphaGrid, input.sab, complete_betagrid, complete_sab );
    std::swap( input.betaGrid, complete_betagrid );
    std::swap( input.sab, complete_sab );
    input.knltype = ScatKnlData::KnlType::SCALED_SAB;
  }

  if ( input.knltype == ScatKnlData::KnlType::SCALED_SAB ) {

    //Unscale, i.e. calculate S=S_scaled*exp(-beta/2):

    for ( auto&& beta : enumerate(input.betaGrid) ) {
      auto slice = sliceSABAtBetaIdxMutable( input.sab, input.alphaGrid.size(), beta.idx );
      const double exparg = -0.5*beta.val;
      if ( exparg < 700.0 ) {
        //Straightforward, when beta is not below -1400, we can safely and
        //accurately calculate exp(-beta/2):
        const double expfact = std::exp( exparg );
        for ( double& sabEntry : slice )
          sabEntry *= expfact;
      } else {
        //The factor of exp(-beta/2) overflows at double precision. Fortunately,
        //this normally happens when S_scaled is itself a very small number,
        //thus balancing out the extreme factor (otherwise the physics would
        //probably not make much sense anyway). The numerically safe way to
        //calculate S is in this case to carry out the cancellation in
        //log-space, i.e. by using the more expensive rewrite:
        //
        // S = S_scaled*exp(-beta/2) = exp(-beta/2+log(S_scaled))
        //
        //Of course, when S_scaled=0, S is trivially also 0.
        //
        //Since log(eps), where eps is the smallest positive floating point
        //number in double precision, is roughly -744, it follows that any
        //scaled sab grid specified must have strictly S_scaled=0 for all points
        //at beta less than approximately -2900. Accurate specification of
        //S(alpha,beta) at large negative beta-values thus requires unscaled
        //S(alpha,beta) values to be specified directly (assuming of course that
        //whoever provided the table were able to circumvent numerical issues in
        //their calculations).
        for ( double& sabEntry : slice ) {
          if (sabEntry==0.0)
            continue;//0 is still 0 when multiplied with exp(-beta/2)
          double expargcombined = exparg + std::log(sabEntry);
          if ( expargcombined < 700.0 ) {
            sabEntry = std::exp(expargcombined);
          } else {
            NCRYSTAL_THROW2(BadInput,"Problems unscaling of S(alpha,beta), at point where S_scaled="
                            <<sabEntry<<" and beta="<<beta.val<<", since it requires evaluation"
                            " of exp("<<expargcombined<<") which is infinity at double precision."
                            " Most likely this indicates a problem with the input data.");
          }
        }
      }
    }//beta idx
    input.knltype = ScatKnlData::KnlType::SAB;
  }

  if ( input.knltype == ScatKnlData::KnlType::SQW ) {
    //Convert  S(q,w) -> S(alpha,beta).
    //Postponed for later release.
    NCRYSTAL_THROW2(LogicError,"Support for kernels in S(q,w) format is planned, but not yet implemented.");
#if 0
    //Leftover notes from previous attempt:

    //Still not working/verified (and we also need an SQW_SCALED version...):

    //First thing to note is that w=omega in the input is actually W=hbar*omega
    //(the unit is eV). Because of that, the conversion is (the beta and S
    //conversion would have had a factor of hbar if omega was actually omega
    //rather than hbar*omega):
    //
    // alpha = hbar^2/2*M*kT * qval^2
    // beta = -1/kT * W
    // The S-values themselves needs to be scaled with 1/kT.

    //other miscellaneous leftover notes:
    //        #beta = -B*w (if "w" is actually hbar*w, then B=1/kY)
    //        #alpha = A*q^2
    //        #dw = dbeta/B
    //        #int(S(q,w),dw) = (1/B) * int(S(alpha,beta)dbeta)
    //        #SO: If I get 1 for int(S(alpha,beta)dbeta), it must mean that S_ab = S_qw * B

    const double kT = constant_boltzmann*input.temperature;//unit is eV
    const double invkT = 1.0/kT;
    const double M = input.elementMassAMU*constant_dalton2eVc2;//Unit is eV/(Aa/s)^2 = eV*s^2/Aa^2
    const double hbar = constant_hbar;//eV*s
    const double c_qsq2a = 1e40*hbar*hbar/ (2.0 * M * kT);//Unit is Aa^2 [NOTICE: !!!!!!! 1e40 here is just an unvalidated try!]
    for (auto& e : input.alphaGrid)
      e = c_qsq2a * (e * e);
    for (auto& e : input.betaGrid)
      e *= -invkT;
    //reverse order, since conversion factor was negative.:
    std::reverse(std::begin(input.betaGrid), std::end(input.betaGrid));
    for (auto& e : input.sab)
      e *= invkT;
#endif
    input.knltype = ScatKnlData::KnlType::SAB;
  }

  nc_assert_always( input.knltype == ScatKnlData::KnlType::SAB );

  ///////////////////////////////////////////////////////////////////////////////
  // Thicken beta grid to reduce numerical artifacts in subsequent processing. //
  // This is obviously just a workaround and not the correct way to solve this //
  // (TODO: Remove once processing handles sparse grids better)!               //
  ///////////////////////////////////////////////////////////////////////////////

  //FIXME: This should probably be removed in the near future!

  int nminbeta_int = ncgetenv_int("SAB_BETATHICKENING_MINNBETA",500);
  nc_assert_always( nminbeta_int >=0 && nminbeta_int < 20000 );
  const unsigned nminbeta = static_cast<unsigned>(nminbeta_int);

  if ( !input.betaGridOptimised && input.betaGrid.size() < nminbeta ) {
    const unsigned nextra = nminbeta / input.betaGrid.size();//Number of extra beta points
                                                             //to insert between each point
    nc_assert_always( nextra >= 1 );
    auto newNBeta = 1 + ( input.betaGrid.size() -1 ) * ( 1 + nextra );
    VectD newb;
    newb.reserve( newNBeta );
    VectD newS;
    const std::size_t nalpha = input.alphaGrid.size();
    newS.reserve( newNBeta * nalpha );
    auto itB = input.betaGrid.begin();
    auto itBLast = std::prev(input.betaGrid.end());
    std::size_t ibeta(0);
    for ( ; itB != itBLast; ++itB, ++ibeta ) {
      auto itBNext = std::next(itB);
      auto alphaSlice = sliceSABAtBetaIdx( input.sab, nalpha, ibeta );
      auto alphaSliceNext = sliceSABAtBetaIdx( input.sab, nalpha, ibeta+1 );
      //First copy over existing row:
      newb.push_back(*itB);
      std::copy(alphaSlice.begin(), alphaSlice.end(), std::back_inserter(newS));
      //Then insert the new rows, one by one:
      const double dBeta = (*itBNext - *itB)/(nextra+1);
      for ( auto iextra : ncrange(nextra) ) {
        double beta = *itB + (iextra+1)*dBeta;
        newb.push_back(beta);
        for ( auto ialpha : ncrange(nalpha) ) {
          double S = alphaSlice.at(ialpha);
          double SNext = alphaSliceNext.at(ialpha);
          double kkk = ( beta - *itB) / (*itBNext - *itB);
          //newS.push_back( S + (SNext - S) * kkk  );
          newS.push_back( S * ( 1.0 - kkk) + kkk * SNext );
        }
      }
    }
    //Add last slide as well:
    auto alphaSliceLast = sliceSABAtBetaIdx( input.sab, nalpha, ibeta );
    std::copy(alphaSliceLast.begin(), alphaSliceLast.end(), std::back_inserter(newS));
    newb.push_back(*itBLast);
    nc_assert_always( newS.size() == newNBeta * nalpha );
    nc_assert_always( newb.size() == newNBeta );

    //Apply:
    std::swap(newS,input.sab);
    std::swap(newb,input.betaGrid);
  }

  ///////////////////////////////////////////
  // Transfer to SABData object and return //
  ///////////////////////////////////////////

  SABData out{ std::move(input.alphaGrid),
               std::move(input.betaGrid),
               std::move(input.sab),
               input.temperature,
               input.boundXS,
               input.elementMassAMU,
               input.suggestedEmax };

#ifndef NDEBUG
  validateScatKnlData(out);
#endif

  return out;
}

void NC::SABUtils::expandBetaAndSABToAllBetas( Span<const double> halfbetagrid,
                                               Span<const double> alphagrid,
                                               Span<const double> sab_for_halfbetagrid,
                                               VectD& complete_betagrid,
                                               VectD& complete_sab )
{
  //Prepare and check:
  complete_betagrid.clear(  );
  complete_sab.clear();
  const std::size_t nalpha = alphagrid.size();
  const std::size_t nbeta_old = halfbetagrid.size();
  const std::size_t nbeta_positive = nbeta_old - 1;
  const std::size_t nbeta_new = nbeta_positive*2 + 1;
  nc_assert_always(!halfbetagrid.empty());
  nc_assert_always(halfbetagrid.front()==0.0);
  nc_assert_always( nbeta_old * nalpha == static_cast<std::size_t>(sab_for_halfbetagrid.size()) );

  //Step 1. Create complete beta-grid:
  // -> the negative values and zero:
  complete_betagrid.reserve( nbeta_new );
  for (auto it = halfbetagrid.rbegin(); it != halfbetagrid.rend(); ++it)
    complete_betagrid.emplace_back( - (*it) );

  // -> avoid signed negative zero, for aesthetic reasons :
  nc_assert(complete_betagrid.back()==0.0);
  complete_betagrid.back() = 0.0;
  // -> add the positive values (the +1 below avoids adding the zero again):
  for (auto e: Span<const double>(halfbetagrid.begin()+1,halfbetagrid.end()) )
    complete_betagrid.emplace_back(e);
  nc_assert_always( complete_betagrid.size() == nbeta_new );

  //Step 2. Expand the sab kernel to complete beta range

  //Step 2.1: fill out the parts for beta<0 using S(alpha,-beta) := S(alpha,beta):
  complete_sab.reserve( nbeta_new * alphagrid.size() );
  complete_sab.resize( alphagrid.size() * nbeta_positive, 0. );
  const auto srcbegin = sab_for_halfbetagrid.begin();
  auto targetrowbegin = complete_sab.begin();

  for ( std::size_t i = 0; i < nbeta_positive; ++i )  {
    auto srcrowbegin = srcbegin + (nbeta_positive-i)*nalpha;
    nc_assert( targetrowbegin + (nalpha-1) < complete_sab.end() );
    std::copy( srcrowbegin, srcrowbegin + nalpha, targetrowbegin );
    targetrowbegin += nalpha;
  }

  //Step 2.2: copy over the original values for beta>=0:
  complete_sab.insert( complete_sab.end(),
                       sab_for_halfbetagrid.begin(),
                       sab_for_halfbetagrid.end() );

  //Validate:
  nc_assert_always( complete_betagrid.size() == nbeta_new );
  nc_assert_always( complete_sab.size() == nalpha * nbeta_new );

}
