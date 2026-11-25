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

#include "NCrystal/internal/cif/NCCIFSites.hh"
#include "NCrystal/internal/cif/NCCIFSetting.hh"
#include "NCrystal/internal/cif/NCCIFUtils.hh"
#include "NCrystal/internal/utils/NCAtomUtils.hh"
#include "NCrystal/internal/utils/NCFMA.hh"
#include <unordered_map>

namespace NC = NCrystal;
namespace NCC = NCrystal::CIF;

namespace NCRYSTAL_NAMESPACE {
  namespace CIF {
    namespace {

      using Uij = UijUtils::Uij;
      using Groups = std::vector<std::vector<std::size_t>>;//(site indices)

      //Value i of an optional column, or a StrView without value if the
      //column is absent or the value is "?" or ".":
      StrView colValue( const ValueColumn& c, std::size_t i )
      {
        if ( c.empty() || isNull( c[i] ) )
          return StrView();
        return c[i].text;
      }

      Optional<CIFNumber> number( StrView s )
      {
        return s.has_value() ? parseCIFNumber( s ) : NullOpt;
      }

      //v*num/den, deterministic on all platforms and correctly rounded except
      //possibly very near a tie (exact product and division remainder, cf.
      //ncfma):
      double mulDiv( double v, double num, double den )
      {
        const double p = v * num;
        const double pe = ncfma( v, num, -p );
        const double q = p / den;
        const double r = ncfma( -q, den, p );
        return q + ( r + pe ) / den;
      }

      //Ustd with Ueq=ueq rescaled to Ueq=target (unchanged if equal, or if ueq
      //is not positive, leaving the error for later):
      Uij rescaled( const Uij& u, double ueq, double target )
      {
        if ( !( ueq > 0.0 ) || ueq == target )
          return u;
        Uij res;
        for ( auto i : ncrange( 6 ) )
          res[i] = mulDiv( u[i], target, ueq );
        return res;
      }

      //Periodic distance between fractional positions (cf. merging below):
      double cellDistance( const Vector& p1, const Vector& p2 )
      {
        double sq = 0.0;
        for ( auto i : ncrange( 3 ) ) {
          const double e = ncabs( ( p1[i] - std::floor( p1[i] ) )
                                  - ( p2[i] - std::floor( p2[i] ) ) );
          const double m = ncmin( e, 1.0 - e );
          sq = ncfma( m, m, sq );
        }
        return std::sqrt( sq );
      }

      //A site with positive occupancy, expanded:
      struct Site {
        std::string label, qlabel;//(qlabel: unique, for messages)
        std::uint32_t line;//(of the site row)
        std::string marker;
        Vector pos;//as given (wrapped)
        double occupancy;//(at most 1)
        Optional<PairDD> uiso;//(value and rounding error)
        Optional<Uij> uij, uijerr;
        SGSiteOrbit orbit;
        Optional<double> finalUiso;
        Optional<Uij> finalUij;
        std::vector<Uij> finalUijs;
      };

      //A site for messages, e.g. 'CIF site "Al1" in "x.cif" (data_a, line 7)':
      std::string siteRef( BlockView b, const Site& s )
      {
        return "CIF site \"" + s.label + "\" in "
          + describeLocation( b, s.line );
      }

      //Several sites for messages, e.g. 'CIF sites "A" (line 7) and "B" (line
      //9) in "x.cif" (data_a)' (without the last part if !withLocation):
      std::string sitesRef( BlockView b, const std::vector<const Site*>& v,
                            bool withLocation = true )
      {
        if ( v.size() == 1 && withLocation )
          return siteRef( b, *v.front() );
        std::string res = ( v.size() == 1 ? "CIF site " : "CIF sites " );
        for ( auto i : ncrange( v.size() ) )
          res += ( i == 0 ? "" : ( i + 1 == v.size() ? " and " : ", " ) )
            + ( "\"" + v[i]->label + "\" (line " + std::to_string( v[i]->line )
                + ")" );
        return withLocation ? res + " in " + describeLocation( b ) : res;
      }

      //Factor converting B to U (B = 8 pi^2 U):
      constexpr double bToU = 1.0 / ( 8.0 * kPi * kPi );

      //The six texts of the anisotropic displacements of each site (rows
      //matched to sites by label, in order, as some files reuse labels):
      std::vector<std::array<StrView,6>> anisoTexts( const StructureItems& s )
      {
        std::vector<std::array<StrView,6>> res( s.sites.label.size() );
        if ( !s.aniso.has_value() )
          return res;
        const auto& a = s.aniso.value();
        std::vector<bool> taken( res.size(), false );
        auto itSkip = a.skippedRows.begin();
        for ( auto r : ncrange( a.label.size() ) ) {
          if ( itSkip != a.skippedRows.end() && *itSkip == r ) {
            ++itSkip;
            continue;
          }
          for ( auto i : ncrange( res.size() ) ) {
            if ( taken[i] || s.sites.label[i].text != a.label[r].text )
              continue;
            taken[i] = true;
            for ( auto j : ncrange( 6 ) )
              res[i][j] = a.components[j][r].text;
            break;
          }
        }
        return res;
      }

      //Element marker of each site (e.g. "Al" or "H2"), from its type symbol,
      //or else from its label checked against the formula (with a warning):
      void setMarkers( BlockView b, std::vector<Site>& sites,
                       const std::vector<StrView>& typeSymbols,
                       StrView formula, const CIFWarnFct& warnfct )
      {
        auto marker = []( const CIFElement& e ) {
          return elementZToName( e.Z ) + ( e.A ? std::to_string( e.A )
                                           : std::string() );
        };
        std::vector<std::pair<std::string,std::string>> examples;
        for ( auto i : ncrange( sites.size() ) ) {
          const auto ts = vectAt( typeSymbols, i );
          auto e = elementFromCIFTypeSymbol( ts.has_value() && !ts.empty()
                                             ? ts : StrView::make("?") );
          if ( e.has_value() ) {
            sites[i].marker = marker( e.value() );
            continue;
          }
          const auto& lbl = sites[i].label;
          e = elementFromCIFLabel( lbl, formula.has_value() ? formula
                                   : StrView::make("") );
          if ( !e.has_value() ) {
            if ( formula.has_value() )
              NCRYSTAL_THROW2( BadInput, siteRef( b, sites[i] ) << " has no"
                               " usable _atom_site_type_symbol, and its label"
                               " does not unambiguously start with an element"
                               " of the chemical formula (\"" << formula
                               << "\")" );
            NCRYSTAL_THROW2( BadInput, siteRef( b, sites[i] ) << " has no"
                             " usable _atom_site_type_symbol, and without a"
                             " _chemical_formula_sum its label must be an"
                             " element name followed by nothing but digits"
                             " and parentheses (like \"Al\", \"Fe2\" or"
                             " \"O(1)\")" );
          }
          sites[i].marker = marker( e.value() );
          bool known = false;
          for ( auto& ex : examples )
            known = known || ( ex.first == sites[i].marker );
          if ( !known )
            examples.emplace_back( sites[i].marker, lbl );
        }
        if ( examples.empty() )
          return;
        std::string ex;
        const std::size_t nex = ( examples.size() > 4 ? 4 : examples.size() );
        for ( auto i : ncrange( nex ) )
          ex += ( i ? ", " : "" ) + examples[i].first + " from "
            + examples[i].second;
        if ( examples.size() > 4 )
          ex += ", ...";
        warnfct( "Some CIF sites in " + describeLocation( b ) + " have no"
                 " _atom_site_type_symbol, so their elements were taken from"
                 " their labels (" + ex + ")" );
      }

      //Uiso of a site with a tensor, from its Uiso (if given) and the Ueq of
      //its tensor (both with rounding errors). Differences explained by the
      //rounding errors plus tolUij are accepted:
      double resolveUisoUeq( const std::string& site,
                             const Optional<PairDD>& uiso, const PairDD& ueq,
                             double tolUij, bool permissive,
                             const CIFWarnFct& warnfct )
      {
        if ( !uiso.has_value() )
          return ueq.first;
        const auto& u = uiso.value();
        const double allow = u.second + ueq.second
          + tolUij * ncmax( ueq.first, u.first );
        if ( ncabs( ueq.first - u.first ) > allow ) {
          //%g, since the values can be computed (Ueq, or Uiso from Biso):
          const std::string msg
            = ( site + " has Uiso=" + dbl2shortstr( u.first, "%g" ).to_string()
                + " which differs from the Ueq="
                + dbl2shortstr( ueq.first, "%g" ).to_string()
                + " of its anisotropic displacements by more than the"
                " rounding errors of the printed values and tol_uij allow" );
          if ( !permissive )
            NCRYSTAL_THROW2( BadInput, msg << " (try permissive mode to use"
                             " the Uiso value, or a larger tol_uij)" );
          warnfct( msg + " (using the Uiso value due to permissive mode)" );
        }
        return u.first;
      }

      //Total number of atoms in the cell stated in the data (0 if not, or
      //not usable):
      double expectedAtomCount( const StructureItems& s )
      {
        StableSum n;
        const auto& col = s.atomTypeNumberInCell;
        for ( auto i : ncrange( col.size() ) ) {
          const auto& v = col[i];
          if ( isNull( v ) )
            return 0.0;
          const auto x = parseCIFNumber( v.text );
          if ( !x.has_value() )
            return 0.0;
          n.add( x.value().value );
        }
        return n.sum();
      }

      //Grid of positions (fractional, wrapped) for finding those within 1e-4
      //(cf. cellDistance) of a given one:
      class PositionGrid final {
      public:
        static constexpr int nbins = 1024;//(bin width above 1e-4)
        void add( const Vector& p, std::size_t idx )
        {
          m_bins[ key( bin( p[0] ), bin( p[1] ), bin( p[2] ) ) ]
            .push_back( idx );
        }
        //Indices of positions in the bins around p (including more distant
        //ones):
        template<class Fct>
        void forNear( const Vector& p, Fct&& f ) const
        {
          const int b[3] = { bin( p[0] ), bin( p[1] ), bin( p[2] ) };
          for ( int dx = -1; dx <= 1; ++dx )
            for ( int dy = -1; dy <= 1; ++dy )
              for ( int dz = -1; dz <= 1; ++dz ) {
                auto it = m_bins.find( key( b[0] + dx, b[1] + dy,
                                            b[2] + dz ) );
                if ( it != m_bins.end() )
                  for ( auto idx : it->second )
                    f( idx );
              }
        }
      private:
        static int bin( double x )
        {
          const int i = static_cast<int>( ( x - std::floor( x ) ) * nbins );
          return i < 0 ? 0 : ( i >= nbins ? nbins - 1 : i );
        }
        static std::uint64_t key( int x, int y, int z )
        {
          auto w = []( int i ) {
            return static_cast<std::uint64_t>( ( i % nbins + nbins ) % nbins );
          };
          return ( w( x ) * nbins + w( y ) ) * nbins + w( z );
        }
        std::unordered_map<std::uint64_t,std::vector<std::size_t>> m_bins;
      };

      //Sites at the same positions merged into atoms, with the composition
      //and total occupancy, checked against the expected number of positions
      //(if not 0). The sites of each atom are added to groups:
      std::vector<CIFAtom> mergeSites( BlockView b,
                                       const std::vector<Site>& sites,
                                       double expectedCount,
                                       const SGSymmetry& sym,
                                       const CellParameters& cell,
                                       Groups& groups )
      {
        //(positions within this distance are the same, cf. NCInfoBuilder):
        constexpr double posTolerance = 0.0001;
        PositionGrid grid;
        std::vector<std::pair<std::size_t,std::size_t>> allpos;
        for ( auto i : ncrange( sites.size() ) )
          for ( auto j : ncrange( sites[i].orbit.positions.size() ) ) {
            grid.add( sites[i].orbit.positions[j], allpos.size() );
            allpos.emplace_back( i, j );
          }
        std::vector<bool> done( sites.size(), false );
        std::vector<CIFAtom> res;
        //(atoms of several sites, with their first site):
        std::vector<std::pair<std::size_t,std::size_t>> reexpand;
        for ( auto i : ncrange( sites.size() ) ) {
          if ( done[i] )
            continue;
          done[i] = true;
          const auto& s1 = sites[i];
          const auto& p1 = s1.orbit.positions.front();
          std::vector<std::size_t> group;
          grid.forNear( p1, [&]( std::size_t idx ) {
            const auto j = allpos[idx].first;
            if ( done[j] || sites[j].orbit.positions.size()
                 != s1.orbit.positions.size() )
              return;
            if ( cellDistance( sites[j].orbit.positions[allpos[idx].second],
                               p1 ) < posTolerance )
              group.push_back( j );
          } );
          std::sort( group.begin(), group.end() );
          group.erase( std::unique( group.begin(), group.end() ),
                       group.end() );
          group.insert( group.begin(), i );
          for ( auto j : group ) {
            const auto& s2 = sites[j];
            if ( s1.finalUiso != s2.finalUiso || s1.finalUij != s2.finalUij )
              NCRYSTAL_THROW2( BadInput, sitesRef( b, { &s1, &s2 } ) << " are"
                               " at the same positions but have different"
                               " displacement parameters (Uiso or tensors),"
                               " which is not supported" );
            done[j] = true;
          }
          std::vector<const Site*> gsites;
          for ( auto j : group )
            gsites.push_back( &sites[j] );
          CIFAtom a;
          StableSum occsum;
          for ( auto j : group ) {
            a.labels.push_back( sites[j].label );
            a.composition.emplace_back( sites[j].occupancy, sites[j].marker );
            occsum.add( sites[j].occupancy );
          }
          double occupancy = occsum.sum();
          if ( ncabs( occupancy - 1.0 ) <= 1e-12 )
            occupancy = 1.0;//numerical noise
          nc_assert( occupancy > 0.0 );
          if ( occupancy != 1.0 ) {
            for ( auto& c : a.composition )
              c.first /= occupancy;
            if ( occupancy > 1.0 + 1e-6 ) {
              std::string m = ( sitesRef( b, gsites ) + " are at the same"
                                " positions with a total occupancy above 1 ("
                                + dbl2shortstr( occupancy, "%g" ).to_string()
                                + ")" );
              if ( ncabs( std::floor( occupancy + 0.5 ) - occupancy ) < 1e-10 )
                m += ( ", probably the same site listed twice (e.g. at"
                       " symmetry-equivalent coordinates)" );
              NCRYSTAL_THROW( BadInput, m );
            }
            occupancy = ncmin( 1.0, occupancy );
          }
          std::sort( a.composition.begin(), a.composition.end() );
          a.occupancy = occupancy;
          a.uiso = s1.finalUiso;
          a.uij = s1.finalUij;
          a.positions = s1.orbit.positions;
          a.uijs = s1.finalUijs;
          a.splitPositions = ( s1.orbit.nearCopies
                               ? s1.orbit.nearCopies + 1 : 0 );
          a.splitDistance = s1.orbit.nearestCopyDistance;
          if ( group.size() > 1 )
            reexpand.emplace_back( res.size(), i );
          res.push_back( std::move( a ) );
          groups.push_back( std::move( group ) );
        }
        //Atoms of several sites are expanded again with the total occupancy
        //(cf. the split site rule of expandSiteToOrbit):
        for ( auto& e : reexpand ) {
          auto& a = res[e.first];
          const auto o = expandSiteToOrbit( sym, cell, sites[e.second].pos,
                                            a.occupancy,
                                            SGSiteTolerances::lenient() );
          a.splitPositions = ( o.nearCopies ? o.nearCopies + 1 : 0 );
          a.splitDistance = o.nearestCopyDistance;
        }
        if ( expectedCount ) {
          std::size_t n = 0;
          for ( auto& a : res )
            n += a.positions.size();
          const double nd = static_cast<double>( n );
          if ( ncabs( nd - expectedCount )
               / ncmax( 1e-300, ncabs( nd ) + ncabs( expectedCount ) ) > 0.01 )
            NCRYSTAL_THROW2( BadInput, "The unit cell has " << n << " atoms in"
                             " space group setting "
                             << sym.spaceGroup().toString() << ", but the"
                             " CIF data in " << describeLocation( b )
                             << " states " << fmtg( expectedCount )
                             << " (_atom_type_number_in_cell), so the space"
                             " group or its setting is probably wrong" );
        }
        return res;
      }

      //Setting with more operations (in the same coordinates) whose symmetry
      //the cell and the atoms also have, within tolSymprec (Aa, applied
      //conservatively as fractional deviation times the longest cell axis):
      //the highest order, the default setting among equal orders, else the
      //first in table order. NullOpt if none:
      Optional<SpaceGroup> higherSymmetry( const Document& doc,
                                           const StructureItems& s,
                                           const SpaceGroup& sg,
                                           const CellParameters& cell,
                                           const std::vector<CIFAtom>& atoms,
                                           double tolSymprec )
      {
        //Positions grouped by species (composition and occupancy):
        std::vector<std::vector<Vector>> species;
        std::vector<std::size_t> speciesOf;
        for ( auto i : ncrange( atoms.size() ) ) {
          std::size_t isp = 0;
          while ( isp < i && !( atoms[isp].composition
                                == atoms[i].composition
                                && atoms[isp].occupancy
                                == atoms[i].occupancy ) )
            ++isp;
          if ( isp == i ) {
            speciesOf.push_back( species.size() );
            species.emplace_back();
          } else {
            speciesOf.push_back( speciesOf[isp] );
          }
          auto& v = species[speciesOf.back()];
          v.insert( v.end(), atoms[i].positions.begin(),
                    atoms[i].positions.end() );
        }
        const double maxlen = ncmax( cell.a, ncmax( cell.b, cell.c ) );
        Optional<SpaceGroup> best;
        unsigned bestOrder = 0;
        for ( auto& cand : settingsWithMoreOperations( sg ) ) {
          const unsigned order = SGSymmetry::get( cand ).order();
          const bool betterDefault = ( cand.isDefaultSetting()
                                       && best.has_value()
                                       && !best.value().isDefaultSetting() );
          if ( order < bestOrder || ( order == bestOrder && !betterDefault ) )
            continue;
          if ( !snapCell( doc, s, cand ).cell.has_value() )
            continue;
          if ( !( checkStructureSymmetry( cand, species ).maxDeviation
                  * maxlen <= tolSymprec ) )
            continue;
          best = cand;
          bestOrder = order;
        }
        return best;
      }

      //The atoms regrouped into the orbits of the higher symmetry setting
      //(with the sites of each in groups), or a description of the
      //inconsistent displacement parameters (in which case the returned
      //atoms are empty):
      struct Regrouped {
        std::vector<CIFAtom> atoms;
        Groups groups;
        std::string inconsistency;
      };
      Regrouped regroup( BlockView b, const std::vector<Site>& sites,
                         const std::vector<CIFAtom>& atoms,
                         const Groups& groups, const SGSymmetry& sym,
                         const CellParameters& cell,
                         const UijUtils::CellMetric& metric, double tolUij )
      {
        Regrouped res;
        std::vector<bool> done( atoms.size(), false );
        for ( auto i : ncrange( atoms.size() ) ) {
          if ( done[i] )
            continue;
          //Positions of the atoms of the same species:
          std::vector<std::pair<std::size_t,std::size_t>> refs;
          std::vector<Vector> pos;
          for ( auto j : ncrange( atoms.size() ) ) {
            if ( done[j] || atoms[j].composition != atoms[i].composition
                 || atoms[j].occupancy != atoms[i].occupancy )
              continue;
            done[j] = true;
            for ( auto k : ncrange( atoms[j].positions.size() ) ) {
              refs.emplace_back( j, k );
              pos.push_back( atoms[j].positions[k] );
            }
          }
          for ( auto& orbit : groupPositionsIntoOrbits(
                  sym, cell, pos, atoms[i].occupancy,
                  SGSiteTolerances::lenient() ) ) {
            //The old atoms in the orbit (in order):
            std::vector<std::size_t> olds;
            for ( auto m : orbit )
              if ( std::find( olds.begin(), olds.end(), refs[m].first )
                   == olds.end() )
                olds.push_back( refs[m].first );
            std::sort( olds.begin(), olds.end() );
            const auto& a0 = atoms[olds.front()];
            CIFAtom a;
            Groups::value_type g;
            std::vector<const Site*> gsites;
            for ( auto io : olds ) {
              a.labels.insert( a.labels.end(), atoms[io].labels.begin(),
                               atoms[io].labels.end() );
              g.insert( g.end(), groups[io].begin(), groups[io].end() );
            }
            for ( auto j : g )
              gsites.push_back( &sites[j] );
            a.composition = a0.composition;
            a.occupancy = a0.occupancy;
            const Vector& rep = pos[orbit.front()];
            const auto o = expandSiteToOrbit( sym, cell, rep, a.occupancy,
                                              SGSiteTolerances::lenient() );
            a.positions = o.positions;
            a.splitPositions = ( o.nearCopies ? o.nearCopies + 1 : 0 );
            a.splitDistance = o.nearestCopyDistance;
            //Displacements, which must agree within tolUij:
            const auto sref = sitesRef( b, gsites, false );
            VectD uisos;
            bool anyUiso = false, anyUij = false, allUij = true;
            for ( auto io : olds ) {
              anyUiso = anyUiso || atoms[io].uiso.has_value();
              if ( atoms[io].uiso.has_value() )
                uisos.push_back( atoms[io].uiso.value() );
              anyUij = anyUij || atoms[io].uij.has_value();
              allUij = allUij && atoms[io].uij.has_value();
            }
            if ( anyUiso && uisos.size() != olds.size() ) {
              res.inconsistency = "only some of " + sref + " have Uiso values";
              return res;
            }
            if ( anyUij && !allUij ) {
              res.inconsistency = ( "only some of " + sref + " have"
                                    " anisotropic displacements" );
              return res;
            }
            if ( !uisos.empty() ) {
              const auto mm = std::minmax_element( uisos.begin(),
                                                   uisos.end() );
              if ( *mm.second - *mm.first > tolUij * ncabs( *mm.second ) ) {
                res.inconsistency = sref + " have different Uiso values";
                return res;
              }
              StableSum sum;
              for ( auto u : uisos )
                sum.add( u );
              a.uiso = ( *mm.first == *mm.second ? *mm.first
                         : sum.sum() / static_cast<double>( uisos.size() ) );
            }
            if ( allUij && anyUij ) {
              std::vector<Vector> opos;
              std::vector<Uij> ou;
              for ( auto m : orbit ) {
                opos.push_back( pos[m] );
                ou.push_back( vectAt( atoms[refs[m].first].uijs,
                                      refs[m].second ) );
              }
              auto ur = UijUtils::tryUstdOfOrbit( sym, metric, opos, ou,
                                                  sref, tolUij );
              if ( !ur.uij.has_value() ) {
                //(e.g. "the uij values of the site ... do not obey the site
                //symmetry"):
                res.inconsistency = ur.errorMsg.value();
                auto& c = res.inconsistency.front();
                if ( c >= 'A' && c <= 'Z' )
                  c = static_cast<char>( c - 'A' + 'a' );
                return res;
              }
              Uij u = ur.uij.value();
              const double ueq = UijUtils::calcUeqFromUstar(
                UijUtils::convertUstdToUstar( u, metric ), metric );
              u = rescaled( u, ueq, a.uiso.value() );
              a.uij = u;
              for ( auto& p : a.positions )
                a.uijs.push_back( UijUtils::mapUstdToPosition(
                                    sym, metric, a.positions.front(), p, u ) );
            }
            res.atoms.push_back( std::move( a ) );
            res.groups.push_back( std::move( g ) );
          }
        }
        //(in the order of the first site of each):
        std::vector<std::size_t> order( res.atoms.size() );
        for ( auto i : ncrange( order.size() ) )
          order[i] = i;
        std::sort( order.begin(), order.end(),
                   [&res]( std::size_t x, std::size_t y )
                   { return res.groups[x].front() < res.groups[y].front(); } );
        Regrouped sorted;
        for ( auto i : order ) {
          sorted.atoms.push_back( std::move( res.atoms[i] ) );
          sorted.groups.push_back( std::move( res.groups[i] ) );
        }
        return sorted;
      }
    }
  }
}

NCC::CIFAtoms NCC::getAtoms( const Document& doc, const StructureItems& s,
                             const SpaceGroup& sg, CIFPermissive perm,
                             double tolUij, double tolSymprec,
                             const CIFWarnFct& warnfct )
{
#ifndef NDEBUG
  nc_assert_always( bool(warnfct) );
  nc_assert_always( tolUij > 0.0 && std::isfinite( tolUij ) );
  nc_assert_always( tolSymprec >= 0.0 && std::isfinite( tolSymprec ) );
#endif
  const bool permissive = ( perm == CIFPermissive::YES );
  const auto b = s.block( doc );
  const auto& cols = s.sites;
  const auto aniso = anisoTexts( s );
  //(tensors given as B are converted to U, with their rounding errors):
  const double anisoFactor = ( s.aniso.has_value()
                               && s.aniso.value().form
                               == StructureItems::AnisoForm::B ? bToU : 1.0 );

  //The sites with positive occupancy:
  std::vector<Site> sites;
  std::vector<StrView> typeSymbols;
  for ( auto i : ncrange( cols.label.size() ) ) {
    Site site;
    site.label = cols.label[i].text.to_string();
    site.line = cols.label[i].line;
    const auto occtext = colValue( cols.occupancy, i );
    PairDD occ{ 1.0, 0.0 };
    if ( occtext.has_value() ) {
      const auto x = parseCIFNumber( occtext );
      if ( !x.has_value() )
        NCRYSTAL_THROW2( BadInput, siteRef( b, site ) << " has an invalid"
                         " occupancy (\"" << occtext << "\")" );
      occ = { x.value().value, x.value().roundingError };
    }
    if ( !( occ.first > 0.0 ) )
      continue;
    if ( occ.first > 1.0 + 1e-6 )
      NCRYSTAL_THROW2( BadInput, siteRef( b, site ) << " has an occupancy"
                       " above 1 (" << fmt( occ.first ) << ")" );
    double xyz[3];
    const ValueColumn* xyzcols[3] = { &cols.fractX, &cols.fractY,
                                      &cols.fractZ };
    for ( auto j : ncrange( 3 ) ) {
      const auto x = number( colValue( *xyzcols[j], i ) );
      if ( !x.has_value() || !std::isfinite( x.value().value ) )
        NCRYSTAL_THROW2( BadInput, siteRef( b, site ) << " has missing or"
                         " invalid fractional coordinates" );
      const double r = x.value().value - std::floor( x.value().value );
      xyz[j] = ( r >= 1.0 || r == 0.0 ) ? 0.0 : r;
    }
    site.pos = Vector( xyz[0], xyz[1], xyz[2] );
    if ( vectAt( aniso, i )[0].has_value() ) {
      Uij u, e;
      bool allzero = true;
      for ( auto j : ncrange( 6 ) ) {
        const auto x = parseCIFNumber( vectAt( aniso, i )[j] );
        if ( !x.has_value() )
          NCRYSTAL_THROW2( BadInput, siteRef( b, site ) << " has invalid"
                           " anisotropic displacement values" );
        u[j] = x.value().value * anisoFactor;
        e[j] = x.value().roundingError * anisoFactor;
        allzero = allzero && ( u[j] == 0.0 );
      }
      if ( !allzero ) {//(all zero: a placeholder for no tensor)
        site.uij = u;
        site.uijerr = e;
      }
    }
    auto u = number( colValue( cols.uIsoOrEquiv, i ) );
    if ( u.has_value() && u.value().value > 0.0 ) {
      site.uiso = PairDD{ u.value().value, u.value().roundingError };
    } else {
      const auto bi = number( colValue( cols.bIsoOrEquiv, i ) );
      if ( bi.has_value() && bi.value().value > 0.0 ) {
        site.uiso = PairDD{ bi.value().value * bToU,
                            bi.value().roundingError * bToU };
      }
    }
    site.occupancy = ncmin( occ.first, 1.0 );
    //(unique labels without commas, for messages):
    std::string q = site.label.empty() ? std::string("?") : site.label;
    std::replace( q.begin(), q.end(), ',', ';' );
    for ( bool clash = true; clash; ) {
      clash = false;
      for ( auto& e : sites )
        clash = clash || ( e.qlabel == q );
      if ( clash )
        q += '\'';
    }
    site.qlabel = std::move( q );
    typeSymbols.push_back( colValue( cols.typeSymbol, i ) );
    sites.push_back( std::move( site ) );
  }
  if ( sites.empty() )
    NCRYSTAL_THROW2( BadInput, "The CIF data in " << describeLocation( b )
                     << " has no atom sites with positive occupancy" );
  setMarkers( b, sites, typeSymbols,
              s.formulaSum.has_value() && !isNull( s.formulaSum.value() )
              ? s.formulaSum.value().text : StrView(), warnfct );

  //Cell (with rounding errors, except for angles printed as 90 or 120
  //without uncertainty, which are exact):
  double cv[6], ce[6];
  for ( auto j : ncrange( 6 ) ) {
    const auto x = parseCIFNumber( s.cell[j].text );
    if ( !x.has_value() || !( x.value().value > 0.0 ) )
      NCRYSTAL_THROW2( BadInput, "Invalid cell parameter in CIF data in "
                       << describeLocation( b, s.cell[j].line ) << ": \""
                       << s.cell[j].text << "\"" );
    cv[j] = x.value().value;
    const bool exact = ( j >= 3 && ( cv[j] == 90.0 || cv[j] == 120.0 )
                         && !s.cell[j].text.contains( '(' ) );
    ce[j] = exact ? 0.0 : x.value().roundingError;
  }
  const CellParameters cell{ cv[0], cv[1], cv[2], cv[3], cv[4], cv[5] };
  const CellParameters cellerr{ ce[0], ce[1], ce[2], ce[3], ce[4], ce[5] };
  const auto metric = UijUtils::cellMetric( cell.a, cell.b, cell.c,
                                            cell.alpha, cell.beta,
                                            cell.gamma );
  const SGSymmetry& sym = SGSymmetry::get( sg );

  //Expand all sites (with their tensors snapped to the site symmetry):
  std::vector<Optional<Uij>> snapped;
  std::vector<PairDD> ueqs;
  for ( auto& site : sites ) {
    site.orbit = expandSiteToOrbit( sym, cell, site.pos, site.occupancy,
                                    SGSiteTolerances::lenient() );
    snapped.emplace_back();
    ueqs.emplace_back( 0.0, 0.0 );
    if ( !site.uij.has_value() )
      continue;
    std::ostringstream descr;
    descr << '"' << site.qlabel << "\" (" << fmt( site.pos.x(), "%.6g" )
          << ", " << fmt( site.pos.y(), "%.6g" ) << ", "
          << fmt( site.pos.z(), "%.6g" ) << ")";
    snapped.back() = UijUtils::snapUstdToSiteSymmetry(
      sym, metric, site.orbit.positions.front(), site.uij.value(),
      descr.str(), 0.05, site.uijerr.value() );
    ueqs.back() = UijUtils::calcUeqWithRoundingError( site.uij.value(), cell,
                                                      site.uijerr.value(),
                                                      cellerr );
  }

  //A single magnitude of the displacements of each site (the uiso, with the
  //tensors rescaled to it):
  for ( auto i : ncrange( sites.size() ) ) {
    auto& site = sites[i];
    if ( site.uiso.has_value() )
      site.finalUiso = site.uiso.value().first;
    if ( !snapped[i].has_value() )
      continue;
    const auto& ueq = ueqs[i];
    const double uiso = resolveUisoUeq( siteRef( b, site ), site.uiso, ueq,
                                        tolUij, permissive, warnfct );
    site.finalUiso = uiso;
    site.finalUij = rescaled( snapped[i].value(), ueq.first, uiso );
    for ( auto& p : site.orbit.positions )
      site.finalUijs.push_back(
        rescaled( UijUtils::mapUstdToPosition( sym, metric,
                                               site.orbit.positions.front(),
                                               p, snapped[i].value() ),
                  ueq.first, uiso ) );
  }

  Groups groups;
  auto atoms = mergeSites( b, sites, expectedAtomCount( s ), sym, cell,
                           groups );

  //A higher symmetry of the atoms in the same coordinates is adopted, unless
  //their displacements contradict it:
  SpaceGroup finalSG = sg;
  const auto hs = ( tolSymprec > 0.0
                    ? higherSymmetry( doc, s, sg, cell, atoms, tolSymprec )
                    : NullOpt );
  if ( hs.has_value() ) {
    auto rg = regroup( b, sites, atoms, groups,
                             SGSymmetry::get( hs.value() ), cell, metric,
                             tolUij );
    const std::string given = sg.toString();
    std::string msg = ( "The atoms of the CIF data in " + describeLocation( b )
                        + " have the higher symmetry of space group setting "
                        + hs.value().toString() + " (in the same coordinates"
                        " as the given setting " + given + ")" );
    if ( !rg.inconsistency.empty() ) {
      msg += ", but " + rg.inconsistency;
      if ( !permissive )
        NCRYSTAL_THROW2( BadInput, msg << " (try permissive mode to keep"
                         " setting " << given << ")" );
      warnfct( msg + " (setting " + given + " kept due to permissive"
               " mode)" );
    } else {
      warnfct( msg + ", which is used" );
      atoms = std::move( rg.atoms );
      groups = std::move( rg.groups );
      finalSG = hs.value();
    }
  }

  //Split sites:
  for ( auto ia : ncrange( atoms.size() ) ) {
    const auto& a = atoms[ia];
    if ( !a.splitPositions )
      continue;
    std::vector<const Site*> gsites;
    for ( auto j : groups[ia] )
      gsites.push_back( &sites[j] );
    const bool several = ( gsites.size() > 1 );
    std::ostringstream msg;
    msg << sitesRef( b, gsites ) << ( several ? " are" : " is" ) << " a split"
        " site: " << a.splitPositions << " of " << ( several ? "their" : "its" )
        << " symmetry-equivalent positions are "
        << fmt( a.splitDistance, "%.2g" ) << " Aa apart (occupancy "
        << fmt( a.occupancy, "%g" ) << "). NCrystal treats all vacancies as"
        " uncorrelated, so incoherent elastic scattering is overestimated for"
        " this material";
    if ( !permissive )
      NCRYSTAL_THROW2( BadInput, msg.str() << " (try permissive mode to load"
                       " it anyway)" );
    warnfct( msg.str() + " (accepted due to permissive mode)" );
  }
  return { finalSG, std::move( atoms ) };
}
