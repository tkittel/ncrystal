#ifndef NCrystal_IsoXSUtils_hh
#define NCrystal_IsoXSUtils_hh

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

#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  //Utilities for IsoXSDescr (NCTypes.hh). Its Tabulation t describes a cross
  //section xs(e) (energy e in eV, xs in barn) as:
  //
  //  e < E[0]:           xs = tail[0]/sqrt(e) + tail[1]
  //  E[i] <= e < E[i+1]: e*xs = P_i(u) + sqrt(e)*Q_i(u) + r[i]/e
  //  e >= E.back():      xs = tail[3] + tail[2]/e + tail[4]/e^2
  //
  //with u=(e-E[i])/(E[i+1]-E[i]), and P_i and Q_i cubics in u with
  //coefficients p[4i..4i+3] and q[4i..4i+3] (lowest order first). Empty p, q
  //or r means zero. E must be positive and strictly increasing, with at least 2
  //entries. The r[i] terms allow a sum of tabulations to be a tabulation.
  //
  //IsoXSDescr::breakpts holds the energies of known kinks and discontinuities
  //of the cross section, sorted.

  //Evaluate the tabulated cross section at e (finding the interval with a
  //binary search, so meant for occasional use):
  double evalIsoXS( const IsoXSDescr::Tabulation&, double e );

}

#endif
