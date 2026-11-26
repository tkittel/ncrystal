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

namespace NC = NCrystal;

//Used in Length::as_wavelength (checked here rather than in the header, since
//it could fail for user code built with x87 floating point on 32bit x86):
static_assert( NC::Length::angstrom == 1e-10, "" );

namespace NCRYSTAL_NAMESPACE {
  //As good a place as any for all of these SmallVector asserts (since we do not
  //have NCSmallVector.cc):
  static_assert(SmallVector_IC<int,1>::mode==SVMode::FASTACCESS,"");
  static_assert(SmallVector_IC<int,1,SVMode::FASTACCESS>::mode==SVMode::FASTACCESS,"");
  static_assert(SmallVector_IC<int,1,SVMode::LOWFOOTPRINT>::mode==SVMode::LOWFOOTPRINT,"");
  static_assert(SmallVector<int,1>::mode==SVMode::FASTACCESS,"");
  static_assert(SmallVector<int,1,SVMode::FASTACCESS>::mode==SVMode::FASTACCESS,"");
  static_assert(SmallVector<int,1,SVMode::LOWFOOTPRINT>::mode==SVMode::LOWFOOTPRINT,"");
  static_assert(std::is_copy_constructible<SmallVector_IC<int,1,SVMode::FASTACCESS>>::value,"");
  static_assert(std::is_copy_constructible<SmallVector_IC<int,1,SVMode::LOWFOOTPRINT>>::value,"");
  static_assert(std::is_copy_assignable<SmallVector_IC<int,1,SVMode::FASTACCESS>>::value,"");
  static_assert(std::is_copy_assignable<SmallVector_IC<int,1,SVMode::LOWFOOTPRINT>>::value,"");
  static_assert(!std::is_copy_constructible<SmallVector<int,1,SVMode::FASTACCESS>>::value,"");
  static_assert(!std::is_copy_constructible<SmallVector<int,1,SVMode::LOWFOOTPRINT>>::value,"");
  static_assert(!std::is_copy_assignable<SmallVector<int,1,SVMode::FASTACCESS>>::value,"");
  static_assert(!std::is_copy_assignable<SmallVector<int,1,SVMode::LOWFOOTPRINT>>::value,"");
  static_assert(std::is_move_constructible<SmallVector_IC<int,1,SVMode::FASTACCESS>>::value,"");
  static_assert(std::is_move_constructible<SmallVector_IC<int,1,SVMode::LOWFOOTPRINT>>::value,"");
  static_assert(std::is_move_assignable<SmallVector_IC<int,1,SVMode::FASTACCESS>>::value,"");
  static_assert(std::is_move_assignable<SmallVector_IC<int,1,SVMode::LOWFOOTPRINT>>::value,"");
  static_assert(std::is_move_constructible<SmallVector<int,1,SVMode::FASTACCESS>>::value,"");
  static_assert(std::is_move_constructible<SmallVector<int,1,SVMode::LOWFOOTPRINT>>::value,"");
  static_assert(std::is_move_assignable<SmallVector<int,1,SVMode::FASTACCESS>>::value,"");
  static_assert(std::is_move_assignable<SmallVector<int,1,SVMode::LOWFOOTPRINT>>::value,"");
}

NC::DataSourceName::DataSourceName()
  : m_str( [](){ static auto s_def = makeSO<std::string>(); return s_def; }() )
{
  static_assert( EnergyDomain::null().isNull(), "" );
}

std::ostream& NC::operator<<( std::ostream& os, const OrientDir& od )
{
  using Basis = OrientDir::CrysDir::Basis;
  const auto b = od.crystal.basis;
  const auto& v = od.crystal.v;
  os << ( b == Basis::Cartesian ? "@crys:"
          : ( b == Basis::HKL ? "@crys_hkl:" : "@crys_uvw:" ) )
     << fmt(v[0]) << ',' << fmt(v[1]) << ',' << fmt(v[2])
     << "@lab:" << fmt(od.lab[0]) << ',' << fmt(od.lab[1]) << ','
     << fmt(od.lab[2]);
  return os;
}

std::ostream& NC::operator<<( std::ostream& os, const LCAxisState& a )
{
  using Form = LCAxisState::Form;
  if ( a.form == Form::Primary )
    return os << "@dir1";
  os << ( a.form == Form::Cartesian ? ""
          : ( a.form == Form::HKL ? "@crys_hkl:"
              : ( a.form == Form::UVW ? "@crys_uvw:" : "@lab:" ) ) );
  return os << fmt(a.v[0]) << ',' << fmt(a.v[1]) << ',' << fmt(a.v[2]);
}

std::ostream& NC::operator<< (std::ostream& os, const DensityState& ds)
{
  if ( ds.type == DensityState::Type::SCALEFACTOR ) {
    os << fmt(ds.value)<<"x";
  } else if ( ds.type == DensityState::Type::DENSITY ) {
    os << fmt(ds.value)<<"gcm3";
  } else {
    nc_assert( ds.type == DensityState::Type::NUMBERDENSITY );
    os << fmt(ds.value)<<"perAa3";
  }
  return os;
}

