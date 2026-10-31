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

#include "NCrystal/interfaces/NCProcImpl.hh"
namespace NC = NCrystal;

namespace {

  class MyIsoScat final : public NC::ProcImpl::ScatterIsotropicMat {
  public:
    const char* name() const noexcept final { return "MyIsoScat"; }
    NC::CrossSect crossSectionIsotropic(NC::CachePtr&, NC::NeutronEnergy ) const final
    {
      return NC::CrossSect{1.0};
    }
    NC::ScatterOutcomeIsotropic sampleScatterIsotropic(NC::CachePtr&, NC::RNG&, NC::NeutronEnergy ) const final
    {
      return { NC::NeutronEnergy{0.0}, NC::CosineScatAngle{0.5} };
    }
  };

  //Isotropic process with a configurable domain, keeping the default
  //isoXSDescr():
  class MyIsoScatDomain final : public NC::ProcImpl::ScatterIsotropicMat {
  public:
    MyIsoScatDomain( NC::EnergyDomain d ) : m_d( d ) {}
    const char* name() const noexcept final { return "MyIsoScatDomain"; }
    NC::EnergyDomain domain() const noexcept final { return m_d; }
    NC::CrossSect crossSectionIsotropic(NC::CachePtr&, NC::NeutronEnergy ) const final
    {
      return NC::CrossSect{1.0};
    }
    NC::ScatterOutcomeIsotropic sampleScatterIsotropic(NC::CachePtr&, NC::RNG&, NC::NeutronEnergy ekin ) const final
    {
      return { ekin, NC::CosineScatAngle{0.5} };
    }
  private:
    NC::EnergyDomain m_d;
  };

  //Isotropic process overriding isoXSDescr() with a table (xs = 2 barn):
  class MyIsoScatTabulated final : public NC::ProcImpl::ScatterIsotropicMat {
  public:
    const char* name() const noexcept final { return "MyIsoScatTabulated"; }
    NC::CrossSect crossSectionIsotropic(NC::CachePtr&, NC::NeutronEnergy ) const final
    {
      return NC::CrossSect{2.0};
    }
    NC::ScatterOutcomeIsotropic sampleScatterIsotropic(NC::CachePtr&, NC::RNG&, NC::NeutronEnergy ekin ) const final
    {
      return { ekin, NC::CosineScatAngle{0.5} };
    }
    NC::IsoXSDescr isoXSDescr() const override
    {
      NC::IsoXSDescr d;
      d.breakpts = { 0.5 };
      NC::IsoXSDescr::Tabulation t;
      t.E = { 1.0, 2.0 };
      t.p = { 2.0, 2.0, 0.0, 0.0 };//e*xs = 2*e = 2 + 2*u on [1,2]
      t.tail[1] = 2.0;
      t.tail[3] = 2.0;
      d.table = std::move( t );
      return d;
    }
  };

  class MyAnisoScat final : public NC::ProcImpl::ScatterAnisotropicMat {
  public:
    const char* name() const noexcept final { return "MyAnisoScat"; }

    NC::CrossSect crossSection(NC::CachePtr&, NC::NeutronEnergy, const NC::NeutronDirection& ) const final
    {
      return NC::CrossSect{1.0};
    }

    NC::ScatterOutcome sampleScatter(NC::CachePtr&, NC::RNG&, NC::NeutronEnergy, const NC::NeutronDirection& ) const final
    {
      return { NC::NeutronEnergy(0.0), NC::NeutronDirection{0,0,1} };
    }

  };

  class MyIsoAbs final : public NC::ProcImpl::AbsorptionIsotropicMat {
  public:
    const char* name() const noexcept final { return "MyIsoAbs"; }
    NC::CrossSect crossSectionIsotropic(NC::CachePtr&, NC::NeutronEnergy ) const final
    {
      return NC::CrossSect{1.0};
    }
  };

}

int main() {

  //The testing is really just compile-time at this point.
  (void)std::make_shared<MyIsoScat>();
  (void)std::make_shared<MyAnisoScat>();
  (void)std::make_shared<MyIsoAbs>();
  {
    NC::ProcImpl::ProcComposition::ComponentList cl;
    cl.emplace_back( 0.5, NC::makeSO<MyIsoScat>() );
    cl.emplace_back( 2.0, NC::makeSO<MyIsoScat>() );
    (void)std::make_shared<NC::ProcImpl::ProcComposition>( std::move(cl) );
  }
  (void)std::make_shared<NC::ProcImpl::NullScatter>();
  (void)std::make_shared<NC::ProcImpl::NullAbsorption>();

  //Default isoXSDescr(): breakpts at the finite, non-zero domain bounds, no table:
  {
    using NC::NeutronEnergy;
    auto check = []( NC::EnergyDomain d, const NC::VectD& expected )
    {
      MyIsoScatDomain p( d );
      const NC::ProcImpl::ScatterIsotropicMat& base = p;
      NC::IsoXSDescr x = base.isoXSDescr();
      nc_assert_always( x.breakpts == expected );
      nc_assert_always( !x.table.has_value() );
    };
    check( { NeutronEnergy{0.0}, NeutronEnergy{NC::kInfinity} }, {} );
    check( { NeutronEnergy{1e-3}, NeutronEnergy{2.0} }, { 1e-3, 2.0 } );
    check( { NeutronEnergy{0.0}, NeutronEnergy{2.0} }, { 2.0 } );
    check( { NeutronEnergy{1e-3}, NeutronEnergy{NC::kInfinity} }, { 1e-3 } );
    check( NC::EnergyDomain::null(), {} );
    nc_assert_always( MyIsoScat().isoXSDescr().breakpts.empty() );
  }
  //An override is returned intact through the base class:
  {
    MyIsoScatTabulated p;
    const NC::ProcImpl::ScatterIsotropicMat& base = p;
    NC::IsoXSDescr x = base.isoXSDescr();
    nc_assert_always( ( x.breakpts == NC::VectD{ 0.5 } ) );
    nc_assert_always( x.table.has_value() );
    const auto& t = x.table.value();
    nc_assert_always( ( t.E == NC::VectD{ 1.0, 2.0 } ) && t.p.size() == 4 );
    nc_assert_always( t.tail[1] == 2.0 && t.tail[3] == 2.0 );
  }

}
