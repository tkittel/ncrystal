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

#include "NCrystal/internal/cfgutils/NCCfgVars.hh"
#include <iostream>
#include <sstream>
#include "NCrystal/internal/cfgutils/NCCfgManip.hh"

#define REQUIRE(x) nc_assert_always(x)

namespace NC = NCrystal;

namespace {

  std::string sanitise_whitespace( std::string s ) {
    std::string s2 = std::move(s);
    NC::strreplace( s2, "\r", "\\r" );
    NC::strreplace( s2, "\n", "\\n" );
    return s2;
  }

  void test_varids() {

    //unsigned enumval_p = 0;
    for (auto& p : { "bad","density","phasechoice","atomdb","dcutoff","dcutoffup","infofactory","temp",
                    "absnfactory", "bkgd", "bragg", "coh_elas", "elas", "incoh_elas", "inelas",
                    "vdoslux", "scatfactory", "dir1", "dir2", "dirtol",
                    "lcaxis", "lcmode", "mos", "mosprec", "rot", "sccutoff"} ) {
      auto pstr = NC::StrView(p);
      auto vid = NC::Cfg::varIdFromName( pstr );
      auto str2 = vid.has_value() ? NC::Cfg::varName(vid.value()) : "invalid";
      std::cout<<p<< "  ->  ";
      if ( !vid.has_value() ) {
        std::cout<<"NOTFOUND"<<std::endl;
      } else {
        std::cout<<NC::enumAsInt(vid.value())<< "  ->  "<<str2<<std::endl;
        nc_assert_always( pstr == str2 );
        //nc_assert_always( enumval_p++ == NC::enumAsInt(vid.value()) );
      }
      //nc_assert_always( enumval_p == 25 );
    }
  }

  void test_parsing() {
    using NC::Cfg::CfgManip;

    auto doTest_impl = [](const char*cstr,bool expectbad)
    {
      std::cout<<"\n------> Parsing \""<<sanitise_whitespace(cstr)<<"\":\n";
      try {
        auto l = NC::Cfg::CfgData();
        auto toplvl = CfgManip::applyStrCfg( l, cstr );
        std::cout<<"  => Parsed "<<CfgManip::size(l)<<" entries (+"<<toplvl.size()<<" toplvl entries)\n";
        std::cout<<"  => cfgstr: "<<l<<"\n";
        std::cout<<"  => .get_temp(): "<<CfgManip::get_temp(l)<<std::endl;
        std::cout<<"  => .set_temp(500)"<<std::endl;
        CfgManip::set_temp(l,NC::Temperature{500});
        std::cout<<"  => cfgstr: "<<l<<"\n";
        std::cout<<"  => .get_temp(): "<<CfgManip::get_temp(l)<<std::endl;


      } catch ( NC::Error::BadInput& e ) {
        if (expectbad) {
          std::cout<<"  => Got expected ERROR: NC::"<<e.getTypeName()<< ": "<<e.what() << std::endl;
          return;
        }
        throw;
      }
      if ( expectbad )
        NCRYSTAL_THROW(LogicError,"Did not fail as expected!!");
    };
    auto doTest = [&doTest_impl](const char*cstr) { return doTest_impl(cstr,false); };
    auto doTestBad = [&doTest_impl](const char*cstr) { return doTest_impl(cstr,true); };

    doTestBad("vdoslux=34");
    doTest("atomdb=Al is Al27");
    doTest("atomdb=Al is :: Al27");
    doTestBad("atomdb2=Al is Al27");
    doTest("; ; ;\t ;atomdb=Al is Al27@B is B10");
    doTest("; ; \r;\t ;atomdb=Al is Al27@B is B10");
    doTest("atomdb = ;");
    doTest("   atomdb =  ; dcutoff  = 0.5   ");
    doTestBad("vdoslux=2;atomdb={Al is Al27@B is B10}");
    doTestBad("vdoslux=2; ; ; ;  ;    ;;;;atomdb=  Al is Al27@B is B10 and lala ; ; ;vdoslux\n\t=34");
    doTest("vdoslux=2; ; ; ;  ;    ;;;;atomdb=  Al is Al27@B is B10 ; ; ;vdoslux\n\t=3");
    doTestBad("vdoslux=2;atomdb=Al is Al27@B is B10=Al is Al27@B is B10");
    doTestBad("vdoslux=2;atomdb==Al is Al27@B is B10");
    doTestBad("vdoslux=2;atomdb=Al is Al27@B is B10  =");
    //TODO: we could add a lot more tests here...
  }

  void test_parse_and_fill() {
    using NC::Cfg::CfgManip;
    auto doTest_impl = [](const char*cstr,bool expectbad)
    {
      std::cout<<"\n------> Parsing \""<<sanitise_whitespace(cstr)<<"\":\n";
      try {
        auto l = NC::Cfg::CfgData();
        auto toplvl = CfgManip::applyStrCfg(l,cstr);
        std::cout<<"  => Parsed "<<CfgManip::size(l)<<" entries (+"<<toplvl.size()<<" toplvl entries)\n";
        std::cout<<"  => cfgstr: "<<l<<"\n";
        std::cout<<"  => .get_temp(): "<<CfgManip::get_temp(l)<<std::endl;
        std::cout<<"  => .set_temp(500)"<<std::endl;
        CfgManip::set_temp(l,NC::Temperature{500});
        std::cout<<"  => cfgstr: "<<l<<"\n";
        std::cout<<"  => .get_temp(): "<<CfgManip::get_temp(l)<<std::endl;
      } catch ( NC::Error::BadInput& e ) {
        if (expectbad) {
          std::cout<<"  => Got expected ERROR: NC::"<<e.getTypeName()<< ": "<<e.what() << std::endl;
          return;
        }
        throw;
      }
      if ( expectbad )
        NCRYSTAL_THROW(LogicError,"Did not fail as expected!!");
    };
    auto doTest = [&doTest_impl](const char*cstr) { return doTest_impl(cstr,false); };
    auto doTestBad = [&doTest_impl](const char*cstr) { return doTest_impl(cstr,true); };
    doTest("temp=300.0");
    doTestBad("temp=100.0fahrenheit");
    doTest("temp=100.0000000000000000000000000000000000");
    doTest("temp=0.123456789012345");
    doTest("temp=0.1234567890123456");
    doTest("temp=0.12345678901234567");
    doTest("temp=0.123456789012345678");
    doTest("temp=0.1234567890123456789");
    doTest("temp=0.1234567890123456789000000000000000000000000000000000000000000000000000010");
    doTest("temp=100.123456789123456789123456789");
    //    std::terminate();

    doTestBad("temp=100.123456789123456789123456789e300");
    doTestBad("temp=-2");
    doTest("temp=-1");
    doTestBad("temp=-0");
    doTest("temp=100.0F");
    doTest("temp=300.0K");
    doTest("temp=30.0C");
    doTest("temp=-20.0C");
    doTest("temp= 5 K ");
    doTestBad("temp=-5 K");
    doTestBad("temp=-infC");
    doTestBad("temp=infC");
    doTestBad("temp=inf");
    doTestBad("temp=+inf");
    doTestBad("temp=-inf");
    doTestBad("temp=nanK");
    doTestBad("temp=nan");
    doTest("temp=-400F");
    doTestBad("dcutoff=inf");
    doTestBad("dcutoff=+inf");

    doTest("dcutoff=0");
    doTest("dcutoff=-1");//backwards compat
    doTest("dcutoff=0.4");
    doTest("dcutoff=0.4Aa");
    doTest("dcutoff=0.04nm");
    doTestBad("dcutoff=-1.1");
    doTestBad("dcutoff=-2");
    doTestBad("dcutoff=-0.1");
    doTestBad("dcutoff=1e-4Aa");
    doTest("dcutoff=1e-3Aa");
    //Normalisation of numbers followed by a unit (shorter of original text
    //and generated text):
    doTest("dcutoff=0.0400nm");
    doTest("dcutoff= 1.50 Aa ");
    doTest("dcutoff=1.0e-1nm");
    doTest("temp=293.150K");
    doTest("temp=+20.0 C");
    doTest("mos=30.000deg");
    doTest("mos=1e-3deg");
    doTest("mos=0.50 arcmin");
    doTest("mos=1.0000000000000000001deg");
    doTest("dcutoff=1e5");
    doTestBad("dcutoff=1.1e5");


    doTest("; ; ;\t ;atomdb=Al is Al27@B is B10");
    doTest("; ; \r;\t ;atomdb=Al is Al27@B is B10");
    doTest("atomdb = ;");
    doTest("   atomdb =  ; dcutoff  = 0.5   ");
    doTest("vdoslux=2; ; ; ;  ;    ;;;;atomdb=  Al is Al27@B is B10  ; ; ;vdoslux\n\t=3");
    doTestBad("vdoslux=2; ; ; ;  ;    ;;;;atomdb=  Al is Al27@B is B10 and ala ; ; ;vdoslux\n\t=34");
    doTestBad("vdoslux=2; ; ; ;  ;    ;;;;atomdb=  Al is Al27@B is B10 ; ; ;vdoslux\n\t=34");

    doTest("vdoslux=2; ; ; ;  ;    ;;;;atomdb=  Al is Al27@B is B10 ; ; ;vdoslux=3;temp=20.85C;dcutoff=0.6Aa");
    doTest("vdoslux=2; ; ; ;  ;    ;;;;atomdb=  Al is Al27@B is B10 ; ; ;vdoslux=3;temp=20.85C;dcutoff=0.6Aa;elas=1");
    doTest("vdoslux=2; ; ; ;  ;    ;;;;atomdb=  Al is Al27@B is B10 ; ; ;vdoslux=3;temp=20.85C;dcutoff=0.6Aa;elas=1;bragg=0");

    //Test smallvector buffer overflow:
    static_assert(NC::Cfg::VarBufVector::nsmall==7,"");
    doTest("temp=20;dcutoff=1;dcutoffup=2;sccutoff=0.5;sans=0;incoh_elas=0;inelas=freegas;absnfactory=stdabs;atomdb=H:is:H2;scatfactory=stdscat");

    doTest("dir1=@crys_uvw:1,2,3@lab:0,0,1");
    doTest(" dir1 = @ crys_uvw : 1 , 2e0 , 3 @ lab : 0 , \t0 , 1\n ");
    doTestBad("dir1=@crys_xyz:1,2,3@lab:0,0,1");
    doTestBad("dir1=@crys_uvw:0,0,0@lab:0,0,1");
    doTestBad("dir1=@crys_uvw:1,2@lab:0,0,1");
    doTestBad("dir1=@crys:1e200,0,0@lab:0,0,1");
    doTestBad("dir1=@crys:0,0,1@lab:1e155,1e155,0");

    doTest("dir2=arbitrary");
    doTest(" dir2 = \tarbitrary ");
    doTestBad("dir2=arbitrarily");
    doTestBad("dir2=@crys:arbitrary@lab:0,0,1");
    doTestBad("dir1=arbitrary");

    doTest("lcaxis=0,0,1");
    doTest("lcaxis=@crys:0,0,1");
    doTest(" lcaxis = @ crys_hkl : 0 , 0 ,\t2e0 ");
    doTest("lcaxis=@crys_uvw:1,0,0.5");
    doTest("lcaxis=@lab:0,-1,1e-3");
    doTestBad("lcaxis=@crys_xyz:0,0,1");
    doTestBad("lcaxis=@lab:0,0,0");
    doTestBad("lcaxis=@lab:0,1");
    doTestBad("lcaxis=@lab0,1,0");
    doTestBad("lcaxis=lab:0,1,0");
    doTestBad("lcaxis=@:0,1,0");
    doTestBad("lcaxis=@crys_hkl:0,0,inf");
    doTestBad("lcaxis=@crys_hkl:0,0,nan");
    doTestBad("lcaxis=1e200,1e200,0");
    doTestBad("lcaxis=@lab:1e155,0,1e155");
    doTestBad("lcaxis=@crys_hkl:0,0,1@lab:0,1,0");
    doTest("lcaxis=0.50,0,1.0");
    doTest("lcaxis=@crys_uvw:0.50,-0,1e0");
    doTest("lcaxis=@dir1;dir1=@crys:0,0,1@lab:0,1,0");
    doTest(" lcaxis = @ dir1 ;dir1=@crys:0,0,1@lab:0,1,0");
    doTest("dir1avg=@crys_hkl:0,0,1@lab:0,1,0");
    doTest("fibre= @crys_uvw:1,1,0@lab:1,0,0 ");
    doTestBad("lcaxis=@dir2");
    doTestBad("lcaxis=@dir1:0,0,1");
    doTestBad("lcaxis=dir1");
    doTestBad("dir1avg=0,0,1");
    doTestBad("fibre=arbitrary");

    doTest("rot=0,1,0:30deg");
    doTest(" rot = 0 , 1.0 ,\t0 : 30 deg ");
    doTest("rot=0,0,2:1.5707963267948966");
    doTest("rot=1,1,1:0.5rad");
    doTest("rot=1,-1e-3,1:-1e-3 arcsec");
    doTest("rot=0,1,0:720arcmin");
    doTest("rot=0,1,0:123456789.123456789123deg");
    doTest("rot=0,1,0:0");
    doTest("rot=0,1,0:0.50");
    doTest("rot=0,1,0:0.500 rad");
    doTest("rot=0.50,1,0:+2.0e1deg");
    doTestBad("rot=0,1,0");
    doTestBad("rot=0,1:30deg");
    doTestBad("rot=0,1,0,0:30deg");
    doTestBad("rot=0,0,0:30deg");
    doTestBad("rot=0,1,0:30degs");
    doTestBad("rot=0,1,0:");
    doTestBad("rot=:1");
    doTestBad("rot=a,1,0:3");
    doTestBad("rot=1,0,0:1:2");
    doTestBad("rot=0,1,0:inf");
    doTestBad("rot=inf,1,0:1");
    doTestBad("rot=0,1e200,0:30deg");
    doTestBad("rot=1e155,1e155,0:30deg");
    doTestBad("rot=0,1,0:nan");
    doTestBad("rot=@lab:0,1,0:30deg");

    {
      NC::Cfg::CfgData cd{};
      std::cout<<"Test C++ interface:"<<std::endl;
      std::cout<<"  => cfgstr: "<<cd<<"\n";
      std::cout<<"  => .get_infofactory(): "<<CfgManip::get_infofactory(cd)<<std::endl;
      CfgManip::set_infofactory(cd,NC::StrView::make("lala")); std::cout<<"  => .set_infofactory(StrView(\"lala\")) : "<<CfgManip::get_infofactory(cd)<<std::endl;
      CfgManip::set_infofactory(cd,"bla_blu-bla"); std::cout<<"  => .set_infofactory(const char*(\"bla_blu-bla\")) : "<<CfgManip::get_infofactory(cd)<<std::endl;
      CfgManip::set_infofactory_stdstr(cd,std::string("wuhu")); std::cout<<"  => .set_infofactory(std::string(\"wuhu\")) : "<<CfgManip::get_infofactory(cd)<<std::endl;
      std::cout<<"  => .set_mos(1arcmin)"<<std::endl;
      CfgManip::set_mos( cd, NC::MosaicityFWHM{ 1*NC::kArcMin } );
      std::cout<<"  => .get_mos() : "<<CfgManip::get_mos(cd)<<std::endl;

      using LCForm = NC::LCAxisState::Form;
      std::cout<<"  => .set_lcaxis( {0,0,1} )"<<std::endl;
      CfgManip::set_lcaxis( cd, { LCForm::Cartesian, {0,0,1} } );
      std::cout<<"  => .get_lcaxis() : "<<CfgManip::get_lcaxis(cd).value()<<std::endl;
      CfgManip::set_lcaxis( cd, { LCForm::Cartesian, {1,1,1} } );
      std::cout<<"  => .get_lcaxis() : "<<CfgManip::get_lcaxis(cd).value()<<std::endl;
      std::cout<<"  => cfgstr: "<<cd<<"\n";
      (void)CfgManip::applyStrCfg(cd,"mos=1arcmin");
      std::cout<<"  => .applyStrCfg(\"mos=1arcmin\") : "<<std::endl;
      std::cout<<"  => cfgstr: "<<cd<<"\n";
      (void)CfgManip::applyStrCfg(cd,"mos= 1 \tarcmin ");
      std::cout<<"  => .applyStrCfg(\"mos= 1 \\tarcmin \") : "<<std::endl;
      std::cout<<"  => cfgstr: "<<cd<<"\n";

      auto applyAndPrint = [&cd](const char* cfgstr)
      {
        (void)CfgManip::applyStrCfg(cd,cfgstr);
        std::cout<<"  => .applyStrCfg(\""<<cfgstr<<"\") : "<<std::endl;
        std::cout<<"  => cfgstr: "<<cd<<"\n";
      };
      applyAndPrint("dir1=@crys:0,0,1@lab:1,2,3");
      applyAndPrint("dir1=@crys_hkl:4,5,6@lab:1,2,3");
      applyAndPrint(" dir1 = @ crys_hkl : 4 , 5 , 6 @ lab : 1 , \t2 , 3\n ");
      applyAndPrint("dir1=@crys_hkl:4,5,6@lab:1,2,3");
      applyAndPrint("dir1=@crys_uvw:4,5,6@lab:1,2,3");
      applyAndPrint("rot= 0,0,2 : 0.25 deg");
      applyAndPrint("dir2=@crys_hkl:4,5,6@lab:1,2,3");
      applyAndPrint("dir2=arbitrary");
      {
        const auto r = CfgManip::get_rot(cd).value();
        std::cout<<"  => .get_rot() : "<<r.first<<" "<<NC::fmt(r.second)
                 <<std::endl;
      }
      std::cout<<"  => JSON: ";
      CfgManip::streamJSON(cd,std::cout);
      std::cout<<"\n";
      applyAndPrint("lcaxis=@crys_hkl:0,0,2");
      std::cout<<"  => JSON: ";
      CfgManip::streamJSON(cd,std::cout);
      std::cout<<"\n";
    }
  }


  //ValVector is not used by any parameter in NCrystal itself, but is
  //available for plugins (e.g. for directions). Test it via test-only
  //variable definitions, normalised (the default) or not:
  struct testvar_vec final : public NC::Cfg::ValVector<testvar_vec> {
    static constexpr auto name = "testvec";
  };
  struct testvar_vecraw final : public NC::Cfg::ValVector<testvar_vecraw> {
    static constexpr auto name = "testvecraw";
    static constexpr bool auto_normalise = false;
  };

  template<class TVar>
  void testValVectorStr( const char * str )
  {
    using namespace NC::Cfg;
    std::cout<<"ValVector \""<<TVar::name<<"\" from \""<<str<<"\": ";
    try {
      const auto buf = TVar::from_str( VarId::lcaxis, str );
      const auto v = TVar::get_val( buf );
      std::ostringstream ss;
      TVar::stream_val( ss, buf );
      std::cout<<ss.str()<<" JSON=";
      TVar::asJSONObject( std::cout, buf );
      std::cout<<std::endl;
      //Round trip, and same as set_val:
      const std::string text = ss.str();
      const auto buf2 = TVar::from_str( VarId::lcaxis, text );
      REQUIRE( TVar::cmp( buf, buf2 ) == 0 );
      REQUIRE( TVar::cmp( buf, TVar::set_val( VarId::lcaxis, v ) ) == 0 );
      if ( TVar::auto_normalise )
        REQUIRE( v.isUnitVector() );
    } catch ( NC::Error::BadInput& e ) {
      std::cout<<"BadInput: "<<e.what()<<std::endl;
    }
  }

  void testValVector()
  {
    using namespace NC::Cfg;
    for ( auto s : { "1,2,2", " 1 , 2 ,\t2 ", "0,0,-5", "1e-200,0,0", "1,2",
                     "1,2,3,4", "a,b,c", "1,,2", "0,0,0", "inf,0,0", "nan,0,0",
                     "1e200,1e200,0" } ) {
      testValVectorStr<testvar_vec>( s );
      testValVectorStr<testvar_vecraw>( s );
    }
    //Comparisons are lexicographic in the components:
    const auto a = testvar_vecraw::set_val( VarId::lcaxis, { 1.0, 2.0, 3.0 } );
    const auto b = testvar_vecraw::set_val( VarId::lcaxis, { 1.0, 2.5, 0.0 } );
    REQUIRE( testvar_vecraw::cmp( a, b ) < 0 && testvar_vecraw::cmp( b, a ) > 0 );
    REQUIRE( testvar_vecraw::cmp( a, a ) == 0 );
  }

  void test_cfgdata() {

    static constexpr auto name = "vdoslux";
    static constexpr auto var_id = NC::Cfg::constexpr_varIdFromName(name);
    static constexpr auto var_name = NC::Cfg::varName( var_id );
    static_assert(var_id==NC::Cfg::VarId::vdoslux,"");
    static_assert(NC::constexpr_cstrequal(var_name,name),"");
    std::cout<<"\""<<name<<"\" -> "<<NC::enumAsInt(var_id)<<" -> \""<<var_name<<"\""<<std::endl;

    nc_assert_always( !( NC::StrView("temp") < NC::StrView("dcutoff") ) ) ;
    nc_assert_always( NC::StrView("coh") < NC::StrView("dcutoff") ) ;
    nc_assert_always( NC::StrView("coh_elas") < NC::StrView("dcutoff") ) ;
    nc_assert_always( NC::StrView("dcutoff") < NC::StrView("dcutoffup") ) ;
    nc_assert_always( NC::StrView("dcutoffup") < NC::StrView("incoh_elas") ) ;

    // make_varinfo<vardef_coh_elas>(),
    // make_varinfo<vardef_dcutoff>(),
    // make_varinfo<vardef_dcutoffup>(),
    // make_varinfo<vardef_incoh_elas>(),
    // make_varinfo<vardef_infofactory>(),
    // make_varinfo<vardef_temp>(),
    // make_varinfo<vardef_vdoslux>()

    for ( const auto& varinfo : NC::Cfg::varlist ) {
      //std::cout<<"\""<<varinfo.name()<<"\""<<std::endl;
      auto opt_varid = NC::Cfg::varIdFromName( varinfo.nameSV() );
      nc_assert_always(opt_varid.has_value());
      auto varid = opt_varid.value();
      //inline constexpr unsigned constexpr_varName2Id( const char * name, unsigned detail_idx = 0 )
      nc_assert_always( &NC::Cfg::varInfo(varid) == &varinfo );
      std::cout<<"\""<<varinfo.name()<<"\" -> "<<NC::enumAsInt(varid)<<" -> \""<<NC::Cfg::varInfo(varid).name()<<"\""<<std::endl;
    }

    {
      auto opt_varid = NC::Cfg::varIdFromName( "temp" );
      nc_assert_always(opt_varid.has_value());
      const auto & varinfo = NC::Cfg::varInfo(opt_varid.value());
      auto varid = opt_varid.value();
      auto strval = "120F";
      auto buf_temp = varinfo.from_str( varid, strval );
      std::cout<<" setting \""<<varinfo.name()<<"\" to \""<<strval
               <<"\" -> "<<NC::Cfg::vardef_temp::get_val(buf_temp)<<" -> \"";
      varinfo.stream( std::cout, buf_temp );
      std::cout<<"\"\n";
    }
  }

  void test_vdoslux_values() {
    //Hidden vdoslux values: 1000..1005 (legacy) are kept, 2000..2006 are
    //normalised to 0..6 (so they give identical cfg strings and results):
    using NC::Cfg::CfgManip;
    for ( auto v : { "3", "2003", "2000", "2006", "1003", "1000", "1005",
                     "7", "1006", "2007", "999", "-1" } ) {
      std::cout<<"vdoslux="<<v<<": ";
      try {
        auto l = NC::Cfg::CfgData();
        const std::string cfgstr = std::string("vdoslux=") + v;
        auto toplvl = CfgManip::applyStrCfg( l, cfgstr );
        nc_assert_always( toplvl.empty() );
        std::cout<<"get_vdoslux()="<<CfgManip::get_vdoslux(l)
                 <<" cfgstr=\""<<l<<"\""<<std::endl;
      } catch ( NC::Error::BadInput& e ) {
        std::cout<<"BadInput: "<<e.what()<<std::endl;
      }
    }
  }
}

int main() {
  testValVector();
  test_cfgdata();
  //  std::cout<<sizeof(B)<<std::endl;
  //  std::cout<<sizeof(NC::SmallVector<B,8,NC::SVMode::FASTACCESS>)<<std::endl;
  test_varids();
  test_parsing();
  test_parse_and_fill();
  test_vdoslux_values();
  NC::Cfg::dumpCfgVarList( std::cout, NC::Cfg::CfgVarListMode::TXT_FULL );
  std::cout<<std::endl;
  NC::Cfg::dumpCfgVarList( std::cout, NC::Cfg::CfgVarListMode::TXT_SHORT );
  std::cout<<std::endl;
  NC::Cfg::dumpCfgVarList( std::cout, NC::Cfg::CfgVarListMode::JSON );
  std::cout<<std::endl;
}
