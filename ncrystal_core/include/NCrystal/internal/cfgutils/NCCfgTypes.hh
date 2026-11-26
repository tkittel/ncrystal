#ifndef NCrystal_CfgTypes_hh
#define NCrystal_CfgTypes_hh

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

#include "NCrystal/core/NCImmutBuf.hh"
#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/core/NCVariant.hh"

namespace NCRYSTAL_NAMESPACE {

  //////////////////////////////////////////////////////////////
  // Infrastructure needed for the variables in NCCfgVars.hh  //
  //////////////////////////////////////////////////////////////

  namespace Cfg {

    using VarId = detail::VarId;
    using VarBuf = detail::VarBuf;
    using VarBufVector = detail::VarBufVector;

    static constexpr StrView forbidden_chars_non_multiphase = StrView::make("\"'|(){}[]*&><#$");
    static constexpr StrView forbidden_chars_value_strreps = StrView::make("\"'|(){}[]*&><;=#$");//+ "=;"
    static constexpr StrView forbidden_chars_multiphase = StrView::make("\"'|(){}[]#$");//*&<> are allowed only in certain locations of multiphase strings

    enum class VarGroupId { Info, ScatterBase, ScatterExtra, Absorption };

    //A few utilities for doubles:
    inline double sanitiseDblValue(double v,const char * varname) {
      if ( std::isnan(v) )
        NCRYSTAL_THROW2(BadInput,"NAN (not-a-number) value provided"
                        " for parameter \""<<varname<<"\"");
      return ( v == 0.0 ? 0.0 : v );
    }
    using ValDbl_ShortStrOrigRep = ShortStr<VarBuf::buffer_local_size
                                            - sizeof(double)>;
    struct units_notavail {
      static constexpr auto actual_unit = nullptr;
      static void listAvailableUnits(std::ostream&){}
      static Optional<std::pair<double,ValDbl_ShortStrOrigRep>> parse(StrView){ return NullOpt; }
    };

    //Normalised text of a number in a cfg value: the shorter of its original
    //text and the generated one (dbl2shortstr of its value), preferring the
    //latter if equally long. The result refers to one of the two inputs:
    inline StrView shorterNumStr( StrView orig, const ShortStrDbl& generated )
    {
      return ( generated.size() <= orig.size() ? generated.to_view() : orig );
    }

    //Helper util for parsing doubles with units:
    struct ParsedDblValue {
      double value;
      StrView unit;
      //Normalised str rep (number normalised with shorterNumStr, followed by
      //the unit), kept in a dedicated buffer since internal whitespace must
      //be edited away, e.g. "1.0 Aa" -> "1Aa".
      ValDbl_ShortStrOrigRep origStrRep;//empty if did not fit.
    };
    Optional<ParsedDblValue> unitSplit(StrView sv);

    void standardInputStrSanityCheck( const char * parname, StrView strrep );

    //Throws BadInput about a syntax error in the value of a parameter:
    [[noreturn]] void throwSyntaxError( const char * parname, StrView value );

    //Numbers in cfg values which must be finite (BadInput exception
    //otherwise). Returned with sanitiseDblValue applied:
    double checkFiniteDbl( double, const char * parname );

    //Vectors in cfg values. The components and the squared length must be
    //finite and the vector non-null (BadInput exception otherwise).
    //Components are returned with
    //sanitiseDblValue applied:
    Vector checkVector( const Vector&, const char * parname );
    //Parses "x,y,z" (spaces allowed), returning NullOpt in case of syntax
    //errors. The checked vector is returned along with its normalised text,
    //where each component is normalised with shorterNumStr:
    Optional<std::pair<Vector,std::string>> parseVector( StrView,
                                                         const char * parname );
    //Normalised text of a vector given only by its value:
    std::string vectorStr( const Vector& );

    //Bases of crystal directions, as named in cfg-strings ("crys", "crys_hkl"
    //and "crys_uvw") and in JSON ("cartesian", "hkl" and "uvw"):
    const char * crysBasisName( OrientDir::CrysDir::Basis );
    const char * crysBasisJSONName( OrientDir::CrysDir::Basis );
    Optional<OrientDir::CrysDir::Basis> crysBasisFromName( StrView );

    //Storage of structured values as N doubles, a tag (char), and the
    //normalised text (null-terminated). Comparisons are by tag (in
    //descending order), then by the doubles, then by the text:
    template<std::size_t N>
    struct DblsTagText {
      std::array<double,N> vals;
      unsigned char tag;
      StrView text;//refers to the VarBuf
    };
    template<std::size_t N>
    VarBuf encodeDblsTagText( VarId, const std::array<double,N>&,
                              unsigned char tag, StrView text );
    template<std::size_t N>
    DblsTagText<N> decodeDblsTagText( const VarBuf& );
    template<std::size_t N>
    int cmpDblsTagText( const VarBuf&, const VarBuf& );

    template <class Derived, class TValueType>
    class ValBase {
    public:
      using units = units_notavail;
      using value_type = TValueType;

      constexpr static bool has_default_value()
      {
        static_assert ( std::is_same<decltype(Derived::default_value()),value_type>::value
                        || std::is_same<decltype(Derived::default_value()),NullOptType>::value, "" );
        return std::is_same<decltype(Derived::default_value()),value_type>::value;
      }

      static VarBuf default_value_in_varbuf()
      {
        nc_assert(has_default_value());
        struct detail2 {
          static VarBuf default_varbuf( NullOptType ) noexcept { return VarBuf(NullOpt); }//not actually called
          static VarBuf default_varbuf( const value_type& val ) noexcept { return Derived::set_val( static_cast<VarId>(0), val ); }//dummy varid
        };
        return detail2::default_varbuf( Derived::default_value() );
      }

      static void stream_default_value( std::ostream& os )
      {
        //This does not need to be super effiecient, we only need it to generate
        //documentation. It is more important to be consistent, and for that
        //reason we first place the default value in a VarBuf and then stream it.
        if (!has_default_value())
          return;
        Derived::stream_val(os,default_value_in_varbuf());
      }
      struct detail_dvon {
        //Next fct won't actually get called
        static constexpr value_type dvon( NullOptType )
        {
          nc_assert_always(false);
          return *(value_type*)(nullptr);// cppcheck-suppress nullPointer
        }
        static constexpr value_type dvon( const value_type& v ) { return v; }
      };
      constexpr static value_type default_value_or_nonsense()
      {
        return detail_dvon::dvon( Derived::default_value() );
      }

      static void stream_default_value_json( std::ostream& os )
      {
        //Like stream_default_value, but with proper type for json. Outputs json
        //"null" if there is no default value.
        if ( has_default_value() )
          Derived::asJSONObject( os, default_value_in_varbuf() );
        else
          streamJSON(os,json_null_t{});
      }

    };

    template <class Derived>
    class ValDbl : public ValBase<Derived,double> {
      // As ever, floating point numbers are a bit weird and annoying in
      // several ways, even for IEEE 754 doubles: Contain values that prevent
      // proper sorting (NaN), bit patterns that are not unique for a given
      // value (-0 and +0), string representations that are not unique for a
      // given value ("0.1" and "0.99999999999999999999999"). To prevent bugs
      // and weird user experiences, we adopt the following strategy. First of
      // all, we always store the actual floating point (double) value in the
      // first sizeof(double) chars of the buffer. The following chars
      // optionally provides a null-terminated c-string with the string
      // representation of the value. However, this might be absent
      // (i.e. start with a null char) in case the string representation does
      // not fit in the local buffer, in which case a string representation
      // will be created on-the-fly from the double value when needed.
      //
      // A few take-aways:
      //
      // 1) To avoid weird sorting and other issues: Forbid NAN, and map -0 to 0.
      // 2) Always store the double itself in first 8 bytes
      // 3) Store str-rep in last bytes if and only if it fits there (and when
      //    setting from string, prefer the original version only if it is shorter).
      // 4) When sorting, first sort the double value, then presence and
      //    contents of str-reps.
      // 5) When streaming and strrep is missing, generate it with the same
      //    method call as during setting.

      static constexpr auto buf_maxstrlen = VarBuf::buffer_local_size - sizeof(double) - 1;//-1 for null char
      static_assert(VarBuf::buffer_local_size>=detail::varbuf_calc::ValDbl_buf_minsize,"");
      static_assert(VarBuf::buffer_local_size>=sizeof(double),"");
      static_assert(alignof(VarBuf)%alignof(double)==0,"");
      static VarBuf detail_set_val_with_str( VarId varid, double value, StrView valstr )
      {
        //Assumes value already sanitised, validated, and the shortest str rep was found.
        char tmp[VarBuf::buffer_local_size];
        std::memcpy(&tmp[0],(void*)&value,sizeof(double));
        std::size_t nused = sizeof(double) + 1;//at least the double and a null char
        if ( valstr.size() <= buf_maxstrlen ) {
          std::memcpy( &tmp[sizeof(double)], valstr.data_begin(), valstr.size() );
          tmp[sizeof(double)+valstr.size()] = '\0';
          nused += valstr.size();//also the characters
        } else {
          //don't store strrep, indicate with a single null char
          tmp[sizeof(double)] = '\0';
        }
        return VarBuf( &tmp[0], nused, varid );
      }
    public:
      using Base = ValDbl;
      static_assert( ValDbl_ShortStrOrigRep::bufsize == buf_maxstrlen + 1,"");

      static constexpr auto value_type_descr = "floating point number";

      static double value_validate( double value ) { return value; }

      static int cmp( const VarBuf& a, const VarBuf& b )
      {
        //First just sort on values:
        double valA = get_val( a );
        double valB = get_val( b );
        //NB: sanitiseDblValue guarantees no NaNs:
        nc_assert( !std::isnan(valA) );
        nc_assert( !std::isnan(valB) );
        if ( valA != valB )
          return ( valA < valB ? -1 : 1 );
        //valA = valB. To consider the full state (and thus avoid weird issues),
        //we now also consider the str-rep's if any:
        const char * cstr_A = ( a.dataAssertLocal() + sizeof(double) );
        const char * cstr_B = ( b.dataAssertLocal() + sizeof(double) );
        if ( bool(cstr_A) != bool(cstr_B) )
          return bool(cstr_A) ? -1 : 1;//those with strreps before those without
        if (!cstr_A)
          return 0;//both strreps are absent
        //Both strreps are there:
        return std::strcmp(cstr_A,cstr_B);
      }

      static double get_val( const VarBuf& buf )
      {
        double v;
        std::memcpy( &v, buf.dataAssertLocal(), sizeof(v) );
        return v;
      }

      static VarBuf set_val( VarId varid, double value ) {
        value = Derived::value_validate( sanitiseDblValue(value,Derived::name) );
        return detail_set_val_with_str( varid, value, dbl2shortstr(value).to_view() );
      }
      static void stream_val( std::ostream& os, const VarBuf& buf )
      {
        const char * buf_data = buf.dataAssertLocal();
        const char * cstr = ( buf_data + sizeof(double) );
        if ( *cstr ) {
          //str rep (including null termination) was short enough to keep cached:
          os << cstr;
        } else {
          //str rep was not short enough to cache, generate on demand (passing
          //through same method as the cached would have been, for
          //consistency):
          double v;
          std::memcpy( &v, buf_data, sizeof(v) );
          os << dbl2shortstr( v );
        }
      }

      static void asJSONObject( std::ostream& os, const VarBuf& b )
      {
        streamJSON(os,get_val(b));
      }

      static VarBuf from_str( VarId varid, StrView orig_str ) {
        standardInputStrSanityCheck(Derived::name,orig_str);
        auto opt_val = Derived::units::parse(orig_str);//this allows derived class to handle units and ranges
        if (!opt_val.has_value())
          NCRYSTAL_THROW2(BadInput,"Syntax error - invalid value \""
                          <<orig_str<<"\" provided for parameter \""<<Derived::name<<"\"");
        double value = opt_val.value().first;
        const auto& orig_str_cleaned = opt_val.value().second;
        value = Derived::value_validate( sanitiseDblValue( value, Derived::name ) );
        auto valstr = dbl2shortstr( value );
        return detail_set_val_with_str( varid,
                                        value,
                                        ( ( orig_str_cleaned.empty() || valstr.size() <= orig_str_cleaned.size() )
                                          ? valstr.to_view()
                                          : orig_str_cleaned.to_view() ) );
      }
    };

    template <class Derived>
    class ValVector : public ValBase<Derived,Vector> {
    public:
      //Simply stores a Vector directly in the buffer.
      static_assert(VarBuf::buffer_local_size>=sizeof(Vector),"");
      static_assert(alignof(VarBuf)%alignof(Vector)==0,"");
      static_assert(std::is_trivially_copyable<Vector>::value,"");
      static_assert(std::is_trivially_destructible<Vector>::value,"");

      static constexpr bool auto_normalise = true;//will normalise all vectors to unit length unless overridden to false.
      using Base = ValVector;
      static constexpr auto value_type_descr = "vector (3D)";
      static Vector get_val( const VarBuf& buf )
      {
        Vector v;
        std::memcpy( &v, buf.dataAssertLocal(), sizeof(v) );
        return v;
      }
      static int cmp( const VarBuf& a, const VarBuf& b )
      {
        return get_val( a ).lexCmp( get_val( b ) );
      }
      static VarBuf set_val( VarId varid, const Vector& value ) {
        static_assert(std::is_trivially_destructible<Vector>::value,"");
        static_assert(std::is_trivially_copyable<Vector>::value,"");
        alignas(Vector) char tmp[sizeof(Vector)];
        Vector & v = *(new (&tmp[0]) Vector(no_init));
        v[0] = sanitiseDblValue( value[0], Derived::name );
        v[1] = sanitiseDblValue( value[1], Derived::name );
        v[2] = sanitiseDblValue( value[2], Derived::name );
        if ( Derived::auto_normalise ) {
          double m2 = v.mag2();
          if ( !m2 )
            NCRYSTAL_THROW2(BadInput,"Null vector provided for parameter \""<<Derived::name<<"\"");
          v.normalise();//slower but safer than "v /= std::sqrt(m2)"
        }
        //Extra sanity checks:
        for ( auto& e : v )
          e = sanitiseDblValue( e, Derived::name );
        Derived::extraChecks( v );
        return VarBuf( &tmp[0], sizeof(tmp), varid );
      }

      static void stream_val( std::ostream& os, const VarBuf& buf )
      {
        const Vector v = get_val(buf);
        os << fmt( v[0] ) << ',' << fmt( v[1] ) << ',' << fmt( v[2] );
      }

      static void asJSONObject( std::ostream& os, const VarBuf& b )
      {
        streamJSON(os,get_val(b));
      }

      static VarBuf from_str( VarId varid, StrView sv_in ) {
        standardInputStrSanityCheck(Derived::name,sv_in);
        Vector v(no_init);
        auto sv = sv_in;
        for ( auto idx : ncrange(3) ) {
          auto i = sv.find(',');
          Optional<double> opt_dblval;
          if ( ( i == StrView::npos ) == ( idx == 2 ) ) {
            opt_dblval = sv.substr(0,i).trimmed().toDbl();
          }
          if ( !opt_dblval.has_value() )
            NCRYSTAL_THROW2(BadInput,"Syntax error - invalid value \""
                            <<sv_in<<"\" provided for parameter \""<<Derived::name<<"\"");
          v[idx] = opt_dblval.value();
          sv = sv.substr( i < StrView::npos ? i+1 : StrView::npos );//+1 to skip past comma
        }
        return set_val( varid, v );
      }

      static void extraChecks( const Vector& ) {}
    };

    template <class Derived>
    class ValOrientDir : public ValBase<Derived,OrientDir> {
    public:
      //Stored as the crystal and lab directions (6 doubles), the basis of the
      //crystal direction (tag) and the normalised text.
      using Base = ValOrientDir;
      static constexpr auto value_type_descr = "crystal axis orientation";
      using Basis = OrientDir::CrysDir::Basis;

      static OrientDir get_val( const VarBuf& buf )
      {
        const auto d = decodeDblsTagText<6>( buf );
        OrientDir od;
        od.crystal.basis = static_cast<Basis>( d.tag );
        od.crystal.v = { d.vals[0], d.vals[1], d.vals[2] };
        od.lab = { d.vals[3], d.vals[4], d.vals[5] };
        return od;
      }

      static int cmp( const VarBuf& a, const VarBuf& b )
      {
        return cmpDblsTagText<6>( a, b );
      }

      static VarBuf set_val( VarId varid, const OrientDir& od )
      {
        const auto& v = od.crystal.v;
        const auto c = checkVector( Vector( v[0], v[1], v[2] ), Derived::name );
        const auto l = checkVector( od.lab.as<Vector>(), Derived::name );
        return detail_encode( varid, od.crystal.basis, c, l,
                              vectorStr( c ), vectorStr( l ) );
      }

      static void stream_val( std::ostream& os, const VarBuf& buf )
      {
        os << decodeDblsTagText<6>( buf ).text;
      }

      static void asJSONObject( std::ostream& os, const VarBuf& buf )
      {
        const auto od = get_val( buf );
        os << "{\"basis\":";
        streamJSON( os, crysBasisJSONName( od.crystal.basis ) );
        os << ",\"crystal\":";
        streamJSON( os, od.crystal.v );
        os << ",\"lab\":";
        streamJSON( os, od.lab );
        os << '}';
      }

      static VarBuf from_str( VarId varid, StrView sv_in )
      {
        //"@crys:c1,c2,c3@lab:l1,l2,l3", with "@crys_hkl:" or "@crys_uvw:" as
        //alternatives to "@crys:", and spaces allowed around all delimiters
        //or at string ends.
        standardInputStrSanityCheck( Derived::name, sv_in );
        auto throw_error = [&sv_in]() { throwSyntaxError( Derived::name, sv_in ); };
        const StrView sv = sv_in.trimmed();
        const auto i0 = ( sv.startswith('@')
                          ? sv.substr(1).find('@') : StrView::npos );
        if ( i0 == StrView::npos )
          throw_error();
        const auto i = i0 + 1;
        auto splitName = [&throw_error]( StrView part,
                                         StrView& name ) -> StrView
        {
          const auto j = part.find(':');
          if ( j == StrView::npos )
            throw_error();
          name = part.substr(0,j).trimmed();
          return part.substr(j+1);
        };
        StrView cname, lname;
        const auto c = parseVector( splitName( sv.substr(1,i-1), cname ),
                                    Derived::name );
        const auto l = parseVector( splitName( sv.substr(i+1), lname ),
                                    Derived::name );
        const auto basis = crysBasisFromName( cname );
        if ( !c.has_value() || !l.has_value() || !basis.has_value()
             || lname != StrView::make("lab") )
          throw_error();
        return detail_encode( varid, basis.value(), c.value().first,
                              l.value().first, c.value().second,
                              l.value().second );
      }

    private:
      static VarBuf detail_encode( VarId varid, Basis basis, const Vector& c,
                                   const Vector& l, const std::string& ctext,
                                   const std::string& ltext )
      {
        std::string text = std::string("@") + crysBasisName( basis ) + ':'
          + ctext + "@lab:" + ltext;
        return encodeDblsTagText<6>( varid,
                                     { c[0], c[1], c[2], l[0], l[1], l[2] },
                                     static_cast<unsigned char>( basis ),
                                     text );
      }
    };

    template <class Derived>
    class ValStr : public ValBase<Derived,StrView> {
    public:
      using Base = ValStr;
      static constexpr auto value_type_descr = "string";//NB: Don't change this without follow up
                                                        //changes elsewhere (used to recognise type
                                                        //needing quotation marks in printouts)

      //Encode value in buffer - including null char.
      static StrView get_val( const VarBuf& buf )
      {
        return StrView( buf.data() );
      }

      static int cmp( const VarBuf& a, const VarBuf& b )
      {
        StrView valA = get_val( a );
        StrView valB = get_val( b );
        return valA == valB ? 0 : ( valA < valB ? -1 : 1 );
      }

      static VarBuf set_val( VarId varid, StrView value ) {
        return actual_set_val( varid, value );
      }

      static VarBuf set_val( VarId varid, const std::string& value ) {
        return actual_set_val( varid, StrView(value) );
      }

      static VarBuf set_val( VarId varid, const char * value ) {
        return actual_set_val( varid, StrView(value) );
      }

      static VarBuf from_str( VarId varid, StrView sv_in ) {
        return actual_set_val( varid, sv_in );
      }

      static void stream_val( std::ostream& os, const VarBuf& buf )
      {
        os << buf.data();
      }

      static void asJSONObject( std::ostream& os, const VarBuf& b )
      {
        streamJSON(os,get_val(b));
      }

    protected:
      static Variant<StrView,std::string> str2val( StrView sv ) { return sv; }
    private:
      static VarBuf actual_set_val( VarId varid, StrView sv_in ) {
        standardInputStrSanityCheck(Derived::name,sv_in);
        //NB: Not sure why "auto opt_value = ..." in the next line does not work!!
        Variant<StrView,std::string> opt_value = Derived::str2val(sv_in);
        static_assert(std::is_same<Variant<StrView,std::string>,decltype(opt_value)>::value,"");
        if (opt_value.empty())
          NCRYSTAL_THROW2(BadInput,"Syntax error - invalid value \""
                          <<sv_in<<"\" provided for parameter \""<<Derived::name<<"\"");

        if ( opt_value.has_value<std::string>() ) {
          const auto& value = opt_value.get<std::string>();
          return VarBuf( &value[0], value.size()+1, varid );//+1 for null terminator!
        } else {
          nc_assert(opt_value.has_value<StrView>());
          const auto& value = opt_value.get<StrView>();
          //Must copy in order to null terminate (NB: we could possibly avoid
          //this by a bit of care and an additional constructor on VarBuf):
          SmallVector<char,256> tmp;
          tmp.setByCopy( value.begin(),value.end() );
          tmp.push_back('\0');
          return VarBuf( &tmp[0], tmp.size(), varid );
        }
      }
    };

    template <class Derived>
    class ValBool : public ValBase<Derived,bool> {
      //ValBuf is just holding a single char with the value 0 or 1.
    public:
      using Base = ValBool;
      static constexpr auto value_type_descr = "boolean";

      static int cmp( const VarBuf& a, const VarBuf& b )
      {
        bool valA = get_val( a );
        bool valB = get_val( b );
        return valA == valB ? 0 : ( valA ? -1 : 1 );
      }

      static bool get_val( const VarBuf& buf )
      {
        static_assert(VarBuf::buffer_local_size>=1,"");
        return 0 != buf.dataAssertLocal()[0];
      }
      static VarBuf set_val( VarId varid, bool value ) {
        char c(value?1:0);
        return VarBuf( (const char*)&c, 1, varid );
      }
      static void stream_val( std::ostream& os, const VarBuf& buf )
      {
        os << (get_val(buf)?"1":"0");
      }

      static void asJSONObject( std::ostream& os, const VarBuf& b )
      {
        streamJSON(os,get_val(b));
      }

      static VarBuf from_str( VarId varid, StrView sv ) {
        standardInputStrSanityCheck(Derived::name,sv);
        if ( sv == StrView::make("true") || sv == StrView::make("1") ) {
          return set_val( varid, true );
        } else if ( sv == StrView::make("false") || sv == StrView::make("0") ) {
          return set_val( varid, false );
        } else {
          NCRYSTAL_THROW2(BadInput,"Could not convert \""<<sv
                          <<"\" to boolean value (should be \"true\", \"1\", \"false\" or \"0\")");
        }
      }
    };

    template <class Derived>
    class ValInt : public ValBase<Derived,int64_t> {
      //ValBuf is just holding the bytes of a int64_t value directly (we don't
      //have the same value<->string issues as for a floating point number).
    public:
      using Base = ValInt;
      static constexpr auto value_type_descr = "integer";

      static int cmp( const VarBuf& a, const VarBuf& b )
      {
        int64_t valA = get_val( a );
        int64_t valB = get_val( b );
        return valA == valB ? 0 : ( valA < valB ? -1 : 1 );
      }

      static int64_t get_val( const VarBuf& buf )
      {
        static_assert(VarBuf::buffer_local_size>=sizeof(int64_t),"");
        static_assert(alignof(VarBuf)%alignof(int64_t)==0,"");
        int64_t v;
        std::memcpy( &v, buf.dataAssertLocal(), sizeof(v) );
        return v;
      }
      static VarBuf set_val( VarId varid, int64_t value ) {
        value = Derived::value_validate(value);
        return VarBuf( (const char*)&value, sizeof(int64_t), varid );
      }

      static int64_t value_validate( int64_t value ) { return value; }

      static void stream_val( std::ostream& os, const VarBuf& buf )
      {
        os << get_val(buf);
      }

      static void asJSONObject( std::ostream& os, const VarBuf& b )
      {
        streamJSON(os,get_val(b));
      }

      static Optional<int64_t> str2val( StrView sv ) { return sv.toInt(); }

      static VarBuf from_str( VarId varid, StrView str ) {
        standardInputStrSanityCheck(Derived::name,str);
        Optional<int64_t> opt_val = Derived::str2val(str);
        if (!opt_val.has_value())
          NCRYSTAL_THROW2(BadInput,"Syntax error - invalid value \""
                          <<str<<"\" provided for parameter \""<<Derived::name<<"\"");
        return set_val( varid, opt_val.value() );
      }
    };

    using VarFromStrFct = VarBuf(*)(VarId, StrView);
    using StreamVarFct  = void(*)(std::ostream&, const VarBuf&);
    using BufCmpFct = int(*)(const VarBuf&,const VarBuf&);
    using BufToJSONFct = void(*)(std::ostream&, const VarBuf&);

    class VarInfo final {
    public:

      using StreamFct = void (*)( std::ostream& );

      constexpr VarInfo( VarFromStrFct fs,
                         StreamVarFct sv,
                         BufCmpFct bcf,
                         VarGroupId grId,
                         const char * thename,
                         const char * thedescr,
                         StreamFct streamDefaultFct,
                         StreamFct streamDefaultFctJSON,
                         StreamFct streamUnitsFct,
                         const char * actualUnitName,
                         const char * valueTypeDescr,
                         BufToJSONFct jsonfct )
        : m_fromStrFct(fs),
          m_streamFct(sv),
          m_bcf(bcf),
          m_groupId(grId),
          m_namecstr(thename),
          m_name(StrView::constexpr_t(),thename),
          m_description(thedescr),
          m_sdf(streamDefaultFct),
          m_sdfjson(streamDefaultFctJSON),
          m_suf(streamUnitsFct),
          m_aun(actualUnitName),
          m_vtd(valueTypeDescr),
          m_jsonfct(jsonfct)
      {
      }
      VarBuf from_str( VarId varid, StrView sv ) const { return m_fromStrFct(varid,sv); }
      void stream( std::ostream& os, const VarBuf& buf) const { m_streamFct(os,buf); }
      void streamAsJSON( std::ostream& os, const VarBuf& buf ) const { nc_assert(m_jsonfct); m_jsonfct(os,buf); }
      int bufCmp( const VarBuf& a, const VarBuf& b ) const { return m_bcf(a,b); }
      constexpr VarGroupId groupId() const noexcept { return m_groupId; }
      constexpr const char * name() const noexcept { return m_namecstr; }
      constexpr const StrView& nameSV() const noexcept { return m_name; }
      constexpr const char * description() const noexcept { return m_description; }
      constexpr bool hasDefaultValue() const noexcept { return m_sdf != nullptr; }
      void streamDefaultValue( std::ostream& os ) const { nc_assert(m_sdf); m_sdf(os); }
      void streamDefaultValueJSON( std::ostream& os ) const { nc_assert(m_sdfjson); m_sdfjson(os); }
      constexpr bool hasUnits() const noexcept { return m_suf != nullptr; }
      void streamUnitsDescr( std::ostream& os ) const { nc_assert(m_suf); m_suf(os); }
      const char * actualUnitName() const { return m_aun; }
      const char * valueTypeDescr() const { return m_vtd; }

    private:
      VarFromStrFct m_fromStrFct;
      StreamVarFct m_streamFct;
      BufCmpFct m_bcf;
      VarGroupId m_groupId;
      const char * m_namecstr;
      StrView m_name;//NB: Could potentially keep name in local buffer, for
                     //better search locality (ShortStr<32>, but requires
                     //constexpr ShortStr constructor)?
      const char * m_description;
      StreamFct m_sdf;
      StreamFct m_sdfjson;
      StreamFct m_suf;
      const char * m_aun;
      const char * m_vtd;
      BufToJSONFct m_jsonfct;
    };

    template<class TVarDef>
    inline constexpr VarInfo make_varinfo()
    {
      static_assert( constexpr_strlen(TVarDef::name) >= 3,"" );
      static_assert( constexpr_strlen(TVarDef::name) <= 17,"" );//"absorptionfactory" is 17
      return VarInfo( TVarDef::from_str,
                      TVarDef::stream_val,
                      TVarDef::cmp,
                      TVarDef::group,
                      TVarDef::name,
                      TVarDef::description,
                      (TVarDef::has_default_value() ? TVarDef::stream_default_value : nullptr),
                      (TVarDef::has_default_value() ? TVarDef::stream_default_value_json : nullptr),
                      ( TVarDef::units::actual_unit!=nullptr ? TVarDef::units::listAvailableUnits : VarInfo::StreamFct(nullptr) ),
                      TVarDef::units::actual_unit,
                      TVarDef::value_type_descr,
                      TVarDef::asJSONObject );
    }

    struct units_temperature {
      static constexpr auto actual_unit = "K";
      static void listAvailableUnits(std::ostream&);
      static Optional<std::pair<double,ValDbl_ShortStrOrigRep>> parse(StrView sv);
    };

    struct units_length {
      static constexpr auto actual_unit = "Aa";
      static void listAvailableUnitsNoDefault( std::ostream& );
      static void listAvailableUnits( std::ostream& );
      static Optional<std::pair<double,ValDbl_ShortStrOrigRep>>
      parse( StrView sv );
    };

    struct units_angle {
      static constexpr auto actual_unit = "rad";
      static constexpr auto units_descr = "rad [default], deg, arcmin, arcsec";
      static void listAvailableUnits(std::ostream&);
      static Optional<std::pair<double,ValDbl_ShortStrOrigRep>> parse(StrView sv);
    };

    struct units_purenumberonly {
      static constexpr const char* actual_unit = nullptr;
      static void listAvailableUnits(std::ostream&){}
      static Optional<std::pair<double,ValDbl_ShortStrOrigRep>> parse(StrView sv);
    };

    template <class Derived>
    class ValRot : public ValBase<Derived,std::pair<LabAxis,double>> {
    public:
      //Rotation around a lab axis, stored as the axis and the angle in rad (4
      //doubles) and the normalised text (which keeps the angle unit).
      using Base = ValRot;
      using value_type = std::pair<LabAxis,double>;
      static constexpr auto value_type_descr = "rotation";

      static value_type get_val( const VarBuf& buf )
      {
        const auto d = decodeDblsTagText<4>( buf );
        return { LabAxis{ d.vals[0], d.vals[1], d.vals[2] }, d.vals[3] };
      }

      static int cmp( const VarBuf& a, const VarBuf& b )
      {
        return cmpDblsTagText<4>( a, b );
      }

      static VarBuf set_val( VarId varid, const value_type& r )
      {
        const auto a = checkVector( r.first.as<Vector>(), Derived::name );
        const double angle = checkFiniteDbl( r.second, Derived::name );
        return detail_encode( varid, a, angle, vectorStr( a ) + ':'
                              + dbl2shortstr( angle ).to_view().to_string() );
      }

      static void stream_val( std::ostream& os, const VarBuf& buf )
      {
        os << decodeDblsTagText<4>( buf ).text;
      }

      static void asJSONObject( std::ostream& os, const VarBuf& buf )
      {
        const auto r = get_val( buf );
        os << "{\"axis\":";
        streamJSON( os, r.first );
        os << ",\"angle\":";
        streamJSON( os, r.second );
        os << '}';
      }

      static VarBuf from_str( VarId varid, StrView sv_in )
      {
        //"a1,a2,a3:angle", with optional angle unit and spaces allowed around
        //all delimiters or at string ends.
        standardInputStrSanityCheck( Derived::name, sv_in );
        auto throw_error = [&sv_in]() { throwSyntaxError( Derived::name, sv_in ); };
        const auto i = sv_in.find(':');
        if ( i == StrView::npos )
          throw_error();
        const auto a = parseVector( sv_in.substr(0,i), Derived::name );
        const auto ang = units_angle::parse( sv_in.substr(i+1) );
        if ( !a.has_value() || !ang.has_value() )
          throw_error();
        const double angle = checkFiniteDbl( ang.value().first, Derived::name );
        //Angle text as given (normalised by units_angle::parse), but with a
        //plain number shortened as for ValDbl:
        const auto generated = dbl2shortstr( angle );
        StrView atext = ang.value().second.to_view().trimmed();
        if ( atext.empty() || atext.toDbl().has_value() )
          atext = shorterNumStr( atext.empty() ? generated.to_view() : atext,
                                 generated );
        return detail_encode( varid, a.value().first, angle,
                              a.value().second + ':' + atext.to_string() );
      }

    private:
      static VarBuf detail_encode( VarId varid, const Vector& a, double angle,
                                   const std::string& text )
      {
        return encodeDblsTagText<4>( varid, { a[0], a[1], a[2], angle },
                                     0, text );
      }
    };

    namespace detail {
      template<class T>
      inline constexpr unsigned constexpr_name2Idx_helper( const T& list, const char * name, unsigned i)
      {
        return constexpr_cstrequal(name,list[i].name()) ? i : constexpr_name2Idx_helper(list, name,i+1);
      }
    }

    template<class T>
    inline constexpr unsigned constexpr_name2Idx( const T& list, const char * name )
    {
      return detail::constexpr_name2Idx_helper<T>(list,name,0);
    }

  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace Cfg {

    template<std::size_t N>
    inline VarBuf encodeDblsTagText( VarId varid,
                                     const std::array<double,N>& vals,
                                     unsigned char tag, StrView text )
    {
      std::string buf( N*sizeof(double) + 1 + text.size() + 1, '\0' );
      std::memcpy( &buf[0], vals.data(), N*sizeof(double) );
      buf[N*sizeof(double)] = static_cast<char>( tag );
      std::memcpy( &buf[N*sizeof(double)+1], text.data_begin(), text.size() );
      return VarBuf( buf.data(), buf.size(), varid );
    }

    template<std::size_t N>
    inline DblsTagText<N> decodeDblsTagText( const VarBuf& buf )
    {
      const char * data = buf.data();
      DblsTagText<N> res;
      std::memcpy( res.vals.data(), data, N*sizeof(double) );
      res.tag = static_cast<unsigned char>( data[N*sizeof(double)] );
      res.text = StrView( data + N*sizeof(double) + 1 );
      return res;
    }

    template<std::size_t N>
    inline int cmpDblsTagText( const VarBuf& a, const VarBuf& b )
    {
      const auto da = decodeDblsTagText<N>( a );
      const auto db = decodeDblsTagText<N>( b );
      if ( da.tag != db.tag )
        return da.tag > db.tag ? -1 : 1;
      for ( auto i : ncrange(N) )
        if ( da.vals[i] != db.vals[i] )
          return da.vals[i] < db.vals[i] ? -1 : 1;
      return da.text == db.text ? 0 : ( da.text < db.text ? -1 : 1 );
    }

  }
}

#endif
