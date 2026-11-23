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

#include "NCrystal/internal/cif/NCCIFParse.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC = NCrystal;
namespace NCC = NCrystal::CIF;

namespace NCRYSTAL_NAMESPACE {
  namespace CIF {
    namespace {

      //The data is null-terminated and null-free (cf. RawStrData), so '\0'
      //marks its end. CIF 1.1 character set: isSimpleASCII (syntax,
      //paragraph 22). Fast paths for the bytes allowed inside words (33-126)
      //and inside quoted values, text fields and comments (32-126):
      inline bool isNonBlank( char c )
      {
        return unsigned( static_cast<unsigned char>( c ) ) - 33u < 94u;
      }
      inline bool isPrintable( char c )
      {
        return unsigned( static_cast<unsigned char>( c ) ) - 32u < 95u;
      }
      inline bool isLineEnd( char c ) { return c == '\n' || c == '\r'; }
      inline bool isTokenEnd( char c )
      {
        return c == '\0' || isWhiteSpace( c );
      }

      std::string hexByte( char c )//e.g. "0xEF"
      {
        constexpr char hex[] = "0123456789ABCDEF";
        const auto u = static_cast<unsigned char>( c );
        return std::string( "0x" ) + hex[u / 16] + hex[u % 16];
      }

      inline char lowerASCII( char c )
      {
        return ( c >= 'A' && c <= 'Z' ) ? char( c - 'A' + 'a' ) : c;
      }

      int ciCompare( StrView a, StrView b )
      {
        const std::size_t n = std::min( a.size(), b.size() );
        for ( std::size_t i = 0; i < n; ++i ) {
          auto ca = static_cast<unsigned char>( lowerASCII( a[i] ) );
          auto cb = static_cast<unsigned char>( lowerASCII( b[i] ) );
          if ( ca != cb )
            return ca < cb ? -1 : 1;
        }
        return a.size() == b.size() ? 0 : ( a.size() < b.size() ? -1 : 1 );
      }

      bool ciStartsWith( StrView s, StrView prefix )
      {
        return ( s.size() >= prefix.size()
                 && ciCompare( s.substr( 0, prefix.size() ), prefix ) == 0 );
      }

      //Sort indices of named objects case-insensitively (ties in index
      //order), and return the pair of equal names whose second has the
      //lowest index, if any:
      template<class TGetName>
      Optional<std::pair<std::uint32_t,std::uint32_t>>
      sortAndFindDuplicate( std::uint32_t* b, std::uint32_t* e,
                            TGetName getName )
      {
        std::sort( b, e, [&getName]( std::uint32_t i, std::uint32_t j )
        {
          const int c = ciCompare( getName( i ), getName( j ) );
          return c == 0 ? i < j : c < 0;
        } );
        Optional<std::pair<std::uint32_t,std::uint32_t>> res;
        for ( auto it = b; it != e && std::next( it ) != e; ++it ) {
          const auto i = *it;
          const auto j = *std::next( it );
          if ( ciCompare( getName( i ), getName( j ) ) == 0
               && ( !res.has_value() || j < res.value().second ) )
            res.emplace( i, j );
        }
        return res;
      }

      //Start of a value for error messages (first line, at most 40 chars,
      //bytes outside printable ASCII as e.g. "\xEF"):
      std::string shortText( StrView s )
      {
        std::size_t n = 0;
        while ( n < s.size() && n < 40 && !isLineEnd( s[n] ) )
          ++n;
        std::string res;
        for ( char c : s.substr( 0, n ) ) {
          if ( isPrintable( c ) )
            res += c;
          else
            res += "\\x" + hexByte( c ).substr( 2 );
        }
        return n < s.size() ? res + "..." : res;
      }

      std::string normaliseLineEndings( StrView s )
      {
        std::string res;
        res.reserve( s.size() );
        for ( std::size_t i = 0; i < s.size(); ++i ) {
          if ( s[i] != '\r' ) {
            res += s[i];
          } else {
            res += '\n';
            if ( i + 1 < s.size() && s[i + 1] == '\n' )
              ++i;
          }
        }
        return res;
      }

      constexpr std::uint32_t s_maxU32
        = std::numeric_limits<std::uint32_t>::max();

      //NB: MP note fallback, cf. Document::parse:
      constexpr auto s_mpTag
        = StrView::make("_[local]_ncrystal_mp_database_code");
      constexpr auto s_mpNote
        = StrView::make("# Note from NCrystal.cifutils: Data from"
                        " materialsproject.org / mp-");
    }

    class Document::Parser {
    public:
      Parser( shared_obj<const TextData>, CIFParseLenient,
              Optional<Span<const StrView>> );
      Document parse();

    private:
      enum class TokType { End, Tag, Value, Loop, Data, Save, Global, Stop };
      struct Token {
        TokType type;
        StrView text;//tag, value, or block name for Data
        ValueKind kind;
        std::uint32_t line;
        bool hasCR;//text field with CR line endings
      };

      Document m_doc;
      const DataSourceName& m_dsn;
      Optional<Span<const StrView>> m_filterPrefixes;
      const char* m_p;
      const char* m_lineStart;
      std::uint32_t m_line = 1;
      bool m_strict;

      [[noreturn]] void throwBadInput( std::uint32_t line,
                                       const std::string& reason );
      //For text containing bytes outside the CIF 1.1 character set:
      [[noreturn]] void throwBadByte( std::uint32_t line, StrView text );
      void checkLineLength( const char* eol );
      //Advance to the next line ending or the end of the data, returning
      //whether bytes outside the CIF 1.1 character set were passed:
      bool skipToLineEnd();
      void consumeEOL();
      void skipComment();
      void skipWhiteSpace();
      Token nextToken();
      Token lexQuoted();
      Token lexTextField();
      Token lexWord();
      bool storeValues( StrView tag ) const;
      void addValue( const Token& );
      void addTag( const Token&, std::uint32_t item, std::uint32_t column );
      std::uint32_t u32( std::size_t, std::uint32_t line );
      //Parse the loop starting at the loop_ token t, leaving t at the first
      //token after the loop:
      void parseLoop( Token& t );
      void addMPNoteFallback();
      //Finish the last block: tag lookup table, and duplicate tag check:
      void endBlock();
      [[noreturn]] void throwUnexpectedToken( const Token&,
                                              const char* context );
    };
  }
}

NCC::Document::Document( shared_obj<const TextData> td )
  : m_data( std::move( td ) )
{
}

NCC::Document::Parser::Parser( shared_obj<const TextData> td,
                               CIFParseLenient lenient,
                               Optional<Span<const StrView>> filterPrefixes )
  : m_doc( std::move( td ) ),
    m_dsn( m_doc.m_data->dataSourceName() ),
    m_filterPrefixes( filterPrefixes ),
    m_p( m_doc.m_data->rawData().begin() ),
    m_lineStart( m_p ),
    m_strict( lenient == CIFParseLenient::NO )
{
  nc_assert_always( *m_doc.m_data->rawData().end() == '\0' );
}

void NCC::Document::Parser::throwBadInput( std::uint32_t line,
                                           const std::string& reason )
{
  NCRYSTAL_THROW2( BadInput, "Invalid CIF data in \"" << m_dsn
                   << "\" (line " << line << "): " << reason );
}

void NCC::Document::Parser::throwBadByte( std::uint32_t line, StrView text )
{
  for ( char c : text ) {
    if ( !isSimpleASCII( c ) ) {
      if ( m_strict )
        throwBadInput( line, "byte " + hexByte( c ) + " is not in the CIF 1.1"
                       " character set (printable ASCII, tab and line"
                       " endings)" );
      throwBadInput( line, "byte " + hexByte( c ) + " is not allowed outside"
                     " values and comments" );
    }
  }
  nc_assert_always( false );
  throwBadInput( line, "unexpected error" );
}

bool NCC::Document::Parser::skipToLineEnd()
{
  bool bad = false;
  while ( true ) {
    const char c = *m_p;
    if ( isPrintable( c ) || c == '\t' ) {
      ++m_p;
    } else if ( isLineEnd( c ) || c == '\0' ) {
      return bad;
    } else {
      bad = true;
      ++m_p;
    }
  }
}

void NCC::Document::Parser::checkLineLength( const char* eol )
{
  //CIF 1.1 syntax, paragraph 28 (not counting the line ending):
  if ( m_strict && eol - m_lineStart > 2048 )
    throwBadInput( m_line, "line longer than 2048 characters" );
}

void NCC::Document::Parser::consumeEOL()
{
  nc_assert( *m_p == '\n' || *m_p == '\r' );
  checkLineLength( m_p );
  //LF, CRLF and CR (CIF 1.1 syntax, paragraphs 25 and 42):
  if ( *m_p == '\r' && m_p[1] == '\n' )
    ++m_p;
  ++m_p;
  m_lineStart = m_p;
  if ( m_line == s_maxU32 )
    throwBadInput( m_line, "too many lines" );
  ++m_line;
}

void NCC::Document::Parser::skipComment()
{
  nc_assert( *m_p == '#' );
  const char* b = m_p;
  if ( skipToLineEnd() && m_strict )
    throwBadByte( m_line, StrView( b, m_p - b ) );
}

void NCC::Document::Parser::skipWhiteSpace()
{
  //Only called where a token has ended (always at white space) or at the
  //start of the data, so a "#" here starts a comment (CIF 1.1 syntax,
  //paragraph 45):
  while ( true ) {
    const char c = *m_p;
    if ( c == ' ' || c == '\t' )
      ++m_p;
    else if ( c == '\n' || c == '\r' )
      consumeEOL();
    else if ( c == '#' )
      skipComment();
    else
      return;
  }
}

NCC::Document::Parser::Token NCC::Document::Parser::nextToken()
{
  skipWhiteSpace();
  const char c = *m_p;
  if ( c == '\0' ) {
    checkLineLength( m_p );
    return { TokType::End, StrView(), ValueKind::Unquoted, m_line, false };
  }
  //CIF 1.1 syntax, paragraphs 17 and 57 (";" only delimits text fields at
  //the start of a line, elsewhere it can start an unquoted value):
  if ( c == ';' && m_p == m_lineStart )
    return lexTextField();
  if ( c == '\'' || c == '"' )
    return lexQuoted();
  return lexWord();
}

NCC::Document::Parser::Token NCC::Document::Parser::lexQuoted()
{
  //CIF 1.1 syntax, paragraphs 14, 15 and 56: ends at the first matching
  //quote followed by white space, on the same line:
  const char q = *m_p++;
  const char* b = m_p;
  bool bad = false;
  while ( true ) {
    const char c = *m_p;
    if ( c == q && isTokenEnd( m_p[1] ) )
      break;
    if ( isPrintable( c ) || c == '\t' ) {
      ++m_p;
    } else if ( isLineEnd( c ) || c == '\0' ) {
      throwBadInput( m_line, q == '\'' ? "unterminated single-quoted value"
                     : "unterminated double-quoted value" );
    } else {
      bad = true;
      ++m_p;
    }
  }
  StrView text( b, m_p - b );
  ++m_p;
  if ( bad && m_strict )
    throwBadByte( m_line, text );
  return { TokType::Value, text,
           q == '\'' ? ValueKind::SingleQuoted : ValueKind::DoubleQuoted,
           m_line, false };
}

NCC::Document::Parser::Token NCC::Document::Parser::lexTextField()
{
  //CIF 1.1 syntax, paragraphs 17, 18 and 56: ends at the next line starting
  //with ";", which must be followed by white space. The line ending before
  //it is not part of the value:
  const std::uint32_t line0 = m_line;
  const char* b = ++m_p;
  const bool skipFirstEOL = !m_strict && isLineEnd( *m_p );
  bool first = true;
  bool bad = false;
  bool hasCR = false;
  while ( true ) {
    bad = skipToLineEnd() || bad;
    if ( *m_p == '\0' )
      throwBadInput( line0, "unterminated text field" );
    const char* eol = m_p;
    consumeEOL();
    if ( first && skipFirstEOL )
      b = m_p;
    first = false;
    hasCR = hasCR || ( eol >= b && *eol == '\r' );
    if ( *m_p == ';' ) {
      StrView text( b, eol >= b ? eol - b : 0 );
      ++m_p;
      if ( !isTokenEnd( *m_p ) )
        throwBadInput( m_line, "text field terminator \";\" not followed by"
                       " white space" );
      if ( bad && m_strict )
        throwBadByte( line0, text );
      return { TokType::Value, text, ValueKind::TextField, line0, hasCR };
    }
  }
}

NCC::Document::Parser::Token NCC::Document::Parser::lexWord()
{
  const char* b = m_p;
  bool extra = false;
  while ( true ) {
    const char c = *m_p;
    if ( isNonBlank( c ) ) {
      ++m_p;
    } else if ( isTokenEnd( c ) ) {
      break;
    } else {
      extra = true;
      ++m_p;
    }
  }
  const StrView w( b, m_p - b );
  Token t{ TokType::Value, w, ValueKind::Unquoted, m_line, false };
  const char c = w[0];
  if ( c == '_' ) {
    if ( extra )
      throwBadByte( m_line, w );
    //CIF 1.1 syntax, paragraphs 29 and 57:
    if ( w.size() == 1 )
      throwBadInput( m_line, "tag \"_\" without a name" );
    if ( m_strict && w.size() > 75 )
      throwBadInput( m_line, "tag longer than 75 characters" );
    t.type = TokType::Tag;
    return t;
  }
  //CIF 1.1 syntax, paragraphs 19, 32 and 50:
  if ( c == '[' || c == ']' || c == '$' )
    throwBadInput( m_line, std::string( "unquoted value starting with"
                                        " reserved character \"" )
                   + c + "\"" );
  //Reserved words (CIF 1.1 syntax, paragraphs 11 and 57; "data_" and
  //"save_" as prefixes):
  if ( ciStartsWith( w, "data_" ) ) {
    if ( extra )
      throwBadByte( m_line, w );
    t.type = TokType::Data;
    t.text = w.substr( 5 );
    if ( t.text.empty() )
      throwBadInput( m_line, "empty block name" );
    //CIF 1.1 syntax, paragraph 30:
    if ( m_strict && t.text.size() > 75 )
      throwBadInput( m_line, "block name longer than 75 characters" );
    return t;
  }
  if ( ciStartsWith( w, "save_" ) )
    t.type = TokType::Save;
  else if ( ciCompare( w, "loop_" ) == 0 )
    t.type = TokType::Loop;
  else if ( ciCompare( w, "global_" ) == 0 )
    t.type = TokType::Global;
  else if ( ciCompare( w, "stop_" ) == 0 )
    t.type = TokType::Stop;
  if ( t.type != TokType::Value )
    return t;
  if ( m_strict && extra )
    throwBadByte( m_line, w );
  //CIF 1.1 syntax, paragraph 57:
  if ( w.size() == 1 && ( c == '?' || c == '.' ) )
    t.kind = ( c == '?' ? ValueKind::Unknown : ValueKind::Inapplicable );
  return t;
}

bool NCC::Document::Parser::storeValues( StrView tag ) const
{
  if ( !m_filterPrefixes.has_value() )
    return true;
  for ( auto& p : m_filterPrefixes.value() )
    if ( ciStartsWith( tag, p ) )
      return true;
  return false;
}

std::uint32_t NCC::Document::Parser::u32( std::size_t n, std::uint32_t line )
{
  if ( n > s_maxU32 )
    throwBadInput( line, "too much data (count exceeds 32 bit range)" );
  return static_cast<std::uint32_t>( n );
}

void NCC::Document::Parser::addValue( const Token& t )
{
  nc_assert( t.type == TokType::Value );
  if ( !t.hasCR ) {
    m_doc.m_values.push_back( Value{ t.text, t.kind, t.line } );
    return;
  }
  const auto& s = *m_doc.m_normalised.insert( normaliseLineEndings( t.text ) )
    .first;
  m_doc.m_values.push_back( Value{ StrView( s ), t.kind, t.line } );
}

void NCC::Document::Parser::addTag( const Token& t, std::uint32_t item,
                                    std::uint32_t column )
{
  u32( m_doc.m_tags.size() + 1, t.line );
  m_doc.m_tags.push_back( Tag{ t.text, item, column, t.line } );
}

void NCC::Document::Parser::throwUnexpectedToken( const Token& t,
                                                  const char* context )
{
  switch ( t.type ) {
  case TokType::Save:
    //Allowed for data-file-only parsers (CIF 1.1 syntax, paragraph 5):
    throwBadInput( t.line, "save frames are not supported" );
  case TokType::Global:
    throwBadInput( t.line, "global_ is not supported (CIF 1.1 syntax,"
                   " paragraph 33)" );
  case TokType::Stop:
    throwBadInput( t.line, "stop_ is not supported (CIF 1.1 syntax,"
                   " paragraph 11)" );
  default:
    break;
  }
  const std::string what
    = ( t.type == TokType::Tag ? "tag " + t.text.to_string()
        : ( t.type == TokType::Loop ? std::string( "loop_" )
            : "value \"" + shortText( t.text ) + "\"" ) );
  throwBadInput( t.line, what + " " + context );
}

void NCC::Document::Parser::parseLoop( Token& t )
{
  nc_assert( t.type == TokType::Loop );
  const std::uint32_t loopLine = t.line;
  auto& d = m_doc;
  const auto item = u32( d.m_items.size(), loopLine );
  const auto tagBegin = u32( d.m_tags.size(), loopLine );
  bool doStore = false;
  t = nextToken();
  for ( ; t.type == TokType::Tag; t = nextToken() ) {
    addTag( t, item, static_cast<std::uint32_t>( d.m_tags.size() - tagBegin ) );
    doStore = doStore || storeValues( t.text );
  }
  const auto nTags = static_cast<std::uint32_t>( d.m_tags.size() - tagBegin );
  if ( nTags == 0 )
    throwBadInput( loopLine, "loop_ without tags" );
  const auto valueBegin = u32( d.m_values.size(), loopLine );
  std::uint64_t n = 0;
  if ( doStore ) {
    for ( ; t.type == TokType::Value; t = nextToken(), ++n )
      addValue( t );
  } else {
    for ( ; t.type == TokType::Value; t = nextToken() )
      ++n;
  }
  if ( n == 0 )
    throwBadInput( loopLine, "loop without values" );
  //CIF 1.1 syntax, paragraph 63:
  if ( n % nTags )
    throwBadInput( loopLine, "loop has " + std::to_string( n ) + " values,"
                   " which is not a multiple of its "
                   + std::to_string( nTags ) + " tags" );
  const auto nInFile = u32( n, loopLine );
  u32( d.m_values.size(), loopLine );
  d.m_items.push_back( Item{ tagBegin, nTags, valueBegin,
                             doStore ? nInFile : 0, nInFile, true, doStore } );
}

void NCC::Document::Parser::addMPNoteFallback()
{
  //NB: MP note fallback, cf. Document::parse. Called at the end of the data,
  //before ending the last block:
  auto& d = m_doc;
  const auto& b = d.m_blocks.back();
  for ( auto i : ncrange( b.tagBegin, u32( d.m_tags.size(), m_line ) ) )
    if ( ciCompare( vectAt( d.m_tags, i ).name, s_mpTag ) == 0 )
      return;
  const auto& raw = d.m_data->rawData();
  const StrView all( raw.begin(), raw.end() - raw.begin() );
  const std::size_t n = std::min<std::size_t>( all.size(), 500 );
  const StrView tail = all.substr( all.size() - n );
  const auto pos = tail.find( s_mpNote );
  if ( pos == StrView::npos )
    return;
  const auto ib = pos + s_mpNote.size() - 3;//at "mp-"
  auto ie = ib + 3;
  while ( ie < tail.size() && tail[ie] >= '0' && tail[ie] <= '9' )
    ++ie;
  if ( ie == ib + 3 || ( ie < tail.size() && !isWhiteSpace( tail[ie] ) ) )
    return;
  //Line of the note (m_line is the last line):
  std::uint32_t line = m_line;
  char prev = '\0';
  for ( char c : tail.substr( pos ) ) {
    if ( c == '\r' || ( c == '\n' && prev != '\r' ) )
      --line;
    prev = c;
  }
  const auto item = u32( d.m_items.size(), line );
  const auto valueBegin = u32( d.m_values.size(), line );
  addTag( Token{ TokType::Tag, s_mpTag, ValueKind::Unquoted, line, false },
          item, 0 );
  const bool k = storeValues( s_mpTag );
  if ( k )
    d.m_values.push_back( Value{ tail.substr( ib, ie - ib ),
                                 ValueKind::Unquoted, line } );
  d.m_items.push_back( Item{ u32( d.m_tags.size() - 1, line ), 1, valueBegin,
                             k ? 1u : 0u, 1, false, k } );
}

void NCC::Document::Parser::endBlock()
{
  auto& d = m_doc;
  auto& b = d.m_blocks.back();
  b.nItems = static_cast<std::uint32_t>( d.m_items.size() - b.itemBegin );
  b.nTags = static_cast<std::uint32_t>( d.m_tags.size() - b.tagBegin );
  b.lookupBegin = static_cast<std::uint32_t>( d.m_lookup.size() );
  for ( auto i : ncrange( b.tagBegin, b.tagBegin + b.nTags ) )
    d.m_lookup.push_back( i );
  //CIF 1.1 syntax, paragraph 7:
  auto dup = sortAndFindDuplicate( d.m_lookup.data() + b.lookupBegin,
                                   d.m_lookup.data() + d.m_lookup.size(),
                                   [&d]( std::uint32_t i )
                                   { return vectAt( d.m_tags, i ).name; } );
  if ( dup.has_value() ) {
    const auto& t1 = vectAt( d.m_tags, dup.value().first );
    const auto& t2 = vectAt( d.m_tags, dup.value().second );
    throwBadInput( t2.line, "duplicate tag " + t2.name.to_string()
                   + " in block (first in line " + std::to_string( t1.line )
                   + ")" );
  }
}

NCC::Document NCC::Document::Parser::parse()
{
  auto& d = m_doc;
  d.m_values.reserve( 1024 );
  d.m_tags.reserve( 1024 );
  d.m_items.reserve( 1024 );
  d.m_blocks.reserve( 16 );
  d.m_lookup.reserve( 1024 );

  if ( !m_strict && std::strncmp( m_p, "\xEF\xBB\xBF", 3 ) == 0 )
    m_lineStart = ( m_p += 3 );

  //CIF 2.0 starts with this version comment (CIF 1.1 uses "#\#CIF_1.1",
  //paragraph 34):
  if ( std::strncmp( m_p, "#\\#CIF_2.0", 10 ) == 0 )
    throwBadInput( 1, "CIF 2.0 data is not supported" );

  Token t = nextToken();
  if ( t.type == TokType::End )
    throwBadInput( t.line, "no data block found" );
  if ( t.type != TokType::Data )
    throwUnexpectedToken( t, "before the first data block" );

  while ( t.type == TokType::Data ) {
    d.m_blocks.push_back( Block{ t.text, t.line,
                                 u32( d.m_items.size(), t.line ), 0,
                                 u32( d.m_tags.size(), t.line ), 0, 0 } );
    u32( d.m_blocks.size(), t.line );
    t = nextToken();
    while ( true ) {
      if ( t.type == TokType::Tag ) {
        const Token tag = t;
        const auto item = u32( d.m_items.size(), t.line );
        const auto tagIndex = u32( d.m_tags.size(), t.line );
        const auto valueBegin = u32( d.m_values.size(), t.line );
        addTag( tag, item, 0 );
        t = nextToken();
        if ( t.type == TokType::Tag || t.type == TokType::End )
          throwBadInput( tag.line, "tag " + tag.text.to_string()
                         + " without a value" );
        if ( t.type != TokType::Value ) {
          //CIF 1.1 syntax, paragraphs 8, 11 and 57:
          const std::string w = ( t.type == TokType::Data
                                  ? "data_" + shortText( t.text )
                                  : shortText( t.text ) );
          throwBadInput( t.line, "reserved word " + w + " used as an unquoted"
                         " value (of tag " + tag.text.to_string() + ")" );
        }
        const bool k = storeValues( tag.text );
        if ( k )
          addValue( t );
        d.m_items.push_back( Item{ tagIndex, 1, valueBegin, k ? 1u : 0u, 1,
                                   false, k } );
        t = nextToken();
      } else if ( t.type == TokType::Loop ) {
        parseLoop( t );
      } else if ( t.type == TokType::Data || t.type == TokType::End ) {
        break;
      } else {
        throwUnexpectedToken( t, "where a tag is expected" );
      }
    }
    if ( t.type == TokType::End )
      addMPNoteFallback();
    endBlock();
  }
  nc_assert_always( t.type == TokType::End );

  //CIF 1.1 syntax, paragraph 6:
  d.m_blockLookup.reserve( d.m_blocks.size() );
  for ( auto i : ncrange( d.m_blocks.size() ) )
    d.m_blockLookup.push_back( static_cast<std::uint32_t>( i ) );
  auto dup = sortAndFindDuplicate( d.m_blockLookup.data(),
                                   d.m_blockLookup.data()
                                   + d.m_blockLookup.size(),
                                   [&d]( std::uint32_t i )
                                   {
                                     return vectAt( d.m_blocks, i ).name;
                                   } );
  if ( dup.has_value() ) {
    const auto& b1 = vectAt( d.m_blocks, dup.value().first );
    const auto& b2 = vectAt( d.m_blocks, dup.value().second );
    throwBadInput( b2.line, "duplicate block name data_"
                   + b2.name.to_string() + " (first in line "
                   + std::to_string( b1.line ) + ")" );
  }

  d.m_values.shrink_to_fit();
  d.m_tags.shrink_to_fit();
  d.m_items.shrink_to_fit();
  d.m_blocks.shrink_to_fit();
  d.m_lookup.shrink_to_fit();
  return std::move( m_doc );
}

NCC::Document NCC::Document::parse( shared_obj<const TextData> td,
                                    CIFParseLenient lenient,
                                    Optional<Span<const StrView>> prefixes )
{
  return Document::Parser( std::move( td ), lenient, prefixes ).parse();
}

NC::Optional<std::uint32_t> NCC::BlockView::findTag( StrView tag ) const
{
  const auto& d = *m_doc;
  const auto& b = vectAt( d.m_blocks, m_block );
  auto lb = d.m_lookup.begin() + b.lookupBegin;
  auto le = lb + b.nTags;
  auto it = std::lower_bound( lb, le, tag, [&d]( std::uint32_t i, StrView key )
  {
    return ciCompare( vectAt( d.m_tags, i ).name, key ) < 0;
  } );
  if ( it == le || ciCompare( vectAt( d.m_tags, *it ).name, tag ) != 0 )
    return NullOpt;
  return *it;
}

void NCC::BlockView::checkLoaded( std::uint32_t itag ) const
{
  const auto& t = vectAt( m_doc->m_tags, itag );
  if ( !vectAt( m_doc->m_items, t.item ).loaded )
    NCRYSTAL_THROW2( LogicError, "Values of CIF tag " << t.name
                     << " in block data_" << name()
                     << " were not loaded (filtered out)" );
}

NCC::ValueColumn NCC::BlockView::values( StrView tag ) const
{
  auto itag = findTag( tag );
  if ( !itag.has_value() )
    return {};
  checkLoaded( itag.value() );
  const auto& t = vectAt( m_doc->m_tags, itag.value() );
  const auto& it = vectAt( m_doc->m_items, t.item );
  return { m_doc->m_values.data() + it.valueBegin + t.column,
           it.nValues / it.nTags, it.nTags, it.isLoop };
}

NC::Optional<NCC::LoopView> NCC::BlockView::loopOf( StrView tag ) const
{
  auto itag = findTag( tag );
  if ( !itag.has_value() )
    return NullOpt;
  checkLoaded( itag.value() );
  const auto item = vectAt( m_doc->m_tags, itag.value() ).item;
  if ( !vectAt( m_doc->m_items, item ).isLoop )
    return NullOpt;
  return LoopView( m_doc, item );
}

NC::Optional<std::size_t> NCC::LoopView::column( StrView tag ) const
{
  const auto& it = vectAt( m_doc->m_items, m_item );
  for ( auto i : ncrange( it.nTags ) )
    if ( ciCompare( vectAt( m_doc->m_tags, it.tagBegin + i ).name,
                    tag ) == 0 )
      return std::size_t( i );
  return NullOpt;
}

NC::Optional<NCC::BlockView> NCC::Document::findBlock( StrView name ) const
{
  const auto le = m_blockLookup.end();
  auto it = std::lower_bound( m_blockLookup.begin(), le, name,
                              [this]( std::uint32_t i, StrView key )
  {
    return ciCompare( vectAt( m_blocks, i ).name, key ) < 0;
  } );
  if ( it == le || ciCompare( vectAt( m_blocks, *it ).name, name ) != 0 )
    return NullOpt;
  return BlockView( this, *it );
}
