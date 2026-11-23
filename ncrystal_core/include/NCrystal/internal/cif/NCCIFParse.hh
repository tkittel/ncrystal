#ifndef NCrystal_CIFParse_hh
#define NCrystal_CIFParse_hh

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

#include "NCrystal/text/NCTextData.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCrystal/internal/utils/NCSpan.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace CIF {

    //Parser for the syntax of CIF 1.1 data files (data blocks with single
    //data items and loops, no save frames or global_ blocks).

    //Whether Document::parse accepts some deviations from the specification
    //(line and name lengths, non-ASCII bytes in values and comments):
    enum class CIFParseLenient { YES, NO };

    //A data value. Text is without quotes or text field delimiters (text
    //fields have their line endings normalised to "\n", and in lenient mode
    //do not start with the line ending after the opening ";" if nothing else
    //is on that line, cf. Document::parse). Unknown and Inapplicable are the
    //unquoted "?" and "." (text is then "?" or "."):
    enum class ValueKind : std::uint8_t { Unquoted, SingleQuoted,
                                          DoubleQuoted, TextField,
                                          Unknown, Inapplicable };
    struct Value {
      StrView text;
      ValueKind kind;
      std::uint32_t line;//line number in the data (1-based)
    };

    //True for the unquoted "?" and ".":
    inline bool isNull( const Value& v )
    {
      return v.kind == ValueKind::Unknown || v.kind == ValueKind::Inapplicable;
    }

    class Document;
    class BlockView;
    class LoopView;

    //View (like Span) of the values of one tag: a single value, or the
    //column of its loop (every nTags'th value). Does not own the values,
    //which stay in the Document:
    class ValueColumn {
    public:
      std::size_t size() const { return m_size; }
      bool empty() const { return m_size == 0; }
      bool inLoop() const { return m_inLoop; }
      const Value& operator[]( std::size_t i ) const;//nc_assert on range
      ValueColumn() = default;//empty
    private:
      friend class BlockView;
      friend class LoopView;
      ValueColumn( const Value*, std::uint32_t size,
                   std::uint32_t stride, bool inLoop );
      const Value* m_first = nullptr;//into Document::m_values
      std::uint32_t m_size = 0;
      std::uint32_t m_stride = 1;//nTags for loop columns, 1 otherwise
      bool m_inLoop = false;
    };

    //View of a loop in the Document: its tags and its values row by row:
    class LoopView {
    public:
      std::size_t nTags() const;
      std::size_t nRows() const;
      StrView tag( std::size_t itag ) const;
      const Value& value( std::size_t irow, std::size_t itag ) const;
      //Column of a tag in this loop, or NullOpt if not one of its tags:
      Optional<std::size_t> column( StrView tag ) const;
      //All values of the tag in column itag:
      ValueColumn values( std::size_t itag ) const;
    private:
      friend class BlockView;
      LoopView( const Document* d, std::uint32_t i ) : m_doc(d), m_item(i) {}
      const Document* m_doc;
      std::uint32_t m_item;
    };

    //View of a data block (data_<name>) in the Document:
    class BlockView {
    public:
      StrView name() const;
      std::uint32_t line() const;//line of the data_ header

      //Whether the tag is present in the block (loaded or filtered out):
      bool has( StrView tag ) const;

      //Values of the tag: empty if absent, one value for a single item, or
      //the column of its loop. Throws LogicError if the tag is present but
      //its values were not loaded (cf. filterValuesByPrefixes in
      //Document::parse):
      ValueColumn values( StrView tag ) const;

      //The loop containing the tag, or NullOpt if absent or not in a loop.
      //Throws LogicError if present but not loaded:
      Optional<LoopView> loopOf( StrView tag ) const;

      //All tags in the block, in file order (as written in the data):
      std::size_t nTags() const;
      StrView tag( std::size_t i ) const;
    private:
      friend class Document;
      BlockView( const Document* d, std::uint32_t i ) : m_doc(d), m_block(i) {}
      //Index in Document::m_tags of the tag, or NullOpt:
      Optional<std::uint32_t> findTag( StrView ) const;
      //Throws LogicError if values of the tag (index in m_tags) not loaded:
      void checkLoaded( std::uint32_t itag ) const;
      const Document* m_doc;
      std::uint32_t m_block;
    };

    //A parsed CIF 1.1 document. Keeps the TextData alive, and all StrView's
    //point into it (or into an internal store of normalised text fields),
    //so they stay valid as long as the document exists, also after moves.
    //Views (BlockView etc.) refer to the Document object itself, so must be
    //obtained after any moves:
    class Document final : private MoveOnly {
    public:
      //Parse CIF 1.1 data. Throws BadInput, with the data source name and
      //line number, on any syntax error. If filterValuesByPrefixes is given,
      //values are only stored for tags starting with one of the prefixes
      //(case-insensitive, e.g. "_cell_"); a loop is stored if any of its tags
      //match. Values of all other tags are only counted (an empty list
      //stores no values at all). Tags are always stored, and all syntax
      //checks are done in all cases. With CIFParseLenient::YES, besides the
      //deviations named above, a leading UTF-8 byte-order mark (EF BB BF) is
      //skipped, and a text field with nothing after the opening ";" on its
      //line does not include that line's ending (";\nfoo\n;" gives "foo",
      //strictly it is "\nfoo").
      //
      //NB: Fallback for CIF data downloaded from materialsproject.org by
      //older NCrystal versions, which only noted the id in a final comment
      //("# Note from NCrystal.cifutils: Data from materialsproject.org /
      //mp-<digits>", searched in the last 500 bytes): if the last block lacks
      //the tag _[local]_ncrystal_mp_database_code, it is added with the value
      //"mp-<digits>" from the note.
      static Document parse( shared_obj<const TextData>,
                             CIFParseLenient = CIFParseLenient::YES,
                             Optional<Span<const StrView>>
                             filterValuesByPrefixes = NullOpt );

      std::size_t nBlocks() const { return m_blocks.size(); }
      BlockView block( std::size_t i ) const;
      Optional<BlockView> findBlock( StrView name ) const;//case-insensitive
      const DataSourceName& dataSourceName() const;
      //Unique per parsed document (kept when moved):
      UniqueIDValue getUniqueID() const noexcept;
      Document( Document&& ) = default;
      Document& operator=( Document&& ) = default;
    private:
      friend class BlockView;
      friend class LoopView;
      class Parser;
      Document( shared_obj<const TextData> );
      struct Item {
        std::uint32_t tagBegin, nTags;//range in m_tags
        std::uint32_t valueBegin, nValues;//range in m_values (if loaded)
        std::uint32_t nValuesInFile;//not loaded items have no m_values
        bool isLoop, loaded;
      };
      struct Tag {
        StrView name;
        std::uint32_t item;//index into m_items
        std::uint32_t column;//position within the item's tags
        std::uint32_t line;
      };
      struct Block {
        StrView name;
        std::uint32_t line;
        std::uint32_t itemBegin, nItems;//range in m_items
        std::uint32_t tagBegin, nTags;//range in m_tags (file order)
        std::uint32_t lookupBegin;//range in m_lookup has nTags entries
      };
      shared_obj<const TextData> m_data;
      UniqueID m_uid;
      std::set<std::string> m_normalised;//text fields with changed endings
      std::vector<Value> m_values;
      std::vector<Tag> m_tags;
      std::vector<Item> m_items;
      std::vector<Block> m_blocks;
      std::vector<std::uint32_t> m_lookup;//per block: indices into m_tags,
                                          //sorted by tag (case-insensitive)
      std::vector<std::uint32_t> m_blockLookup;//indices into m_blocks, sorted
                                               //by name (case-insensitive)
    };
  }
}

////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace CIF {
    inline ValueColumn::ValueColumn( const Value* f, std::uint32_t n,
                                     std::uint32_t stride, bool inLoop )
      : m_first(f), m_size(n), m_stride(stride), m_inLoop(inLoop)
    {
    }
    inline const Value& ValueColumn::operator[]( std::size_t i ) const
    {
      nc_assert( i < m_size );
      return m_first[ i * m_stride ];
    }
    inline std::size_t LoopView::nTags() const
    {
      return vectAt( m_doc->m_items, m_item ).nTags;
    }
    inline std::size_t LoopView::nRows() const
    {
      const auto& it = vectAt( m_doc->m_items, m_item );
      return it.nValues / it.nTags;
    }
    inline StrView LoopView::tag( std::size_t itag ) const
    {
      const auto& it = vectAt( m_doc->m_items, m_item );
      nc_assert( itag < it.nTags );
      return vectAt( m_doc->m_tags, it.tagBegin + itag ).name;
    }
    inline const Value& LoopView::value( std::size_t irow,
                                         std::size_t itag ) const
    {
      const auto& it = vectAt( m_doc->m_items, m_item );
      nc_assert( itag < it.nTags && irow < it.nValues / it.nTags );
      return vectAt( m_doc->m_values, it.valueBegin + irow * it.nTags + itag );
    }
    inline ValueColumn LoopView::values( std::size_t itag ) const
    {
      const auto& it = vectAt( m_doc->m_items, m_item );
      nc_assert( itag < it.nTags );
      return { m_doc->m_values.data() + it.valueBegin + itag,
               it.nValues / it.nTags, it.nTags, true };
    }
    inline StrView BlockView::name() const
    {
      return vectAt( m_doc->m_blocks, m_block ).name;
    }
    inline std::uint32_t BlockView::line() const
    {
      return vectAt( m_doc->m_blocks, m_block ).line;
    }
    inline bool BlockView::has( StrView tag ) const
    {
      return findTag( tag ).has_value();
    }
    inline std::size_t BlockView::nTags() const
    {
      return vectAt( m_doc->m_blocks, m_block ).nTags;
    }
    inline StrView BlockView::tag( std::size_t i ) const
    {
      const auto& b = vectAt( m_doc->m_blocks, m_block );
      nc_assert( i < b.nTags );
      return vectAt( m_doc->m_tags, b.tagBegin + i ).name;
    }
    inline BlockView Document::block( std::size_t i ) const
    {
      nc_assert( i < m_blocks.size() );
      return { this, static_cast<std::uint32_t>( i ) };
    }
    inline const DataSourceName& Document::dataSourceName() const
    {
      return m_data->dataSourceName();
    }
    inline UniqueIDValue Document::getUniqueID() const noexcept
    {
      return m_uid.getUniqueID();
    }
  }
}

#endif
