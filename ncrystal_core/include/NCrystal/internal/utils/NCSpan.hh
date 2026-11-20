#ifndef NCrystal_Span_hh
#define NCrystal_Span_hh

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

#include "NCrystal/core/NCDefs.hh"
#include <type_traits>
#include <iterator>

namespace NCRYSTAL_NAMESPACE {

  template <class T>
  class Span final {

    // Span class, similar to C++20's std::span, but without support for static
    // extents. Has more error checking in debug builds than in the std::span
    // (e.g. operator[] will check for out-of-range in debug builds).

  public:
    typedef T element_type;
    typedef typename std::remove_cv< T >::type value_type;
    typedef T &       reference;
    typedef T *       pointer;
    typedef T const * const_pointer;
    typedef T const & const_reference;
    typedef pointer        iterator;
    typedef const_pointer  const_iterator;
    typedef std::ptrdiff_t difference_type;
    typedef std::size_t    size_type;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    //Default constructed, empty:
    constexpr Span() noexcept : m_begin(nullptr), m_end(nullptr) {}

    //From "iterators":
    constexpr Span(pointer bb, pointer ee) ncnoexceptndebug : m_begin(bb), m_end(ee) { nc_assert(m_end>=m_begin); }

    //From other spans:
    template<class TOther> constexpr Span( Span<TOther>&& o ) noexcept : m_begin(o.begin()), m_end(o.end()) {}
    template<class TOther> constexpr Span( const Span<TOther>& o ) noexcept : m_begin(o.begin()), m_end(o.end()) {}
    template<class TOther> constexpr Span( Span<TOther>& o ) noexcept : m_begin(o.begin()), m_end(o.end()) {}

    //From C-style arrays:
    template<size_type N> constexpr Span( value_type (&carray)[N] ) noexcept : m_begin(&carray[0]), m_end(&carray[0]+N) {}
    template<size_type N> constexpr Span( const value_type (&carray)[N] ) noexcept : m_begin(&carray[0]), m_end(&carray[0]+N) {}

    //From containers like std::vector or std::array with data() and size():
    template<class TAny> constexpr Span( TAny&& v ) noexcept
    {
      size_type vsize = static_cast<size_type>(v.size());
      if (vsize>0) {
        m_begin = v.data();
        m_end = m_begin + vsize;
      } else {
        m_begin = m_end = nullptr;
      }
    }

    //Assignment:
    template<class TAny>
    constexpr Span& operator=(TAny&& t) noexcept
    {
      Span<T> sp(t);
      m_begin = sp.m_begin;
      m_end = sp.m_end;
      return *this;
    }

    //Destructor:
    ~Span() = default;

    //Iterators:
    constexpr iterator begin() noexcept { return m_begin; }
    constexpr iterator end() noexcept { return m_end; }
    constexpr const_iterator begin() const noexcept { return m_begin; }
    constexpr const_iterator end() const noexcept { return m_end; }
    constexpr const_iterator cbegin() const noexcept { return m_begin; }
    constexpr const_iterator cend() const noexcept { return m_end; }
    constexpr reverse_iterator rbegin() noexcept { return reverse_iterator(end()); }
    constexpr reverse_iterator rend() noexcept { return reverse_iterator(begin()); }
    constexpr const_reverse_iterator rbegin() const noexcept { return reverse_iterator(end()); }
    constexpr const_reverse_iterator rend() const noexcept { return reverse_iterator(begin()); }

    //Accessors:
    constexpr size_type size() const noexcept { return static_cast<size_type>(m_end-m_begin); }
    constexpr size_type size_bytes() const noexcept { return size() * sizeof(element_type); }
    constexpr bool empty() const noexcept { return m_begin == m_end; }
    constexpr pointer data() noexcept { return m_begin; }
    constexpr const_pointer data() const noexcept { return m_begin; }

    constexpr reference operator[](size_type i) ncnoexceptndebug { nc_assert(i<size()); return *(m_begin+i); }
    constexpr const_reference operator[](size_type i) const ncnoexceptndebug  { nc_assert(i<size()); return *(m_begin+i); }
    constexpr reference at(size_type i) { nc_assert_always(i<size()); return *(m_begin+i); }
    constexpr const_reference at(size_type i) const { nc_assert_always(i<size()); return *(m_begin+i); }
    constexpr reference front() const ncnoexceptndebug { nc_assert(!empty()&&m_begin); return *m_begin; }
    constexpr reference back() const ncnoexceptndebug { nc_assert(!empty()&&m_begin); return *std::prev(m_end); }

    //Subspans:
    constexpr Span subspan(size_type offset, size_type count = std::numeric_limits<size_type>::max() ) const ncnoexceptndebug {
      auto n = size();
      if (offset>=n)
        return Span();//empty
      auto mbo = m_begin+offset;
      if ( count >= n - offset )//like this, since count+offset might overflow
        return Span(mbo,m_end);
      //return Span(m_begin+offset,m_begin+std::min<size_type>(n,offset+count));
      nc_assert(mbo+count<m_end);
      return Span(mbo,mbo+count);
    }

    constexpr Span first(size_type count) const ncnoexceptndebug {
      auto n = size();
      if (n==0||count==0)
        return Span();
      return Span(m_begin,m_begin+std::min<size_type>(n,count) );
    }

    constexpr Span last(size_type count) const ncnoexceptndebug {
      auto n = size();
      if ( count >= n )
        return Span(m_begin,m_end);
      if (count==0)
        return Span();
      return Span( m_begin+(n-count), m_end );
    }

  private:
    pointer m_begin;
    pointer m_end;
  };


}

#endif
