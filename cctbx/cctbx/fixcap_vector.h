// derived from boost::array
// TODO: complete std::vector interface
// XXX copyright

#ifndef FIXCAP_VECTOR_H
#define FIXCAP_VECTOR_H

#include <cstddef>
#include <stdexcept>
#include <iterator>
#include <algorithm>

// FIXES for broken compilers
#include <boost/config.hpp>

namespace cctbx {

    template<class T, std::size_t N>
    class fixcap_vector {
      public:
        T elems[N];    // fixed-size array of elements of type T

      public:
        // type definitions
        typedef T              value_type;
        typedef T*             iterator;
        typedef const T*       const_iterator;
        typedef T&             reference;
        typedef const T&       const_reference;
        typedef std::size_t    size_type;
        typedef std::ptrdiff_t difference_type;

        fixcap_vector() : m_size(0) {}
        fixcap_vector(std::size_t n) : m_size(n) {
          std::fill_n(begin(),n,value_type());
        }
        fixcap_vector(const T& value) : m_size(0) { assign(value); }

        // iterator support
        iterator begin() { return elems; }
        const_iterator begin() const { return elems; }
        iterator end() { return elems+size(); }
        const_iterator end() const { return elems+size(); }

        // reverse iterator support
#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) && !defined(BOOST_MSVC_STD_ITERATOR)
        typedef std::reverse_iterator<iterator> reverse_iterator;
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
#else
        // workaround for broken reverse_iterator implementations
        typedef std::reverse_iterator<iterator,T> reverse_iterator;
        typedef std::reverse_iterator<const_iterator,T> const_reverse_iterator;
#endif

        reverse_iterator rbegin() { return reverse_iterator(end()); }
        const_reverse_iterator rbegin() const {
            return const_reverse_iterator(end());
        }
        reverse_iterator rend() { return reverse_iterator(begin()); }
        const_reverse_iterator rend() const {
            return const_reverse_iterator(begin());
        }

        // operator[]
        reference operator[](size_type i) { return elems[i]; }
        const_reference operator[](size_type i) const { return elems[i]; }

        // at() with range check
        reference at(size_type i) { rangecheck(i); return elems[i]; }
        const_reference at(size_type i) const { rangecheck(i); return elems[i]; }

        // front() and back()
        reference front() { return elems[0]; }
        const_reference front() const { return elems[0]; }
        reference back() { return elems[size()-1]; }
        const_reference back() const { return elems[size()-1]; }

        size_type size() const { return m_size; }
        bool empty() const { if (size() == 0) return true; return false; }
        static size_type max_size() { return size(); }
        enum { fixcap_size = N };

        // swap (note: linear complexity)
        void swap (fixcap_vector<T,N>& y) {
            std::swap_ranges(begin(),end(),y.begin());
        }

        // direct access to data
        const T* data() const { return elems; }

        // assignment with type conversion
        template <typename T2>
        fixcap_vector<T,N>& operator= (const fixcap_vector<T2,N>& rhs) {
            std::copy(rhs.begin(),rhs.end(), begin());
            return *this;
        }

        // assign one value to all elements
        void assign (const T& value)
        {
            std::fill_n(begin(),N,value);
        }

        void push_back (const T& x)
        {
          if (m_size == fixcap_size)
            { throw std::range_error("fixcap_vector"); }
          elems[m_size] = x;
          m_size++;
        }

      private:
        // check range (may be private because it is fixcap)
        void rangecheck (size_type i) {
            if (i >= size()) { throw std::range_error("fixcap_vector"); }
        }

      private:
        size_type m_size;
    };

    // comparisons
    template<class T, std::size_t N>
    bool operator== (const fixcap_vector<T,N>& x, const fixcap_vector<T,N>& y) {
        return std::equal(x.begin(), x.end(), y.begin());
    }
    template<class T, std::size_t N>
    bool operator< (const fixcap_vector<T,N>& x, const fixcap_vector<T,N>& y) {
        return std::lexicographical_compare(x.begin(),x.end(),y.begin(),y.end());
    }
    template<class T, std::size_t N>
    bool operator!= (const fixcap_vector<T,N>& x, const fixcap_vector<T,N>& y) {
        return !(x==y);
    }
    template<class T, std::size_t N>
    bool operator> (const fixcap_vector<T,N>& x, const fixcap_vector<T,N>& y) {
        return y<x;
    }
    template<class T, std::size_t N>
    bool operator<= (const fixcap_vector<T,N>& x, const fixcap_vector<T,N>& y) {
        return !(y<x);
    }
    template<class T, std::size_t N>
    bool operator>= (const fixcap_vector<T,N>& x, const fixcap_vector<T,N>& y) {
        return !(x<y);
    }

    // global swap()
    template<class T, std::size_t N>
    inline void swap (fixcap_vector<T,N>& x, fixcap_vector<T,N>& y) {
        x.swap(y);
    }

} // namespace cctbx

#endif // FIXCAP_VECTOR_H
