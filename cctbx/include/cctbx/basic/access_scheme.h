/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (NK Sauter)
 */

#ifndef CCTBX_ACSCHEME_H
#define CCTBX_ACSCHEME_H

#include <cctbx/carray.h>
#include <cctbx/basic/shared_storage.h>
#include <cctbx/vector/reductions.h>
#include <cctbx/loops.h>


namespace cctbx {

  template <std::size_t N>
  struct IndexConvert {
    template <typename ExtendArrayType, typename IndexArrayType>
    std::size_t operator()(const ExtendArrayType& stride, 
                           const ExtendArrayType& base,  
                           const IndexArrayType& i) {
      return IndexConvert<N-1>()(stride,base, i) + base[N-1]+i[N-1]*stride[N-1];
    }
  };

  template<>
  struct IndexConvert<1> {
    template <typename ExtendArrayType, typename IndexArrayType>
    std::size_t operator()(const ExtendArrayType& stride, 
                           const ExtendArrayType& base,  
                           const IndexArrayType& i) {
      return base[0]+i[0]*stride[0];
    }
  };

  class Range {
    public:
      explicit Range (std::size_t begin) :
        m_begin(begin),m_end(begin),m_stride(1){}
      explicit Range (std::size_t begin,std::size_t end,std::size_t stride = 1) 
       : m_begin(begin),m_end(end),m_stride(stride){
        assert (stride>0);} //Otherwise meaningless
      static Range all() {return Range();}
      const std::size_t& begin() const {return m_begin;}
      const std::size_t& end() const {return m_end;}
      const std::size_t& stride() const {return m_stride;}
      bool isAll() {return m_stride==0;}
      std::size_t extent() const {
        assert (m_stride!=0); //internal representation of Range::all()
        assert ((m_end-m_begin)%m_stride == 0); //integral # of strides in range
        return 1+(m_end-m_begin)/m_stride;
      }

    protected:
      std::size_t m_begin;
      std::size_t m_end;
      std::size_t m_stride;
    private:
      Range():m_begin(0),m_end(0),m_stride(0){}//constructor for Range::all()

  };

  template <std::size_t D, typename Index1dType = IndexConvert<D> >
  class AccessScheme 
  {
    public:
      typedef carray<int, D> internal_tuple_type;
      typedef carray<int, D> IndexTuple;
      typedef AccessScheme<D,Index1dType> access_type;
      typedef cctbx::nested_loop<internal_tuple_type> loop_type;

      AccessScheme() {};
      AccessScheme(const internal_tuple_type& N) {
        std::copy(N.begin(), N.end(), extent.begin());
        for (int i = 0, istride=1; i<D; ++i) { 
          base[i] = 0; stride[D-1-i] = istride; istride*=extent[D-1-i];}
      }
      AccessScheme(std::size_t n0) {
        extent = internal_tuple_type(n0);
        for (int i = 0, istride=1; i<D; ++i) { 
          base[i] = 0; stride[D-1-i] = istride; istride*=extent[D-1-i];}
      }
      AccessScheme(std::size_t n0, std::size_t n1) {
        extent = internal_tuple_type(n0,n1);
        for (int i = 0, istride=1; i<D; ++i) { 
          base[i] = 0; stride[D-1-i] = istride; istride*=extent[D-1-i];}
      }
      AccessScheme(std::size_t n0, std::size_t n1, std::size_t n2) {
        extent = internal_tuple_type(n0,n1,n2);
        for (int i = 0, istride=1; i<D; ++i) { 
          base[i] = 0; stride[D-1-i] = istride; istride*=extent[D-1-i];}
      }

      std::size_t size1d() const { return cctbx::vector::product(extent); }
      std::size_t length(const std::size_t& i) const 
        {return extent[i]; }
      const internal_tuple_type& shape() const {return extent;}
  
      //Convenience functions for slicing
      //The operator functions return a new AccessScheme, while slice() is a
      // helper function that slices a given dimension of an AccessScheme
      access_type operator()(Range r){
        access_type result(*this); //copy constructor
        result.slice(0,r);
        return result;
      }
      access_type operator()(Range r0, Range r1) {
        access_type result(*this); //copy constructor
        result.slice(0,r0);
        result.slice(1,r1);
        return result;
      }
      access_type operator()(Range r0, Range r1, Range r2){
        access_type result(*this); //copy constructor
        result.slice(0,r0);
        result.slice(1,r1);
        result.slice(2,r2);
        return result;
      }
      void slice(std::size_t slice, Range r){
        if (r.isAll()) return;
        base[slice]+=r.begin()*stride[slice];
        stride[slice]*=r.stride();
        extent[slice]=r.extent();
      }

      bool isContiguous() {
        for (int i = 0, istride=1; i<D; ++i) {
          if (base[i]!=0) return false;
          if (stride[D-1-i]!=istride) return false;
          istride*=extent[D-1-i];
        }
        return true;
      }

      template <typename IndexTupleType>
      std::size_t operator()(const IndexTupleType& I) const {
        return Index1dType()(stride, base, I);
      }

      template <typename IndexTupleType>
      bool is_valid_index(const IndexTupleType& I) const {
        if (I.size() != extent.size()) return false;
        for(std::size_t j=0;j<extent.size();j++) {
          std::size_t i = I[j];
          if (i >= extent[j]) return false;
        }
        return true;
      }

      // Iteration over the access_scheme
      loop_type loop() const {
        internal_tuple_type zeroes;
        for (int i = 0; i<D; ++i) { zeroes[i]=0; }
        return loop_type(zeroes, extent);
      }
        
    protected:
      internal_tuple_type extent;
      internal_tuple_type stride;
      internal_tuple_type base;
  };
}

#endif // CCTBX_ACSCHEME_H
