// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

/*! \file
    Toolbox for Miller indices.
 */

#ifndef CCTBX_MILLER_H
#define CCTBX_MILLER_H

#if defined(__GNUC__) && __GNUC__ < 3
# include <iostream>
#else
# include <ostream>
#endif
#include <vector>
#include <map>
#include <algorithm>
#include <cctbx/fixes/cstdlib>
#include <cctbx/error.h>
#include <cctbx/coordinates.h>
#include <cctbx/indexed_value.h>
#include <cctbx/array_family/tiny_types.h>
#include <cctbx/array_family/shared.h>
#include <cctbx/math/utils.h>

namespace cctbx {
  //! %Miller index namespace.
  namespace Miller {

    //! Enumeration for symbolic subscripting (e.g. MillerIndex[H]).
    enum {H, K, L};

    //! Miller index class.
    class Index : public af::int3 {
      public:
        //! @name Constructors.
        //@{
        Index() {
          for(std::size_t i=0;i<3;i++) elems[i] = 0;
        }
        Index(const af::int3& v) {
          for(std::size_t i=0;i<3;i++) elems[i] = v[i];
        }
        explicit Index(const int* hkl) {
          for(std::size_t i=0;i<3;i++) elems[i] = hkl[i];
        }
        Index(int h, int k, int l) {
          elems[0] = h; elems[1] = k; elems[2] = l;
        }
        //@}

        //! @name Convenience methods.
        //@{
        bool is000() const {
          return !(elems[0] || elems[1] || elems[2]);
        }
        Index operator-() const {
          return Index(-elems[0], -elems[1], -elems[2]);
        }
        //@}

        //! @name Definition of sort order for human-readable listings.
        //@{
        /*! This comparison is computationally more expensive than
            the Miller::hashCompare below.
         */
        bool operator<(const Index& m2) const
        {
          const int P[3] = {2, 0, 1};
          std::size_t i;
          for(i=0;i<3;i++) {
            if (elems[P[i]] >= 0 && m2[P[i]] <  0) return true;
            if (elems[P[i]] <  0 && m2[P[i]] >= 0) return false;
          }
          for(i=0;i<3;i++) {
            if (math::abs(elems[P[i]]) < math::abs(m2[P[i]])) return true;
            if (math::abs(elems[P[i]]) > math::abs(m2[P[i]])) return false;
          }
          return false;
        }
        bool operator>(const Index& m2) const {
          return !(*this < m2);
        }
        //@}
    };

    //! Multiplication of Miller indices and fractional coordiantes.
    template <class FloatType>
    inline FloatType
    operator*(const Index& lhs, const fractional<FloatType>& rhs) {
      FloatType result = 0.;
      for(std::size_t i=0;i<3;i++) result += lhs[i] * rhs[i];
      return result;
    }

    /*! \brief Definition of fast comparison for use in,
        e.g., std::map<Miller::Index>.
     */
    class hashCompare {
      public:
        //! This fast comparison function is implemented as operator().
        bool operator()(const Index& m1,const Index& m2) const {
          for(std::size_t i=0;i<3;i++) {
            if (m1[i] < m2[i]) return true;
            if (m1[i] > m2[i]) return false;
          }
          return false;
        }
    };

    //! ostream output operator for class Miller::Index.
    inline std::ostream& operator<<(std::ostream& os, const Miller::Index& MIx)
    {
      os << "H=" << MIx[0] << " K=" << MIx[1] << " L=" << MIx[2];
      return os;
    }

    /*! \brief Determine min(H[i]) and max(H[i])+1, i=1..3,
         for an array of Miller indices.
     */
    struct index_span : af::tiny<af::tiny<int, 2>, 3>
    {
      typedef af::tiny<int, 2> min_end;
      typedef af::tiny<min_end, 3> base_class;

      index_span() {}

      index_span(af::shared<Index> index_array)
      {
        this->fill(min_end(0, 0));
        if (index_array.size()) {
          for(std::size_t j=0;j<3;j++) {
            (*this)[j] = min_end().fill(index_array[0][j]);
          }
        }
        for(std::size_t i=1;i<index_array.size();i++) {
          for(std::size_t j=0;j<3;j++) {
            math::update_min((*this)[j][0], index_array[i][j]);
            math::update_max((*this)[j][1], index_array[i][j]);
          }
        }
        for(std::size_t j=0;j<3;j++) (*this)[j][1]++;
      }

      af::int3 min() const
      {
        af::int3 result;
        for(std::size_t j=0;j<3;j++) result[j] = (*this)[j][0];
        return result;
      }

      af::int3 max() const
      {
        af::int3 result;
        for(std::size_t j=0;j<3;j++) result[j] = (*this)[j][1] - 1;
        return result;
      }

      af::int3 abs_range() const
      {
        af::int3 result;
        std::size_t j;
        for(j=0;j<3;j++) {
          result[j] = math::abs((*this)[j][0]);
          math::update_max(result[j], math::abs((*this)[j][1]-1));
        }
        for(j=0;j<3;j++) result[j] += 1;
        return result;
      }

      af::int3 map_grid() const
      {
        af::int3 result = abs_range();
        for(std::size_t j=0;j<3;j++) {
          result[j] = (result[j] - 1) * 2 + 1;
        }
        return result;
      }

      bool is_in_domain(Index const& h) const
      {
        for(std::size_t j=0;j<3;j++) {
          if (!((*this)[j][0] <= h[j] && h[j] < (*this)[j][1])) return false;
        }
        return true;
      }

      std::size_t pack(Index const& h) const
      {
        return ((h[0] - (*this)[0][0]) * range_((*this)[1])
              + (h[1] - (*this)[1][0])) * range_((*this)[2])
              + (h[2] - (*this)[2][0]);
      }

      private:
        static int range_(min_end const& span)
        {
          return span[1] - span[0];
        }
    };

    class join_sets
    {
      public:
        join_sets() {}

        join_sets(const af::shared<Index>& a1,
                  const af::shared<Index>& a2)
        {
          typedef std::map<Index, std::size_t> lookup_map_type;
          lookup_map_type lookup_map;
          std::size_t i;
          for(i=0;i<a2.size();i++) lookup_map[a2[i]] = i;
          std::vector<bool> a2_flags(a2.size(), false);
          for(i=0;i<a1.size();i++) {
            lookup_map_type::const_iterator l = lookup_map.find(a1[i]);
            if (l == lookup_map.end()) {
              singles_[0].push_back(i);
            }
            else {
              pairs_.push_back(af::tiny<std::size_t, 2>(i, l->second));
              a2_flags[l->second] = true;
            }
          }
          for(i=0;i<a2.size();i++) {
            if (!a2_flags[i]) singles_[1].push_back(i);
          }
        }

        af::shared<af::tiny<std::size_t, 2> > pairs() const {
          return pairs_;
        }

        af::shared<std::size_t> singles(std::size_t i) const {
          if (i) return singles_[1];
          return singles_[0];
        }
      protected:
        af::shared<af::tiny<std::size_t, 2> > pairs_;
        af::shared<std::size_t> singles_[2];
    };

    template <typename DataType,
              typename SortCmpFunctor>
    void
    inplace_sort(
      af::shared<Index> miller_indices,
      af::shared<DataType> data,
      SortCmpFunctor sort_op)
    {
      cctbx_assert(miller_indices.size() == data.size());
      typedef indexed_value<Index, DataType, SortCmpFunctor> ivalue_type;
      af::shared<ivalue_type> ivalues;
      ivalues.reserve(miller_indices.size());
      for(std::size_t i=0;i<miller_indices.size();i++) {
        ivalues.push_back(ivalue_type(miller_indices[i], data[i]));
      }
      std::sort(ivalues.begin(), ivalues.end());
      for(std::size_t i=0;i<miller_indices.size();i++) {
        miller_indices[i] = ivalues[i].index;
        data[i] = ivalues[i].value;
      }
    }

  } // namespace Miller
} // namespace cctbx

#endif // CCTBX_MILLER_H
