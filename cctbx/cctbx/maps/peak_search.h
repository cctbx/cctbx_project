// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MAPS_PEAK_SEARCH_H
#define CCTBX_MAPS_PEAK_SEARCH_H

#include <algorithm>
#include <functional>
#include <cctbx/shared_storage.h>

namespace cctbx { namespace maps {

  //! Peak search in a 3d-map covering the unit cell.
  /*! Dependencies due to symmetry (space group, Euclidean normalizer,
      or similar) are encoded in a 3d-map of tags.
      <p>
      Requirements:
      <ul>
        <li>physical dimensions of %maps are equal to generic dimensions
        <li>dimensions of data and tags are equal
      </ul>
      <pre>
      tags:
        on input: < 0: independent grid point
                  >= 0: dependent grid point,
                        flag is 1d-index of corresponding
                        independent grid point
        on output: -1: independent, not a peak
                   -2: independent, peak
                   tags for dependent grid points are unchanged

      level = 1: compare to the 6 nearest neighbors
            = 2: also compare to the 12 second-nearest neighbors
            > 2: also compare to the 8 third-nearest neighbors</pre>
   */
  template <typename DataVecRefNdType,
            typename TagsVecRefNdType>
  void
  peak_search_unit_cell(const DataVecRefNdType& data,
                        TagsVecRefNdType& tags,
                        int level)
  {
    typedef typename DataVecRefNdType::value_type data_value_type;
    typedef typename TagsVecRefNdType::value_type tags_value_type;

    const data_value_type* pdata = data.begin();
    tags_value_type* ptags = tags.begin();
    int Nx = data.dim()[0];
    int Ny = data.dim()[1];
    int Nz = data.dim()[2];
    int       nk =           Nz;
    int    nj_nk =      Ny * Nz;
    int ni_nj_nk = Nx * Ny * Nz;

    int im, jm, km;
    int i0, j0, k0;
    int ip, jp, kp;
    int ibreak, jbreak, kbreak;

    // reset tags for independent grid points
    for (int i = 0; i < ni_nj_nk; i++) {
      if (ptags[i] < 0) ptags[i] = -1;
    }

    const data_value_type* pd = pdata;
    tags_value_type* pf = ptags;
    im = ni_nj_nk - nj_nk;
    i0 = 0;
    ip = nj_nk;
    ibreak = ni_nj_nk;
    while (ip < ibreak) {
      jm = nj_nk - nk;
      j0 = 0;
      jp = nk;
      jbreak = nj_nk;
      while (jp < jbreak) {
        km = nk - 1;
        k0 = 0;
        kp = 1;
        kbreak = nk;
        while (kp < kbreak) {
          tags_value_type* indep_pf;
          if (*pf < 0) indep_pf = pf;
          else         indep_pf = &ptags[*pf];
          while (! (*indep_pf < -1)) {
            if (level >= 1) {
              /* m00 0m0 00m
                 p00 0p0 00p
               */
              if (*pd < pdata[im + j0 + k0]) break;
              if (*pd < pdata[ip + j0 + k0]) break;
              if (*pd < pdata[i0 + jm + k0]) break;
              if (*pd < pdata[i0 + jp + k0]) break;
              if (*pd < pdata[i0 + j0 + km]) break;
              if (*pd < pdata[i0 + j0 + kp]) break;
              if (level >= 2) {
                /* mm0 m0m 0mm mp0 m0p 0mp
                   pp0 p0p 0pp pm0 p0m 0pm
                 */
                if (*pd < pdata[im + jm + k0]) break;
                if (*pd < pdata[ip + jp + k0]) break;
                if (*pd < pdata[im + j0 + km]) break;
                if (*pd < pdata[ip + j0 + kp]) break;
                if (*pd < pdata[i0 + jm + km]) break;
                if (*pd < pdata[i0 + jp + kp]) break;
                if (*pd < pdata[im + jp + k0]) break;
                if (*pd < pdata[ip + jm + k0]) break;
                if (*pd < pdata[im + j0 + kp]) break;
                if (*pd < pdata[ip + j0 + km]) break;
                if (*pd < pdata[i0 + jm + kp]) break;
                if (*pd < pdata[i0 + jp + km]) break;
                if (level >= 3) {
                  /* mmm mmp mpm mpp
                     ppp ppm pmp pmm
                   */
                  if (*pd < pdata[im + jm + km]) break;
                  if (*pd < pdata[ip + jp + kp]) break;
                  if (*pd < pdata[im + jm + kp]) break;
                  if (*pd < pdata[ip + jp + km]) break;
                  if (*pd < pdata[im + jp + km]) break;
                  if (*pd < pdata[ip + jm + kp]) break;
                  if (*pd < pdata[im + jp + kp]) break;
                  if (*pd < pdata[ip + jm + km]) break;
                }
              }
            }
            *indep_pf = -2;
            break;
          }
          pd++;
          pf++;
          km = k0;
          k0 = kp;
          kp++;
          if (kp == nk) { kp = 0; kbreak = 1; }
        }
        jm = j0;
        j0 = jp;
        jp += nk;
        if (jp == nj_nk) { jp = 0; jbreak = nk; }
      }
      im = i0;
      i0 = ip;
      ip += nj_nk;
      if (ip == ni_nj_nk) { ip = 0; ibreak = nj_nk; }
    }
  }

  template <typename IndexType,
            typename ValueType,
            typename SortCmpFunctor = std::less<ValueType> >
  struct indexed_value
  {
    typedef IndexType index_type;
    typedef ValueType value_type;

    indexed_value() {}
    indexed_value(const index_type& i, const value_type& v)
      : index(i), value(v)
    {}

    bool
    operator<(
    const indexed_value<IndexType, ValueType, SortCmpFunctor>& rhs) const {
      return SortCmpFunctor()(this->value, rhs.value);
    }

    index_type index;
    value_type value;
  };

  template <typename DataVecRefNdType,
            typename TagsVecRefNdType>
  shared_storage<
    indexed_value<
      typename DataVecRefNdType::dimension_type::index_tuple_type,
      typename DataVecRefNdType::value_type,
      std::greater<typename DataVecRefNdType::value_type> >
  >
  collect_peaks(const DataVecRefNdType& data,
                const TagsVecRefNdType& tags,
                const typename DataVecRefNdType::value_type& cutoff,
                bool use_cutoff = true)
  {
    typedef typename
    DataVecRefNdType::dimension_type::index_tuple_type index_tuple_type;
    typedef
      indexed_value<
        index_tuple_type,
        typename DataVecRefNdType::value_type,
        std::greater<typename DataVecRefNdType::value_type>
      >
      iv_type;
    std::vector<iv_type> result;
    nested_loop<index_tuple_type> loop(data.dim());
    for (const index_tuple_type& pivot = loop(); !loop.over(); loop.incr()) {
      if (tags(pivot) == -2 && (!use_cutoff || data(pivot) >= cutoff)) {
        result.push_back(iv_type(pivot, data(pivot)));
      }
    }
    return result;
  }

  template <typename ValueType,
            typename CountType = long>
  struct peak_histogram
  {
    typedef ValueType value_type;
    typedef CountType count_type;

    peak_histogram() {}

    template <typename DataVecRefType,
              typename TagsVecRefType>
    peak_histogram(const DataVecRefType& data,
                   const TagsVecRefType& tags,
                   std::size_t n_slots = 1000)
      : slots(n_slots)
    {
      std::size_t i;
      for(i=0;i<data.size();i++) {
        if (tags[i] == -2) break;
      }
      if (i == data.size()) {
        data_min = 0;
        data_max = 0;
      }
      else {
        data_min = data[i];
        data_max = data[i];
        for(i++;i<data.size();i++) {
          if (tags[i] != -2) continue;
          if (data_min > data[i]) data_min = data[i];
          if (data_max < data[i]) data_max = data[i];
        }
      }
      slot_width = (data_max - data_min) / slots.size();
      std::fill(slots.begin(), slots.end(), 0);
      for(i=0;i<data.size();i++) {
        if (tags[i] != -2) continue;
        value_type d = data[i] - data_min;
        std::size_t i_slot = 0;
        if (d != 0 && d >= slot_width) {
              i_slot = std::size_t(d / slot_width);
          if (i_slot >= slots.size()) i_slot = slots.size() - 1;
        }
        slots[i_slot]++;
      }
    }

    value_type get_cutoff(const count_type& max_peaks,
                          const value_type& tolerance = 1.e-4) const
    {
      count_type cum = 0;
      std::size_t i = slots.size();
      for (; i; i--) {
        cum += slots[i-1];
        if (cum > max_peaks) break;
      }
      return data_min + i * slot_width + slot_width * tolerance;
    }

    value_type data_min;
    value_type data_max;
    value_type slot_width;
    shared_storage<count_type> slots;
  };


  template <typename ValueType,
            typename IndexTupleType = carray<int, 3> >
  class peak_list
    : public shared_storage<
        indexed_value<IndexTupleType, ValueType, std::greater<ValueType> > >
  {
    public:
      typedef ValueType value_type;
      typedef IndexTupleType index_tuple_type;
      typedef
        indexed_value<index_tuple_type, value_type, std::greater<value_type> >
        entry_type;

      peak_list() {}

      template <typename DataVecRefType,
                typename TagsVecRefType>
      peak_list(const DataVecRefType& data,
                TagsVecRefType& tags,
                int peak_search_level = 1,
                std::size_t max_peaks = 0)
      {
        peak_search_unit_cell(data, tags, peak_search_level);
        value_type peak_cutoff = 0;
        bool use_cutoff = false;
        if (max_peaks) {
          peak_histogram<value_type> hist(data, tags);
          peak_cutoff = hist.get_cutoff(max_peaks);
          use_cutoff = true;
        }
        *(static_cast<shared_storage<entry_type>*>(this)) = collect_peaks(
          data, tags, peak_cutoff, use_cutoff);
        std::sort(this->begin(), this->end());
      }
  };

}} // namespace cctbx::maps

#endif // CCTBX_MAPS_PEAK_SEARCH_H
