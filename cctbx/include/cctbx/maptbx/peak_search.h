/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Jan: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MAPTBX_PEAK_SEARCH_H
#define CCTBX_MAPTBX_PEAK_SEARCH_H

#include <scitbx/indexed_value.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <scitbx/array_family/loops.h>
#include <scitbx/array_family/sort.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <cctbx/error.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace maptbx {

  //! Peak search in a 3d-map covering the unit cell.
  /*! Dependencies due to symmetry (space group, Euclidean normalizer,
      or similar) are encoded in a 3d-map of tags.
      <p>
      Requirements:
      <ul>
        <li>physical dimensions of maps are equal to generic dimensions
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
  template <typename DataType,
            typename TagType>
  void
  peak_search_unit_cell(
    af::const_ref<DataType, af::c_grid_padded<3> > const& data,
    af::ref<TagType, af::c_grid<3> > const& tags,
    int level)
  {
    CCTBX_ASSERT(tags.accessor().all_eq(data.accessor().focus()));
    CCTBX_ASSERT(!data.accessor().is_padded()); // XXX
    const DataType* pdata = data.begin();
    TagType* ptags = tags.begin();
    int nx = tags.accessor()[0];
    int ny = tags.accessor()[1];
    int nz = tags.accessor()[2];
    int       nk =           nz;
    int    nj_nk =      ny * nz;
    int ni_nj_nk = nx * ny * nz;

    int im, jm, km;
    int i0, j0, k0;
    int ip, jp, kp;
    int ibreak, jbreak, kbreak;

    // reset tags for independent grid points
    for (int i = 0; i < ni_nj_nk; i++) {
      if (ptags[i] < 0) ptags[i] = -1;
    }

    const DataType* pd = pdata;
    TagType* pf = ptags;
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
          TagType* indep_pf;
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

  template <typename ValueType,
            typename CountType = long>
  struct peak_histogram
  {
    peak_histogram() {}

    template <typename DataType,
              typename TagType>
    peak_histogram(
      af::const_ref<DataType, af::c_grid_padded<3> > const& data,
      af::ref<TagType, af::c_grid<3> > const& tags,
      std::size_t n_slots=1000)
    :
      slots(n_slots)
    {
      CCTBX_ASSERT(data.accessor().focus().all_eq(tags.accessor()));
      CCTBX_ASSERT(n_slots > 0);
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
      slots.fill(0);
      for(i=0;i<data.size();i++) {
        if (tags[i] != -2) continue;
        ValueType d = data[i] - data_min;
        std::size_t i_slot = 0;
        if (d != 0 && d >= slot_width) {
              i_slot = std::size_t(d / slot_width);
          if (i_slot >= slots.size()) i_slot = slots.size() - 1;
        }
        slots[i_slot]++;
      }
    }

    ValueType
    get_cutoff(CountType const& max_peaks,
               ValueType const& tolerance=1.e-4) const
    {
      CountType cum = 0;
      std::size_t i = slots.size();
      for (; i; i--) {
        cum += slots[i-1];
        if (cum > max_peaks) break;
      }
      return data_min + i * slot_width + slot_width * tolerance;
    }

    ValueType data_min;
    ValueType data_max;
    ValueType slot_width;
    af::shared<CountType> slots;
  };

  template <typename GridIndexType = af::tiny<long, 3>,
            typename SiteType = scitbx::vec3<double>,
            typename ValueType = double>
  class peak_list
  {
    public:
      peak_list() {}

      template <typename DataType,
                typename TagType>
      peak_list(
        af::const_ref<DataType, af::c_grid_padded<3> > const& data,
        af::ref<TagType, af::c_grid<3> > const& tags,
        int peak_search_level=1,
        std::size_t max_peaks=0,
        bool interpolate=true)
      :
        gridding_(data.accessor().focus())
      {
        peak_search_unit_cell(data, tags, peak_search_level);
        ValueType peak_cutoff = 0;
        bool use_cutoff = false;
        if (max_peaks) {
          peak_histogram<ValueType> hist(data, tags);
          peak_cutoff = hist.get_cutoff(max_peaks);
          use_cutoff = true;
        }
        process_peaks(data, tags, peak_cutoff, use_cutoff, interpolate);
      }

      template <typename DataType,
                typename TagType>
      peak_list(
        af::const_ref<DataType, af::c_grid_padded<3> > const& data,
        af::ref<TagType, af::c_grid<3> > const& tags,
        int peak_search_level,
        ValueType peak_cutoff,
        std::size_t max_peaks,
        bool interpolate=true)
      :
        gridding_(data.accessor().focus())
      {
        peak_search_unit_cell(data, tags, peak_search_level);
        if (max_peaks) {
          peak_histogram<ValueType> hist(data, tags);
          peak_cutoff = std::max(peak_cutoff, hist.get_cutoff(max_peaks));
        }
        process_peaks(data, tags, peak_cutoff, true, interpolate);
      }

      GridIndexType const&
      gridding() const { return gridding_; }

      std::size_t
      size() const
      {
        CCTBX_ASSERT(grid_heights().size() == grid_indices().size());
        CCTBX_ASSERT(sites().size() == grid_indices().size());
        CCTBX_ASSERT(heights().size() == grid_indices().size());
        return grid_indices().size();
      }

      /*! Python: grid_indices(i)
       */
      af::shared<GridIndexType>
      grid_indices() const { return grid_indices_; }

      af::shared<ValueType>
      grid_heights() const { return grid_heights_; }

      af::shared<SiteType>
      sites() const { return sites_; }

      af::shared<ValueType>
      heights() const { return heights_; }

    protected:
      GridIndexType gridding_;
      af::shared<GridIndexType> grid_indices_;
      af::shared<ValueType> grid_heights_;
      af::shared<SiteType> sites_;
      af::shared<ValueType> heights_;

      template <typename DataType,
                typename TagType>
      void
      process_peaks(
        af::const_ref<DataType, af::c_grid_padded<3> > const& data,
        af::const_ref<TagType, af::c_grid<3> > const& tags,
        ValueType const& cutoff,
        bool use_cutoff,
        bool interpolate)
      {
        collect_peaks(data, tags, cutoff, use_cutoff);
        if (interpolate) interpolate_sites_and_heights();
        else copy_sites_and_heights();
        sort();
      }

      template <typename DataType,
                typename TagType>
      void
      collect_peaks(
        af::const_ref<DataType, af::c_grid_padded<3> > const& data,
        af::const_ref<TagType, af::c_grid<3> > const& tags,
        ValueType const& cutoff,
        bool use_cutoff)
      {
        typedef typename af::c_grid<3>::index_type index_type;
        af::nested_loop<index_type> loop(data.accessor().focus());
        for (index_type const& pivot = loop(); !loop.over(); loop.incr()) {
          if (tags(pivot) == -2 && (!use_cutoff || data(pivot) >= cutoff)) {
            grid_indices_.push_back(pivot);
            grid_heights_.push_back(data(pivot));
          }
        }
      }

      void
      copy_sites_and_heights()
      {
        af::const_ref<GridIndexType> gi = grid_indices_.const_ref();
        af::tiny<ValueType, 3> gr(gridding());
        sites_.reserve(gi.size());
        for(std::size_t i=0;i<gi.size();i++) {
          sites_.push_back(af::tiny<ValueType, 3>(gi[i]) / gr);
        }
        heights_.assign(grid_heights_);
      }

      void
      interpolate_sites_and_heights()
      {
        throw CCTBX_NOT_IMPLEMENTED();
      }

      void
      sort()
      {
        af::shared<std::size_t>
          perm = af::sort_permutation(heights_.const_ref(), true);
        af::const_ref<std::size_t> p = perm.const_ref();
        grid_indices_ = af::shuffle(grid_indices_.const_ref(), p);
        grid_heights_ = af::shuffle(grid_heights_.const_ref(), p);
        sites_ = af::shuffle(sites_.const_ref(), p);
        heights_ = af::shuffle(heights_.const_ref(),p);
      }
  };

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_PEAK_SEARCH_H
