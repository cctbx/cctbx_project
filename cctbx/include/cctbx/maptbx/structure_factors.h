/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Copyright (c) 2002 Airlie McCoy.

   Revision history:
     2002 Nov: Modified fragments of cctbx/sftbx/sfmap.h (rwgk)
     2002 May: Created based on phaser/src/MapFFT.cc by Airlie McCoy (rwgk)
 */

#ifndef CCTBX_MAPTBX_STRUCTURE_FACTORS_H
#define CCTBX_MAPTBX_STRUCTURE_FACTORS_H

#include <cctbx/sgtbx/reciprocal_space_asu.h>
#include <cctbx/miller/sym_equiv.h>
#include <cctbx/maptbx/utils.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <scitbx/array_family/versa.h>

namespace cctbx { namespace maptbx { namespace structure_factors {

  template <typename FloatType = double>
  class to_map
  {
    public:
      template <typename OtherFloatType>
      to_map(
        sgtbx::space_group const& space_group,
        bool anomalous_flag,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<std::complex<OtherFloatType> > const& structure_factors,
        af::flex_grid<> const& flex_grid,
        bool conjugate_flag)
      :
        complex_map_(af::c_grid_padded<3>(flex_grid))
      {
        typename af::c_grid_padded<3>::index_type const&
          n_complex = complex_map_.accessor().focus();
        for(std::size_t i=0;i<miller_indices.size();i++) {
          miller::sym_equiv_indices h_eq(space_group, miller_indices[i]);
          for(int e=0;e<h_eq.multiplicity(anomalous_flag);e++) {
            miller::sym_equiv_index h_eq_e = h_eq(e);
            miller::index<> h = h_eq_e.h();
            if (conjugate_flag) h = -h;
            if (!anomalous_flag && h[2] < 0) continue;
            miller::index<> ih = h_as_ih_array(anomalous_flag, h, n_complex);
            CCTBX_ASSERT(ih.const_ref().all_ge(0));
            complex_map_(ih) = h_eq_e.complex_eq(structure_factors[i]);
          }
        }
      }

      af::versa<std::complex<FloatType>, af::c_grid_padded<3> >
      complex_map() const { return complex_map_; }

    protected:
      af::versa<std::complex<FloatType>, af::c_grid_padded<3> > complex_map_;
  };

  template <typename FloatType = double>
  class to_under_sampled_map
  {
    public:
      template <typename OtherFloatType>
      to_under_sampled_map(
        sgtbx::space_group const& space_group,
        bool anomalous_flag,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<std::complex<OtherFloatType> > const& structure_factors,
        af::int3 const& n_real,
        af::flex_grid<> const& map_grid,
        bool conjugate_flag)
      :
        complex_map_(af::c_grid_padded<3>(map_grid))
      {
        CCTBX_ASSERT(map_grid.all().all_ge(map_grid.focus()));
        af::int3 map_grid_focus(complex_map_.accessor().focus());
        std::size_t count_n_real_not_equal_map_grid_focus = 0;
        for(std::size_t i=0;i<3;i++) {
          CCTBX_ASSERT(   map_grid_focus[i] == n_real[i]
                       || map_grid_focus[i] == n_real[i]/2+1);
          if (map_grid_focus[i] != n_real[i]) {
            count_n_real_not_equal_map_grid_focus++;
          }
        }
        CCTBX_ASSERT(count_n_real_not_equal_map_grid_focus <= 1);
        CCTBX_ASSERT(   anomalous_flag == false
                     || count_n_real_not_equal_map_grid_focus == 0);
        for(std::size_t i=0;i<miller_indices.size();i++) {
          miller::sym_equiv_indices h_eq(space_group, miller_indices[i]);
          for(int e=0;e<h_eq.multiplicity(anomalous_flag);e++) {
            miller::sym_equiv_index h_eq_e = h_eq(e);
            miller::index<> h = h_eq_e.h();
            if (conjugate_flag) h = -h;
            af::int3 ih = h_as_ih_array_under_sampled(h, n_real);
            if (!ih.all_lt(map_grid_focus)) continue;
            complex_map_(ih) += h_eq_e.complex_eq(structure_factors[i]);
          }
        }
      }

      af::versa<std::complex<FloatType>, af::c_grid_padded<3> >
      complex_map() const { return complex_map_; }

    protected:
      af::versa<std::complex<FloatType>, af::c_grid_padded<3> > complex_map_;
  };

  template <typename FloatType = double>
  class from_map
  {
    public:
      template <typename OtherFloatType>
      from_map(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group_type const& sg_type,
        bool anomalous_flag,
        double d_min,
        af::const_ref<std::complex<OtherFloatType>,
                      af::c_grid_padded<3> > const& complex_map,
        bool conjugate_flag)
      {
        CCTBX_ASSERT(d_min > 0);
        double d_star_sq_max = 1. / (d_min * d_min);
        typedef typename af::c_grid_padded<3>::index_type index_type;
        index_type n_complex = complex_map.accessor().focus();
        sgtbx::reciprocal_space::asu asu(sg_type);
        sgtbx::space_group const& space_group = sg_type.group();
        miller::index<> h;
        miller::index<> mh;
        uctbx::incremental_d_star_sq<double> incr_d_star_sq(unit_cell);
        index_type ih;
        ih[2] = 1;
        for(ih[0]=0;ih[0]<n_complex[0];ih[0]++) {
          h[0] = ih_as_h(ih[0], n_complex[0]);
          mh[0] = -h[0];
          incr_d_star_sq.update0(h[0]);
        for(ih[1]=0;ih[1]<n_complex[1];ih[1]++,ih[2]=0) {
          h[1] = ih_as_h(ih[1], n_complex[1]);
          mh[1] = -h[1];
          incr_d_star_sq.update1(h[1]);
        for(;ih[2]<n_complex[2];ih[2]++) {
          if (anomalous_flag) {
            h[2] = ih_as_h(ih[2], n_complex[2]);
          }
          else {
            h[2] = ih[2];
          }
          if (incr_d_star_sq.get(h[2]) > d_star_sq_max) continue;
          CCTBX_ASSERT(ih[0]*2 != n_complex[0]);
          CCTBX_ASSERT(ih[1]*2 != n_complex[1]);
          if (anomalous_flag) {
            CCTBX_ASSERT(ih[2]*2 != n_complex[2]);
          }
          mh[2] = -h[2];
          int asu_which = asu.which(h, mh);
          if (asu_which == 0) continue;
          sgtbx::phase_info phase_info(space_group, h, false);
          if (phase_info.is_sys_absent()) continue;
          bool f_conj = false;
          if (!anomalous_flag) {
            if (asu_which > 0) {
              miller_indices_.push_back(h);
              f_conj = conjugate_flag;
            }
            else {
              if (h[2] == 0) continue;
              miller_indices_.push_back(mh);
              f_conj = !conjugate_flag;
            }
          }
          else {
            if (   ((asu_which < 0) != conjugate_flag)
                && phase_info.is_centric()) {
              continue;
            }
            if (conjugate_flag) miller_indices_.push_back(mh);
            else                miller_indices_.push_back(h);
          }
          if (f_conj) data_.push_back(std::conj(complex_map(ih)));
          else        data_.push_back(complex_map(ih));
        }}}
      }

      template <typename OtherFloatType>
      from_map(
        bool anomalous_flag,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<std::complex<OtherFloatType>,
                      af::c_grid_padded<3> > const& complex_map,
        bool conjugate_flag)
      {
        typename af::c_grid_padded<3>::index_type const&
          n_complex = complex_map.accessor().focus();
        data_.reserve(miller_indices.size());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          miller::index<> h = miller_indices[i];
          bool f_conj = conjugate_flag;
          if (!anomalous_flag) {
            if (h[2] < 0) {
              h = -h;
              f_conj = !f_conj;
            }
          }
          else {
            if (f_conj) {
              h = -h;
              f_conj = false;
            }
          }
          miller::index<> ih = h_as_ih_array(anomalous_flag, h, n_complex);
          CCTBX_ASSERT(ih.const_ref().all_ge(0));
          if (!f_conj) data_.push_back(complex_map(ih));
          else         data_.push_back(std::conj(complex_map(ih)));
        }
      }

      af::shared<miller::index<> >
      miller_indices() const { return miller_indices_; }

      af::shared<std::complex<FloatType> >
      data() const { return data_; }

    protected:
      af::shared<miller::index<> > miller_indices_;
      af::shared<std::complex<FloatType> > data_;
  };
}}} // namespace cctbx::maptbx::structure_factors

#endif // CCTBX_MAPTBX_STRUCTURE_FACTORS_H
