#ifndef CCTBX_MAPTBX_STRUCTURE_FACTORS_H
#define CCTBX_MAPTBX_STRUCTURE_FACTORS_H

#include <cctbx/sgtbx/reciprocal_space_asu.h>
#include <cctbx/miller/sym_equiv.h>
#include <cctbx/maptbx/utils.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <scitbx/array_family/versa.h>

namespace cctbx { namespace maptbx { namespace structure_factors {

  //! Copies a structure factor array to a 3-dimensional complex map.
  /*! Reference: David A. Langs (2002), J. Appl. Cryst. 35, 505.
   */
  template <typename FloatType=double>
  class to_map
  {
    public:
      to_map() {}

      template <typename OtherFloatType>
      to_map(
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
            af::int3 ih = h_as_ih_mod_array(h, n_real);
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

  //! Extracts structure factors from a 3-dimensional complex map.
  template <typename FloatType=double>
  class from_map
  {
    public:
      from_map() {}

      //! Extracts Miller indices and data up to a given resolution.
      template <typename OtherFloatType>
      from_map(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group_type const& sg_type,
        bool anomalous_flag,
        double d_min,
        af::const_ref<std::complex<OtherFloatType>,
                      af::c_grid_padded<3> > const& complex_map,
        bool conjugate_flag,
        bool discard_indices_affected_by_aliasing=false)
      :
        n_indices_affected_by_aliasing_(0)
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
          if (discard_indices_affected_by_aliasing) {
            bool discard = false;
            if      (ih[0]*2 == n_complex[0]) discard = true;
            else if (ih[1]*2 == n_complex[1]) discard = true;
            else if (anomalous_flag) {
              if (ih[2]*2 == n_complex[2]) discard = true;
            }
            if (discard) {
              n_indices_affected_by_aliasing_++;
              continue;
            }
          }
          else {
            CCTBX_ASSERT(ih[0]*2 != n_complex[0]);
            CCTBX_ASSERT(ih[1]*2 != n_complex[1]);
            if (anomalous_flag) {
              CCTBX_ASSERT(ih[2]*2 != n_complex[2]);
            }
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

      //! Extracts data for given Miller indices.
      template <typename OtherFloatType>
      from_map(
        bool anomalous_flag,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<std::complex<OtherFloatType>,
                      af::c_grid_padded<3> > const& complex_map,
        bool conjugate_flag,
        bool allow_miller_indices_outside_map=false)
      :
        n_indices_affected_by_aliasing_(0)
      {
        af::int3 map_grid_focus = complex_map.accessor().focus();
        data_.reserve(miller_indices.size());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          array_access aa(
            anomalous_flag, map_grid_focus, conjugate_flag, miller_indices[i]);
          if (aa.ih.all_ge(0)) {
            if (!aa.f_conj) data_.push_back(complex_map(aa.ih));
            else            data_.push_back(std::conj(complex_map(aa.ih)));
          }
          else if (allow_miller_indices_outside_map) {
            outside_map_.push_back(data_.size());
            data_.push_back(0);
          }
          else {
            throw_error_not_in_map();
          }
        }
      }

      //! Symmetry sum of data for given Miller indices.
      template <typename OtherFloatType>
      from_map(
        sgtbx::space_group const& space_group,
        bool anomalous_flag,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<std::complex<OtherFloatType>,
                      af::c_grid_padded<3> > const& complex_map,
        bool conjugate_flag)
      :
        n_indices_affected_by_aliasing_(0)
      {
        typedef FloatType f_t;
        typedef std::complex<FloatType> c_t;
        af::int3 map_grid_focus = complex_map.accessor().focus();
        data_.reserve(miller_indices.size());
        bool sum_bijvoet_pairs = (anomalous_flag && space_group.is_centric());
        f_t two_pi_t_den = scitbx::constants::two_pi / space_group.t_den();
        c_t shift_inv_t(0,0);
        for(std::size_t i=0;i<miller_indices.size();i++) {
          miller::index<> h = miller_indices[i];
          if (space_group.is_centric()) {
            f_t phi = two_pi_t_den * f_t(h * space_group.inv_t());
            shift_inv_t = c_t(std::cos(phi), std::sin(phi));
          }
          c_t f(0,0);
          for(std::size_t i_smx=0;i_smx<space_group.n_smx();i_smx++) {
            sgtbx::rt_mx const& s = space_group.smx(i_smx);
            miller::index<> hr = h * s.r();
            array_access aa(
              anomalous_flag, map_grid_focus, conjugate_flag, hr);
            if (!aa.ih.all_ge(0)) throw_error_not_in_map();
            f_t phi = two_pi_t_den * f_t(h * s.t());
            c_t shift(std::cos(phi), std::sin(phi));
            if (!aa.f_conj) f += complex_map(aa.ih) * shift;
            else            f += std::conj(complex_map(aa.ih)) * shift;
            if (sum_bijvoet_pairs) {
              aa = array_access(
                anomalous_flag, map_grid_focus, conjugate_flag, -hr);
              if (!aa.ih.all_ge(0)) throw_error_not_in_map();
              CCTBX_ASSERT(!aa.f_conj);
              f += complex_map(aa.ih) * std::conj(shift) * shift_inv_t;
            }
          }
          if (!anomalous_flag && space_group.is_centric()) {
            f += std::conj(f) * shift_inv_t;
          }
          f *= space_group.n_ltr();
          data_.push_back(f);
        }
      }

      af::shared<miller::index<> >
      miller_indices() const { return miller_indices_; }

      af::shared<std::complex<FloatType> >
      data() const { return data_; }

      std::size_t
      n_indices_affected_by_aliasing() const
      {
        return n_indices_affected_by_aliasing_;
      }

      af::shared<std::size_t>
      outside_map() const { return outside_map_; }

    protected:
      af::shared<miller::index<> > miller_indices_;
      af::shared<std::complex<FloatType> > data_;
      std::size_t n_indices_affected_by_aliasing_;
      af::shared<std::size_t> outside_map_;

      static
      void
      throw_error_not_in_map()
      {
        throw error("Miller index not in structure factor map.");
      }

      struct array_access
      {
        array_access(
          bool anomalous_flag,
          af::int3 const& map_grid_focus,
          bool conjugate_flag,
          miller::index<> h)
        {
          f_conj = conjugate_flag;
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
          ih = h_as_ih_exact_array(anomalous_flag, h, map_grid_focus);
        }

        bool f_conj;
        af::int3 ih;
      };
  };

}}} // namespace cctbx::maptbx::structure_factors

#endif // CCTBX_MAPTBX_STRUCTURE_FACTORS_H
