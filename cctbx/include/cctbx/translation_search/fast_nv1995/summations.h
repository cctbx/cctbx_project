/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     Oct 2002: Modified copy of phenix/fast_translation/summations.h (rwgk)
     Dec 2001: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_TRANSLATION_SEARCH_FAST_NV1995_SUMMATIONS_H
#define CCTBX_TRANSLATION_SEARCH_FAST_NV1995_SUMMATIONS_H

// Navaza, J. & Vernoslova, E. (1995). Acta Cryst. A51, 445-449.

#include <cctbx/sgtbx/space_group.h>
#include <cctbx/sgtbx/miller_ops.h>
#include <cctbx/maptbx/utils.h>
#include <scitbx/array_family/versa_plain.h>

namespace cctbx { namespace translation_search { namespace fast_nv1995_detail {

  class hermitian_accessor
  {
    public:
      typedef miller::index<> index_type;
      typedef miller::index<>::value_type index_value_type;

      struct index_1d_flag_conj
      {
        long index_1d;
        bool is_conj;
      };

      hermitian_accessor() {}

      hermitian_accessor(bool anomalous_flag,
                         miller::index<> const& n,
                         bool two_n_minus_one_flag)
      :
        anomalous_flag_(anomalous_flag),
        range_(two_n_minus_one_flag ? miller::index<>(2*n-1) : n)
      {
        if (!anomalous_flag_) range_[2] = n[2];
      }

      std::size_t
      size_1d() const { return range_.product(); }

      index_value_type
      operator[](std::size_t i) const { return range_[i]; }

      bool
      anomalous_flag() const { return anomalous_flag_; }

      index_1d_flag_conj
      operator()(miller::index<> h) const
      {
        miller::index<> i3d;
        index_1d_flag_conj result;
        result.index_1d = -1;
        result.is_conj = false;
        if (!anomalous_flag_) {
          if (h[2] < 0) {
            h = -h;
            result.is_conj = true;
          }
          for(std::size_t i=0;i<2;i++) {
            i3d[i] = maptbx::h_as_ih_exact(h[i], range_[i], false);
          }
          i3d[2] = maptbx::h_as_ih_exact(h[2], range_[2], true);
        }
        else {
          for(std::size_t i=0;i<3;i++) {
            i3d[i] = maptbx::h_as_ih_exact(h[i], range_[i], false);
          }
        }
        if (i3d.min() < 0) return result;
        // Manually optimized for best performance.
        result.index_1d = (i3d[0] * range_[1] + i3d[1]) * range_[2] + i3d[2];
        return result;
      }

    protected:
      bool anomalous_flag_;
      miller::index<> range_;
  };

  // Class for accessing the 3d array of f_calc.
  template <typename FloatType>
  class f_calc_map
  {
    public:
      typedef std::complex<FloatType> complex_type;
      typedef af::versa_plain<complex_type, hermitian_accessor>
                data_array_type;

      f_calc_map() {}

      f_calc_map(bool anomalous_flag,
                 miller::index<> const& abs_range)
      :
        data_(hermitian_accessor(anomalous_flag, abs_range, true))
      {}

      bool
      anomalous_flag() const { return data_.accessor().anomalous_flag(); }

      void import(af::const_ref<miller::index<> > const& miller_indices,
                  af::const_ref<complex_type> const& f_calc)
      {
        CCTBX_ASSERT(miller_indices.size() == f_calc.size());
        for(std::size_t i=0;i<f_calc.size();i++) {
          set(miller_indices[i], f_calc[i]);
        }
      }

      complex_type
      operator()(miller::index<> const& h) const
      {
        hermitian_accessor::index_1d_flag_conj ic = data_.accessor()(h);
        if (ic.index_1d < 0) return complex_type(0);
        if (ic.is_conj) return std::conj(data_[ic.index_1d]);
        return data_[ic.index_1d];
      }

    protected:
      void
      set(miller::index<> const& h, complex_type const& val)
      {
        hermitian_accessor::index_1d_flag_conj ic = data_.accessor()(h);
        CCTBX_ASSERT(ic.index_1d >= 0);
        if (ic.is_conj) data_[ic.index_1d] = std::conj(val);
        else            data_[ic.index_1d] =           val;
        if (anomalous_flag() || h[2] != 0) return;
        ic = data_.accessor()(-h);
        CCTBX_ASSERT(ic.index_1d >= 0);
        if (ic.is_conj) data_[ic.index_1d] =           val;
        else            data_[ic.index_1d] = std::conj(val);
      }

      data_array_type data_;
  };

  // Class for accumulating the results of the summations.
  template <typename FloatType>
  class summation_accumulator
  {
    private:
      typedef std::complex<FloatType> complex_type;
      typedef typename miller::index<>::value_type mi_v_t;

    public:
      summation_accumulator() {}

      summation_accumulator(complex_type* begin, miller::index<> const& n)
      : begin_(begin), accessor_(true, n, false)
      {}

      // Adds cf if -h is in the positive halfspace.
      // Ignores cf otherwise.
      // Manually optimized for best performance.
      void
      minus(miller::index<> const& h, complex_type const& cf)
      {
        if (0 > -h[2] || -h[2] >= accessor_[2]) return;
        mi_v_t h1 = -h[1];
        if (h1 < 0) h1 += accessor_[1];
        mi_v_t h0 = -h[0];
        if (h0 < 0) h0 += accessor_[0];
        begin_[(h0 * accessor_[1] + h1) * accessor_[2] - h[2]] += cf;
      }

      // Adds conj(cf) if +h is in the positive halfspace.
      // Adds      cf  if -h is in the positive halfspace.
      // Manually optimized for best performance.
      void
      plus_minus(miller::index<> const& h, complex_type const& cf)
      {
        if (!(0 > h[2] || h[2] >= accessor_[2])) {
          mi_v_t h1 = h[1];
          if (h1 < 0) h1 += accessor_[1];
          mi_v_t h0 = h[0];
          if (h0 < 0) h0 += accessor_[0];
          begin_[(h0 * accessor_[1] + h1) * accessor_[2] + h[2]]
            += std::conj(cf);
        }
        if (!(0 > -h[2] || -h[2] >= accessor_[2])) {
          mi_v_t h1 = -h[1];
          if (h1 < 0) h1 += accessor_[1];
          mi_v_t h0 = -h[0];
          if (h0 < 0) h0 += accessor_[0];
          begin_[(h0 * accessor_[1] + h1) * accessor_[2] - h[2]]
            += cf;
        }
      }

    protected:
      complex_type* begin_;
      hermitian_accessor accessor_;
  };

  // Precomputes f~(hs) (Navaza & Vernoslova (1995), p.447, following Eq. (14))
  // ftil = f_calc(hr) * exp(2*pi*i*ht)
  template <typename FloatTypeSummation>
  void set_ftilde(sgtbx::space_group const& space_group,
                  f_calc_map<FloatTypeSummation> const& fc_map,
                  miller::index<> const& h,
                  miller::index<>* hs,
                  std::complex<FloatTypeSummation>* fts)
  {
    for(std::size_t i=0;i<space_group.order_z();i++) {
      sgtbx::rt_mx s = space_group(i);
      hs[i] = h * s.r();
      fts[i] = fc_map(hs[i]) * std::polar(1.,
        2 * (h * s.t()) * scitbx::constants::pi / s.t().den());
    }
  }

  // lhs * conj(rhs)
  // With some compilers much faster than lhs * std::conj(rhs).
  template <typename T>
  inline
  std::complex<T>
  mul_conj(std::complex<T> const& lhs, std::complex<T> const& rhs)
  {
    return std::complex<T>(
      lhs.real() * rhs.real() + lhs.imag() * rhs.imag(),
      lhs.imag() * rhs.real() - lhs.real() * rhs.imag());
  }

  // Sum according to Eq. (15) of Navaza & Vernoslova (1995), p. 447.
  // Manually optimized for best performance.
  template <typename FloatTypeSummation>
  void
  summation_eq15(
    sgtbx::space_group const& space_group,
    af::const_ref<miller::index<> > const& miller_indices,
    af::const_ref<std::complex<FloatTypeSummation> > const& f_part,
    f_calc_map<FloatTypeSummation> const& p1_f_calc,
    summation_accumulator<FloatTypeSummation>& sum)
  {
    CCTBX_ASSERT(   f_part.size() == 0
                 || f_part.size() == miller_indices.size());
    typedef miller::index<> mi_t;
    typedef std::complex<FloatTypeSummation> cx_t;
    std::size_t order_z = space_group.order_z();
    bool have_f_part = f_part.size();
    std::vector<mi_t> hs_vector(order_z);
    mi_t* hs = &*hs_vector.begin();
    std::vector<cx_t> fts_vector(order_z);
    cx_t* fts = &*fts_vector.begin();
    cx_t fpi(0);
    cx_t fpi_sq(0);
    cx_t two_fpi_sq_conj_fpi(0);
    cx_t four_conj_fpi_fpi(0);
    cx_t two_fpi(0);
    for (std::size_t ih = 0; ih < miller_indices.size(); ih++) {
      mi_t h = miller_indices[ih];
      double mh = space_group.multiplicity(h, p1_f_calc.anomalous_flag());
      set_ftilde(space_group, p1_f_calc, h, hs, fts);
      if (have_f_part) {
        fpi = f_part[ih];
        fpi_sq = fpi * fpi;
        two_fpi_sq_conj_fpi = 2. * fpi_sq * std::conj(fpi);
        four_conj_fpi_fpi = 4. * std::conj(fpi) * fpi;
        two_fpi = 2. * fpi;
      }
      for (std::size_t is0 = 0; is0 < order_z; is0++) {
        const mi_t& hm0 = hs[is0];
        cx_t mh_ftil0c = mh * std::conj(fts[is0]);
        if (have_f_part) {
          cx_t cf = mh_ftil0c * two_fpi_sq_conj_fpi;
          sum.plus_minus(hm0, cf);
        }
        for (std::size_t is1 = 0; is1 < order_z; is1++) {
          mi_t hm01(hm0 - hs[is1]);
          cx_t mh_ftil0c_ftil1(mh_ftil0c * fts[is1]);
          if (have_f_part) {
            cx_t cf = mh_ftil0c_ftil1 * four_conj_fpi_fpi;
            sum.minus(hm01, cf);
            mi_t hm0p1(hm0 + hs[is1]);
            cf = mul_conj(mh_ftil0c, fts[is1]) * fpi_sq;
            sum.plus_minus(hm0p1, cf);
          }
          for (std::size_t is2 = 0; is2 < order_z; is2++) {
            mi_t hm01p2(hm01 + hs[is2]);
            cx_t mh_ftil0c_ftil1_ftil2c(
              mul_conj(mh_ftil0c_ftil1, fts[is2]));
            if (have_f_part) {
              cx_t cf = mh_ftil0c_ftil1_ftil2c * two_fpi;
              sum.plus_minus(hm01p2, cf);
            }
            for (std::size_t is3 = 0; is3 < order_z; is3++) {
              mi_t hm01p23(hm01p2 - hs[is3]);
              cx_t cf = mh_ftil0c_ftil1_ftil2c * fts[is3];
              sum.minus(hm01p23, cf);
            }
          }
        }
      }
    }
  }

  // Sum according to Eq. (14) of Navaza & Vernoslova (1995), p. 447.
  template <typename FloatTypeIobs,
            typename FloatTypeSummation>
  void
  summation_eq14(
    sgtbx::space_group const& space_group,
    af::const_ref<miller::index<> > const& miller_indices,
    af::const_ref<FloatTypeIobs> const& d_i_obs,
    af::const_ref<std::complex<FloatTypeSummation> > const& f_part,
    f_calc_map<FloatTypeSummation> const& p1_f_calc,
    summation_accumulator<FloatTypeSummation>& sum)
  {
    CCTBX_ASSERT(   d_i_obs.size() == 0
                 || d_i_obs.size() == miller_indices.size());
    CCTBX_ASSERT(   f_part.size() == 0
                 || f_part.size() == miller_indices.size());
    typedef miller::index<> mi_t;
    typedef std::complex<FloatTypeSummation> cx_t;
    std::size_t order_z = space_group.order_z();
    bool have_f_part = f_part.size();
    std::vector<mi_t> hs_vector(order_z);
    mi_t* hs = &*hs_vector.begin();
    std::vector<cx_t> fts_vector(order_z);
    cx_t* fts = &*fts_vector.begin();
    for (std::size_t ih = 0; ih < miller_indices.size(); ih++) {
      mi_t h = miller_indices[ih];
      double mh_di = space_group.multiplicity(h, p1_f_calc.anomalous_flag());
      if (d_i_obs.size()) mh_di *= d_i_obs[ih];
      set_ftilde(space_group, p1_f_calc, h, hs, fts);
      for (std::size_t is0 = 0; is0 < order_z; is0++) {
        const mi_t& hm0 = hs[is0];
        cx_t mh_di_ftil0c = mh_di * std::conj(fts[is0]);
        if (have_f_part) {
          cx_t cf = mh_di_ftil0c * f_part[ih];
          sum.plus_minus(hm0, cf);
        }
        for (std::size_t is1 = 0; is1 < order_z; is1++) {
          mi_t hm01(hm0 - hs[is1]);
          cx_t cf = mh_di_ftil0c * fts[is1];
          sum.minus(hm01, cf);
        }
      }
    }
  }

}}} // namespace cctbx::translation_search::fast_nv1995_detail

#endif // CCTBX_TRANSLATION_SEARCH_FAST_NV1995_SUMMATIONS_H
