#ifndef CCTBX_TRANSLATION_SEARCH_FAST_TERMS_H
#define CCTBX_TRANSLATION_SEARCH_FAST_TERMS_H

#include <cctbx/translation_search/fast_nv1995/summations.h>
#include <cctbx/maptbx/copy.h>
#include <scitbx/fftpack/real_to_complex_3d.h>

namespace cctbx { namespace translation_search {

  //! Computation of fast translation function terms.
  /*! Navaza, J. & Vernoslova, E. (1995). Acta Cryst. A51, 445-449.
   */
  template <typename FloatType=double>
  class fast_terms
  {
    public:
      //! Initialization of intermediate arrays.
      /*! Optionally, multiple summations can reuse the intermediate arrays.
       */
      fast_terms(
        af::int3 const& gridding,
        bool anomalous_flag,
        af::const_ref<miller::index<> > const& miller_indices_p1_f_calc,
        af::const_ref<std::complex<FloatType> > const& p1_f_calc)
      :
        rfft_(gridding),
        fc_map_(miller_indices_p1_f_calc, p1_f_calc, anomalous_flag)
      {}

      //! Summations according to Navaza & Vernoslova (1995).
      /*! squared_flag == false: summation according to equation 14.<br>
          The array m corresponds to mH delta_I_obs_H in equation 14.
          <p>
          squared_flag == true: summation according to equation 15.<br>
          The array m corresponds to mH in equation 15.
          <p>
          The summation is only over the space_group.order_p()
          symmetry operations. Compared to the summation over
          all space_group.order_z() symmetry operations, the
          results are smaller by a factor space_group.n_ltr()**p,
          with p=2 for squared_flag=false,
          and p=4 for squared_flag=true.
          <p>
          Target functions other than the correlation coefficient
          may be evaluated by passing different factors in the
          array m.
       */
      fast_terms&
      summation(
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices_f_obs,
        af::const_ref<FloatType> const& m,
        af::const_ref<std::complex<FloatType> > const& f_part,
        bool squared_flag)
      {
        accu_set_to_zero();
        if (squared_flag) {
          fast_nv1995_detail::summation_eq15(
            space_group, miller_indices_f_obs, m, f_part,
            fc_map_, accu_);
        }
        else {
          fast_nv1995_detail::summation_eq14(
            space_group, miller_indices_f_obs, m, f_part,
            fc_map_, accu_);
        }
        return *this;
      }

      /*! \brief In-place fast Fourier transformation of the result
          of the summation().
       */
      fast_terms&
      fft()
      {
        rfft_.backward(accu_mem_.ref());
        return *this;
      }

      //! Reference to final result.
      /*! Not available in Python.
       */
      af::const_ref<FloatType, af::c_grid_padded<3> >
      accu_real_const_ref() const
      {
        return af::const_ref<FloatType, af::c_grid_padded<3> >(
          reinterpret_cast<const FloatType*>(accu_mem_.begin()),
          af::c_grid_padded<3>(rfft_.m_real(), rfft_.n_real()));
      }

      //! Copy of final result.
      af::versa<FloatType, af::c_grid<3> >
      accu_real_copy() const
      {
        af::versa<FloatType, af::c_grid<3> > target_map(
          af::c_grid<3>(rfft_.n_real()));
        maptbx::copy(accu_real_const_ref(), target_map.ref());
        return target_map;
      }

    protected:
      void
      accu_set_to_zero()
      {
        if (accu_mem_.size() == 0) {
          accu_mem_.resize(af::c_grid<3>(rfft_.n_complex()));
          accu_ = fast_nv1995_detail::summation_accumulator<FloatType>(
            accu_mem_.begin(),
            miller::index<>(rfft_.n_real()),
            miller::index<>(rfft_.n_complex()));
        }
        else {
          accu_mem_.fill(0);
        }
      }

      scitbx::fftpack::real_to_complex_3d<FloatType> rfft_;
      miller::f_calc_map<FloatType> fc_map_;
      af::versa<std::complex<FloatType>, af::c_grid<3> > accu_mem_;
      fast_nv1995_detail::summation_accumulator<FloatType> accu_;
  };

}} // namespace cctbx::translation_search

#endif // CCTBX_TRANSLATION_SEARCH_FAST_TERMS_H
