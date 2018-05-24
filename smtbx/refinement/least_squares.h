#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_H

/// Crystallographic least-squares

#include <scitbx/lstbx/normal_equations.h>
#include <scitbx/sparse/matrix.h>
#include <scitbx/array_family/ref_reductions.h>

#include <cctbx/xray/extinction.h>
#include <cctbx/xray/observations.h>

#include <smtbx/error.h>
#include <smtbx/structure_factors/direct/standard_xray.h>

#include <algorithm>


namespace smtbx { namespace refinement { namespace least_squares {

  namespace lstbx = scitbx::lstbx;

  /** \brief Build normal equations for the given data, model, weighting
       and constraints. Optionally builds the design matrix.

      The constraints is performed with a reparametrisation whose Jacobian
      transpose is passed as an argument.
   */
  template <typename FloatType, bool build_design_matrix>
  struct build_design_matrix_and_normal_equations
  {
    //! Default constructor. Some data members are not initialized!
    build_design_matrix_and_normal_equations() {}

    template <class NormalEquations,
              template<typename> class WeightingScheme,
              class OneMillerIndexFcalc>
    build_design_matrix_and_normal_equations(
      NormalEquations &normal_equations,
      cctbx::xray::observations<FloatType> const &reflections,
      af::const_ref<std::complex<FloatType> > const &f_mask,
      WeightingScheme<FloatType> const &weighting_scheme,
      boost::optional<FloatType> scale_factor,
      OneMillerIndexFcalc &f_calc_function,
      scitbx::sparse::matrix<FloatType> const
        &jacobian_transpose_matching_grad_fc,
      cctbx::xray::extinction_correction<FloatType> const &exti,
      bool objective_only)
    :
      f_calc_(reflections.size()),
      observables_(reflections.size()),
      weights_(reflections.size()),
      design_matrix_(af::c_grid<2>(build_design_matrix ? reflections.size() : 0,
        build_design_matrix ? jacobian_transpose_matching_grad_fc.n_rows() : 0))

    {
      // Accumulate equations Fo(h) ~ Fc(h)
      SMTBX_ASSERT(!(reflections.has_twin_components() && f_mask.size()));
      SMTBX_ASSERT((!f_mask.size() || f_mask.size() == reflections.size()))
                  (f_mask.size())(reflections.size());
      bool compute_grad = !objective_only;
      reflections.update_prime_fraction();
      af::shared<FloatType> gradients;
      if (compute_grad) {
        gradients.resize(jacobian_transpose_matching_grad_fc.n_rows());
      }
      for (int i_h=0; i_h<reflections.size(); ++i_h) {
        miller::index<> const &h = reflections.index(i_h);
        if (f_mask.size()) {
          f_calc_function.compute(h, f_mask[i_h], compute_grad);
        }
        else {
          f_calc_function.compute(h, boost::none, compute_grad);
        }
        if (compute_grad) {
          gradients =
            jacobian_transpose_matching_grad_fc*f_calc_function.grad_observable;
        }
        // sort out twinning
        FloatType observable =
          process_twinning(reflections, i_h, f_calc_function,
            gradients, jacobian_transpose_matching_grad_fc, compute_grad);
        // extinction correction
        af::tiny<FloatType,2> exti_k = exti.compute(h, observable, compute_grad);
        observable *= exti_k[0];
        f_calc_[i_h] = f_calc_function.f_calc*std::sqrt(exti_k[0]);
        observables_[i_h] = observable;

        FloatType weight = weighting_scheme(reflections.fo_sq(i_h),
          reflections.sig(i_h), observable, scale_factor);
        weights_[i_h] = weight;
        if (objective_only) {
          normal_equations.add_residual(observable,
            reflections.fo_sq(i_h), weight);
        }
        else {
          if (exti.grad_value()) {
            int grad_index = exti.get_grad_index();
            SMTBX_ASSERT(!(grad_index < 0 || grad_index >= gradients.size()));
            gradients[grad_index] += exti_k[1];
          }
          normal_equations.add_equation(observable,
            gradients.ref(), reflections.fo_sq(i_h), weight);
        }
        if (build_design_matrix) {
          for (int i_g = 0; i_g < gradients.size(); i_g++) {
            design_matrix_(i_h, i_g) = gradients[i_g];
          }
        }
      }
      normal_equations.finalise(objective_only);
    }

    af::shared<std::complex<FloatType> > f_calc() { return f_calc_; }

    af::shared<FloatType> observables() { return observables_; }

    af::shared<FloatType> weights() { return weights_; }

  protected:
    af::versa<FloatType, af::c_grid<2> > design_matrix() { return design_matrix_; }

    template<class OneMillerIndexFcalc>
    FloatType process_twinning(
      cctbx::xray::observations<FloatType> const &reflections,
      int i_h,
      OneMillerIndexFcalc &f_calc_function,
      af::shared<FloatType> &gradients,
      scitbx::sparse::matrix<FloatType> const
        &jacobian_transpose_matching_grad_fc,
      bool compute_grad)
    {
      FloatType obs = f_calc_function.observable;
      if (reflections.has_twin_components()) {
        typename cctbx::xray::observations<FloatType>::iterator_ itr =
          reflections.iterator(i_h);
        scitbx::af::shared<cctbx::xray::twin_fraction<FloatType> const*>
          to_update;
        FloatType identity_part = 0,
          obs_scale = reflections.scale(i_h);
        obs *= obs_scale;
        if (compute_grad) {
          gradients *= obs_scale;
          if (itr.measured_fraction != 0 && itr.measured_fraction->grad) {
            gradients[itr.measured_fraction->grad_index] +=
              f_calc_function.observable;
            to_update.push_back(itr.measured_fraction);
          }
        }
        while (itr.has_next()) {
          typename cctbx::xray::observations<FloatType>::index_twin_component
            twc = itr.next();
          f_calc_function.compute(twc.h, boost::none, compute_grad);
          obs += twc.scale()*f_calc_function.observable;
          if (compute_grad) {
            af::shared<FloatType> tmp_gradients =
              jacobian_transpose_matching_grad_fc*f_calc_function.grad_observable;
            gradients += twc.scale()*tmp_gradients;
            if (twc.fraction != 0 && twc.fraction->grad) {
              SMTBX_ASSERT(!(twc.fraction->grad_index < 0 ||
                twc.fraction->grad_index >= gradients.size()));
              gradients[twc.fraction->grad_index] += f_calc_function.observable;
              to_update.push_back(twc.fraction);
            }
            else if (twc.fraction == 0) {
              identity_part = f_calc_function.observable;
            }
          }
        }
        if (identity_part != 0) {
          for (size_t i=0; i < to_update.size(); i++) {
            gradients[to_update[i]->grad_index] -= identity_part;
          }
        }
      }
      return obs;
    }

  protected:
    af::shared<std::complex<FloatType> > f_calc_;
    af::shared<FloatType> observables_;
    af::shared<FloatType> weights_;
    af::versa<FloatType, af::c_grid<2> > design_matrix_;
  };

  /** \brief Build normal equations for the given data, model, weighting
  and constraints.

  The constraints is performed with a reparametrisation whose Jacobian
  transpose is passed as an argument.
  */
  template <typename FloatType>
  struct build_normal_equations
    : public build_design_matrix_and_normal_equations<FloatType, false>
  {
    //! Default constructor. Some data members are not initialized!
    build_normal_equations()
      :
      build_design_matrix_and_normal_equations<FloatType, false>()
    {}

    template <class NormalEquations,
      template<typename> class WeightingScheme,
      class OneMillerIndexFcalc>
    build_normal_equations(
      NormalEquations &normal_equations,
      cctbx::xray::observations<FloatType> const &reflections,
      af::const_ref<std::complex<FloatType> > const &f_mask,
      WeightingScheme<FloatType> const &weighting_scheme,
      boost::optional<FloatType> scale_factor,
      OneMillerIndexFcalc &f_calc_function,
      scitbx::sparse::matrix<FloatType> const
        &jacobian_transpose_matching_grad_fc,
      cctbx::xray::extinction_correction<FloatType> const &exti,
      bool objective_only = false)
      :
      build_design_matrix_and_normal_equations<FloatType, false>(
        normal_equations,
        reflections, f_mask, weighting_scheme, scale_factor, f_calc_function,
        jacobian_transpose_matching_grad_fc, exti,
        objective_only)
    {}
  };

  /** \brief Build normal equations for the given data, model, weighting
  and constraints and the buld the design matrix

  The constraints is performed with a reparametrisation whose Jacobian
  transpose is passed as an argument.
  */
  template <typename FloatType>
  struct build_design_matrix
    : public build_design_matrix_and_normal_equations<FloatType, true>
  {
    //! Default constructor. Some data members are not initialized!
    build_design_matrix()
      :
      build_design_matrix_and_normal_equations<FloatType, true>()
    {}

    template <class NormalEquations,
      template<typename> class WeightingScheme,
      class OneMillerIndexFcalc>
    build_design_matrix(
      NormalEquations &normal_equations,
      cctbx::xray::observations<FloatType> const &reflections,
      af::const_ref<std::complex<FloatType> > const &f_mask,
      WeightingScheme<FloatType> const &weighting_scheme,
      boost::optional<FloatType> scale_factor,
      OneMillerIndexFcalc &f_calc_function,
      scitbx::sparse::matrix<FloatType> const
      &jacobian_transpose_matching_grad_fc,
      cctbx::xray::extinction_correction<FloatType> const &exti,
      bool objective_only = false)
      :
      build_design_matrix_and_normal_equations<FloatType, true>(
        normal_equations,
        reflections, f_mask, weighting_scheme, scale_factor, f_calc_function,
        jacobian_transpose_matching_grad_fc, exti,
        objective_only)
    {}

    af::versa<FloatType, af::c_grid<2> > design_matrix() {
      return build_design_matrix_and_normal_equations<FloatType, true>::design_matrix_;
    }
  };

  template <typename FloatType>
  struct f_calc_function_result
  {
    f_calc_function_result(
      FloatType const &observable_,
      std::complex<FloatType> const &f_calc_,
      af::shared<FloatType> const &grad_observable_)
      :
    observable(observable_),
    f_calc(f_calc_),
    grad_observable(grad_observable_)
    {}

    f_calc_function_result(
      FloatType const &observable_,
      std::complex<FloatType> const &f_calc_)
      :
    observable(observable_),
    f_calc(f_calc_),
    grad_observable()
    {}

    FloatType const observable;
    std::complex<FloatType> const f_calc;
    af::shared<FloatType> const grad_observable;
  };

  /*  A thin wrapper around OneMillerIndexFcalc to enable caching of the
      results for symmetry related indices.
   */
  template <typename FloatType, class OneMillerIndexFcalc>
  struct f_calc_function_with_cache
  {
    f_calc_function_with_cache(
      OneMillerIndexFcalc &f_calc_function_, bool use_cache_=false)
      :
    f_calc_function(f_calc_function_),
    use_cache(use_cache_),
    observable(),
    grad_observable(),
    f_calc(),
    length_sq(0),
    cache()
    {};

    void compute(
      miller::index<> const &h,
      boost::optional<std::complex<FloatType> > const &f_mask=boost::none,
      bool compute_grad=true)
    {
      if (!use_cache) {
        f_calc_function.compute(h, f_mask, compute_grad);
        observable = f_calc_function.observable;
        grad_observable = f_calc_function.grad_observable;
        f_calc = f_calc_function.f_calc;
      }
      else {
        FloatType h_length_sq = h.length_sq();
        if (h_length_sq != length_sq) {
          cache.clear();
          length_sq = h_length_sq;
        }
        typename cache_t::iterator iter = cache.find(h);
        if (iter == cache.end()) {
          f_calc_function.linearise(h, f_mask);
          observable = f_calc_function.observable;
          grad_observable = f_calc_function.grad_observable;
          f_calc = f_calc_function.f_calc;
          cache.insert(
            std::pair<miller::index<>, f_calc_function_result<FloatType> >(
              h, f_calc_function_result<FloatType>(
                  observable,
                  f_calc_function.f_calc,
                  grad_observable.array().deep_copy())));
        }
        else {
          observable = iter->second.observable;
          f_calc = iter->second.f_calc;
          grad_observable =
            af::ref_owning_shared<FloatType>(iter->second.grad_observable);
        }
      }
    }

    void compute(miller::index<> const &h,
                 bool compute_grad=true)
    {
      compute(h, /*f_mask=*/ boost::none, compute_grad);
    }

    typedef
      std::map<miller::index<>, f_calc_function_result<FloatType> > cache_t;

    OneMillerIndexFcalc &f_calc_function;
    bool use_cache;
    FloatType observable;
    af::ref_owning_shared<FloatType> grad_observable;
    std::complex<FloatType> f_calc;
    FloatType length_sq;
    cache_t cache;
  };


  /// A specialisation of build_normal_equations with caching of
  /// the computation of Fc(h) and its derivatives
  /** This avoids performing the same computation twice for the same h.
      At the moment, the caching works only for (pseudo-)merohedral twins,
      and this code is not used in production code.
   */
  template <typename FloatType>
  struct build_normal_equations_with_caching
  {
    //! Default constructor. Some data members are not initialized!
    build_normal_equations_with_caching() {}

    template <class NormalEquations,
              template<typename> class WeightingScheme,
              class OneMillerIndexFcalc>
    build_normal_equations_with_caching(
      NormalEquations &normal_equations,
      cctbx::xray::observations<FloatType> const &reflections,
      af::const_ref<std::complex<FloatType> > const &f_mask,
      WeightingScheme<FloatType> const &weighting_scheme,
      boost::optional<FloatType> scale_factor,
      OneMillerIndexFcalc &f_calc_function,
      scitbx::sparse::matrix<FloatType> const
        &jacobian_transpose_matching_grad_fc,
      cctbx::xray::extinction_correction<FloatType> const &exti,
      bool objective_only=false)
    :
    build_normal_equations_with_caching(
      normal_equations,
      reflections,
      f_mask,
      weighting_scheme,
      scale_factor,
      f_calc_function_with_cache<FloatType, OneMillerIndexFcalc>(
        f_calc_function, reflections.has_twin_components()),
      jacobian_transpose_matching_grad_fc,
      exti,
      objective_only)
    {}
  };

}}}


#endif // GUARD
