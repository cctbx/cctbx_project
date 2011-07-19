#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_H

/// Crystallographic least-squares

#include <scitbx/lstbx/normal_equations.h>
#include <scitbx/math/accumulators.h>
#include <scitbx/sparse/matrix.h>
#include <scitbx/array_family/ref_reductions.h>

#include <cctbx/sgtbx/seminvariant.h>
#include <cctbx/sgtbx/rot_mx.h>
#include <cctbx/xray/extinction.h>
#include <cctbx/xray/observations.h>

#include <smtbx/refinement/constraints/scatterer_parameters.h>
#include <smtbx/structure_factors/direct/standard_xray.h>

#include <algorithm>
#include <iterator>


namespace smtbx { namespace refinement { namespace least_squares {

  namespace lstbx = scitbx::lstbx;

  /** \brief Restraints to prevent the structure to freely move as a whole
   along the continuous shift directions if any.

   The implementation has a key limitation: the continuous shift directions
   must lie along unit cell axes if there are atoms on special positions on
   axes parallel to continuous shift directions.
   It will therefore work for all space-groups in standard settings
   but it will not work in, e.g., a setting of R3 with
   the 3-fold axis along the direction (1,1,1), as in the space group
   Hall: R 3 (-y+z, x+z, -x+y+z), with an atom at the position (0.5, 0.5, 0.5).
   As far as we know, this is a limitation of the programs ShelXL
   and Crystals too. In order to keep it simple, an exception is thrown if
   a continuous shift directions is found to be non-trivial, whether there are
   atoms on the wrong kind of special positions or not.

   The weight of the restraints is dynamically adjusted to the normal matrix
   they are added to, as a relative weight times the smallest
   relevant diagonal element.
   */
  template <typename FloatType>
  class floating_origin_restraints
  {
  public:
    typedef FloatType scalar_t;

    scalar_t floating_origin_restraint_relative_weight;

    /// Singular directions in the space of crystallographic parameters
    af::small<af::shared<scalar_t>, 3> singular_directions;

    af::small<scitbx::vec3<scalar_t>, 3> origin_shifts;

    /// A mask pinpointing which parameters span the singular space
    af::shared<bool> singular_space;

  public:
    /// Construct with the given relative weight
    /** It precomputes the singular vectors of the normal matrix
        for the crystallographic parameters.
     */
    floating_origin_restraints(
      sgtbx::space_group const &space_group,
      af::const_ref<constraints::scatterer_parameters> const &params,
      scalar_t floating_origin_restraint_relative_weight)

      : floating_origin_restraint_relative_weight(
          floating_origin_restraint_relative_weight)
    {
      // Floating origin restraints: pre-compute singular directions
      sgtbx::structure_seminvariants seminvariants(space_group);
      SMTBX_ASSERT(seminvariants.continuous_shifts_are_principal());
      af::small<sgtbx::ss_vec_mod, 3> const
      &vm = seminvariants.vectors_and_moduli();
      for (int i=0; i<vm.size(); ++i) {
        if (vm[i].m != 0) continue;
        // allowed continuous origin shift:
        scitbx::vec3<scalar_t> v(vm[i].v[0], vm[i].v[1], vm[i].v[2]);
        origin_shifts.push_back(v);

        // Fill singular_dir with the singular direction corresponding to v
        af::shared<scalar_t> singular_dir(af::reserve(5*params.size()));
        std::back_insert_iterator< af::shared<scalar_t> > g(singular_dir);
        for (std::size_t i=0; i<params.size(); ++i) {
          BOOST_FOREACH (constraints::asu_parameter const *p,
                         params[i].ordered())
          {
            if (p == params[i].site) std::copy(v.begin(), v.end(), g);
            else std::fill_n(g, p->size(), scalar_t(0));
          }
        }

        // Save it
        singular_directions.push_back(singular_dir);
      }

      // Pinpoint coordinates participating in the singularity
      if (singular_directions.size()) {
        int n = singular_directions[0].size();
        singular_space = af::shared<bool>(n, false);
        for (int i=0; i<n; ++i) {
          for (int j=0; j<singular_directions.size(); ++j) {
            singular_space[i] = singular_space[i] || singular_directions[j][i];
          }
        }
      }
    }

    /// Add floating origin restraints to the given normal equations
    /** The Jacobian is that of the reparametrisation of crystallographic
        parameters as function of the independent parameters that the
        normal equations are built with
     */
    void add_to(lstbx::normal_equations::linear_ls<scalar_t> &normal_eqns,
                scitbx::sparse::matrix<FloatType> const
                &jacobian_transpose_matching_grad_fc)
    {
      if (!floating_origin_restraint_relative_weight) return;
      if (!singular_space.size()) return;
      af::ref<scalar_t, af::packed_u_accessor>
      a = normal_eqns.normal_matrix().ref();
      scitbx::math::accumulator::min_max_accumulator<scalar_t> acc(a(0,0));
      for (int i=1; i<a.n_rows(); ++i) {
        if (singular_space[i]) acc(a(i,i));
      }
      scalar_t w = floating_origin_restraint_relative_weight * acc.max();
      for (int i=0; i<singular_directions.size(); ++i) {
        af::shared<scalar_t> reparametrised_singular_dir =
        jacobian_transpose_matching_grad_fc*singular_directions[i].const_ref();
        normal_eqns.add_equation(0, reparametrised_singular_dir.const_ref(), w);
      }
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

  /** \brief Build normal equations for the given data, model, weighting
       and constraints.

      The constraints is performed with a reparametrisation whose Jacobian
      transpose is passed as an argument.
   */
  template <typename FloatType>
  struct build_normal_equations
  {
    //! Default constructor. Some data members are not initialized!
    build_normal_equations() {}

    template <template<typename> class NormalEquations,
              template<typename> class WeightingScheme,
              class OneMillerIndexFcalc>
    build_normal_equations(
      NormalEquations<FloatType> &normal_equations,
      cctbx::xray::observations<FloatType> const &reflections,
      af::const_ref<std::complex<FloatType> > const &f_mask,
      WeightingScheme<FloatType> const &weighting_scheme,
      boost::optional<FloatType> scale_factor,
      OneMillerIndexFcalc &f_calc_function,
      scitbx::sparse::matrix<FloatType> const
        &jacobian_transpose_matching_grad_fc,
      cctbx::xray::extinction_correction<FloatType> &exti,
      bool objective_only=false)
    :
      f_calc_(reflections.size()),
      observables_(reflections.size()),
      weights_(reflections.size())
    {
      // Accumulate equations Fo(h) ~ Fc(h)
      SMTBX_ASSERT(!(reflections.has_twin_components() && f_mask.size()));
      SMTBX_ASSERT((!f_mask.size() || f_mask.size() == reflections.size()) ||
                  (f_mask.size()==reflections.size()));
      SMTBX_ASSERT(scale_factor || weighting_scheme.f_calc_independent());
      const bool use_cache = false; //reflections.has_twin_components();
      f_calc_function_with_cache<FloatType, OneMillerIndexFcalc>
        f_calc_func(f_calc_function, use_cache);
      bool compute_grad = !objective_only;
      reflections.update_prime_fraction();
      /* Quite hackish but the assert above makes it safe providing
         the weighting scheme class plays ball.
       */
      FloatType scale_factor_ = scale_factor ? *scale_factor : 0;
      af::shared<FloatType> gradients =
        compute_grad ? af::shared<FloatType> (
                         jacobian_transpose_matching_grad_fc.n_rows(),
                         af::init_functor_null<FloatType>())
                      :
                        af::shared<FloatType>();
      for (int i_h=0; i_h<reflections.size(); ++i_h) {
        miller::index<> const &h = reflections.index(i_h);
        if (f_mask.size()) {
          f_calc_func.compute(h, f_mask[i_h], compute_grad);
        }
        else {
          f_calc_func.compute(h, compute_grad);
        }
        if (compute_grad) {
          gradients =
            jacobian_transpose_matching_grad_fc*f_calc_func.grad_observable;
        }
        // sort out twinning
        FloatType observable =
          process_twinning(reflections, i_h, f_calc_func,
            gradients, jacobian_transpose_matching_grad_fc, compute_grad);
        // extinction correction
        af::tiny<FloatType,2> exti_k = exti.compute(h, observable, compute_grad);
        observable *= exti_k[0];
        f_calc_[i_h] = f_calc_func.f_calc*std::sqrt(exti_k[0]);
        observables_[i_h] = observable;

        FloatType weight = weighting_scheme(reflections.fo_sq(i_h),
          reflections.sig(i_h), observable, scale_factor_);
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
      }
      normal_equations.finalise(objective_only);
    }

    af::shared<std::complex<FloatType> > f_calc() { return f_calc_; }

    af::shared<FloatType> observables() { return observables_; }

    af::shared<FloatType> weights() { return weights_; }

  protected:
    template <class OneMillerIndexFcalc>
    FloatType process_twinning(
      cctbx::xray::observations<FloatType> const &reflections,
      int h_i,
      f_calc_function_with_cache<FloatType, OneMillerIndexFcalc> &f_calc_func,
      af::shared<FloatType> &gradients,
      scitbx::sparse::matrix<FloatType> const
        &jacobian_transpose_matching_grad_fc,
      bool compute_grad)
    {
      FloatType obs = f_calc_func.observable;
      if (reflections.has_twin_components()) {
        typename cctbx::xray::observations<FloatType>::iterator_ itr =
          reflections.iterator(h_i);
        const FloatType identity_part = obs;
        obs *= reflections.scale(h_i);
        if (compute_grad) {
          gradients *= reflections.scale(h_i);
        }
        while (itr.has_next()) {
          typename cctbx::xray::observations<FloatType>::index_twin_component twc =
            itr.next();
          f_calc_func.compute(twc.h, compute_grad);
          obs += twc.fraction->value*f_calc_func.observable;
          if (compute_grad) {
            af::shared<FloatType> tmp_gradients =
              jacobian_transpose_matching_grad_fc*f_calc_func.grad_observable;
            gradients += twc.fraction->value*tmp_gradients;
            if (twc.fraction->grad) {
              SMTBX_ASSERT(!(twc.fraction->grad_index < 0 ||
                twc.fraction->grad_index >= gradients.size()));
              gradients[twc.fraction->grad_index] +=
                f_calc_func.observable - identity_part;
            }
          }
        }
      }
      return obs;
    }

  private:
    af::shared<std::complex<FloatType> > f_calc_;
    af::shared<FloatType> observables_;
    af::shared<FloatType> weights_;
  };

}}}


#endif // GUARD
