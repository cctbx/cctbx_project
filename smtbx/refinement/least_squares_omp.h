/*
The following struct is optimized for cache efficient parallel 
calculations using OpenMP making sure there is no strides
in adressing memory by scheduling with chunksize 1
*/
template < class NormalEquations,
  template<typename> class WeightingScheme,
  class OneMillerIndexFcalc>
struct accumulate_reflection_chunk_omp {
  boost::scoped_ptr<smtbx::error> exception_;
  boost::shared_ptr<NormalEquations> normal_equations_ptr;
  NormalEquations& normal_equations;
  cctbx::xray::observations<FloatType> const& reflections;
  af::const_ref<std::complex<FloatType> > const& f_mask;
  WeightingScheme<FloatType> const& weighting_scheme;
  boost::optional<FloatType> scale_factor;
  boost::shared_ptr<OneMillerIndexFcalc> f_calc_function_ptr;
  OneMillerIndexFcalc& f_calc_function;
  scitbx::sparse::matrix<FloatType> const
      & jacobian_transpose_matching_grad_fc;
  cctbx::xray::extinction_correction<FloatType> const& exti;
  bool objective_only, compute_grad;
  af::ref<std::complex<FloatType> > f_calc;
  af::ref<FloatType> observables;
  af::ref<FloatType> weights;
  af::versa<FloatType, af::c_grid<2> >& design_matrix;
  accumulate_reflection_chunk_omp(
    boost::shared_ptr<NormalEquations> const& normal_equations_ptr,
    cctbx::xray::observations<FloatType> const& reflections,
    af::const_ref<std::complex<FloatType> > const& f_mask,
    WeightingScheme<FloatType> const& weighting_scheme,
    boost::optional<FloatType> scale_factor,
    boost::shared_ptr<OneMillerIndexFcalc> f_calc_function_ptr,
    scitbx::sparse::matrix<FloatType> const
    & jacobian_transpose_matching_grad_fc,
    cctbx::xray::extinction_correction<FloatType> const& exti,
    bool objective_only,
    af::ref<std::complex<FloatType> > f_calc,
    af::ref<FloatType> observables,
    af::ref<FloatType> weights,
    af::versa<FloatType, af::c_grid<2> >& design_matrix)
    : normal_equations_ptr(normal_equations_ptr), normal_equations(*normal_equations_ptr),
    reflections(reflections), f_mask(f_mask), weighting_scheme(weighting_scheme),
    scale_factor(scale_factor),
    f_calc_function_ptr(f_calc_function_ptr), f_calc_function(*f_calc_function_ptr),
    jacobian_transpose_matching_grad_fc(jacobian_transpose_matching_grad_fc),
    exti(exti),
    objective_only(objective_only), compute_grad(!objective_only),
    f_calc(f_calc), observables(observables), weights(weights),
    design_matrix(design_matrix)
  {}

  void operator()() {
    try {
      const int n = reflections.size();
      const int n_rows = jacobian_transpose_matching_grad_fc.n_rows();
      const int threads = get_available_threads();
      std::vector <std::vector<FloatType> > gradients;
      boost::ptr_vector<boost::shared_ptr<OneMillerIndexFcalc> > f_calc_threads;
      f_calc_threads.resize(threads);
      for (int i = 0; i < threads; i++) {
        f_calc_threads[i] = f_calc_function.fork();
      }
      if (compute_grad) {
        gradients.resize(n);
      }
#pragma omp parallel num_threads(threads)
      {
        // Make a gradient vector for each thread
        af::shared<FloatType> gradient(n_rows);
        const int thread = omp_get_thread_num();
#pragma omp for schedule(static,1)
        for (int i_h = 0; i_h < n; ++i_h) {
          miller::index<> const& h = reflections.index(i_h);
          if (f_mask.size()) {
            f_calc_threads[thread]->compute(h, f_mask[i_h], compute_grad);
          }
          else {
            f_calc_threads[thread]->compute(h, boost::none, compute_grad);
          }
          f_calc[i_h] = f_calc_threads[thread]->f_calc;
          FloatType observable;
          //skip hoarding memory if Gradients are not needed.
          if (compute_grad) {
            gradient = jacobian_transpose_matching_grad_fc * 
                        f_calc_threads[thread]->grad_observable;
            // sort out twinning
            observable = process_twinning_with_grads(i_h, 
                                          gradient, f_calc_threads[thread]);
            gradients[i_h].resize(n_rows);
            for (int g = 0; g < gradient.size(); g++) {
              gradients[i_h][g] = gradient[g];
            }
          }
          else {
            // sort out twinning
            observable = process_twinning(i_h, f_calc_threads[thread]);
          }
          af::tiny<FloatType, 2> exti_k = exti.compute(h, 
                                                  observable, compute_grad);
          observable *= exti_k[0];
          f_calc[i_h] *= std::sqrt(exti_k[0]);
          observables[i_h] = observable;

          FloatType weight = weighting_scheme(reflections.fo_sq(i_h),
              reflections.sig(i_h), observable, scale_factor);
          weights[i_h] = weight;
          if (exti.grad_value()) {
            int grad_index = exti.get_grad_index();
            SMTBX_ASSERT(!(grad_index < 0 
              || grad_index >= gradients.size()));
            gradients[i_h][grad_index] += exti_k[1];
          }
          if (build_design_matrix) {
            for (int i_g = 0; i_g < gradients[i_h].size(); i_g++) {
              design_matrix(i_h, i_g) = gradients[i_h][i_g];
            }
          }
        }
      }
      if (objective_only) {
        if (weights.size()) {
          normal_equations.add_residuals_omp(n, observables,
              reflections.data().ref(), weights);
        }
        else {
          normal_equations.add_residuals_omp(n, observables,
              reflections.data().ref());
        }
      }
      else {
        normal_equations.add_equations_omp(n, n_rows, threads,
            observables, gradients, reflections.data().ref(), weights);
      }
    }
    catch (smtbx::error const& e) {
      exception_.reset(new smtbx::error(e));
    }
    catch (std::exception const& e) {
      exception_.reset(new smtbx::error(e.what()));
    }
  }

  FloatType process_twinning_with_grads(int i_h, 
    af::shared<FloatType>& gradients, 
    boost::shared_ptr<OneMillerIndexFcalc> f_calc_thread) {
    typedef typename cctbx::xray::observations<FloatType>::iterator itr_t;
    typedef typename cctbx::xray::twin_fraction<FloatType> twf_t;
    typedef typename cctbx::xray::observations<FloatType>::index_twin_component twc_t;
    FloatType obs = f_calc_thread->observable;
    if (reflections.has_twin_components()) {
      itr_t itr = reflections.iterate(i_h);
      FloatType measured_part = obs,
          identity_part = 0,
          obs_scale = reflections.scale(i_h);
      obs *= obs_scale;
      const twf_t* measured_fraction = reflections.fraction(i_h);
      if (compute_grad) {
        gradients *= obs_scale;
        if (measured_fraction == 0) {
          identity_part = measured_part;
        }
      }
      std::size_t twc_cnt = 0;
      while (itr.has_next()) {
        twc_t twc = itr.next();
        f_calc_thread->compute(twc.h, boost::none, compute_grad);
        obs += twc.scale() * f_calc_thread -> observable;
        if (compute_grad) {
          af::shared<FloatType> tmp_gradients =
              jacobian_transpose_matching_grad_fc * f_calc_thread->grad_observable;
          gradients += twc.scale() * tmp_gradients;
          if (twc.fraction != 0) {
            if (twc.fraction->grad) {
              SMTBX_ASSERT(!(twc.fraction->grad_index < 0 ||
                  twc.fraction->grad_index >= gradients.size()));
              gradients[twc.fraction->grad_index] += f_calc_thread->observable;
            }
          }
          else {
            identity_part += f_calc_thread->observable;
          }
          twc_cnt++;
        }
      }
      if (compute_grad) {
        // consider multiple reflections with the 'prime' scale
        itr.reset();
        while (itr.has_next()) {
          twc_t twc = itr.next();
          if (twc.fraction != 0 && twc.fraction->grad) {
            SMTBX_ASSERT(!(twc.fraction->grad_index < 0 ||
              twc.fraction->grad_index >= gradients.size()));
            gradients[twc.fraction->grad_index] -= identity_part;
          }
        }
        if (twc_cnt != 0 && measured_fraction != 0 && measured_fraction->grad) {
          SMTBX_ASSERT(!(measured_fraction->grad_index < 0 ||
              measured_fraction->grad_index >= gradients.size()));
          gradients[measured_fraction->grad_index] +=
              measured_part - identity_part;
        }
      }
    }
    return obs;
  }
  FloatType process_twinning(int i_h, 
    boost::shared_ptr<OneMillerIndexFcalc> f_calc_thread) {
    typedef typename cctbx::xray::observations<FloatType>::iterator itr_t;
    typedef typename cctbx::xray::twin_fraction<FloatType> twf_t;
    typedef typename cctbx::xray::observations<FloatType>::index_twin_component twc_t;
    FloatType obs = f_calc_thread->observable;
    if (reflections.has_twin_components()) {
      itr_t itr = reflections.iterate(i_h);
      FloatType measured_part = obs,
          identity_part = 0,
          obs_scale = reflections.scale(i_h);
      obs *= obs_scale;
      const twf_t* measured_fraction = reflections.fraction(i_h);
      std::size_t twc_cnt = 0;
      while (itr.has_next()) {
        twc_t twc = itr.next();
        f_calc_thread->compute(twc.h, boost::none, compute_grad);
        obs += twc.scale() * f_calc_thread->observable;
      }
    }
    return obs;
  }
};