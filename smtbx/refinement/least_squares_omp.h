/*
The following struct is optimized for cache efficient parallel 
calculations using OpenMP making sure there is no strides
in adressing memory by scheduling with chunksize 1
*/
template<class NormalEquations,
  template<typename> class WeightingScheme>
struct accumulate_reflection_chunk_omp {
  typedef f_calc_function_base<FloatType>
    f_calc_function_base_t;

  boost::scoped_ptr<smtbx::error> exception_;
  boost::shared_ptr<NormalEquations> normal_equations_ptr;
  NormalEquations& normal_equations;
  cctbx::xray::observations<FloatType> const& reflections;
  af::const_ref<std::complex<FloatType> > const& f_mask;
  twinning_processor<FloatType> const& twp;
  WeightingScheme<FloatType> const& weighting_scheme;
  boost::optional<FloatType> scale_factor;
  boost::shared_ptr<f_calc_function_base_t> f_calc_function_ptr;
  f_calc_function_base_t& f_calc_function;
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
    twinning_processor<FloatType> const& twp,
    WeightingScheme<FloatType> const& weighting_scheme,
    boost::optional<FloatType> scale_factor,
    boost::shared_ptr<f_calc_function_base_t> f_calc_function_ptr,
    scitbx::sparse::matrix<FloatType> const&
      jacobian_transpose_matching_grad_fc,
    cctbx::xray::extinction_correction<FloatType> const& exti,
    bool objective_only,
    af::ref<std::complex<FloatType> > f_calc,
    af::ref<FloatType> observables,
    af::ref<FloatType> weights,
    af::versa<FloatType, af::c_grid<2> >& design_matrix)
    : normal_equations_ptr(normal_equations_ptr), normal_equations(*normal_equations_ptr),
    reflections(reflections), f_mask(f_mask), twp(twp),
    weighting_scheme(weighting_scheme),
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
      const int threads = parent_t::get_available_threads();
      std::vector<FloatType> gradients;
      boost::ptr_vector<boost::shared_ptr<f_calc_function_base_t> > f_calc_threads;
      f_calc_threads.resize(threads);
      for (int i = 0; i < threads; i++) {
        f_calc_threads[i] = f_calc_function.fork();
      }
      if (compute_grad && !build_design_matrix) {
        gradients.resize(n*n_rows);
      }
      twinning_processor<FloatType> twp(reflections, f_mask, compute_grad,
        jacobian_transpose_matching_grad_fc);
#pragma omp parallel num_threads(threads)
      {
        /* Make a gradient vector for each thread.
        Must pre-allocate or Jt.G causes a crash
        */
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
          f_calc[i_h] = f_calc_threads[thread]->get_f_calc();
          //skip hoarding memory if Gradients are not needed.
          if (compute_grad) {
            if (f_calc_threads[thread]->raw_gradients()) {
              gradient = jacobian_transpose_matching_grad_fc *
                f_calc_threads[thread]->get_grad_observable();
            }
            else {
              gradient = af::shared<FloatType>(
                f_calc_threads[thread]->get_grad_observable().begin(),
                f_calc_threads[thread]->get_grad_observable().end());
            }
          }
          // sort out twinning
          FloatType observable = twp.process(
            i_h, *f_calc_threads[thread], gradient);
          // extinction correction
          af::tiny<FloatType, 2> exti_k = exti.compute(h, observable, compute_grad);
          observable *= exti_k[0];
          f_calc[i_h] *= std::sqrt(exti_k[0]);
          observables[i_h] = observable;

          FloatType weight = weighting_scheme(reflections.fo_sq(i_h),
            reflections.sig(i_h), observable, scale_factor);
          weights[i_h] = weight;
          if (objective_only) {
            normal_equations.add_residual(observable,
              reflections.fo_sq(i_h), weight);
          }
          else {
            if (exti.grad_value()) {
              FloatType exti_der = (exti_k[0] + pow(exti_k[0], 3)) / 2;
              int grad_index = exti.get_grad_index();
              SMTBX_ASSERT(!(grad_index < 0 || grad_index >= gradients.size()));
              gradient[grad_index] += exti_k[1] * exti_der;
            }
            if (!build_design_matrix) {
              memcpy(&gradients[i_h * n_rows], gradient.begin(), sizeof(FloatType) * n_rows);
            }
          }
          if (build_design_matrix) {
            memcpy(&design_matrix(i_h, 0), gradient.begin(),
              gradient.size() * sizeof(FloatType));
          }
        }
      }
      if (objective_only) {
        if (weights.size()) {
          normal_equations.add_residuals_omp(n, observables,
              reflections.data().const_ref(), weights);
        }
        else {
          normal_equations.add_residuals_omp(n, observables,
              reflections.data().const_ref());
        }
      }
      else if (!build_design_matrix) {
        normal_equations.add_equations_omp(n, n_rows, threads,
            observables, gradients, reflections.data().const_ref(), weights);
      }
    }
    catch (smtbx::error const& e) {
      exception_.reset(new smtbx::error(e));
    }
    catch (std::exception const& e) {
      exception_.reset(new smtbx::error(e.what()));
    }
  }

};