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
  int max_memory;

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
    af::versa<FloatType, af::c_grid<2> >& design_matrix,
    int &max_memory)
    : normal_equations_ptr(normal_equations_ptr), normal_equations(*normal_equations_ptr),
    reflections(reflections), f_mask(f_mask), twp(twp),
    weighting_scheme(weighting_scheme),
    scale_factor(scale_factor),
    f_calc_function_ptr(f_calc_function_ptr), f_calc_function(*f_calc_function_ptr),
    jacobian_transpose_matching_grad_fc(jacobian_transpose_matching_grad_fc),
    exti(exti),
    objective_only(objective_only), compute_grad(!objective_only),
    f_calc(f_calc), observables(observables), weights(weights),
    design_matrix(design_matrix), max_memory(max_memory)
  {}

  void operator()() {
    try {
      const int n = reflections.size();
      const int n_rows = jacobian_transpose_matching_grad_fc.n_rows();
      const int threads = 2*parent_t::get_available_threads();
      //const int threads = 2;
      bool error_flag = false;
      std::string error_string;
      std::vector<FloatType> gradients, matrix, yo_dot_grad_yc_, yc_dot_grad_yc_;
      boost::ptr_vector<boost::shared_ptr<f_calc_function_base_t> > f_calc_threads;
      f_calc_threads.resize(threads);
      for (int i = 0; i < threads; i++) {
        f_calc_threads[i] = f_calc_function.fork();
      }
      const FloatType temp_memory = threads * ((n_rows * (n_rows + 1) / 2) + 3 * n_rows) * sizeof(FloatType) / 1048576.0;
      const FloatType mem_per_size = n_rows * sizeof(FloatType) / 1048576.0;
      int size = n;
      FloatType req_mem = temp_memory + size * mem_per_size;
      while (req_mem > max_memory) {
        size -= threads;
        req_mem -= mem_per_size * threads;
        if (size >= n) {
          break;
        }
      }
      matrix.resize(threads * (n_rows * (n_rows + 1) / 2));
      yo_dot_grad_yc_.resize(threads * n_rows);
      yc_dot_grad_yc_.resize(threads * n_rows);
      if (compute_grad && !build_design_matrix && size < n) {
        gradients.resize(size * n_rows);
      }
      for (int i_h = 0; i_h * size < n; i_h++) {
        const int start = i_h * size;
        //check whetehr last chunk is smaller
        if (start + size >= n) {
          size = n - start;
          if (compute_grad && !build_design_matrix) {
            gradients.resize(size * n_rows);
          }
        }
#pragma omp parallel num_threads(threads)
        {
#pragma omp for nowait
          for (int dummy = 0; dummy < 3; dummy++) {
            if (dummy == 0) {
              std::fill(matrix.begin(), matrix.end(), 0);
            }
            else if (dummy == 1) {
              std::fill(yo_dot_grad_yc_.begin(), yo_dot_grad_yc_.end(), 0);
            }
            else if (dummy == 2) {
              std::fill(yc_dot_grad_yc_.begin(), yc_dot_grad_yc_.end(), 0);
            }
          }
          /* Make a gradient vector for each thread
             Must pre-allocate or Jt.G causes a crash
          */
          af::shared<FloatType> gradient(n_rows);
          const int thread = omp_get_thread_num();
#pragma omp for
          for (int run = 0; run < size; run++) {
            if (error_flag) {
              continue;
            }
            const int refl_i = start + run;
            miller::index<> const& h = reflections.index(refl_i);
            const twin_fraction<FloatType>* fraction = reflections.fraction(i_h);
            try {
              if (f_mask.size()) {
                f_calc_threads[thread]->compute(h, f_mask[refl_i], fraction, compute_grad);
              }
              else {
                f_calc_threads[thread]->compute(h, boost::none, fraction, compute_grad);
              }
            }
            catch (smtbx::error const& e) {
              error_flag = true;;
              error_string = e.what();
              continue;
            }
            catch (std::exception const& e) {
              error_flag = true;
              error_string = e.what();
              continue;
            }
            f_calc[refl_i] = f_calc_threads[thread]->get_f_calc();
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
            f_calc[refl_i] *= std::sqrt(exti_k[0]);
            observables[refl_i] = observable;

            FloatType weight = weighting_scheme(reflections.fo_sq(refl_i),
              reflections.sig(refl_i), observable, scale_factor);
            weights[refl_i] = weight;
            if (!objective_only) {
              if (exti.grad_value()) {
                FloatType exti_der = (exti_k[0] + pow(exti_k[0], 3)) / 2;
                int grad_index = exti.get_grad_index();
                SMTBX_ASSERT(!(grad_index < 0 || grad_index >= gradient.size()));
                gradient[grad_index] += exti_k[1] * exti_der;
              }
              if (!build_design_matrix) {
                memcpy(&gradients[run * n_rows], gradient.begin(), sizeof(FloatType) * n_rows);
              }
            }
            if (build_design_matrix) {
              memcpy(&design_matrix(refl_i, 0), gradient.begin(),
                gradient.size() * sizeof(FloatType));
            }
          }
        }
        if (!error_flag) {
          if (objective_only) {
            if (weights.size()) {
              normal_equations.add_residuals_omp(size, start, threads, observables,
                reflections.data().const_ref(), weights);
            }
            else {
              normal_equations.add_residuals_omp(size, start, threads, observables,
                reflections.data().const_ref());
            }
          }
          else if (!build_design_matrix) {
            normal_equations.add_equations_omp(n, n_rows, size, start, threads,
              matrix, yo_dot_grad_yc_, yc_dot_grad_yc_,
              observables, gradients, reflections.data().const_ref(), weights);
          }
        }
        if (error_flag) {
          break;
        }
      }
      if (error_flag) {
        exception_.reset(new smtbx::error(error_string));
      }
      //Cleaning up memory in parallel using a faster way than simply clear()
      //Compare http://www.gotw.ca/gotw/054.htm
#pragma omp parallel for num_threads(threads)
      for (int dummy = 0; dummy < 3; dummy++) {
        if (dummy == 0) {
          matrix.clear(), std::vector<FloatType>(matrix).swap(matrix);
        }
        else if (dummy == 1) {
          yo_dot_grad_yc_.clear(), std::vector<FloatType>(yo_dot_grad_yc_).swap(yo_dot_grad_yc_);
        }
        else if (dummy == 2) {
          yc_dot_grad_yc_.clear(), std::vector<FloatType>(yc_dot_grad_yc_).swap(yc_dot_grad_yc_);
        }
        else if (dummy == 3) {
          gradients.clear(), std::vector<FloatType>(gradients).swap(gradients);
        }
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