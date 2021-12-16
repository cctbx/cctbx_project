#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_H

/// Crystallographic least-squares

#include <scitbx/lstbx/normal_equations.h>
#include <scitbx/sparse/matrix.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/matrix/tensors.h>

#include <cctbx/xray/extinction.h>
#include <cctbx/xray/observations.h>

#include <smtbx/error.h>
#include <smtbx/structure_factors/direct/standard_xray.h>

#include <algorithm>
#include <vector>
#include <omp.h>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/smart_ptr/scoped_ptr.hpp>
#include <boost/thread.hpp>

namespace smtbx { namespace refinement { namespace least_squares {

  namespace lstbx = scitbx::lstbx;

  /** \brief Build normal equations for the given data, model, weighting
  and constraints. Optionally builds the design matrix.

  The constraints is performed with a reparametrisation whose Jacobian
  transpose is passed as an argument.
  */

  template <typename FloatType, bool build_design_matrix>
  struct build_design_matrix_and_normal_equations {
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
      bool objective_only=false,
      bool may_parallelise_=false)
    :
      f_calc_(reflections.size()),
      observables_(reflections.size()),
      weights_(reflections.size()),
      design_matrix_(af::c_grid<2>(build_design_matrix ? reflections.size() : 0,
        build_design_matrix ? jacobian_transpose_matching_grad_fc.n_rows() : 0))
    {
      typedef accumulate_reflection_chunk<
                NormalEquations, WeightingScheme, OneMillerIndexFcalc>
              accumulate_reflection_chunk_t;
      typedef boost::shared_ptr<NormalEquations>
              normal_equations_ptr_t;
      typedef boost::shared_ptr<accumulate_reflection_chunk_t>
              accumulate_reflection_chunk_ptr_t;
      typedef boost::shared_ptr<OneMillerIndexFcalc>
              one_miller_index_fcalc_ptr_t;

      // Accumulate equations Fo(h) ~ Fc(h)
      SMTBX_ASSERT((!f_mask.size() || f_mask.size() == reflections.size()))
                  (f_mask.size())(reflections.size());
      reflections.update_prime_fraction();
      
      if(may_parallelise_) {
#if defined(_OPENMP)
        typedef accumulate_reflection_chunk_omp<
            NormalEquations, WeightingScheme, OneMillerIndexFcalc>
            accumulate_reflection_chunk_omp_t;
        accumulate_reflection_chunk_omp_t job(
            normal_equations_ptr_t(&normal_equations, null_deleter()),
            reflections, f_mask, weighting_scheme, scale_factor,
            one_miller_index_fcalc_ptr_t(&f_calc_function, null_deleter()),
            jacobian_transpose_matching_grad_fc,
            exti, objective_only,
            f_calc_.ref(), observables_.ref(), weights_.ref(),
            design_matrix_);
        job();
        if (job.exception_) {
            throw* job.exception_.get();
        }
        normal_equations.finalise(objective_only);
#else 
        //!!
        scitbx::matrix::tensors::initialise<FloatType>();

        const int thread_count = get_available_threads();
        const int equi_chunk_size = reflections.size()/thread_count;
        const int number_of_threads_doing_one_more = reflections.size() % thread_count;
        boost::thread_group pool;
        std::vector<accumulate_reflection_chunk_ptr_t> accumulators;
        std::size_t chunk_end = 0;
        if (fast_linalg::is_initialised()) {
          //When using fast_linalg in parallel mode it is reducing cache misses 
          //when leaving parallelisation to this loop
          openblas_set_num_threads(1);
        }
        for(int thread_idx=0; thread_idx<thread_count; thread_idx++) {
          std::size_t chunk_start = chunk_end;
          chunk_end +=
            thread_idx < number_of_threads_doing_one_more ? equi_chunk_size + 1
                                                          : equi_chunk_size;
          normal_equations_ptr_t chunk_normal_equations(
            new NormalEquations(normal_equations.n_parameters()));
          accumulate_reflection_chunk_ptr_t accumulator(
            new accumulate_reflection_chunk_t(
              chunk_start, chunk_end,
              chunk_normal_equations,
              reflections, f_mask, weighting_scheme, scale_factor,
              one_miller_index_fcalc_ptr_t(f_calc_function.fork()),
              jacobian_transpose_matching_grad_fc,
              exti, objective_only,
              f_calc_.ref(), observables_.ref(), weights_.ref(),
              design_matrix_));
          accumulators.push_back(accumulator);
          pool.create_thread(boost::ref(*accumulator));
        }
        pool.join_all();
        for(int thread_idx=0; thread_idx<thread_count; thread_idx++) {
          if (accumulators[thread_idx]->exception_) {
            throw *accumulators[thread_idx]->exception_.get();
          }
          normal_equations += accumulators[thread_idx]->normal_equations;
        }
        normal_equations.finalise(objective_only);
        if (fast_linalg::is_initialised()) {
            openblas_set_num_threads(thread_count);
        }
#endif
      }
      else {
        accumulate_reflection_chunk_t job(
          0, reflections.size(),
          normal_equations_ptr_t(&normal_equations, null_deleter()),
          reflections, f_mask, weighting_scheme, scale_factor,
          one_miller_index_fcalc_ptr_t(&f_calc_function, null_deleter()),
          jacobian_transpose_matching_grad_fc,
          exti, objective_only,
          f_calc_.ref(), observables_.ref(), weights_.ref(),
          design_matrix_);
        job();
        if (job.exception_) {
          throw *job.exception_.get();
        }
        normal_equations.finalise(objective_only);
      }

    }

    af::shared<std::complex<FloatType> > f_calc() { return f_calc_; }

    af::shared<FloatType> observables() { return observables_; }

    af::shared<FloatType> weights() { return weights_; }

    static int get_available_threads() {
      int &available = available_threads_var();
      if (available == -1) {
        available = std::max(1,
          static_cast<int>(boost::thread::physical_concurrency()));
      }
      return available;
    }

    static void set_available_threads(int thread_count) {
      // limit to the logical cores count
      available_threads_var() =
        std::max(1, std::min(
          static_cast<int>(boost::thread::hardware_concurrency()),
          thread_count));
    }

  protected:
    af::versa<FloatType, af::c_grid<2> > design_matrix() { return design_matrix_; }
    /// Accumulate from reflections whose indices are
    /// in the range [begin, end)
    template <class NormalEquations,
      template<typename> class WeightingScheme,
      class OneMillerIndexFcalc>
    struct accumulate_reflection_chunk {
      boost::scoped_ptr<smtbx::error> exception_;
      int begin, end;
      boost::shared_ptr<NormalEquations> normal_equations_ptr;
      NormalEquations &normal_equations;
      cctbx::xray::observations<FloatType> const &reflections;
      af::const_ref<std::complex<FloatType> > const &f_mask;
      WeightingScheme<FloatType> const &weighting_scheme;
      boost::optional<FloatType> scale_factor;
      boost::shared_ptr<OneMillerIndexFcalc> f_calc_function_ptr;
      OneMillerIndexFcalc &f_calc_function;
      scitbx::sparse::matrix<FloatType> const
        &jacobian_transpose_matching_grad_fc;
      cctbx::xray::extinction_correction<FloatType> const &exti;
      bool objective_only, compute_grad;
      af::ref<std::complex<FloatType> > f_calc;
      af::ref<FloatType> observables;
      af::ref<FloatType> weights;
      af::versa<FloatType, af::c_grid<2> > &design_matrix;
      accumulate_reflection_chunk(
        int begin, int end,
        boost::shared_ptr<NormalEquations> const &normal_equations_ptr,
        cctbx::xray::observations<FloatType> const &reflections,
        af::const_ref<std::complex<FloatType> > const &f_mask,
        WeightingScheme<FloatType> const &weighting_scheme,
        boost::optional<FloatType> scale_factor,
        boost::shared_ptr<OneMillerIndexFcalc> const &f_calc_function_ptr,
        scitbx::sparse::matrix<FloatType> const
          &jacobian_transpose_matching_grad_fc,
        cctbx::xray::extinction_correction<FloatType> const &exti,
        bool objective_only,
        af::ref<std::complex<FloatType> > f_calc,
        af::ref<FloatType> observables,
        af::ref<FloatType> weights,
        af::versa<FloatType, af::c_grid<2> > &design_matrix)
      : begin(begin), end(end),
        normal_equations_ptr(normal_equations_ptr), normal_equations(*normal_equations_ptr),
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
          af::shared<FloatType> gradients;
          if (compute_grad) {
            gradients.resize(jacobian_transpose_matching_grad_fc.n_rows());
          }
          for (int i_h = begin; i_h < end; ++i_h) {
            miller::index<> const& h = reflections.index(i_h);
            if (f_mask.size()) {
              f_calc_function.compute(h, f_mask[i_h], compute_grad);
            }
            else {
              f_calc_function.compute(h, boost::none, compute_grad);
            }
            f_calc[i_h] = f_calc_function.f_calc;
            if (compute_grad) {
              gradients =
                jacobian_transpose_matching_grad_fc * f_calc_function.grad_observable;
            }
            // sort out twinning
            FloatType observable =
              process_twinning(i_h, gradients);
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
                int grad_index = exti.get_grad_index();
                SMTBX_ASSERT(!(grad_index < 0 || grad_index >= gradients.size()));
                gradients[grad_index] += exti_k[1];
              }
              normal_equations.add_equation(observable,
                gradients.ref(), reflections.fo_sq(i_h), weight);
            }
            if (build_design_matrix) {
              FloatType exti_der = (exti_k[0] + pow(exti_k[0], 3)) / 2;
              for (int i_g = 0; i_g < gradients.size(); i_g++) {
                design_matrix(i_h, i_g) = gradients[i_g];
                if (exti.grad_value()) {
                  if (exti.get_grad_index() == i_g) {
                    design_matrix(i_h, i_g) *= exti_der;
                  }
                }
              }
            }
          }
        }
        catch (smtbx::error const &e) {
          exception_.reset(new smtbx::error(e));
        }
        catch (std::exception const &e) {
          exception_.reset(new smtbx::error(e.what()));
        }
      }

      FloatType process_twinning(int i_h, af::shared<FloatType> &gradients) {
        typedef typename cctbx::xray::observations<FloatType>::iterator itr_t;
        typedef typename cctbx::xray::twin_fraction<FloatType> twf_t;
        typedef typename cctbx::xray::observations<FloatType>::index_twin_component twc_t;
        FloatType obs = f_calc_function.observable;
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
            f_calc_function.compute(twc.h, boost::none, compute_grad);
            obs += twc.scale() * f_calc_function.observable;
            if (compute_grad) {
              af::shared<FloatType> tmp_gradients =
                jacobian_transpose_matching_grad_fc * f_calc_function.grad_observable;
              gradients += twc.scale() * tmp_gradients;
              if (twc.fraction != 0) {
                if (twc.fraction->grad) {
                  SMTBX_ASSERT(!(twc.fraction->grad_index < 0 ||
                    twc.fraction->grad_index >= gradients.size()));
                  gradients[twc.fraction->grad_index] += f_calc_function.observable;
                }
              }
              else {
                identity_part += f_calc_function.observable;
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
    };

#if defined(_OPENMP)
    //The following struct is optimized for cache efficient parallel 
    //calculations using OpenMP making sure there is no strides
    //in adressing memory by scheduling with chunksize 1
    template <class NormalEquations,
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
          std::vector <std::vector<FloatType>> gradients;
          boost::ptr_vector<boost::shared_ptr<OneMillerIndexFcalc>> f_calc_threads;
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
            af::shared<FloatType> gradient;
            gradient.resize(n_rows);
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
#endif

  private:
    struct null_deleter {
      void operator()(void const *) const {}
    };

    af::shared<std::complex<FloatType> > f_calc_;
    af::shared<FloatType> observables_;
    af::shared<FloatType> weights_;

    static int& available_threads_var() {
      static int available = -1;
      return available;
    }
  protected:
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
        bool objective_only = false,
        bool may_parallelise_ = false)
      :
      build_design_matrix_and_normal_equations<FloatType, false>(
        normal_equations,
        reflections, f_mask, weighting_scheme, scale_factor, f_calc_function,
        jacobian_transpose_matching_grad_fc, exti,
        objective_only, may_parallelise_)
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
        bool objective_only = false,
        bool may_parallelise_ = false)
      :
      build_design_matrix_and_normal_equations<FloatType, true>(
        normal_equations,
        reflections, f_mask, weighting_scheme, scale_factor, f_calc_function,
        jacobian_transpose_matching_grad_fc, exti,
        objective_only, may_parallelise_)
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
