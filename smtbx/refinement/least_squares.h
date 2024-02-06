#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_H

/// Crystallographic least-squares

#include <scitbx/lstbx/normal_equations.h>
#include <scitbx/sparse/matrix.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/matrix/tensors.h>

#include <cctbx/xray/fc_correction.h>
#include <cctbx/xray/observations.h>

#include <smtbx/error.h>
#include <smtbx/structure_factors/direct/standard_xray.h>
#include <smtbx/refinement/least_squares_twinning.h>

#include <algorithm>
#include <vector>
#if defined(_OPENMP)
  #include <omp.h>
#endif
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/smart_ptr/scoped_ptr.hpp>
#include <boost/thread.hpp>

namespace smtbx { namespace refinement { namespace least_squares {

  namespace lstbx = scitbx::lstbx;

  template <typename FloatType>
  class builder_base {
  public:
    virtual ~builder_base() {}
    virtual af::versa<FloatType, af::c_grid<2> > const& design_matrix() const = 0;
    virtual bool has_design_matrix() const = 0;
    virtual cctbx::xray::observations<FloatType> const& reflections() const = 0;
    virtual af::shared<std::complex<FloatType> > const& f_calc() const = 0;
    virtual af::shared<FloatType> const& observables() const = 0;
    virtual af::shared<FloatType> const& weights() const = 0;

    static int get_available_threads() {
      int& available = available_threads_var();
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

    static bool has_openmp() {
#if defined(_OPENMP)
      return true;
#endif
      return false;
    }
  private:
    static int& available_threads_var() {
      static int available = -1;
      return available;
    }
  };


  /** \brief Build normal equations for the given data, model, weighting
  and constraints. Optionally builds the design matrix.

  The constraints is performed with a reparametrisation whose Jacobian
  transpose is passed as an argument.
  */
  template <typename FloatType,
    bool build_design_matrix>
    struct build_design_matrix_and_normal_equations : public builder_base<FloatType> {

    typedef builder_base<FloatType> parent_t;

    typedef f_calc_function_base<FloatType>
      f_calc_function_base_t;

    typedef boost::shared_ptr<f_calc_function_base_t>
      one_miller_index_fcalc_ptr_t;
    typedef boost::shared_ptr<fc_correction<FloatType> >
      fc_correction_ptr_t;

    build_design_matrix_and_normal_equations(
      cctbx::xray::observations<FloatType> const& reflections,
      MaskData<FloatType> const& f_mask_data,
      boost::optional<FloatType> scale_factor,
      f_calc_function_base_t& f_calc_function,
      scitbx::sparse::matrix<FloatType> const&
        jacobian_transpose_matching_grad_fc,
      cctbx::xray::fc_correction<FloatType> const& fc_cr,
      bool objective_only = false,
      bool may_parallelise = false,
      bool use_openmp = false,
      int max_memory = 300)
      :
      reflections_(reflections),
      f_mask_data(f_mask_data),
      scale_factor(scale_factor),
      f_calc_function(f_calc_function),
      jacobian_transpose_matching_grad_fc(jacobian_transpose_matching_grad_fc),
      fc_cr(fc_cr),
      objective_only(objective_only),
      may_parallelise(may_parallelise),
      use_openmp(use_openmp),
      max_memory(max_memory),
      built(false),
      f_calc_(reflections.size()),
      observables_(reflections.size()),
      weights_(reflections.size()),
      design_matrix_(af::c_grid<2>(build_design_matrix ? reflections.size() : 0,
        build_design_matrix ? jacobian_transpose_matching_grad_fc.n_rows() : 0))
    {}

    template<class NormalEquations,
      template<typename> class WeightingScheme>
    build_design_matrix_and_normal_equations(
      NormalEquations& normal_equations,
      cctbx::xray::observations<FloatType> const& reflections,
      MaskData<FloatType> const& f_mask_data,
      WeightingScheme<FloatType> const& weighting_scheme,
      boost::optional<FloatType> scale_factor,
      f_calc_function_base_t &f_calc_function,
      scitbx::sparse::matrix<FloatType> const &
        jacobian_transpose_matching_grad_fc,
      cctbx::xray::fc_correction<FloatType> const& fc_cr,
      bool objective_only = false,
      bool may_parallelise = false,
      bool use_openmp = false,
      int max_memory = 300)
      :
      reflections_(reflections),
      f_mask_data(f_mask_data),
      scale_factor(scale_factor),
      f_calc_function(f_calc_function),
      jacobian_transpose_matching_grad_fc(jacobian_transpose_matching_grad_fc),
      fc_cr(fc_cr),
      objective_only(objective_only),
      may_parallelise(may_parallelise),
      use_openmp(use_openmp),
      max_memory(max_memory),
      built(false),
      f_calc_(reflections.size()),
      observables_(reflections.size()),
      weights_(reflections.size()),
      design_matrix_(af::c_grid<2>(build_design_matrix ? reflections.size() : 0,
        build_design_matrix ? jacobian_transpose_matching_grad_fc.n_rows() : 0))
    {
      build(normal_equations, weighting_scheme);
    }

    template<class NormalEquations,
      template<typename> class WeightingScheme>
    void build(NormalEquations& normal_equations,
      WeightingScheme<FloatType> const& weighting_scheme)
    {
      typedef boost::shared_ptr<NormalEquations>
              normal_equations_ptr_t;
      typedef accumulate_reflection_chunk<
        NormalEquations, WeightingScheme>
        accumulate_reflection_chunk_t;
      typedef boost::shared_ptr<accumulate_reflection_chunk_t>
              accumulate_reflection_chunk_ptr_t;
      if (built) {
        return;
      }
      // Accumulate equations Fo(h) ~ Fc(h)
      reflections_.update_prime_fraction();
      twinning_processor<FloatType> twp(reflections_, f_mask_data, !objective_only,
        jacobian_transpose_matching_grad_fc);
      if (may_parallelise) {
        //!!
        scitbx::matrix::tensors::initialise<FloatType>();
#if defined(_OPENMP)
        if (use_openmp) {
          typedef accumulate_reflection_chunk_omp<
            NormalEquations, WeightingScheme>
            accumulate_reflection_chunk_omp_t;
          /**
           * @brief A pointer to the normal equations object for local refinement.
           */
          normal_equations_ptr_t local_NE(new NormalEquations(normal_equations.n_parameters()));
          accumulate_reflection_chunk_omp_t job(
            local_NE,
            reflections_, f_mask_data, twp, weighting_scheme, scale_factor,
            one_miller_index_fcalc_ptr_t(&f_calc_function, null_deleter()),
            jacobian_transpose_matching_grad_fc,
            fc_cr, objective_only,
            f_calc_.ref(), observables_.ref(), weights_.ref(),
            design_matrix_, max_memory);
          job();
          if (job.exception_) {
            throw* job.exception_.get();
          }
          if (!build_design_matrix) {
            normal_equations = *local_NE;
            normal_equations.finalise(objective_only);
          }
          built = true;
          return;
        }
#endif 
        const int thread_count = parent_t::get_available_threads();
        boost::thread_group pool;
        std::vector<accumulate_reflection_chunk_ptr_t> accumulators;
        Scheduler scheduler(reflections_.size());
        for(int thread_idx=0; thread_idx<thread_count; thread_idx++) {
          normal_equations_ptr_t chunk_normal_equations(
            new NormalEquations(normal_equations.n_parameters()));
          accumulate_reflection_chunk_ptr_t accumulator(
            new accumulate_reflection_chunk_t(
              scheduler,
              chunk_normal_equations,
              reflections_, f_mask_data, twp, weighting_scheme, scale_factor,
              one_miller_index_fcalc_ptr_t(f_calc_function.fork()),
              jacobian_transpose_matching_grad_fc,
              fc_correction_ptr_t(fc_cr.fork()),
              objective_only,
              f_calc_.ref(), observables_.ref(), weights_.ref(),
              design_matrix_));
          accumulators.push_back(accumulator);
          pool.create_thread(boost::ref(*accumulator));
        }
        pool.join_all();
        if (!build_design_matrix) {
          for (int thread_idx = 0; thread_idx < thread_count; thread_idx++) {
            if (accumulators[thread_idx]->exception_) {
              throw* accumulators[thread_idx]->exception_.get();
            }
            normal_equations += accumulators[thread_idx]->normal_equations;
          }
          normal_equations.finalise(objective_only);
        }
      }
      else {
        Scheduler scheduler(reflections_.size());
        accumulate_reflection_chunk_t job(
          scheduler,
          normal_equations_ptr_t(&normal_equations, null_deleter()),
          reflections_, f_mask_data, twp, weighting_scheme, scale_factor,
          one_miller_index_fcalc_ptr_t(f_calc_function.fork()),
          jacobian_transpose_matching_grad_fc,
          fc_correction_ptr_t(fc_cr.fork()),
          objective_only,
          f_calc_.ref(), observables_.ref(), weights_.ref(),
          design_matrix_);
        job();
        if (job.exception_) {
          throw *job.exception_.get();
        }
        if (!build_design_matrix) {
          normal_equations.finalise(objective_only);
        }
      }
      built = true;
    }

    virtual cctbx::xray::observations<FloatType> const& reflections() const {
      return reflections_;
    }

    virtual af::shared<std::complex<FloatType> > const& f_calc() const { return f_calc_; }

    virtual af::shared<FloatType> const& observables() const { return observables_; }

    virtual af::shared<FloatType> const& weights() const { return weights_; }

    virtual bool has_design_matrix() const {
      return build_design_matrix && built;
    }

  protected:
    af::versa<FloatType, af::c_grid<2> > const& design_matrix() const { return design_matrix_; }
#if defined(_OPENMP)
    #include "least_squares_omp.h"
#endif
    struct chunk {
      const int idx, size;
      chunk()
        : idx(0), size(0)
      {}
      chunk(int idx, int size)
        : idx(idx), size(size)
      {}
    };
    struct Scheduler {
      int count, current;
      boost::mutex mtx;
      Scheduler(int count)
        : count(count),
        current(0)
      {}
      chunk next() {
        boost::mutex::scoped_lock lock(mtx);
        int left = count - current;
        if (left == 0) {
          return chunk();
        }
        int sz = std::min(left, 256);
        chunk rv(current, sz);
        current += sz;
        return rv;
      }
      void reset() {
        current = 0;
      }
    };

    /// Accumulate from reflections whose indices are
    /// returned by scheduler
    template<class NormalEquations,
      template<typename> class WeightingScheme>
    struct accumulate_reflection_chunk {
      Scheduler& scheduler;
      boost::scoped_ptr<smtbx::error> exception_;
      boost::shared_ptr<NormalEquations> normal_equations_ptr;
      NormalEquations &normal_equations;
      cctbx::xray::observations<FloatType> const &reflections;
      MaskData<FloatType> const& f_mask_data;
      twinning_processor<FloatType> const& twp;
      WeightingScheme<FloatType> const &weighting_scheme;
      boost::optional<FloatType> scale_factor;
      boost::shared_ptr<f_calc_function_base_t> f_calc_function_ptr;
      f_calc_function_base_t &f_calc_function;
      scitbx::sparse::matrix<FloatType> const
        &jacobian_transpose_matching_grad_fc;
      boost::shared_ptr<cctbx::xray::fc_correction<FloatType> > fc_cr;
      bool objective_only, compute_grad;
      af::ref<std::complex<FloatType> > f_calc;
      af::ref<FloatType> observables;
      af::ref<FloatType> weights;
      af::versa<FloatType, af::c_grid<2> > &design_matrix;
      accumulate_reflection_chunk(
        Scheduler& scheduler,
        boost::shared_ptr<NormalEquations> const& normal_equations_ptr,
        cctbx::xray::observations<FloatType> const &reflections,
        MaskData<FloatType> const& f_mask_data,
        twinning_processor<FloatType> const& twp,
        WeightingScheme<FloatType> const &weighting_scheme,
        boost::optional<FloatType> scale_factor,
        boost::shared_ptr<f_calc_function_base_t> const &f_calc_function_ptr,
        scitbx::sparse::matrix<FloatType> const
          &jacobian_transpose_matching_grad_fc,
        boost::shared_ptr<cctbx::xray::fc_correction<FloatType> > const &fc_cr,
        bool objective_only,
        af::ref<std::complex<FloatType> > f_calc,
        af::ref<FloatType> observables,
        af::ref<FloatType> weights,
        af::versa<FloatType, af::c_grid<2> > &design_matrix)
      : scheduler(scheduler),
        normal_equations_ptr(normal_equations_ptr), normal_equations(*normal_equations_ptr),
        reflections(reflections), f_mask_data(f_mask_data), twp(twp),
        weighting_scheme(weighting_scheme),
        scale_factor(scale_factor),
        f_calc_function_ptr(f_calc_function_ptr), f_calc_function(*f_calc_function_ptr),
        jacobian_transpose_matching_grad_fc(jacobian_transpose_matching_grad_fc),
        fc_cr(fc_cr),
        objective_only(objective_only), compute_grad(!objective_only),
        f_calc(f_calc), observables(observables), weights(weights),
        design_matrix(design_matrix)
      {}

      void operator()() {
        try {
          af::shared<FloatType> gradients;
          int n_params = jacobian_transpose_matching_grad_fc.n_rows();
          if (compute_grad) {
            gradients.resize(n_params);
          }
          while (true) {
            chunk ch = scheduler.next();
            if (ch.size == 0) {
              break;
            }
            for (int i = 0; i < ch.size; i++) {
              int i_h = ch.idx + i;
              miller::index<> const& h = reflections.index(i_h);
              const twin_fraction<FloatType>* fraction = reflections.fraction(i_h);
              if (f_mask_data.size()) {
                f_calc_function.compute(h, f_mask_data.find(h), fraction, compute_grad);
              }
              else {
                f_calc_function.compute(h, boost::none, fraction, compute_grad);
              }
              f_calc[i_h] = f_calc_function.get_f_calc();
              if (compute_grad) {
                if (f_calc_function.raw_gradients()) {
                  gradients =
                    jacobian_transpose_matching_grad_fc * f_calc_function.get_grad_observable();
                }
                else {
                  gradients = af::shared<FloatType>(
                    f_calc_function.get_grad_observable().begin(),
                    f_calc_function.get_grad_observable().end());
                }
              }
              // sort out twinning
              FloatType observable = twp.process(
                i_h, f_calc_function, gradients);
              // Fc correction
              FloatType fc_k = fc_cr->compute(h, observable, compute_grad);
              if (fc_k != 1) {
                observable *= fc_k;
                f_calc[i_h] *= std::sqrt(fc_k);
              }
              observables[i_h] = observable;

              FloatType weight = weighting_scheme(reflections.fo_sq(i_h),
                reflections.sig(i_h), observable, scale_factor);
              weights[i_h] = weight;
              if (objective_only) {
                normal_equations.add_residual(observable,
                  reflections.fo_sq(i_h), weight);
              }
              else {
                if (fc_cr->grad) {
                  int grad_idx = fc_cr->get_grad_index();
                  af::const_ref<FloatType> fc_cr_grads = fc_cr->get_gradients();
                  SMTBX_ASSERT(grad_idx >= 0 &&
                    grad_idx+fc_cr_grads.size() <= gradients.size());
                  FloatType grad_m = fc_cr->get_grad_Fc_multiplier();
                  if (grad_m != 1) {
                    for (int gi = 0; gi < gradients.size(); gi++) {
                     gradients[gi] *= grad_m;
                    }
                  }
                  for (int gi = 0; gi < fc_cr_grads.size(); gi++) {
                    gradients[grad_idx + gi] = fc_cr_grads[gi];
                  }
                }
                if (!build_design_matrix) {
                  normal_equations.add_equation(observable,
                    gradients.ref(), reflections.fo_sq(i_h), weight);
                }
              }
              if (build_design_matrix) {
                memcpy(&design_matrix(i_h, 0), gradients.begin(),
                  gradients.size() * sizeof(FloatType));
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
    };


  private:
    struct null_deleter {
      void operator()(void const *) const {}
    };

  protected:
    cctbx::xray::observations<FloatType> const& reflections_;
    MaskData<FloatType> const& f_mask_data;
    boost::optional<FloatType> scale_factor;
    f_calc_function_base_t& f_calc_function;
    scitbx::sparse::matrix<FloatType> const&
      jacobian_transpose_matching_grad_fc;
    cctbx::xray::fc_correction<FloatType> const& fc_cr;
    bool objective_only,
      may_parallelise,
      use_openmp,
      built;
    int max_memory;

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
    typedef build_design_matrix_and_normal_equations<FloatType, false> parent_t;
    build_normal_equations(
      cctbx::xray::observations<FloatType> const& reflections,
      MaskData<FloatType> const& f_mask_data,
      boost::optional<FloatType> scale_factor,
      f_calc_function_base<FloatType>& f_calc_function,
      scitbx::sparse::matrix<FloatType> const
      & jacobian_transpose_matching_grad_fc,
      cctbx::xray::fc_correction<FloatType> const& fc_cr,
      bool objective_only = false,
      bool may_parallelise = false,
      bool use_openmp = false)
      : parent_t(
        reflections, f_mask_data, scale_factor, f_calc_function,
        jacobian_transpose_matching_grad_fc, fc_cr,
        objective_only, may_parallelise, use_openmp)
    {}

    template<class NormalEquations,
      template<typename> class WeightingScheme>
    build_normal_equations(
       NormalEquations &normal_equations,
       cctbx::xray::observations<FloatType> const &reflections,
       MaskData<FloatType> const& f_mask_data,
       WeightingScheme<FloatType> const &weighting_scheme,
       boost::optional<FloatType> scale_factor,
       f_calc_function_base<FloatType> &f_calc_function,
       scitbx::sparse::matrix<FloatType> const
       &jacobian_transpose_matching_grad_fc,
       cctbx::xray::fc_correction<FloatType> const &fc_cr,
       bool objective_only = false,
       bool may_parallelise = false,
       bool use_openmp = false,
       int max_memory = 300)
       : parent_t(
        normal_equations,
        reflections, f_mask_data, weighting_scheme, scale_factor, f_calc_function,
        jacobian_transpose_matching_grad_fc, fc_cr,
        objective_only, may_parallelise, use_openmp, max_memory)
    {}
     virtual af::versa<FloatType, af::c_grid<2> > const& design_matrix() const {
       SMTBX_NOT_IMPLEMENTED();
       return parent_t::design_matrix_;
     }
  };

  /** \brief Build only thed esign matrix for the given data, model, weighting
  and constraints and the buld the design matrix

  The constraints is performed with a reparametrisation whose Jacobian
  transpose is passed as an argument.
  */
  
  template <typename FloatType>
  struct build_design_matrix
    : public build_design_matrix_and_normal_equations<FloatType, true>
  {
    typedef build_design_matrix_and_normal_equations<FloatType, true> parent_t;
    build_design_matrix(
      cctbx::xray::observations<FloatType> const& reflections,
      MaskData<FloatType> const& f_mask_data,
      boost::optional<FloatType> scale_factor,
      f_calc_function_base<FloatType>& f_calc_function,
      scitbx::sparse::matrix<FloatType> const
      & jacobian_transpose_matching_grad_fc,
      cctbx::xray::fc_correction<FloatType> const& fc_cr,
      bool objective_only = false,
      bool may_parallelise = false,
      bool use_openmp = false)
      : parent_t(
        reflections, f_mask_data, scale_factor, f_calc_function,
        jacobian_transpose_matching_grad_fc, fc_cr,
        objective_only, may_parallelise, use_openmp)
    {}

    template<class NormalEquations,
      template<typename> class WeightingScheme>
    build_design_matrix(
       NormalEquations &normal_equations,
       cctbx::xray::observations<FloatType> const &reflections,
      MaskData<FloatType> const& f_mask_data,
       WeightingScheme<FloatType> const &weighting_scheme,
       boost::optional<FloatType> scale_factor,
       f_calc_function_base<FloatType> &f_calc_function,
       scitbx::sparse::matrix<FloatType> const
       &jacobian_transpose_matching_grad_fc,
       cctbx::xray::fc_correction<FloatType> const &fc_cr,
       bool objective_only = false,
       bool may_parallelise = false,
       bool use_openmp = false,
       int max_memory = 300)
       : parent_t(
        normal_equations,
        reflections, f_mask_data, weighting_scheme, scale_factor, f_calc_function,
        jacobian_transpose_matching_grad_fc, fc_cr,
        objective_only, may_parallelise, use_openmp, max_memory)
    {}

    virtual af::versa<FloatType, af::c_grid<2> > const& design_matrix() const {
      return parent_t::design_matrix_;
    }
 
  };


}}}


#endif // GUARD
