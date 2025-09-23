#pragma once
#include <boost/thread.hpp>
#include <boost/scoped_ptr.hpp>
#include <scitbx/sparse/matrix.h>
#include <smtbx/refinement/least_squares_fc.h>

namespace smtbx {
  namespace ED {
    using namespace refinement::least_squares;

    template <typename FloatType>
    struct build_kin_thread {
      ED_UTIL_TYPEDEFS;
      typedef f_calc_function_base<FloatType> f_calc_function_base_t;
      typedef typename boost::shared_ptr<f_calc_function_base_t> f_calc_function_base_ptr_t;

      build_kin_thread(
        FloatType Fc2Ug,
        const f_calc_function_base_t& f_calc_function,
        const af::shared<miller::index<> >& indices,
        af::shared<complex_t>& Fcs_kin,
        size_t from,
        size_t to)
        : Fc2Ug(Fc2Ug),
        f_calc_function(f_calc_function.fork()),
        indices(indices),
        Fcs_kin(Fcs_kin),
        from(from),
        to(to)
      {}

      void operator ()() {
        try {
          for (size_t i = from; i < to; i++) {
            f_calc_function->compute(indices[i], boost::none, 0, false);
            Fcs_kin[i] = f_calc_function->get_f_calc() * Fc2Ug;
          }
        }
        catch (smtbx::error const& e) {
          exception_.reset(new smtbx::error(e));
        }
        catch (std::exception const& e) {
          exception_.reset(new smtbx::error(e.what()));
        }
      }

      boost::scoped_ptr<smtbx::error> exception_;
      FloatType Fc2Ug;
      f_calc_function_base_ptr_t f_calc_function;
      const af::shared<miller::index<> >& indices;
      af::shared<complex_t>& Fcs_kin;
      size_t from, to;
    };

    template <typename FloatType>
    void build_kin_mt(int thread_n,
      FloatType Fc2Ug,
      f_calc_function_base<FloatType>& f_calc_function,
      const af::shared<miller::index<> >& indices,
      af::shared<std::complex<FloatType> >& Fcs_kin)
    {
      boost::thread_group pool;
      typedef build_kin_thread<FloatType> builder_t;
      typedef typename boost::shared_ptr<builder_t> build_kin_t;
      size_t chunk_sz = indices.size() / thread_n;
      std::vector<build_kin_t> accumulators;
      for (size_t st = 0; st < indices.size(); st += chunk_sz) {
        size_t end = st + chunk_sz;
        build_kin_t pf(
          new builder_t(Fc2Ug,
            f_calc_function,
            indices,
            Fcs_kin,
            st,
            end >= indices.size() ? indices.size() : end)
        );
        accumulators.push_back(pf);
        pool.create_thread(boost::ref(*pf));
      }
      pool.join_all();
      for (size_t i = 0; i < accumulators.size(); i++) {
        if (accumulators[i]->exception_) {
          throw* accumulators[i]->exception_.get();
        }
      }
    }

    template <typename FloatType>
    struct build_kin_thread_ext {
      ED_UTIL_TYPEDEFS;
      typedef f_calc_function_base<FloatType> f_calc_function_base_t;
      typedef typename boost::shared_ptr<f_calc_function_base_t> f_calc_function_base_ptr_t;

      build_kin_thread_ext(
        const scitbx::sparse::matrix<FloatType>& Jt_matching_grad_fc,
        FloatType Fc2Ug,
        const f_calc_function_base_t& f_calc_function,
        const af::shared<miller::index<> >& indices,
        af::shared<complex_t>& Fcs_kin,
        cmat_t& D_kin,
        size_t from,
        size_t to,
        bool compute_grad)
        : Jt_matching_grad_fc(Jt_matching_grad_fc),
        Fc2Ug(Fc2Ug),
        f_calc_function(f_calc_function.fork()),
        indices(indices),
        Fcs_kin(Fcs_kin),
        D_kin(D_kin),
        from(from),
        to(to),
        compute_grad(compute_grad)
      {}

      void operator ()() {
        try {
          for (size_t i = from; i < to; i++) {
            f_calc_function->compute(indices[i], boost::none, 0, compute_grad);
            Fcs_kin[i] = f_calc_function->get_f_calc() * Fc2Ug;
            if (compute_grad && f_calc_function->get_grad_f_calc().size() > 0) {
              size_t col_n = D_kin.accessor().n_columns();
              if (f_calc_function->raw_gradients()) {
                af::shared<complex_t> grads =
                  Jt_matching_grad_fc * f_calc_function->get_grad_f_calc();
                for (size_t j = 0; j < col_n; j++) {
                  D_kin(i, j) = grads[j] * Fc2Ug;
                }
              }
              else {
                af::const_ref<complex_t> grads = f_calc_function->get_grad_f_calc();
                for (size_t j = 0; j < col_n; j++) {
                  D_kin(i, j) = grads[j] * Fc2Ug;
                }
              }
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

      boost::scoped_ptr<smtbx::error> exception_;
      const scitbx::sparse::matrix<FloatType>& Jt_matching_grad_fc;
      FloatType Fc2Ug;
      f_calc_function_base_ptr_t f_calc_function;
      const af::shared<miller::index<> >& indices;
      af::shared<complex_t>& Fcs_kin;
      cmat_t& D_kin;
      size_t from, to;
      bool compute_grad;
    };

    template <typename FloatType>
    void build_kin_mt_ext(int thread_n,
      const scitbx::sparse::matrix<FloatType>& Jt_matching_grad_fc,
      FloatType Fc2Ug,
      f_calc_function_base<FloatType>& f_calc_function,
      const af::shared<miller::index<> >& indices,
      af::shared<std::complex<FloatType> >& Fcs_kin,
      af::versa<std::complex<FloatType>, af::mat_grid>& design_matrix_kin,
      bool compute_grad)
    {
      boost::thread_group pool;
      typedef build_kin_thread_ext<FloatType> builder_t;
      typedef typename boost::shared_ptr<builder_t> build_kin_t;
      size_t chunk_sz = indices.size() / thread_n;
      std::vector<build_kin_t> accumulators;
      for (size_t st = 0; st < indices.size(); st += chunk_sz) {
        size_t end = st + chunk_sz;
        build_kin_t pf(
          new builder_t(Jt_matching_grad_fc,
            Fc2Ug,
            f_calc_function,
            indices,
            Fcs_kin,
            design_matrix_kin,
            st,
            end >= indices.size() ? indices.size() : end,
            compute_grad)
        );
        accumulators.push_back(pf);
        pool.create_thread(boost::ref(*pf));
      }
      pool.join_all();
      for (size_t i = 0; i < accumulators.size(); i++) {
        if (accumulators[i]->exception_) {
          throw* accumulators[i]->exception_.get();
        }
      }
    }
  }
}
