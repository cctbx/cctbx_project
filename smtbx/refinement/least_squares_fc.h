#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_FC_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_FC_H

#include <cctbx/xray/twin_component.h>
#include <cctbx/xray/dispersion_radial.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/sparse/matrix.h>
#include <boost/shared_ptr.hpp>
#include <vector>

namespace smtbx {
  namespace refinement {
    namespace least_squares {
      using namespace cctbx;
      using namespace cctbx::xray;

      /** @brief Three sets of parameters, and the chain rule between them.

      A refinement here never differentiates the observable with respect to the
      quantities it actually refines. It goes through two changes of variable,
      and most of the machinery in this file and the next exists to carry
      gradients back along them.

      Innermost are the *structure factor parameters* c: for each scatterer its
      site, its ADPs, its occupancy, its f' and f''. These are what Fc is a
      function of, and what the one_h linearisation differentiates -- it hands
      back dFc/dc, complex, one entry per c in a fixed layout (see
      standard_xray.h, where the layout is written).

      Next is the *observable*. Fc is complex and the measurement is not; for
      the usual case y = |Fc|^2 the chain rule is

          dy/dc = 2 Re( conj(Fc) dFc/dc )

      which turns a complex gradient into a real one without changing its
      length. get_grad_observable() is the result.

      Outermost are the *refined parameters* x, and they are not the c. Riding
      hydrogens, rigid groups, shared or constrained ADPs, a free variable
      driving several occupancies: each is a statement that some c are functions
      of fewer, or different, x. The reparametrisation knows those functions and
      supplies their Jacobian J = dc/dx. Gradients travel the other way, so what
      is applied per reflection is its transpose:

          dy/dx_k = sum_c (dy/dc)(dc/dx_k),   i.e.  grad_x = J^T grad_c

      J is very sparse -- most c depend on one or two x, and an unconstrained c
      is its own x -- hence the sparse matrix, and hence flattened_jacobian_
      transpose below, which is about making that product cheap rather than
      about the mathematics.

      A functor whose gradients are already in x -- the dynamical electron
      diffraction one, which cannot separate the two changes of variable and so
      does both itself -- says so via raw_gradients() and is not put through
      either step.
      */

      /** @brief The reparametrisation's Jacobian transpose, flattened once so
          that applying it does not walk a sparse structure per reflection.

      scitbx::sparse::matrix is a vector of sparse columns, each column its own
      heap block. Applying it the way it applies itself means, for every
      reflection, stepping through every column object in turn and following
      each to its own storage -- and the columns are many and nearly all of them
      hold one or two entries, so that walk costs more than the arithmetic it
      carries. The structure does not change during a build, only the vector it
      is applied to does, so the walk is done once here and what is left is a
      flat sequential pass.

      The entries are laid down in the order the sparse traversal visits them
      and accumulated in that order, so a column with a repeated row index --
      which is allowed, the columns not being compacted -- still adds up the
      same way, and the result is the one the sparse product gives, bit for bit.

      Kept grouped by column rather than as a flat list of triples so that a
      caller which produces the vector's elements on demand may produce each
      one once, when its column comes up. See apply_converted.
      */
      template <typename FloatType>
      struct flattened_jacobian_transpose {
        struct entry {
          int row;
          FloatType value;
        };
        /// a column with at least one entry, and where its entries end
        struct column {
          int index;
          std::size_t end;
        };

        flattened_jacobian_transpose(
          scitbx::sparse::matrix<FloatType> const &jt)
          : n_rows(jt.n_rows())
        {
          entries.reserve(jt.non_zeroes());
          for (int j = 0; j < jt.n_cols(); j++) {
            typename scitbx::sparse::matrix<FloatType>::const_row_iterator
              p = jt.col(j).begin();
            if (p == jt.col(j).end()) {
              continue;
            }
            for (; p != jt.col(j).end(); ++p) {
              entry e = { static_cast<int>(p.index()), *p };
              entries.push_back(e);
            }
            column c = { j, entries.size() };
            columns.push_back(c);
          }
        }

        /// w = J^T v
        void apply(af::const_ref<FloatType> const &v, FloatType *w) const {
          std::fill(w, w + n_rows, FloatType(0));
          std::size_t at = 0;
          for (std::size_t c = 0; c < columns.size(); c++) {
            FloatType const v_c = v[columns[c].index];
            for (; at < columns[c].end; at++) {
              w[entries[at].row] += entries[at].value*v_c;
            }
          }
        }

        /** @brief w = J^T v, with v[j] handed over by convert(j) rather than
            read out of an array.

        The elements are asked for one per column, in column order, which is
        exactly as many times as an array of them would have been written. So a
        caller whose v is a transformation of something it already holds may
        skip materialising v at all, and the pass which would have written it
        and the pass which would have read it back both disappear.
        */
        template <class Convert>
        void apply_converted(Convert const &convert, FloatType *w) const {
          std::fill(w, w + n_rows, FloatType(0));
          std::size_t at = 0;
          for (std::size_t c = 0; c < columns.size(); c++) {
            FloatType const v_c = convert(columns[c].index);
            for (; at < columns[c].end; at++) {
              w[entries[at].row] += entries[at].value*v_c;
            }
          }
        }

        std::vector<entry> entries;
        std::vector<column> columns;
        int n_rows;
      };

      /* Need inheritance to achive more flexibility */
      template <typename FloatType>
      struct f_calc_function_base {
        typedef std::complex<FloatType> complex_t;
        virtual ~f_calc_function_base() {}

        virtual void compute(
          miller::index<> const& h,
          boost::optional<complex_t > const& f_mask = boost::none,
          twin_fraction<FloatType> const* fraction = 0,
          bool compute_grad = true) = 0;

        void compute(
          miller::index<> const& h,
          twin_fraction<FloatType> const& fraction,
          complex_t const& f_mask,
          bool compute_grad = true)
        {
          compute(h, f_mask, &fraction, compute_grad);
        }

        /// Evaluate the structure factors
        void evaluate(miller::index<> const& h)
        {
          compute(h, boost::none, 0, false);
        }

        /// Evaluate the structure factors
        void evaluate(miller::index<> const& h,
          complex_t const& f_mask)
        {
          compute(h, f_mask, 0, false);
        }

        /// Linearise the structure factors
        void linearise(miller::index<> const& h)
        {
          compute(h, boost::none, 0, true);
        }

        /// Linearise the structure factors
        void linearise(miller::index<> const& h,
          complex_t const& f_mask)
        {
          compute(h, f_mask, 0, true);
        }

        virtual boost::shared_ptr<f_calc_function_base> fork() const = 0;

        virtual FloatType get_observable() const = 0;
        virtual complex_t get_f_calc() const = 0;
        virtual af::const_ref<complex_t > get_grad_f_calc() const = 0;
        virtual af::const_ref<FloatType> get_grad_observable() const = 0;
        /* returns true if grads are for all and not independent only params */
        virtual bool raw_gradients() const { return true; }

        /** @brief Apply the Jacobian to the gradients of the observable
            without ever assembling them.

        The gradients of Fc are complex and the least squares wants gradients of
        the observable, a conversion elementwise in the components. Done in the
        obvious order that is a pass which writes the converted vector and a
        pass which reads it back to apply the Jacobian; done here it is one
        pass, each component converted as its column comes up. The conversion is
        the observable's own, so which observable it is stays where it is known.

        Only meaningful when the functor was told to defer -- otherwise the
        converted vector has already been assembled and there is nothing to
        save. Returns false if this functor has no fused form, and the caller
        should then apply the Jacobian to get_grad_observable() itself.
        */
        virtual bool apply_jacobian_to_grad_observable(
          flattened_jacobian_transpose<FloatType> const &,
          FloatType *) const
        {
          return false;
        }

        /** @brief Stop assembling the gradients of the observable, because the
            caller will take them through apply_jacobian_to_grad_observable.

        Off by default and set once per build, not per reflection. While it is
        on get_grad_observable() is stale and must not be read.
        */
        virtual void set_defer_grad_observable(bool) {}
        /* The radial correction of f' and f'', or null if there is none.

        It is owned by the structure factor functor rather than sitting beside
        the normal equations the way fc_correction does, because its gradients
        are accumulated during the computation of Fc. Whoever reads them has to
        read them from the very functor that computed them, which under
        threading is a fork of the original -- hence going through here rather
        than holding a pointer of one's own.
        */
        virtual cctbx::xray::dispersion_radial_correction<FloatType> const*
        get_dispersion_correction() const
        {
          return 0;
        }
      };

      /** @brief Put the radial f'/f'' correction's gradients in their slots.

      They sit at the tail of the gradient vector next to BASF, EXTI and the
      rest. The sparse Jacobian product which fills the rest of the vector
      leaves them at zero -- they belong to no scatterer -- so this assigns
      rather than adds. The correction has to be the one belonging to this very
      f_calc_function, which under threading is a fork of the original.
      */
      template <typename FloatType>
      void write_dispersion_gradients(
        f_calc_function_base<FloatType> const &f_calc_function,
        af::shared<FloatType> &gradients)
      {
        cctbx::xray::dispersion_radial_correction<FloatType> const *dc =
          f_calc_function.get_dispersion_correction();
        if (dc == 0 || !dc->grad) {
          return;
        }
        af::const_ref<FloatType> dg = dc->get_gradients();
        SMTBX_ASSERT(dc->grad_index >= 0
          && dc->grad_index + dg.size() <= gradients.size());
        for (std::size_t gi = 0; gi < dg.size(); gi++) {
          gradients[dc->grad_index + gi] = dg[gi];
        }
      }

      /* A thin wrapper around the concrete implementation */
      template <typename FloatType,
        class OneMillerIndexFcalc>
      struct f_calc_function_default : public f_calc_function_base<FloatType> {
        typedef std::complex<FloatType> complex_t;
        typedef f_calc_function_base<FloatType> f_calc_function_base_t;

        f_calc_function_default(boost::shared_ptr<OneMillerIndexFcalc> f_calc_function)
          : f_calc_function(f_calc_function)
        {}

        virtual void compute(
          miller::index<> const& h,
          boost::optional<complex_t > const& f_mask = boost::none,
          twin_fraction<FloatType> const* fraction = 0,
          bool compute_grad = true)
        {
          f_calc_function->compute(h, f_mask, compute_grad);
        }

        virtual boost::shared_ptr<f_calc_function_base_t> fork() const {
          return boost::shared_ptr<f_calc_function_base_t>(
            new f_calc_function_default(f_calc_function->fork()));
        }

        virtual FloatType get_observable() const {
          return f_calc_function->get_observable();
        }
        virtual complex_t get_f_calc() const {
          return f_calc_function->f_calc;
        }
        virtual af::const_ref<complex_t> get_grad_f_calc() const {
          return f_calc_function->grad_f_calc.const_ref();
        }
        virtual af::const_ref<FloatType> get_grad_observable() const {
          return f_calc_function->get_grad_observable().const_ref();
        }
        virtual cctbx::xray::dispersion_radial_correction<FloatType> const*
        get_dispersion_correction() const
        {
          return f_calc_function->disp_cr.get();
        }
        /* Forwarded to the linearisation, which is where the observable is
           known and so where the conversion has to happen.
         */
        virtual bool apply_jacobian_to_grad_observable(
          flattened_jacobian_transpose<FloatType> const &jt,
          FloatType *w) const
        {
          return f_calc_function->apply_jacobian_to_grad_observable(jt, w);
        }
        virtual void set_defer_grad_observable(bool defer) {
          f_calc_function->set_defer_grad_observable(defer);
        }

        boost::shared_ptr<OneMillerIndexFcalc> f_calc_function;
      };


      /*  A thin wrapper around concrete implementation to enable caching of
      the results for symmetry related indices.
       */
      template <typename FloatType>
      struct f_calc_function_with_cache : public f_calc_function_base<FloatType>
      {
        typedef std::complex<FloatType> complex_t;
        typedef f_calc_function_base<FloatType> f_calc_function_base_t;
        struct f_calc_function_result {
          f_calc_function_result(
            FloatType const& observable,
            complex_t const& f_calc,
            af::const_ref<complex_t> const& grad_f_calc,
            af::const_ref<FloatType> const& grad_observable)
            :
            observable(observable),
            f_calc(f_calc),
            grad_f_calc(grad_f_calc.begin(), grad_f_calc.end()),
            grad_observable(grad_observable.begin(), grad_observable.end())
          {}

          f_calc_function_result(
            FloatType const& observable,
            complex_t const& f_calc)
            :
            observable(observable),
            f_calc(f_calc),
            grad_observable()
          {}

          FloatType const observable;
          complex_t const f_calc;
          af::shared<complex_t> grad_f_calc;
          af::shared<FloatType> grad_observable;
        };

        f_calc_function_with_cache(
          boost::shared_ptr<f_calc_function_base_t> f_calc_function,
            bool use_cache = false)
          : f_calc_function(f_calc_function),
          use_cache(use_cache),
          length_sq(0)
        {
          /* A cache hit skips compute() altogether, which would leave the
             radial correction's gradient accumulator holding whatever the last
             reflection that missed put there, while the observable came from
             the cache. Caching the correction's gradients alongside the rest
             would fix it; refusing the combination is what is warranted until
             something asks for it, this wrapper being unreachable from Python.
           */
          SMTBX_ASSERT(!use_cache
            || f_calc_function->get_dispersion_correction() == 0);
        }

        virtual void compute(
          miller::index<> const& h,
          boost::optional<complex_t > const& f_mask = boost::none,
          twin_fraction<FloatType> const* fraction = 0,
          bool compute_grad = true)
        {
          if (!use_cache) {
            f_calc_function->compute(h, f_mask, fraction, compute_grad);
            observable = f_calc_function->get_observable();
            grad_f_calc = f_calc_function->get_grad_f_calc();
            grad_observable = f_calc_function->get_grad_observable();
            f_calc = f_calc_function->get_f_calc();
          }
          else {
            FloatType h_length_sq = h.length_sq();
            if (h_length_sq != length_sq) {
              cache.clear();
              length_sq = h_length_sq;
            }
            typename cache_t::iterator iter = cache.find(h);
            if (iter == cache.end()) {
              f_calc_function->compute(h, f_mask, fraction, compute_grad);
              observable = f_calc_function->get_observable();
              grad_f_calc = f_calc_function->get_grad_f_calc();
              grad_observable = f_calc_function->get_grad_observable();
              f_calc = f_calc_function->get_f_calc();
              cache.insert(
                std::pair<miller::index<>, f_calc_function_result>(
                  h, f_calc_function_result(
                    observable,
                    f_calc,
                    grad_f_calc,
                    grad_observable)));
            }
            else {
              observable = iter->second.observable;
              f_calc = iter->second.f_calc;
              grad_f_calc = iter->second.grad_f_calc.const_ref();
              grad_observable = iter->second.grad_observable.const_ref();
            }
          }
        }

        void compute(miller::index<> const& h,
          bool compute_grad = true)
        {
          compute(h, /*f_mask=*/ boost::none, compute_grad);
        }

        virtual boost::shared_ptr<f_calc_function_base_t> fork() const {
          return boost::shared_ptr<f_calc_function_base_t>(
            new f_calc_function_with_cache(f_calc_function->fork(),
              use_cache));
        }

        virtual FloatType get_observable() const {
          return observable;
        }
        virtual complex_t get_f_calc() const {
          return f_calc;
        }
        virtual af::const_ref<complex_t> get_grad_f_calc() const {
          return grad_f_calc;
        }
        virtual af::const_ref<FloatType> get_grad_observable() const {
          return grad_observable;
        }
        virtual cctbx::xray::dispersion_radial_correction<FloatType> const*
        get_dispersion_correction() const
        {
          // safe only because the constructor refuses use_cache with one
          return f_calc_function->get_dispersion_correction();
        }

        typedef std::map<miller::index<>, f_calc_function_result> cache_t;

        boost::shared_ptr<f_calc_function_base_t> f_calc_function;
        FloatType observable;
        af::const_ref<complex_t> grad_f_calc;
        af::const_ref<FloatType> grad_observable;
        complex_t f_calc;
        bool use_cache;
        FloatType length_sq;
        cache_t cache;
      };

    }
  }
}


#endif // GUARD
