#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_H


#include <smtbx/refinement/least_squares_fc.h>
#include <smtbx/import_scitbx_af.h>
#include <scitbx/vec3.h>

namespace smtbx {
  namespace refinement {
    namespace least_squares {

      template <typename FloatType>
      class f_calc_function_ed : public f_calc_function_base<FloatType> {
      public:
        typedef f_calc_function_base<FloatType> f_calc_function_base_t;
        typedef scitbx::vec3<FloatType> cart_t;
        typedef builder_base<FloatType> data_t;
          
        f_calc_function_ed(data_t& data,
          scitbx::mat3<FloatType> const& omat)
          : data(data),
          omat(omat),
          index(-1)
        {
          f_calc = data.f_calc();
          observables = data.observables();
          design_matrix = data.design_matrix();
        }

        virtual void compute(
          miller::index<> const& h,
          boost::optional<std::complex<FloatType> > const& f_mask = boost::none,
          bool compute_grad = true)
        {
          index = -1;
          cctbx::xray::observations<FloatType> const& refs = data.reflections();
          for (int i = 0; i < refs.size(); i++) {
            if (refs.indices()[i] == h) {
              index = i;
              break;
            }
          }
          SMTBX_ASSERT(index >= 0)(index);
        }

        virtual boost::shared_ptr<f_calc_function_base_t> fork() const {
          return boost::shared_ptr<f_calc_function_base_t>(
            new f_calc_function_ed(data, omat));
        }

        virtual FloatType get_observable() const {
          return observables[index];
        }
        virtual std::complex<FloatType> get_f_calc() const {
          return f_calc[index];
        }
        virtual af::const_ref<FloatType> get_grad_observable() const {
          typedef af::versa_plain<FloatType> one_dim_type;
          typedef typename one_dim_type::accessor_type one_dim_accessor_type;
          one_dim_accessor_type a(design_matrix.accessor().n_columns());
          return af::const_ref<FloatType>(&design_matrix(index, 0), a);
        }

        virtual bool raw_gradients() const { return false; }
      private:
        data_t& data;
        scitbx::mat3<FloatType> omat;
        af::shared<std::complex<FloatType> > f_calc;
        af::shared<FloatType> observables;
        af::shared<FloatType> weights;
        af::versa<FloatType, af::c_grid<2> > design_matrix;
        int index;
      };


    }
  }
}


#endif // GUARD
