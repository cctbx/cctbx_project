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
          
        f_calc_function_ed(cctbx::xray::observations<FloatType> const& reflections,
          af::shared<std::complex<FloatType> > const& Fc,
          af::shared<cart_t> const& beams,
          af::versa<FloatType, af::c_grid<2> > const& design_matrix)
          : reflections(reflections),
          Fc(Fc),
          beams(beams),
          design_matrix(design_matrix),
          index(-1)
        {
        }

        virtual void compute(
          miller::index<> const& h,
          boost::optional<std::complex<FloatType> > const& f_mask = boost::none,
          bool compute_grad = true)
        {
          index = -1;
          for (int i = 0; i < reflections.size(); i++) {
            if (reflections.indices()[i] == h) {
              index = i;
              break;
            }
          }
          SMTBX_ASSERT(index >= 0)(index);
        }

        virtual boost::shared_ptr<f_calc_function_base_t> fork() const {
          return boost::shared_ptr<f_calc_function_base_t>(
            new f_calc_function_ed(reflections, Fc, beams, design_matrix));
        }

        virtual FloatType get_observable() const {
          return std::abs(Fc[index]);
        }
        virtual std::complex<FloatType> get_f_calc() const {
          return Fc[index];
        }
        virtual af::const_ref<FloatType> get_grad_observable() const {
          typedef af::versa_plain<FloatType> one_dim_type;
          typedef typename one_dim_type::accessor_type one_dim_accessor_type;
          one_dim_accessor_type a;
          return af::const_ref<FloatType>(&design_matrix[index], a);
        }

      private:
        cctbx::xray::observations<FloatType> const& reflections;
        af::shared<std::complex<FloatType> > const &Fc;
        af::versa<FloatType, af::c_grid<2> > const &design_matrix;
        af::shared<cart_t> const &beams;
        int index;
      };


    }
  }
}


#endif // GUARD
