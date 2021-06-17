#include <cctbx/boost_python/flex_fwd.h>

#include <scitbx/matrix/tensors.h>
#include <cctbx/sgtbx/tensors.h>
#include <scitbx/array_family/versa_matrix.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace cctbx {
  namespace sgtbx {
    namespace boost_python {

      namespace {

        template <class w_t>
        struct tensor_constraints_wrappers {
          static
            af::versa<int, af::c_grid<2> >
            row_echelon_form_as_versa(w_t const& self)
          {
            return af::mat_const_ref_as_versa(self.row_echelon_form());
          }

          static
            af::versa<double, af::c_grid<2> >
            gradient_sum_matrix_as_versa(w_t const& self)
          {
            return af::mat_const_ref_as_versa(self.gradient_sum_matrix());
          }

          static af::versa<int, af::c_grid<2> > get_indices() {
            const std::vector<std::vector<int> > &indices = w_t::get_indices();
            af::versa<int, af::c_grid<2> > rv(af::c_grid<2>(
              indices.size(), indices[0].size()));
            for (size_t i = 0; i < indices.size(); i++) {
              for (size_t j = 0; j < indices[i].size(); j++) {
                rv(i, j) = indices[i][j];
              }
            }
            return rv;
          }

          static void wrap(const char *name) {
            using namespace boost::python;
            typedef return_value_policy<return_by_value> rbv;
            class_<w_t>(name, no_init)
              .def(init<sgtbx::space_group const&, bool>((
                arg("space_group"),
                arg("reciprocal_space"))))
              .def(init<af::shared<rt_mx> const&, std::size_t, bool>((
                arg("symmetry_matrices"),
                arg("i_first_matrix_to_use"),
                arg("reciprocal_space"))))
              .def("row_echelon_form", row_echelon_form_as_versa)
              .add_property("independent_indices",
                make_getter(&w_t::independent_indices, rbv()))
              .def("gradient_sum_matrix", gradient_sum_matrix_as_versa)
              .def("n_independent_params", &w_t::n_independent_params)
              .def("n_dependent_params", &w_t::n_dependent_params)
              .def("independent_params", &w_t::independent_params,
              (arg("all_params")))
              .def("all_params", &w_t::all_params, (arg("independent_params")))
              .def("independent_gradients", &w_t::independent_gradients,
              (arg("all_gradients")))
              .def("independent_curvatures", &w_t::independent_curvatures,
              (arg("all_curvatures")))
              .def("indices", &get_indices, rbv())
              .staticmethod("indices")
              .def("initialise", &w_t::initialise)
              .staticmethod("initialise")
              .def("cleanup", &w_t::cleanup)
              .staticmethod("cleanup")
              ;
          }
        };

      } // namespace <anoymous>

      void wrap_tensor_constraints() {
        using namespace scitbx::matrix::tensors;
        tensor_constraints_wrappers<tensors::constraints<double, tensor_rank_2<double> > >
          ::wrap("rank_2_tensor_constraints");
        tensor_constraints_wrappers<tensors::constraints<double, tensor_rank_3<double> > >
          ::wrap("rank_3_tensor_constraints");
        tensor_constraints_wrappers<tensors::constraints<double, tensor_rank_4<double> > >
          ::wrap("rank_4_tensor_constraints");
      }

    }
  }
} // namespace cctbx::sgtbx::boost_python
