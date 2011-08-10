#include <boost/python/class.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/python/pure_virtual.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/tuple.hpp>

#include <smtbx/refinement/restraints/origin_fixing.h>

namespace smtbx { namespace refinement { namespace restraints {
  namespace boost_python {


  template <typename FloatType>
  struct origin_fixing_wrapper
  {
    typedef origin_fixing<FloatType> wt;
    typedef typename wt::scalar_t scalar_t;

    struct origin_fixing_scaffold : wt, boost::python::wrapper<wt>
    {
      origin_fixing_scaffold(sgtbx::space_group const &space_group)
      : wt(space_group)
      {}

      af::shared<scalar_t>
      weights(lstbx::normal_equations::linear_ls<scalar_t> &normal_eqns,
              scitbx::sparse::matrix<scalar_t> const
              &jacobian_transpose_matching_grad_fc,
              af::shared<constraints::scatterer_parameters> const &params)
      {
        return this->get_override("weights")(normal_eqns,
                                             jacobian_transpose_matching_grad_fc,
                                             params);
      }
    };

    static boost::python::tuple singular_directions(wt const &self) {
      using namespace boost::python;
      af::small<af::shared<scalar_t>, 3> const
      &whole_deltas = self.singular_directions();
      switch (whole_deltas.size()) {
        case 3:
          return make_tuple(whole_deltas[0],
                            whole_deltas[1],
                            whole_deltas[2]);
        case 2:
          return make_tuple(whole_deltas[0],
                            whole_deltas[1]);
        case 1:
          return make_tuple(whole_deltas[0]);
        case 0:
        default:
          return boost::python::tuple();

      }
    }

    static void wrap(char const *name) {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      class_<origin_fixing_scaffold, boost::noncopyable>(name, no_init)
        .def(init<sgtbx::space_group const &>())
        .def("add_to", &wt::add_to,
             (arg("normal_equations"),
              arg("jacobian_transpose_matching_grad_fc"),
              arg("scatterer_parameters")))
        .def("weights", pure_virtual(&wt::weights),
             (arg("normal_equations"),
              arg("jacobian_transpose_matching_grad_fc"),
              arg("scatterer_parameters")))
        .add_property("origin_shifts",
                      make_function(&wt::origin_shifts, rbv))
        .add_property("has_floating_directions",  &wt::has_floating_directions)
        .add_property("singular_directions", singular_directions)
        ;
    }
  };

  void wrap_origin_fixing_restraints() {
    origin_fixing_wrapper<double>::wrap("origin_fixing");
  }


}}}}
