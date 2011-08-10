#include <smtbx/refinement/restraints/origin_fixing.h>

namespace smtbx { namespace refinement { namespace restraints {
  namespace boost_python {

  template <typename FloatType>
  struct floating_origin_restraints_wrapper
  {
    typedef floating_origin_restraints<FloatType> wt;
    typedef typename wt::scalar_t scalar_t;

    static boost::python::tuple singular_directions(wt const &self) {
      using namespace boost::python;
      switch (self.singular_directions.size()) {
        case 3:
          return make_tuple(self.singular_directions[0],
                            self.singular_directions[1],
                            self.singular_directions[2]);
        case 2:
          return make_tuple(self.singular_directions[0],
                            self.singular_directions[1]);
        case 1:
          return make_tuple(self.singular_directions[0]);
        case 0:
        default:
          return boost::python::tuple();

      }
    }

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def(init<sgtbx::space_group const &,
                  af::const_ref<constraints::scatterer_parameters> const &,
                  scalar_t>
             ((arg("space_group"),
               arg("all_scatterer_parameters"),
               arg("floating_origin_restraint_relative_weight"))))
        .def("add_to", &wt::add_to,
             (arg("normal_equations"),
              arg("jacobian_transpose_matching_grad_fc")))
        .add_property("singular_directions", singular_directions)
        ;
    }
  };

  void wrap_origin_fixing_restraints() {
    floating_origin_restraints_wrapper<
      double>::wrap("floating_origin_restraints");

  }


}}}}
