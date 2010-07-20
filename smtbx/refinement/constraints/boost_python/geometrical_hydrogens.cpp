#include <boost/python/class.hpp>
#include <boost/python/implicit.hpp>

#include <smtbx/refinement/constraints/geometrical_hydrogens.h>

namespace smtbx { namespace refinement { namespace constraints {
namespace boost_python {

  struct terminal_linear_ch_site_wrapper
  {
    typedef terminal_linear_ch_site wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt, bases<crystallographic_parameter>,
             std::auto_ptr<wt>,
             boost::noncopyable>("terminal_linear_ch_site", no_init)
        .def(init<site_parameter *,
                  site_parameter *,
                  independent_scalar_parameter *,
                  wt::scatterer_type *>
             ((arg("pivot"), arg("pivot_neighbour"), arg("length"),
               arg("hydrogen"))))
        ;
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  void wrap_geometrical_hydrogens() {
    terminal_linear_ch_site_wrapper::wrap();
  }


}}}}
