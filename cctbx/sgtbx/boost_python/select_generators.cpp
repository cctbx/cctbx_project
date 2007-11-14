#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/list.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_by_value.hpp>

#include <cctbx/sgtbx/select_generators.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

  namespace {

    struct any_generator_set_wrapper
    {
      typedef select_generators::any wt;

      static boost::python::tuple generators(
        tr_vec inv_t, rt_mx const *gen, int n_gen)
      {
        using namespace boost::python;
        list result;
        if (inv_t.is_valid()) {
          result.append(rt_mx(rot_mx(1,-1), inv_t));
        }
        for(int i=0; i < n_gen; i++) result.append(gen[i]);
        return tuple(result);
      }

      static boost::python::tuple z_gen(wt const& self) {
        return generators(self.z_inv_t, self.z_gen, self.n_gen);
      }

      static boost::python::tuple p_gen(wt const& self) {
        return generators(self.p_inv_t, self.p_gen, self.n_gen);
      }

      static void wrap() {
        using namespace boost::python;
        typedef return_value_policy<return_by_value> rbv;
        typedef return_value_policy<copy_const_reference> ccr;
        class_<wt>("any_generator_set", no_init)
          .def(init<space_group const&, int, int>((
              arg_("space_group"),
              arg_("z2p_r_den")=cb_r_den,
              arg_("z2p_t_den")=cb_t_den)))
          .def("z_gen", z_gen)
          .def("p_gen", p_gen)
          .def("set_primitive", &wt::set_primitive)
        ;
      }
    };

  } // anonymous namespace

  void wrap_select_generators()
  {
    any_generator_set_wrapper::wrap();
  }

}}} // namespace cctbx::sgtbx::select_generators::boost_python
