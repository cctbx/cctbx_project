#include <cctbx/boost_python/flex_fwd.h>

#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <cctbx/maptbx/bcr/bcr.h>

#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>

namespace cctbx { namespace maptbx { namespace boost_python {

  boost::python::tuple
  getinitargs_bcr_model(bcr_model<> const& self)
  {
    return boost::python::make_tuple(
      self.site_frac,
      self.b,
      self.c,
      self.r);
  }

namespace {

  void init_module()
  {
    using namespace boost::python;

    typedef return_value_policy<return_by_value> rbv;
    class_<bcr_model<> >("bcr_model")
      .def(init<scitbx::vec3<double> const&,
                af::shared<double> const&,
                af::shared<double> const&,
                af::shared<double> const& >((arg("site_frac"),arg("b"),arg("c"),
                                               arg("r"))))
      .add_property("site_frac", make_getter(&bcr_model<>::site_frac, rbv()))
      .add_property("b",         make_getter(&bcr_model<>::b,         rbv()))
      .add_property("c",         make_getter(&bcr_model<>::c,         rbv()))
      .add_property("r",         make_getter(&bcr_model<>::r,         rbv()))
      .enable_pickling()
      .def("__getinitargs__", getinitargs_bcr_model)
    ;

    {
      class_<image<> >("BCRimage", no_init) // why no_init ?
        .def(init<
          cctbx::uctbx::unit_cell const&,
          boost::python::list const&,
          int const& >(
                    (arg("unit_cell"),
                     arg("bcr_models"),
                     arg("step"))))
        .def("fsc",   &image<>::cc)
        .def("d",     &image<>::d)
        .def("d_inv", &image<>::d_inv)
      ;
    }

  }

} // namespace <anonymous>

}}} // namespace cctbx::maptbx::boost_python

BOOST_PYTHON_MODULE(cctbx_maptbx_bcr_bcr_ext)
{
  cctbx::maptbx::boost_python::init_module();
}
