#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_by_value.hpp>

#include <scitbx/stl/vector_wrapper.h>
#include <cctbx/xray/observations.h>

namespace cctbx { namespace xray { namespace boost_python {
namespace {

  template <typename FloatType>
  struct observations_wrapper {

    typedef observations<FloatType> obst;
    static boost::python::tuple detwin(obst const& self,
      cctbx::sgtbx::space_group const& space_group,
      bool anomalous_flag,
      scitbx::af::const_ref<cctbx::miller::index<> > const& fo_sq_indices,
      scitbx::af::const_ref<FloatType> const& fc_sqs)
    {
      scitbx::af::tiny<scitbx::af::shared<FloatType>, 2> result =
        self.detwin(space_group, anomalous_flag, fo_sq_indices, fc_sqs);
      return boost::python::make_tuple(result[0], result[1]);
    }

    typedef typename obst::index_twin_component itct;
    static twin_fraction<FloatType> get_twin_fractions(itct const& self) {
      return *self.fraction;
    }

    static void wrap() {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;

      class_<obst>("observations", no_init)
        .def(init<scitbx::af::shared<cctbx::miller::index<> > const&,
                  scitbx::af::shared<FloatType> const&,
                  scitbx::af::shared<FloatType> const&,
                  scitbx::af::shared<
                    cctbx::xray::twin_component<FloatType>*> const& >
             ((arg("indices"),
               arg("data"),
               arg("sigmas"),
               arg("merohedral_components"))))
        .def(init<scitbx::af::shared<cctbx::miller::index<> > const&,
                  scitbx::af::shared<FloatType> const&,
                  scitbx::af::shared<FloatType> const&,
                  scitbx::af::shared<int> const&,
                  scitbx::af::shared<cctbx::xray::twin_fraction<FloatType>*> const&,
                  scitbx::af::shared<
                    cctbx::xray::twin_component<FloatType>*> const& >
             ((arg("indices"),
               arg("data"),
               arg("sigmas"),
               arg("scale_indices"),
               arg("twin_fractions"),
               arg("merohedral_components"))))
        .def(init<observations<FloatType> const&,
                  scitbx::af::shared<
                    cctbx::xray::twin_fraction<FloatType>*> const&,
                  scitbx::af::shared<
                    cctbx::xray::twin_component<FloatType>*> const& >
             ((arg("observations"),
               arg("twin_fractions"),
               arg("merohedral_components"))))
        .add_property("indices", &obst::indices)
        .add_property("data", &obst::data)
        .add_property("sigmas", &obst::sigmas)
        .add_property("twin_fractions", &obst::twin_fractions)
        .add_property("merohedral_components", &obst::merohedral_components)
        .add_property("measured_scale_indices", &obst::measured_scale_indices)
        .def("iterator", &obst::iterator)
        .def("detwin", detwin)
        ;

      typedef typename obst::iterator_ itrt;
      class_<itrt>("iterator", no_init)
        .def("has_next", &itrt::has_next)
        .def("next", &itrt::next)
        ;


      typedef typename obst::filter_result frt;
      class_<frt>("filter_result", no_init)
        .def_readonly("omitted_count", &frt::omitted_count)
        .def_readonly("sys_abs_count", &frt::sys_abs_count)
        .add_property("selection", &frt::get_selection)
        ;

      class_<itct>("index_twin_component", no_init)
        .def(init<cctbx::miller::index<> const&,
                  cctbx::xray::twin_fraction<FloatType> const*>
             ((arg("index"),
               arg("fraction"))))
        .add_property("h", make_getter(&itct::h, rbv()))
        .add_property("fraction", &get_twin_fractions)
        ;

      typedef observations<FloatType>::filter ft;
      class_<ft>("filter", no_init)
        .def(init<uctbx::unit_cell const& ,
                  sgtbx::space_group const&,
                  bool,
                  scitbx::af::const_ref<miller::index<> > const&,
                  FloatType, FloatType, FloatType>
             ((arg("unit_cell"),
               arg("space_group"),
               arg("anomalous_flag"),
               arg("omit_indices"),
               arg("res_min_d"),
               arg("res_max_d"),
               arg("min_i_o_sig"))))
      .def("is_to_omit", &ft::is_to_omit)
      ;
    }
  };
} // namespace anonymous

  void wrap_observations() {
    using namespace boost::python;
    observations_wrapper<double>::wrap();
    def("filter_data",
        observations<double>::filter_data,
        (arg("indices"), arg("data"), arg("sigmas"),
         arg("scale_indices"),
         arg("filter")));
  }

}}} //end of cctbx::xray::boost_python
