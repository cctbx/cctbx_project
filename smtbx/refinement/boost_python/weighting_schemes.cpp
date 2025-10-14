#include <smtbx/refinement/weighting_schemes.h>

#include <boost/python/class.hpp>


namespace smtbx { namespace refinement { namespace least_squares {
  using namespace cctbx;
  using namespace cctbx::uctbx;

  namespace boost_python {
  
  template<template<typename FloatType> class WeightingScheme>
  struct weighting_scheme_class
    : boost::python::class_<WeightingScheme<double>,
        boost::python::bases<IWeightingScheme<double> > >
  {
    typedef WeightingScheme<double> wt;
    typedef boost::python::class_<wt,
      boost::python::bases<IWeightingScheme<double> > > base_t;

    static af::shared<double> weights(wt const &weighting_scheme,
                                      af::const_ref<double> const &fo_sq,
                                      af::const_ref<double> const &sigmas,
                                      af::const_ref<double> const &fc_sq,
                                      af::const_ref<miller::index<> > const& indices,
                                      double scale_factor)
    {
      return least_squares::weights(
        weighting_scheme, fo_sq, sigmas, fc_sq, indices, scale_factor);
    }

    weighting_scheme_class(char const *name)
      : base_t(name, boost::python::no_init)
    {
      using namespace boost::python;
      this->def("__call__", &wt::operator(),
          (arg("fo_sq"), arg("sigma"), arg("fc_sq"), arg("h"), arg("scale_factor")));
      this->def("compute", &wt::compute,
          (arg("fo_sq"), arg("sigma"), arg("fc_sq"), arg("h"), arg("scale_factor")));
      this->def("__call__", weights,
          (arg("fo_sq"), arg("sigmas"), arg("fc_sq"), arg("indices"), arg("scale_factor")));
      //this->def("__call__", weights,
      //    (arg("fo_sq"), arg("sigmas"), arg("fc_sq"), arg("scale_factor")));
    }
  };

  struct iweighting_scheme_wrapper {
    static void wrap() {
      using namespace boost::python;
      typedef IWeightingScheme<double> wt;
      class_<wt, boost::noncopyable>("weighting_scheme", no_init)
        .def("compute", &wt::compute)
        ;
    }
  };

  struct mainstream_shelx_weighting_wrapper
  {
    static void wrap() {
      using namespace boost::python;
      typedef weighting_scheme_class<mainstream_shelx_weighting> wt;
      wt("mainstream_shelx_weighting")
        .def(init<optional<double, double> >((arg("a"), arg("b"))))
        .def_readwrite("a", &wt::wt::a)
        .def_readwrite("b", &wt::wt::b)
        ;
    }
  };

  struct new_shelx_weighting_wrapper
  {
    static void wrap() {
      using namespace boost::python;
      typedef weighting_scheme_class<new_shelx_weighting> wt;
      wt("new_shelx_weighting")
        .def(init<optional<double, double> >((arg("a"), arg("b"))))
        .def_readwrite("a", &wt::wt::a)
        .def_readwrite("b", &wt::wt::b)
        ;
    }
  };

  struct stl_weighting_wrapper
  {
    static void wrap() {
      using namespace boost::python;
      typedef weighting_scheme_class<stl_weighting> wt;
      wt("stl_weighting")
        .def(init<unit_cell const&, optional<double> >((arg("unit_cell"), arg("a"))))
        .def_readwrite("a", &wt::wt::a)
        ;
    }
  };

  struct unit_weighting_wrapper
  {
    static void wrap() {
      using namespace boost::python;
      weighting_scheme_class<unit_weighting>("unit_weighting")
        .def(init<>())
        ;
    }
  };

  struct sigma_weighting_wrapper
  {
    static void wrap() {
      using namespace boost::python;
      weighting_scheme_class<sigma_weighting>("sigma_weighting")
        .def(init<>())
        ;
    }
  };

  void wrap_weighting_schemes() {
    iweighting_scheme_wrapper::wrap();
    mainstream_shelx_weighting_wrapper::wrap();
    new_shelx_weighting_wrapper::wrap();
    unit_weighting_wrapper::wrap();
    sigma_weighting_wrapper::wrap();
    stl_weighting_wrapper::wrap();
  }


}}}}
