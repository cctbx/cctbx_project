#include <cctbx/boost_python/flex_fwd.h>

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <scitbx/array_family/boost_python/ref_pickle_double_buffered.h>
#include <boost/python/return_internal_reference.hpp>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  struct to_string : pickle_double_buffered::to_string
  {
    using pickle_double_buffered::to_string::operator<<;

    to_string& operator<<(cctbx::xray::scatterer<> const& val)
    {
      *this << val.label
            << val.scattering_type
            << val.fp
            << val.fdp;
      *this << val.site.const_ref()
            << val.occupancy
            << val.anisotropic_flag
            << val.u_iso;
      *this << val.u_star.const_ref()
            << val.multiplicity()
            << val.weight_without_occupancy();
      return *this;
    }
  };

  struct from_string : pickle_double_buffered::from_string
  {
    from_string(PyObject* str_obj)
    : pickle_double_buffered::from_string(str_obj)
    {}

    using pickle_double_buffered::from_string::operator>>;

    from_string& operator>>(cctbx::xray::scatterer<>& val)
    {
      int multiplicity;
      cctbx::xray::scatterer<>::float_type weight_without_occupancy;
      *this >> val.label
            >> val.scattering_type
            >> val.fp
            >> val.fdp;
      *this >> val.site.ref()
            >> val.occupancy
            >> val.anisotropic_flag
            >> val.u_iso;
      *this >> val.u_star.ref()
            >> multiplicity
            >> weight_without_occupancy;
      val.setstate(multiplicity, weight_without_occupancy);
      return *this;
    }
  };

}}}} // namespace scitbx::af::boost_python::<anonymous>

namespace cctbx { namespace xray { namespace {

  af::shared<scitbx::vec3<double> >
  extract_sites(af::const_ref<scatterer<> > const& scatterers)
  {
    af::shared<scitbx::vec3<double> >
      result(af::reserve(scatterers.size()));
    for(std::size_t i=0;i<scatterers.size();i++) {
      result.push_back(scatterers[i].site);
    }
    return result;
  }

  void
  set_sites(
    af::ref<scatterer<> > const& scatterers,
    af::const_ref<scitbx::vec3<double> > const& sites)
  {
    CCTBX_ASSERT(scatterers.size() == sites.size());
    for(std::size_t i=0;i<scatterers.size();i++) {
      scatterers[i].site = sites[i];
    }
  }

  af::shared<double>
  extract_occupancies(af::const_ref<scatterer<> > const& scatterers)
  {
    af::shared<double>
      result(af::reserve(scatterers.size()));
    for(std::size_t i=0;i<scatterers.size();i++) {
      result.push_back(scatterers[i].occupancy);
    }
    return result;
  }

  void
  set_occupancies(
    af::ref<scatterer<> > const& scatterers,
    af::const_ref<double> const& occupancies)
  {
    CCTBX_ASSERT(scatterers.size() == occupancies.size());
    for(std::size_t i=0;i<scatterers.size();i++) {
      scatterers[i].occupancy = occupancies[i];
    }
  }

  std::size_t
  count_anisotropic(af::const_ref<scatterer<> > const& scatterers)
  {
    std::size_t result = 0;
    for(std::size_t i=0;i<scatterers.size();i++) {
      if (scatterers[i].anisotropic_flag) result++;
    }
    return result;
  }

}}} // namespace cctbx::xray::<anonymous>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_xray_scatterer()
  {
    using namespace cctbx;

    flex_wrapper<cctbx::xray::scatterer<>,
                 boost::python::return_internal_reference<>
                >::plain("xray_scatterer")
      .def_pickle(flex_pickle_double_buffered<
        cctbx::xray::scatterer<>, to_string, from_string>())
      .def("extract_sites", cctbx::xray::extract_sites)
      .def("set_sites", cctbx::xray::set_sites)
      .def("extract_occupancies", cctbx::xray::extract_occupancies)
      .def("set_occupancies", cctbx::xray::set_occupancies)
      .def("count_anisotropic", cctbx::xray::count_anisotropic)
    ;
  }

}}} // namespace scitbx::af::boost_python
