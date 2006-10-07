#include <cctbx/boost_python/flex_fwd.h>

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <scitbx/array_family/boost_python/ref_pickle_double_buffered.h>
#include <boost/python/args.hpp>
#include <boost/python/return_internal_reference.hpp>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  struct to_string : pickle_double_buffered::to_string
  {
    using pickle_double_buffered::to_string::operator<<;

    to_string()
    {
      unsigned version = 2;
      *this << version;
    }

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
      *this << val.flags.bits
            << val.flags.param;
      return *this;
    }
  };

  struct from_string : pickle_double_buffered::from_string
  {
    from_string(PyObject* str_obj)
    : pickle_double_buffered::from_string(str_obj)
    {
      *this >> version;
      CCTBX_ASSERT(version == 2);
    }

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
      *this >> val.flags.bits
            >> val.flags.param;
      val.setstate(multiplicity, weight_without_occupancy);
      return *this;
    }

    unsigned version;
  };

}}}} // namespace scitbx::af::boost_python::<anonymous>

namespace cctbx { namespace xray { namespace {

  af::shared<std::string>
  extract_labels(af::const_ref<scatterer<> > const& self)
  {
    af::shared<std::string> result(af::reserve(self.size()));
    for(std::size_t i=0;i<self.size();i++) {
      result.push_back(self[i].label);
    }
    return result;
  }

  af::shared<std::string>
  extract_scattering_types(af::const_ref<scatterer<> > const& self)
  {
    af::shared<std::string> result(af::reserve(self.size()));
    for(std::size_t i=0;i<self.size();i++) {
      result.push_back(self[i].scattering_type);
    }
    return result;
  }

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

  af::shared<double>
  extract_u_iso(
    af::const_ref<scatterer<> > const& self)
  {
    af::shared<double> result(af::reserve(self.size()));
    for(std::size_t i=0;i<self.size();i++) {
      result.push_back(self[i].u_iso);
      //XXX Feature change: don't change & return. Instead, return whatever is
      //                    there.
      //if (!self[i].flags.use_u_iso()) {
      //  result.push_back(-1);
      //}
      //else {
      //  result.push_back(self[i].u_iso);
      //}
    }
    return result;
  }

  af::shared<scitbx::sym_mat3<double> >
  extract_u_cart_or_u_cart_plus_u_iso(
    af::ref<scatterer<> > const& self,
    uctbx::unit_cell const& unit_cell)
  {
    af::shared<scitbx::sym_mat3<double> > result(af::reserve(self.size()));
    for(std::size_t i=0;i<self.size();i++) {
      scitbx::sym_mat3<double> u_cart_=scitbx::sym_mat3<double>(0,0,0,0,0,0);
      if (self[i].flags.use_u_aniso()) {
        CCTBX_ASSERT(
           self[i].u_star != scitbx::sym_mat3<double>(-1,-1,-1,-1,-1,-1));
        u_cart_ += adptbx::u_star_as_u_cart(unit_cell, self[i].u_star);
      }
      if (self[i].flags.use_u_iso()) {
        if(self[i].u_star != scitbx::sym_mat3<double>(-1,-1,-1,-1,-1,-1)) {
          CCTBX_ASSERT(self[i].u_iso >= 0.);
          u_cart_[0] += self[i].u_iso;
          u_cart_[1] += self[i].u_iso;
          u_cart_[2] += self[i].u_iso;
        }
        else {
          u_cart_ = scitbx::sym_mat3<double>(-1,-1,-1,-1,-1,-1);
        }
      }
      if (!self[i].flags.use_u_iso() && !self[i].flags.use_u_aniso()) {
       u_cart_ = scitbx::sym_mat3<double>(-1,-1,-1,-1,-1,-1);
      }
      result.push_back(u_cart_);
    }
    return result;
  }

  af::shared<double>
  extract_u_iso_or_u_equiv(
    af::ref<scatterer<> > const& self,
    uctbx::unit_cell const& unit_cell)
  {
    af::shared<double> result(af::reserve(self.size()));
    for(std::size_t i=0;i<self.size();i++) {
      double u_iso_= 0.0;
      if (self[i].flags.use_u_aniso()) {
        CCTBX_ASSERT(
           self[i].u_star != scitbx::sym_mat3<double>(-1,-1,-1,-1,-1,-1));
        u_iso_ += adptbx::u_star_as_u_iso(unit_cell, self[i].u_star);
      }
      if (self[i].flags.use_u_iso()) {
        CCTBX_ASSERT(self[i].u_iso >= 0.);
        u_iso_ += self[i].u_iso;
      }
      if (!self[i].flags.use_u_iso() && !self[i].flags.use_u_aniso()) {
       u_iso_ = -1.0;
      }
      result.push_back(u_iso_);
    }
    return result;
  }

  void
  set_u_iso(
    af::ref<scatterer<> > const& self,
    af::const_ref<double> const& u_iso)
  {
    CCTBX_ASSERT(self.size() == u_iso.size());
    for(std::size_t i=0;i<self.size();i++) {
      if (self[i].flags.use_u_iso()) {
        self[i].u_iso = u_iso[i];
      }
    }
  }

  void
  adjust_u_iso(
    af::ref<scatterer<> > const& self)
  {
    double u_mean = 0.0;
    double b_min  = 0.5;
    double b_max  = 500.0;
    double min_factor = 8.0;
    double max_factor = 6.0;
    int counter = 0;
    for(std::size_t i=0;i<self.size();i++) {
      if(self[i].flags.use() && self[i].u_iso != -1.0) {
         if(self[i].u_iso < 0.0) self[i].u_iso = 0.0;
         u_mean = u_mean + self[i].u_iso;
         counter = counter + 1;
      }
    }
    u_mean = u_mean / counter;
    double u_min = std::max(adptbx::b_as_u(b_min), u_mean / min_factor);
    u_mean = 0.0;
    counter = 0;
    for(std::size_t i=0;i<self.size();i++) {
      if(self[i].flags.use() && self[i].u_iso != -1.0) {
         if(self[i].u_iso < u_min) self[i].u_iso = u_min;
         u_mean = u_mean + self[i].u_iso;
         counter = counter + 1;
      }
    }
    u_mean = u_mean / counter;
    double u_max = std::min(adptbx::b_as_u(b_max), u_mean * max_factor);
    for(std::size_t i=0;i<self.size();i++) {
      if(self[i].flags.use() && self[i].u_iso != -1.0) {
         if(self[i].u_iso > u_max) self[i].u_iso = u_max;
      }
    }
  }

  af::shared<scitbx::sym_mat3<double> >
  extract_u_star(af::const_ref<scatterer<> > const& self)
  {
    af::shared<scitbx::sym_mat3<double> > result(af::reserve(self.size()));
    for(std::size_t i=0;i<self.size();i++) {
      if (self[i].flags.use_u_aniso()) {
        result.push_back(self[i].u_star);
      }
      else {
        result.push_back(scitbx::sym_mat3<double>(-1,-1,-1,-1,-1,-1));
      }
    }
    return result;
  }

  void
  set_u_star(
    af::ref<scatterer<> > const& self,
    af::const_ref<scitbx::sym_mat3<double> > const& u_star)
  {
    CCTBX_ASSERT(self.size() == u_star.size());
    for(std::size_t i=0;i<self.size();i++) {
      if (self[i].flags.use_u_aniso()) {
        self[i].u_star = u_star[i];
      }
    }
  }

  af::shared<scitbx::sym_mat3<double> >
  extract_u_cart(
    af::const_ref<scatterer<> > const& self,
    uctbx::unit_cell const& unit_cell)
  {
    af::shared<scitbx::sym_mat3<double> > result(af::reserve(self.size()));
    for(std::size_t i=0;i<self.size();i++) {
      if(self[i].u_star != scitbx::sym_mat3<double>(-1,-1,-1,-1,-1,-1)) {
        result.push_back(adptbx::u_star_as_u_cart(unit_cell, self[i].u_star));
      }
      else {
        result.push_back(scitbx::sym_mat3<double>(-1,-1,-1,-1,-1,-1));
      }
      //XXX Feature change: don't change & return. Instead, return whatever is
      //                    there.
      //if (self[i].flags.use_u_aniso()) {
      //  result.push_back(adptbx::u_star_as_u_cart(unit_cell, self[i].u_star));
      //}
      //else {
      //  result.push_back(scitbx::sym_mat3<double>(-1,-1,-1,-1,-1,-1));
      //}
    }
    return result;
  }

  void
  set_u_cart(
    af::ref<scatterer<> > const& self,
    uctbx::unit_cell const& unit_cell,
    af::const_ref<scitbx::sym_mat3<double> > const& u_cart)
  {
    CCTBX_ASSERT(self.size() == u_cart.size());
    for(std::size_t i=0;i<self.size();i++) {
      if (self[i].flags.use_u_aniso()) {
        self[i].u_star = adptbx::u_cart_as_u_star(unit_cell, u_cart[i]);
      }
    }
  }

  af::shared<bool>
  extract_use_u_iso(
    af::const_ref<scatterer<> > const& self)
  {
    af::shared<bool> result(af::reserve(self.size()));
    for(std::size_t i=0;i<self.size();i++) {
      result.push_back(self[i].flags.use_u_iso());
    }
    return result;
  }

  af::shared<bool>
  extract_use_u_aniso(
    af::const_ref<scatterer<> > const& self)
  {
    af::shared<bool> result(af::reserve(self.size()));
    for(std::size_t i=0;i<self.size();i++) {
      result.push_back(self[i].flags.use_u_aniso());
    }
    return result;
  }

  void
  convert_to_isotropic(
    af::ref<scatterer<> > const& self,
    uctbx::unit_cell const& unit_cell)
  {
    for(std::size_t i=0;i<self.size();i++) {
      self[i].convert_to_isotropic(unit_cell);
    }
  }

  void
  convert_to_isotropic(
    af::ref<scatterer<> > const& self,
    uctbx::unit_cell const& unit_cell,
    af::const_ref<std::size_t> const& selection)
  {
    for(std::size_t j=0;j<selection.size();j++) {
      std::size_t i=selection[j];
      self[i].convert_to_isotropic(unit_cell);
    }
  }

  void
  convert_to_anisotropic(
    af::ref<scatterer<> > const& self,
    uctbx::unit_cell const& unit_cell)
  {
    for(std::size_t i=0;i<self.size();i++) {
      self[i].convert_to_anisotropic(unit_cell);
    }
  }

  void
  convert_to_anisotropic(
    af::ref<scatterer<> > const& self,
    uctbx::unit_cell const& unit_cell,
    af::const_ref<bool> const& selection)
  {
    for(std::size_t i=0;i<self.size();i++) {
      if(selection[i]) {
        self[i].convert_to_anisotropic(unit_cell);
      }
    }
  }

  void
  convert_to_anisotropic(
    af::ref<scatterer<> > const& self,
    uctbx::unit_cell const& unit_cell,
    af::const_ref<std::size_t> const& selection)
  {
    for(std::size_t j=0;j<selection.size();j++) {
      std::size_t i=selection[j];
      self[i].convert_to_anisotropic(unit_cell);
    }
  }

  std::size_t
  count_anisotropic(af::const_ref<scatterer<> > const& self)
  {
    std::size_t result = 0;
    for(std::size_t i=0;i<self.size();i++) {
      if (self[i].flags.use_u_aniso()) result++;
    }
    return result;
  }

  std::size_t
  count_anomalous(af::const_ref<scatterer<> > const& self)
  {
    std::size_t result = 0;
    for(std::size_t i=0;i<self.size();i++) {
      if (self[i].fdp != 0) result++;
    }
    return result;
  }

}}} // namespace cctbx::xray::<anonymous>

namespace scitbx { namespace af { namespace boost_python {

  template <>
  struct flex_default_element<cctbx::xray::scatterer<> >
  {
    static cctbx::xray::scatterer<>
    get()
    {
      return cctbx::xray::scatterer<>(
        "", cctbx::fractional<>(0,0,0), 0, 0, "", 0, 0);
    }
  };

  void wrap_flex_xray_scatterer()
  {
    using namespace cctbx;
    using namespace boost::python;

    flex_wrapper<cctbx::xray::scatterer<>,
                 return_internal_reference<>
                >::plain("xray_scatterer")
      .def_pickle(flex_pickle_double_buffered<
        cctbx::xray::scatterer<>, to_string, from_string>())
      .def("extract_labels", cctbx::xray::extract_labels)
      .def("extract_scattering_types", cctbx::xray::extract_scattering_types)
      .def("extract_sites", cctbx::xray::extract_sites)
      .def("set_sites", cctbx::xray::set_sites,
        (arg_("sites")))
      .def("extract_occupancies", cctbx::xray::extract_occupancies)
      .def("set_occupancies", cctbx::xray::set_occupancies,
        (arg_("occupancies")))
      .def("adjust_u_iso", cctbx::xray::adjust_u_iso)
      .def("extract_u_iso", cctbx::xray::extract_u_iso)
      .def("extract_use_u_iso", cctbx::xray::extract_use_u_iso)
      .def("extract_use_u_aniso", cctbx::xray::extract_use_u_aniso)
      .def("extract_u_iso_or_u_equiv", cctbx::xray::extract_u_iso_or_u_equiv,
        (arg_("unit_cell")))
      .def("extract_u_cart_or_u_cart_plus_u_iso",
                              cctbx::xray::extract_u_cart_or_u_cart_plus_u_iso,
        (arg_("unit_cell")))
      .def("set_u_iso", cctbx::xray::set_u_iso,
        (arg_("u_iso")))
      .def("extract_u_star", cctbx::xray::extract_u_star)
      .def("set_u_star", cctbx::xray::set_u_star,
        (arg_("u_star")))
      .def("extract_u_cart", cctbx::xray::extract_u_cart,
        (arg_("unit_cell")))
      .def("set_u_cart", cctbx::xray::set_u_cart,
        (arg_("unit_cell"), arg_("u_cart")))
      .def("convert_to_isotropic", (void(*)(
              af::ref<cctbx::xray::scatterer<> > const&,
              uctbx::unit_cell const&)) cctbx::xray::convert_to_isotropic,
              (arg_("unit_cell")))
      .def("convert_to_isotropic", (void(*)(
              af::ref<cctbx::xray::scatterer<> > const&,
              uctbx::unit_cell const&,
              af::const_ref<std::size_t> const&)) cctbx::xray::convert_to_isotropic,
              (arg_("unit_cell"),arg_("selection")))
      .def("convert_to_anisotropic", (void(*)(
              af::ref<cctbx::xray::scatterer<> > const&,
              uctbx::unit_cell const&)) cctbx::xray::convert_to_anisotropic,
              (arg_("unit_cell")))
      .def("convert_to_anisotropic", (void(*)(
              af::ref<cctbx::xray::scatterer<> > const&,
              uctbx::unit_cell const&,
              af::const_ref<bool> const&)) cctbx::xray::convert_to_anisotropic,
              (arg_("unit_cell"),arg_("selection")))
      .def("convert_to_anisotropic", (void(*)(
              af::ref<cctbx::xray::scatterer<> > const&,
              uctbx::unit_cell const&,
              af::const_ref<std::size_t> const&)) cctbx::xray::convert_to_anisotropic,
              (arg_("unit_cell"),arg_("selection")))
      .def("count_anisotropic", cctbx::xray::count_anisotropic)
      .def("count_anomalous", cctbx::xray::count_anomalous)
    ;
  }

}}} // namespace scitbx::af::boost_python
