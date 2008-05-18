#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/scattering_type_registry.h>
#include <scitbx/stl/map_as_python_dict.h>
#include <scitbx/stl/vector_as_python_list.h>
#include <boost/python/class.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/args.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_by_value.hpp>

#if defined(__sgi) && !defined(__GNUC__)
// see comments in scitbx/scitbx/source_generators/flex_fwd_h.py
typedef boost::optional<scitbx::math::gaussian::sum<double> > b_o_g_s_d;
#endif

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  template<class FormFactorType>
  struct scattering_type_registry_traits;

  template<>
  struct scattering_type_registry_traits<eltbx::xray_scattering::gaussian>
    : xray::form_factor_traits<eltbx::xray_scattering::gaussian>
  {
    static
    std::string class_name() {
      return "scattering_type_registry";
    }
  };


  template<class FormFactorType>
  struct scattering_type_registry_wrappers
  {
    typedef generic_scattering_type_registry<FormFactorType> w_t;
    typedef scattering_type_registry_traits<FormFactorType> traits;

    static
    boost::python::dict
    type_index_pairs_as_dict(w_t const& self)
    {
      return scitbx::stl::boost_python::map_as_dict(
        self.type_index_pairs);
    }

    static
    boost::python::list
    unique_form_factors_as_list(w_t const& self)
    {
      return scitbx::stl::boost_python::vector_as_list(
        self.unique_form_factors.const_ref());
    }

    static
    boost::python::tuple
    getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(
        type_index_pairs_as_dict(self),
        unique_form_factors_as_list(self),
        self.unique_counts);
    }

    static
    std::auto_ptr<w_t>
    constructor_for_pickle(
      boost::python::dict const& type_index_pairs,
      boost::python::list const& unique_form_factors,
      typename w_t::unique_counts_t const& unique_counts)
    {
      std::auto_ptr<w_t> self(new w_t);
      scitbx::stl::boost_python::update_map_from_dict(
        self->type_index_pairs, type_index_pairs);
      scitbx::stl::boost_python::update_vector_from_list(
        self->unique_form_factors, unique_form_factors);
      self->unique_counts = unique_counts;
      CCTBX_ASSERT(self->unique_form_factors.size() \
                == self->type_index_pairs.size());
      CCTBX_ASSERT(self->unique_counts.size() \
                == self->type_index_pairs.size());
      return self;
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_value_policy<return_by_value> rbv;
      std::string ff_name = traits::form_factor_name();

      std::string class_name = traits::class_name(); // necessary to work
        // around a bug in some flavours of gcc 3.2 and 3.3
      class_<w_t>(class_name.c_str())
        .def("type_index_pairs_as_dict", type_index_pairs_as_dict)
        .def("unique_form_factors_as_list", unique_form_factors_as_list)
        .def((std::string("unique_")
              + ff_name
              + std::string("s_as_list")).c_str(),
             unique_form_factors_as_list)
        .add_property("unique_counts", make_getter(&w_t::unique_counts, rbv()))
        .def("size", &w_t::size)
        .def("has_key", &w_t::has_key, (arg_("scattering_type")))
        .def("process",
          (std::size_t(w_t::*)(
            std::string const&)) &w_t::process, (
              arg_("scattering_type")))
        .def("process",
          (af::shared<std::size_t>(w_t::*)(
            af::const_ref<scatterer<> > const&)) &w_t::process, (
              arg_("scatterers")))
        .def("unique_index", &w_t::unique_index, (arg_("scattering_type")))
        .def("unique_indices",
          (af::shared<std::size_t>(w_t::*)(
            af::const_ref<scatterer<> > const&) const) &w_t::unique_indices, (
              arg_("scatterers")))
        .def("form_factor",
             &w_t::form_factor,
             (arg_("scattering_type")),
             ccr())
        .def("form_factor_not_optional",
             &w_t::form_factor_not_optional,
             (arg_("scattering_type")),
             ccr())
        .def(ff_name.c_str(),
             &w_t::form_factor,
             (arg_("scattering_type")),
             ccr())
        .def((ff_name + std::string("_not_optional")).c_str(),
              &w_t::form_factor_not_optional,
              (arg_("scattering_type")),
              ccr())
        .def("unassigned_types", &w_t::unassigned_types)
        .def("assign",
             &w_t::assign,
             (arg_("scattering_type"),
              arg_(ff_name.c_str())))
        .def("assign_from_table", &w_t::assign_from_table, (
          arg_("table")))
        .def("unique_form_factors_at_d_star_sq",
          &w_t::unique_form_factors_at_d_star_sq, (
            arg_("d_star_sq")))
        .def("dilated_form_factors_at_d_star_sq",
             &w_t::dilated_form_factors_at_d_star_sq,
             (arg_("d_star_sq"),
              arg_("dilation_coeffs"),
              arg_("unique_indices")))
        .enable_pickling()
        .def("__getinitargs__", getinitargs)
        .def("__init__", make_constructor(constructor_for_pickle),
          "constructor for pickle");
      ;
    }
  };

} // namespace <anoymous>

  void wrap_scattering_type_registry()
  {
    scattering_type_registry_wrappers<eltbx::xray_scattering::gaussian>::wrap();
  }

}}} // namespace cctbx::xray::boost_python
