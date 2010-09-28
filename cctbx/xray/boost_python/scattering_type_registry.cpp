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

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  struct scattering_type_registry_wrappers
  {
    typedef scattering_type_registry w_t;

    static
    boost::python::dict
    type_index_pairs_as_dict(w_t const& self)
    {
      return scitbx::stl::boost_python::map_as_dict(
        self.type_index_pairs);
    }

    static
    boost::python::list
    unique_gaussians_as_list(w_t const& self)
    {
      return scitbx::stl::boost_python::vector_as_list(
        self.unique_gaussians.const_ref());
    }

    static
    boost::python::tuple
    getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(
        type_index_pairs_as_dict(self),
        unique_gaussians_as_list(self),
        self.unique_counts);
    }

    static
    std::auto_ptr<w_t>
    constructor_for_pickle(
      boost::python::dict const& type_index_pairs,
      boost::python::list const& unique_gaussians,
      w_t::unique_counts_t const& unique_counts)
    {
      std::auto_ptr<w_t> self(new w_t);
      scitbx::stl::boost_python::update_map_from_dict(
        self->type_index_pairs, type_index_pairs);
      scitbx::stl::boost_python::update_vector_from_list(
        self->unique_gaussians, unique_gaussians);
      self->unique_counts = unique_counts;
      CCTBX_ASSERT(self->unique_gaussians.size() \
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
      class_<w_t>("scattering_type_registry")
        .def("type_index_pairs_as_dict", type_index_pairs_as_dict)
        .def("unique_gaussians_as_list", unique_gaussians_as_list)
        .add_property("unique_counts", make_getter(&w_t::unique_counts, rbv()))
        .def("size", &w_t::size)
        .def("has_key", &w_t::has_key, (arg("scattering_type")))
        .def("process",
          (std::size_t(w_t::*)(
            std::string const&)) &w_t::process, (
              arg("scattering_type")))
        .def("process",
          (af::shared<std::size_t>(w_t::*)(
            af::const_ref<scatterer<> > const&)) &w_t::process, (
              arg("scatterers")))
        .def("unique_index", &w_t::unique_index, (arg("scattering_type")))
        .def("unique_indices",
          (af::shared<std::size_t>(w_t::*)(
            af::const_ref<scatterer<> > const&) const) &w_t::unique_indices, (
              arg("scatterers")))
        .def("occupancy_sums", &w_t::occupancy_sums<xray::scatterer<> >,
             arg("scatterers"))
        .def("unit_cell_occupancy_sums",
             &w_t::unit_cell_occupancy_sums<xray::scatterer<> >,
             arg("scatterers"))
        .def("gaussian", &w_t::gaussian, (arg("scattering_type")), ccr())
        .def("gaussian_not_optional",
          &w_t::gaussian_not_optional,
            (arg("scattering_type")), ccr())
        .def("unassigned_types", &w_t::unassigned_types)
        .def("assign", &w_t::assign, (
          arg("scattering_type"), arg("gaussian")))
        .def("assign_from_table", &w_t::assign_from_table, (
          arg("table")))
        .def("unique_form_factors_at_d_star_sq",
          &w_t::unique_form_factors_at_d_star_sq, (
            arg("d_star_sq")))
        .def("dilated_form_factors_at_d_star_sq",
             &w_t::dilated_form_factors_at_d_star_sq,
             (arg("d_star_sq"), arg("dilation_coeffs"),
              arg("unique_indices")))
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
    scattering_type_registry_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
