#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/dict.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <iotbx/pdb/input.h>

namespace iotbx { namespace pdb {
namespace {

  struct columns_73_76_evaluator_wrappers
  {
    typedef columns_73_76_evaluator w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("columns_73_76_evaluator", no_init)
        .def(init<
          af::const_ref<std::string> const&,
          optional<
            unsigned,
            unsigned> >(
              (arg_("lines"),
               arg_("is_frequent_threshold_atom_records")=1000,
               arg_("is_frequent_threshold_other_records")=100)))
        .def_readonly("finding", &w_t::finding)
        .def_readonly("is_old_style", &w_t::is_old_style)
      ;
    }
  };

  struct input_atom_labels_wrappers
  {
    typedef input_atom_labels w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_internal_reference<> rir;
      class_<w_t>("input_atom_labels", no_init)
        .def("name", &w_t::name)
        .def("resname", &w_t::resname)
        .def("chain", &w_t::chain)
        .def_readonly("resseq", &w_t::resseq)
        .def("icode", &w_t::icode)
        .def("segid", &w_t::segid)
        .def("altloc", &w_t::altloc)
        .def("pdb_format", &w_t::pdb_format)
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t, rir>::wrap(
          "shared_input_atom_labels");
      }
    }
  };

  struct conformer_selections_wrappers
  {
    typedef conformer_selections w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_internal_reference<> rir;
      class_<w_t>("conformer_selections", no_init)
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t, rir>::wrap(
          "shared_conformer_selections");
      }
    }
  };

  struct input_wrappers
  {
    typedef input w_t;

    static
    boost::python::dict
    record_type_counts_as_dict(w_t const& self)
    {
      using namespace boost::python;
      dict result;
      typedef w_t::record_type_counts_t rtct;
      rtct const& rtc = self.record_type_counts();
      for(rtct::const_iterator i=rtc.begin(); i!= rtc.end(); i++) {
        result[i->first.elems] = i->second;
      }
      return result;
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef return_internal_reference<> rir;
      class_<w_t>("input", no_init)
        .def(init<
          const char*,
          af::const_ref<std::string> const&>((
            arg_("source_info"),
            arg_("lines"))))
        //
        .def("source_info", &w_t::source_info, rbv())
        .def("record_type_counts", record_type_counts_as_dict)
        .def("unknown_section", &w_t::unknown_section, rbv())
        .def("title_section", &w_t::title_section, rbv())
        .def("remark_section", &w_t::remark_section, rbv())
        .def("primary_structure_section",
          &w_t::primary_structure_section, rbv())
        .def("heterogen_section", &w_t::heterogen_section, rbv())
        .def("secondary_structure_section",
          &w_t::secondary_structure_section, rbv())
        .def("connectivity_annotation_section",
          &w_t::connectivity_annotation_section, rbv())
        .def("miscellaneous_features_section",
          &w_t::miscellaneous_features_section, rbv())
        .def("crystallographic_section", &w_t::crystallographic_section, rbv())
        .def("input_atom_labels_list", &w_t::input_atom_labels_list, rbv())
        .def("model_numbers", &w_t::model_numbers, rbv())
        .def("model_indices", &w_t::model_indices, rbv())
        .def("ter_indices", &w_t::ter_indices, rbv())
        .def("chain_indices", &w_t::chain_indices, rbv())
        .def("break_indices", &w_t::break_indices, rbv())
        .def("connectivity_section", &w_t::connectivity_section, rbv())
        .def("bookkeeping_section", &w_t::bookkeeping_section, rbv())
        //
        .def("name_selection_cache", &w_t::name_selection_cache, rir())
        .def("altloc_selection_cache", &w_t::altloc_selection_cache, rir())
        .def("resname_selection_cache", &w_t::resname_selection_cache, rir())
        .def("chain_selection_cache", &w_t::chain_selection_cache, rir())
        .def("resseq_selection_cache", &w_t::resseq_selection_cache, rir())
        .def("icode_selection_cache", &w_t::icode_selection_cache, rir())
        .def("segid_selection_cache", &w_t::segid_selection_cache, rir())
        //
        .def("model_numbers_are_unique", &w_t::model_numbers_are_unique)
        .def("model_atom_counts", &w_t::model_atom_counts)
        .def("find_duplicate_atom_labels", &w_t::find_duplicate_atom_labels)
        .def("find_false_blank_altloc", &w_t::find_false_blank_altloc)
        .def("conformer_selections", &w_t::conformer_selections)
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    columns_73_76_evaluator_wrappers::wrap();
    input_atom_labels_wrappers::wrap();
    conformer_selections_wrappers::wrap();
    input_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_input() { wrap_all(); }

}}} // namespace iotbx::pdb::boost_python
