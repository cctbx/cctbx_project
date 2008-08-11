#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/list.hpp>
#include <boost/python/str.hpp>

#include <iotbx/pdb/rna_dna_atom_names.h>
#include <boost/scoped_array.hpp>
#include <vector>

namespace iotbx { namespace pdb { namespace boost_python {

  struct rna_dna_atom_names_info_wrappers
  {
    typedef rna_dna_atom_names::info w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("rna_dna_atom_names_info", no_init)
        .def(init<const char*>((arg_("atom_name"))))
        .def_readonly("reference_name", &w_t::reference_name)
        .def("compatible_residue_names", &w_t::compatible_residue_names)
        .def("is_compatible_with", &w_t::is_compatible_with, (
          arg_("residue_name")))
        .def("is_hydrogen", &w_t::is_hydrogen)
        .def("is_deuterium", &w_t::is_deuterium)
        .def("is_o2prime", &w_t::is_o2prime)
        .def("is_ho2prime", &w_t::is_ho2prime)
        .def("is_h2primeprime", &w_t::is_h2primeprime)
        .def("is_in_phosphate_group", &w_t::is_in_phosphate_group)
        .def("is_op3_or_hop3", &w_t::is_op3_or_hop3)
        .def("is_ho5prime", &w_t::is_ho5prime)
        .def("is_ho3prime", &w_t::is_ho3prime)
        .def("change_h2primeprime_to_ho2prime",
          &w_t::change_h2primeprime_to_ho2prime)
        .def("change_ho5prime_to_hop3", &w_t::change_ho5prime_to_hop3)
        .def("change_to_unknown", &w_t::change_to_unknown)
      ;
    }
  };

  void
  rna_dna_atom_names_interpretation_core(boost::python::object& self)
  {
    bool have_o2prime = false;
    bool have_ho2prime = false;
    std::vector<unsigned> i_atom_h2primeprime;
    bool have_phosphate = false;
    bool have_op3_or_hop3 = false;
    bool have_ho3prime = false;
    std::vector<unsigned> i_atom_ho5prime;
    boost::python::list atom_names(self.attr("atom_names"));
    unsigned n_atoms = static_cast<unsigned>(boost::python::len(atom_names));
    boost::scoped_array<rna_dna_atom_names::info>
      infos_(new rna_dna_atom_names::info[n_atoms]);
    typedef boost::python::extract<const char*> ecc;
    for(unsigned i_atom=0;i_atom<n_atoms;i_atom++) {
      const char* atom_name = ecc(atom_names[i_atom])();
      infos_[i_atom] = rna_dna_atom_names::info(atom_name);
      rna_dna_atom_names::info const& info = infos_[i_atom];
      if (info.is_o2prime()) {
        have_o2prime = true;
      }
      else if (info.is_ho2prime()) {
        have_ho2prime = true;
      }
      else if (info.is_h2primeprime()) {
        i_atom_h2primeprime.push_back(i_atom);
      }
      if (info.is_in_phosphate_group()) {
        have_phosphate = true;
        if (info.is_op3_or_hop3()) {
          have_op3_or_hop3 = true;
        }
      }
      if (info.is_ho5prime()) {
        i_atom_ho5prime.push_back(i_atom);
      }
      if (info.is_ho3prime()) {
        have_ho3prime = true;
      }
    }
    if (have_phosphate) {
      for(unsigned i=0;i<i_atom_ho5prime.size();i++) {
        infos_[i_atom_ho5prime[i]].change_ho5prime_to_hop3();
        have_op3_or_hop3 = true;
      }
    }
    typedef boost::python::str bps;
    bps residue_name_inp(self.attr("residue_name"));
    if (residue_name_inp[0] == "?") {
      if (have_o2prime) {
        self.attr("residue_name") = bps(residue_name_inp[1]);
        if (!have_ho2prime) {
          for(unsigned i=0;i<i_atom_h2primeprime.size();i++) {
            infos_[i_atom_h2primeprime[i]].change_h2primeprime_to_ho2prime();
          }
          have_ho2prime = true;
        }
      }
      else if (i_atom_h2primeprime.size() != 0) {
        self.attr("residue_name") = bps("D" + residue_name_inp[1]);
      }
      else if (have_ho2prime) {
        self.attr("residue_name") = bps(residue_name_inp[1]);
      }
      else {
        self.attr("residue_name") = bps("D" + residue_name_inp[1]);
      }
    }
    const char* residue_name = ecc(self.attr("residue_name"))();
    boost::python::list infos;
    unsigned n_unexpected = 0;
    for(unsigned i_atom=0;i_atom<n_atoms;i_atom++) {
      rna_dna_atom_names::info& info = infos_[i_atom];
      infos.append(info);
      if (info.reference_name == 0 || !info.is_compatible_with(residue_name)) {
        info.change_to_unknown();
        n_unexpected++;
      }
    }
    self.attr("infos") = infos;
    self.attr("have_o2prime") = have_o2prime;
    self.attr("have_ho2prime") = have_ho2prime;
    self.attr("have_phosphate") = have_phosphate;
    self.attr("have_op3_or_hop3") = have_op3_or_hop3;
    self.attr("have_ho3prime") = have_ho3prime;
    self.attr("n_expected") = n_atoms - n_unexpected;
    self.attr("n_unexpected") = n_unexpected;
  }

  void
  wrap_rna_dna_atom_names()
  {
    rna_dna_atom_names_info_wrappers::wrap();
    boost::python::def(
      "rna_dna_atom_names_interpretation_core",
       rna_dna_atom_names_interpretation_core);
  }

}}} // namespace iotbx::pdb::boost_python
