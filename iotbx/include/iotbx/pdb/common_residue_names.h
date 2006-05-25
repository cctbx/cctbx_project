#ifndef IOTBX_PDB_COMMON_RESIDUE_NAMES_H
#define IOTBX_PDB_COMMON_RESIDUE_NAMES_H

#include <iotbx/pdb/small_str.h>
#include <scitbx/array_family/ref.h>
#include <string>
#include <set>

namespace iotbx { namespace pdb { namespace common_residue_names {

  static const char* amino_acid[] = {
    "GLY ",
    "ALA ",
    "VAL ",
    "LEU ",
    "ILE ",
    "MET ",
    "MSE ",
    "PRO ",
    "PHE ",
    "TRP ",
    "SER ",
    "THR ",
    "ASN ",
    "GLN ",
    "TYR ",
    "CYS ",
    "LYS ",
    "ARG ",
    "HIS ",
    "ASP ",
    "GLU ",
    0
  };

  static const char* rna_dna[] = {
    "A   ",
    "C   ",
    "G   ",
    "T   ",
    "U   ",
    "  A ",
    "  C ",
    "  G ",
    "  T ",
    "  U ",
    "+A  ",
    "+C  ",
    "+G  ",
    "+T  ",
    "+U  ",
    " +A ",
    " +C ",
    " +G ",
    " +T ",
    " +U ",
    "ADE ", // CNS dna-rna.top
    "CYT ", // CNS dna-rna.top
    "GUA ", // CNS dna-rna.top
    "THY ", // CNS dna-rna.top
    "URI ", // CNS dna-rna.top
    "CMP ", // CNS dna-rna.top
    0
  };

  static const char* water[] = {
    "HOH ",
    "H2O ",
    "OH2 ",
    "DOD ",
    "D2O ",
    "OD2 ",
    "WAT ",
    "TIP ",
    "TIP3",
    0
  };

  inline
  void
  initialize_set(
    std::set<str4>& name_set,
    const char** names)
  {
    if (name_set.size() != 0) return;
    for(const char** n=names; *n; n++) name_set.insert(str4(*n));
  }

  inline
  void
  initialize_set(
    std::set<str4>& name_set,
    scitbx::af::const_ref<std::string> const& names)
  {
    for(std::size_t i=0;i<names.size();i++) {
      std::string const& name = names[i];
      if (name.size() > 4) {
        throw std::runtime_error(
          "residue name with more than 4 characters: \""+name+"\"");
      }
      name_set.insert(str4(
        name.data(), name.size(), /*i_begin*/ 0, /*pad_with*/ ' '));
    }
  }

  inline
  const std::set<str4>&
  amino_acid_set()
  {
    static std::set<str4> result;
    initialize_set(result, amino_acid);
    return result;
  }

  inline
  const std::set<str4>&
  rna_dna_set()
  {
    static std::set<str4> result;
    initialize_set(result, rna_dna);
    return result;
  }

  inline
  const std::set<str4>&
  water_set()
  {
    static std::set<str4> result;
    initialize_set(result, water);
    return result;
  }

  inline
  std::string const&
  get_class(str4 const& name)
  {
    static const std::set<str4>& aa_set = amino_acid_set();
    static const std::set<str4>& na_set = rna_dna_set();
    static const std::set<str4>& w_set = water_set();
    static const std::string common_amino_acid("common_amino_acid");
    static const std::string common_rna_dna("common_rna_dna");
    static const std::string common_water("common_water");
    static const std::string other("other");
    if (aa_set.find(name) != aa_set.end()) {
      return common_amino_acid;
    }
    if (na_set.find(name) != na_set.end()) {
      return common_rna_dna;
    }
    if (w_set.find(name) != w_set.end()) {
      return common_water;
    }
    return other;
  }

}}} // namespace iotbx::pdb::common_residue_names

#endif // IOTBX_PDB_COMMON_RESIDUE_NAMES_H
