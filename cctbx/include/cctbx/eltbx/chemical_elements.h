#ifndef CCTBX_ELTBX_CHEMICAL_ELEMENTS_H
#define CCTBX_ELTBX_CHEMICAL_ELEMENTS_H

#include <set>
#include <string>

namespace cctbx { namespace eltbx { namespace chemical_elements {

  static const char* proper_caps[] = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
    "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
    "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
    "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", 0
  };

  static const char* proper_upper[] = {
    "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE", "NA", "MG",
    "AL", "SI", "P", "S", "CL", "AR", "K", "CA", "SC", "TI", "V", "CR",
    "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR",
    "KR", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD",
    "AG", "CD", "IN", "SN", "SB", "TE", "I", "XE", "CS", "BA", "LA",
    "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER",
    "TM", "YB", "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU",
    "HG", "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH",
    "PA", "U", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD",
    "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", 0
  };

  template <typename StringType>
  inline
  void
  initialize_set(
    std::set<StringType>& name_set,
    const char** names)
  {
    if (name_set.size() != 0) return;
    for(const char** n=names; *n; n++) name_set.insert(StringType(*n));
  }

  inline
  std::set<std::string> const&
  proper_caps_set()
  {
    static std::set<std::string> result;
    initialize_set(result, proper_caps);
    return result;
  }

  inline
  std::set<std::string> const&
  proper_upper_set()
  {
    static std::set<std::string> result;
    initialize_set(result, proper_upper);
    return result;
  }

  inline
  std::set<std::string> const&
  proper_and_isotopes_upper_set()
  {
    static std::set<std::string> result;
    if (result.size() == 0) {
      initialize_set(result, proper_upper);
      result.insert("D");
      result.insert("T");
    }
    return result;
  }

}}} // cctbx::eltbx::chemical_elements

#endif // CCTBX_ELTBX_CHEMICAL_ELEMENTS_H
