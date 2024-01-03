#ifndef IOTBX_PDB_COMMON_RESIDUE_NAMES_H
#define IOTBX_PDB_COMMON_RESIDUE_NAMES_H

#include <iotbx/pdb/small_str.h>
#include <scitbx/array_family/ref.h>
#include <string>
#include <set>
#include <iotbx/pdb/modified_aa_names.h>
#include <iotbx/pdb/modified_rna_dna_names.h>

namespace iotbx { namespace pdb { namespace common_residue_names {

  static const char* amino_acid[] = {
    "GLY",
    "ALA",
    "VAL",
    "LEU",
    "ILE",
    "MET",
    "MSE",
    "PRO",
    "PHE",
    "TRP",
    "SER",
    "THR",
    "ASN",
    "GLN",
    "TYR",
    "CYS",
    "LYS",
    "ARG",
    "HIS",
    "ASP",
    "GLU",
    0
  };

  static const char* d_amino_acid[] = {
    "DAL", // ALA
    "DAR", // ARG
    "DAS", // ASP
    "DCY", // CYS
    "DGL", // GLU
    "DGN", // GLN
    "DHI", // HIS
    "DIL", // ILE
    "DLE", // LEU
    "DLY", // LYS
    "DPN", // PHE
    "DPR", // PRO
    "DSG", // ASN
    "DSN", // SER
    "DTH", // THR
    "DTR", // TRP
    "DTY", // TYR
    "DVA", // VAL
    "MED", // MET
    0
  };

  static const char* rna_dna[] = {
    "A  ",
    "C  ",
    "G  ",
    "T  ",
    "U  ",
    "  A",
    "  C",
    "  G",
    "  T",
    "  U",
    " A ",
    " C ",
    " G ",
    " T ",
    " U ",
    "+A ",
    "+C ",
    "+G ",
    "+T ",
    "+U ",
    " +A",
    " +C",
    " +G",
    " +T",
    " +U",
    "DA ", // PDB-V3
    "DC ", // PDB-V3
    "DG ", // PDB-V3
    "DT ", // PDB-V3
    " DA", // PDB-V3
    " DC", // PDB-V3
    " DG", // PDB-V3
    " DT", // PDB-V3
    "ADE", // CNS dna-rna.top
    "CYT", // CNS dna-rna.top
    "GUA", // CNS dna-rna.top
    "THY", // CNS dna-rna.top
    "URI", // CNS dna-rna.top
    0
  };

  static const char* ccp4_mon_lib_rna_dna[] = {
    "AD ",
    // "AR ", old RNA names
    "CD ",
    // "CR ",
    "GD ",
    // "GR ",
    "TD ",
    // "UR ",
    " AD",
    // " AR",
    " CD",
    // " CR",
    " GD",
    // " GR",
    " TD",
    // " UR",
    "Ad ",
    // "Ar ",
    "Cd ",
    // "Cr ",
    "Gd ",
    // "Gr ",
    "Td ",
    // "Ur ",
    " Ad",
    // " Ar",
    " Cd",
    // " Cr",
    " Gd",
    // " Gr",
    " Td",
    // " Ur",
    0
  };

  static const char* water[] = {
    "HOH",
    "H2O",
    "OH2",
    "DOD",
    "OD2",
    "WAT",
    0
  };

  static const char* small_molecule[] = {
    "GOL",
    "PO4",
    "SO4",
    0
  };

  static const char* common_saccharide[] = {
    "NAG",
    "NDG",
    "MAN",
    "BMA",
    "FUC",
    "FUL",
    "BGC",
    "GLC",
    0
  };

  /* Survey of PDB as of 2006/05/25.

     The following element names appear as HETATM residue names, but do
     not correspond to isolated atoms or ions:

     Frequency Element HETNAM
     in PDB    symbol
       36        NO    NITROGEN OXIDE
        1        AS    ADENOSINE-5'-THIO-MONOPHOSPHATE
        1        NP    4-HYDROXY-3-NITROPHENYLACETYL-EPSILON-AMINOCAPROIC ACID
        1        SC    CYTIDINE-5'-THIO-MONOPHOSPHATE
        1        PU    PUROMYCIN-N-AMINOPHOSPHONIC ACID
        1        Y     2'-DEOXY-N6-(S)STYRENE OXIDE ADENOSINE MONOPHOSPHATE

     The element names AS, CM, H, I, IN, N, O, P, S appear as residue
     names but not on HETATM records.

     The following element names appear exclusively as HETATM residue
     names for isolated atoms or ions (these are the only names
     included in the element_set below):

     Frequency Element
     in PDB    symbol
      2579       ZN
      2577       CA
      2376       MG
      1581       CL
      1188       NA
       743       MN
       501       K
       386       FE
       319       CU
       312       CD
       215       HG
       214       NI
       165       CO
        62       BR
        47       XE
        41       SR
        33       CS
        31       PT
        25       BA
        23       TL
        23       PB
        19       SM
        19       AU
        15       RB
        14       YB
        13       LI
        12       KR
        10       MO
         7       LU
         7       CR
         6       OS
         6       GD
         4       TB
         4       LA
         4       F
         4       AR
         4       AG
         3       HO
         3       GA
         3       CE
         2       W
         2       SE
         2       RU
         2       RE
         2       PR
         2       IR
         2       EU
         2       AL
         1       V
         1       TE
         1       SB
         1       PD
   */
  static const char* element[] = {
    "ZN ",
    " ZN",
    "CA ",
    " CA",
    "MG ",
    " MG",
    "CL ",
    " CL",
    "NA ",
    " NA",
    "MN ",
    " MN",
    "K  ",
    " K ",
    "  K",
    "FE ",
    " FE",
    "CU ",
    " CU",
    "CD ",
    " CD",
    "HG ",
    " HG",
    "NI ",
    " NI",
    "CO ",
    " CO",
    "BR ",
    " BR",
    "XE ",
    " XE",
    "SR ",
    " SR",
    "CS ",
    " CS",
    "PT ",
    " PT",
    "BA ",
    " BA",
    "TL ",
    " TL",
    "PB ",
    " PB",
    "SM ",
    " SM",
    "AU ",
    " AU",
    "RB ",
    " RB",
    "YB ",
    " YB",
    "LI ",
    " LI",
    "KR ",
    " KR",
    "MO ",
    " MO",
    "LU ",
    " LU",
    "CR ",
    " CR",
    "OS ",
    " OS",
    "GD ",
    " GD",
    "TB ",
    " TB",
    "LA ",
    " LA",
    "F  ",
    " F ",
    "  F",
    "AR ",
    " AR",
    "AG ",
    " AG",
    "HO ",
    " HO",
    "GA ",
    " GA",
    "CE ",
    " CE",
    "W  ",
    " W ",
    "  W",
    "SE ",
    " SE",
    "RU ",
    " RU",
    "RE ",
    " RE",
    "PR ",
    " PR",
    "IR ",
    " IR",
    "EU ",
    " EU",
    "AL ",
    " AL",
    "V  ",
    " V ",
    "  V",
    "TE ",
    " TE",
    "SB ",
    " SB",
    "PD ",
    " PD",
    0
  };

  inline
  void
  initialize_set(
    std::set<std::string>& name_set,
    const char** names)
  {
    if (name_set.size() != 0) return;
    for(const char** n=names; *n; n++) name_set.insert(std::string(*n));
  }

  inline
  const std::set<std::string>&
  amino_acid_set()
  {
    static std::set<std::string> result;
    initialize_set(result, amino_acid);
    return result;
  }

  inline
  const std::set<std::string>&
  d_amino_acid_set()
  {
    static std::set<std::string> result;
    initialize_set(result, d_amino_acid);
    return result;
  }

  inline
  const std::set<std::string>&
  modified_amino_acid_set()
  {
    static std::set<std::string> result;
    initialize_set(result, modified_amino_acid);
    return result;
  }

  inline
  const std::set<std::string>&
  rna_dna_set()
  {
    static std::set<std::string> result;
    initialize_set(result, rna_dna);
    return result;
  }

  inline
  const std::set<std::string>&
  modified_rna_dna_set()
  {
    static std::set<std::string> result;
    initialize_set(result, modified_rna_dna);
    return result;
  }

  inline
  const std::set<std::string>&
  ccp4_mon_lib_rna_dna_set()
  {
    static std::set<std::string> result;
    initialize_set(result, ccp4_mon_lib_rna_dna);
    return result;
  }

  inline
  const std::set<std::string>&
  water_set()
  {
    static std::set<std::string> result;
    initialize_set(result, water);
    return result;
  }

  inline
  const std::set<std::string>&
  small_molecule_set()
  {
    static std::set<std::string> result;
    initialize_set(result, small_molecule);
    return result;
  }

  inline
  const std::set<std::string>&
  saccharide_set()
  {
    static std::set<std::string> result;
    initialize_set(result, common_saccharide);
    return result;
  }

  inline
  const std::set<std::string>&
  element_set()
  {
    static std::set<std::string> result;
    initialize_set(result, element);
    return result;
  }

  inline
  std::string const&
  get_class(std::string const& name, bool consider_ccp4_mon_lib_rna_dna=false)
  {
    static const std::set<std::string>& aa_set = amino_acid_set();
    static const std::set<std::string>& d_aa_set = d_amino_acid_set();
    static const std::set<std::string>& modified_aa_set = modified_amino_acid_set();
    static const std::set<std::string>& na_set = rna_dna_set();
    static const std::set<std::string>& modified_na_set = modified_rna_dna_set();
    static const std::set<std::string>& ml_na_set = ccp4_mon_lib_rna_dna_set();
    static const std::set<std::string>& w_set = water_set();
    static const std::set<std::string>& sm_set = small_molecule_set();
    static const std::set<std::string>& s_set = saccharide_set();
    static const std::set<std::string>& e_set = element_set();
    static const std::string common_amino_acid("common_amino_acid");
    static const std::string d_amino_acid("d_amino_acid");
    static const std::string modified_amino_acid("modified_amino_acid");
    static const std::string common_rna_dna("common_rna_dna");
    static const std::string modified_rna_dna("modified_rna_dna");
    static const std::string ccp4_mon_lib_rna_dna("ccp4_mon_lib_rna_dna");
    static const std::string common_water("common_water");
    static const std::string common_small_molecule("common_small_molecule");
    static const std::string common_saccharide("common_saccharide");
    static const std::string common_element("common_element");
    static const std::string other("other");
    std::string padded(name);
    if (padded.size() < 3) {
      padded.insert(padded.begin(), 3-padded.size(), ' ');
    }
    if (aa_set.find(padded) != aa_set.end()) {
      return common_amino_acid;
    }
    if (d_aa_set.find(padded) != d_aa_set.end()) {
      return d_amino_acid;
    }
    if (modified_aa_set.find(padded) != modified_aa_set.end()) {
      return modified_amino_acid;
    }
    if (na_set.find(padded) != na_set.end()) {
      return common_rna_dna;
    }
    if (modified_na_set.find(padded) != modified_na_set.end()) {
      return modified_rna_dna;
    }
    if (   consider_ccp4_mon_lib_rna_dna
        && ml_na_set.find(padded) != ml_na_set.end()) {
      return ccp4_mon_lib_rna_dna;
    }
    if (w_set.find(padded) != w_set.end()) {
      return common_water;
    }
    if (sm_set.find(padded) != sm_set.end()) {
      return common_small_molecule;
    }
    if (s_set.find(padded) != s_set.end()) {
      return common_saccharide;
    }
    if (e_set.find(padded) != e_set.end()) {
      return common_element;
    }
    return other;
  }

}}} // namespace iotbx::pdb::common_residue_names

#endif // IOTBX_PDB_COMMON_RESIDUE_NAMES_H
