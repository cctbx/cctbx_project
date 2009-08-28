#ifndef IOTBX_PDB_RNA_DNA_ATOM_NAMES_H
#define IOTBX_PDB_RNA_DNA_ATOM_NAMES_H

#include <string>
#include <cctype>

namespace iotbx { namespace pdb {

//! Highly optimized dictionary of RNA/DNA atom name aliases.
namespace rna_dna_atom_names {

  //! Constants.
  namespace info_flags {

    static const unsigned none            = 0x00000000U;
    static const unsigned a               = 0x00000001U;
    static const unsigned c               = 0x00000002U;
    static const unsigned g               = 0x00000004U;
    static const unsigned u               = 0x00000008U;
    static const unsigned da              = 0x00000010U;
    static const unsigned dc              = 0x00000020U;
    static const unsigned dg              = 0x00000040U;
    static const unsigned dt              = 0x00000080U;
    static const unsigned any_bit         = 0x00000100U;
    static const unsigned any             = 0x000001ffU;
    static const unsigned hydrogen        = 0x00000200U;
    static const unsigned deuterium       = 0x00000400U;
    static const unsigned o2prime         = 0x00000800U;
    static const unsigned ho2prime        = 0x00001000U;
    static const unsigned h2primeprime    = 0x00002000U;
    static const unsigned phosphate_group = 0x00004000U;
    static const unsigned op3_or_hop3     = 0x00008000U;
    static const unsigned ho5prime        = 0x00010000U;
    static const unsigned ho3prime        = 0x00020000U;

  } // namespace info_flags

  //! Translates atom name alias to reference name.
  struct info
  {
    const char* reference_name;
    unsigned flags;

    info() {}

    info(const char* atom_name)
    :
      reference_name(0),
      flags(info_flags::none)
    {
      using namespace info_flags;
      if (atom_name == 0) return;
      while (*atom_name && std::isspace(*atom_name)) atom_name++;
      switch (atom_name[0])
      {
        case '1':
          if (atom_name[1] == 'D' || atom_name[1] == 'd') {
            flags |= deuterium;
          }
          else if (atom_name[1] != 'H' && atom_name[1] != 'h') {
            break;
          }
          flags |= hydrogen;
          switch (atom_name[2])
          {
            case '2':
              if (rest_is_whitespace(&atom_name[3])) {
                reference_name = " H21";
                flags |= g | dg;
                return;
              }
              if (atom_name[3] == '\'' || atom_name[3] == '*') {
                if (rest_is_whitespace(&atom_name[4])) {
                  reference_name = " H2'";
                  flags |= any;
                  return;
                }
              }
              break;

            case '4':
              if (!rest_is_whitespace(&atom_name[3])) break;
              reference_name = " H41";
              flags |= c | dc;
              return;

            case '5':
              if (atom_name[3] == '\'' || atom_name[3] == '*') {
                if (rest_is_whitespace(&atom_name[4])) {
                  reference_name = " H5'";
                  flags |= any;
                  return;
                }
                break;
              }
              if (atom_name[3] == 'M' || atom_name[3] == 'm') {
                if (rest_is_whitespace(&atom_name[4])) {
                  reference_name = " H71";
                  flags |= dt;
                  return;
                }
              }
              break;

            case '6':
              if (!rest_is_whitespace(&atom_name[3])) break;
              reference_name = " H61";
              flags |= a | da;
              return;

            default:
              break;
          }
          flags = none;
          break;

        case '2':
          if (atom_name[1] == 'D' || atom_name[1] == 'd') {
            flags |= deuterium;
          }
          else if (atom_name[1] != 'H' && atom_name[1] != 'h') {
            break;
          }
          flags |= hydrogen;
          switch (atom_name[2])
          {
            case '2':
              if (rest_is_whitespace(&atom_name[3])) {
                reference_name = " H22";
                flags |= g | dg;
                return;
              }
              if (atom_name[3] == '\'' || atom_name[3] == '*') {
                if (rest_is_whitespace(&atom_name[4])) {
                  reference_name = "H2''";
                  flags |= da | dc | dg | dt | h2primeprime;
                  return;
                }
              }
              break;

            case '4':
              if (!rest_is_whitespace(&atom_name[3])) break;
              reference_name = " H42";
              flags |= c | dc;
              return;

            case '5':
              if (atom_name[3] == '\'' || atom_name[3] == '*') {
                if (rest_is_whitespace(&atom_name[4])) {
                  reference_name = "H5''";
                  flags |= any;
                  return;
                }
                break;
              }
              if (atom_name[3] == 'M' || atom_name[3] == 'm') {
                if (rest_is_whitespace(&atom_name[4])) {
                  reference_name = " H72";
                  flags |= dt;
                  return;
                }
              }
              break;

            case '6':
              if (!rest_is_whitespace(&atom_name[3])) break;
              reference_name = " H62";
              flags |= a | da;
              return;

            case 'O':
            case 'o':
              if (atom_name[3] == '\'' || atom_name[3] == '*') {
                if (rest_is_whitespace(&atom_name[4])) {
                  reference_name = "HO2'";
                  flags |= a | c | g | u | ho2prime;
                  return;
                }
                break;
              }
              if (atom_name[3] == 'P' || atom_name[3] == 'p') {
                if (rest_is_whitespace(&atom_name[4])) {
                  reference_name = "HOP2";
                  flags |= any | phosphate_group;
                  return;
                }
              }
              break;

            default:
              break;
          }
          flags = none;
          break;

        case '3':
          if (atom_name[1] == 'D' || atom_name[1] == 'd') {
            flags |= deuterium;
          }
          else if (atom_name[1] != 'H' && atom_name[1] != 'h') {
            break;
          }
          flags |= hydrogen;
          switch (atom_name[2])
          {
            case '5':
              if (atom_name[3] == 'M' || atom_name[3] == 'm') {
                if (rest_is_whitespace(&atom_name[4])) {
                  reference_name = " H73";
                  flags |= dt;
                  return;
                }
              }
              break;

            case 'O':
            case 'o':
              if (atom_name[3] == 'P' || atom_name[3] == 'p') {
                if (rest_is_whitespace(&atom_name[4])) {
                  reference_name = "HOP3";
                  flags |= any | phosphate_group | op3_or_hop3;
                  return;
                }
              }
              break;

            default:
              break;
          }
          flags = none;
          break;

        case 'C':
        case 'c':
          switch (atom_name[1])
          {
            case '1':
              if (atom_name[2] == '\'' || atom_name[2] == '*') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " C1'";
                  flags |= any;
                  return;
                }
              }
              break;

            case '2':
              if (rest_is_whitespace(&atom_name[2])) {
                reference_name = " C2 ";
                flags |= any;
                return;
              }
              if (atom_name[2] == '\'' || atom_name[2] == '*') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " C2'";
                  flags |= any;
                  return;
                }
              }
              break;

            case '3':
              if (atom_name[2] == '\'' || atom_name[2] == '*') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " C3'";
                  flags |= any;
                  return;
                }
              }
              break;

            case '4':
              if (rest_is_whitespace(&atom_name[2])) {
                reference_name = " C4 ";
                flags |= any;
                return;
              }
              if (atom_name[2] == '\'' || atom_name[2] == '*') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " C4'";
                  flags |= any;
                  return;
                }
              }
              break;

            case '5':
              if (rest_is_whitespace(&atom_name[2])) {
                reference_name = " C5 ";
                flags |= any;
                return;
              }
              if (atom_name[2] == '\'' || atom_name[2] == '*') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " C5'";
                  flags |= any;
                  return;
                }
                break;
              }
              if (   atom_name[2] == 'M' || atom_name[2] == 'm'
                  || atom_name[2] == 'A' || atom_name[2] == 'a') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " C7 ";
                  flags |= dt;
                  return;
                }
              }
              break;

            case '6':
              if (!rest_is_whitespace(&atom_name[2])) break;
              reference_name = " C6 ";
              flags |= any;
              return;

            case '7':
              if (!rest_is_whitespace(&atom_name[2])) break;
              reference_name = " C7 ";
              flags |= dt;
              return;

            case '8':
              if (!rest_is_whitespace(&atom_name[2])) break;
              reference_name = " C8 ";
              flags |= a | g | da | dg;
              return;

            default:
              break;
          }
          break;

        case 'D':
        case 'd':
          flags |= deuterium;
        case 'H':
        case 'h':
          flags |= hydrogen;
          switch (atom_name[1])
          {
            case '1':
              if (rest_is_whitespace(&atom_name[2])) {
                reference_name = " H1 ";
                flags |= g | dg;
                return;
              }
              if (atom_name[2] == '\'' || atom_name[2] == '*') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " H1'";
                  flags |= any;
                  return;
                }
              }
              break;

            case '2':
              if (rest_is_whitespace(&atom_name[2])) {
                reference_name = " H2 ";
                flags |= a | da;
                return;
              }
              if (atom_name[2] == '\'' || atom_name[2] == '*') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " H2'";
                  flags |= any;
                  return;
                }
                if (atom_name[3] == '1') {
                  if (rest_is_whitespace(&atom_name[4])) {
                    reference_name = " H2'";
                    flags |= any;
                    return;
                  }
                  break;
                }
                if (atom_name[3] == '\'' || atom_name[3] == '2') {
                  if (rest_is_whitespace(&atom_name[4])) {
                    reference_name = "H2''";
                    flags |= da | dc | dg | dt | h2primeprime;
                    return;
                  }
                }
                break;
              }
              if (atom_name[2] == '1') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " H21";
                  flags |= g | dg;
                  return;
                }
                break;
              }
              if (atom_name[2] == '2') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " H22";
                  flags |= g | dg;
                  return;
                }
              }
              break;

            case '3':
              if (rest_is_whitespace(&atom_name[2])) {
                reference_name = " H3 ";
                flags |= u | dt;
                return;
              }
              if (atom_name[2] == '\'' || atom_name[2] == '*') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " H3'";
                  flags |= any;
                  return;
                }
                break;
              }
              if (atom_name[2] == 'T' || atom_name[2] == 't') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = "HO3'";
                  flags |= any | ho3prime;
                  return;
                }
              }
              break;

            case '4':
              if (atom_name[2] == '\'' || atom_name[2] == '*') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " H4'";
                  flags |= any;
                  return;
                }
                break;
              }
              if (atom_name[2] == '1') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " H41";
                  flags |= c | dc;
                  return;
                }
                break;
              }
              if (atom_name[2] == '2') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " H42";
                  flags |= c | dc;
                  return;
                }
              }
              break;

            case '5':
              if (rest_is_whitespace(&atom_name[2])) {
                reference_name = " H5 ";
                flags |= c | u | dc;
                return;
              }
              if (atom_name[2] == '\'') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " H5'";
                  flags |= any;
                  return;
                }
                if (atom_name[3] == '1') {
                  if (rest_is_whitespace(&atom_name[4])) {
                    reference_name = " H5'";
                    flags |= any;
                    return;
                  }
                  break;
                }
                if (atom_name[3] == '2' || atom_name[3] == '\'') {
                  if (rest_is_whitespace(&atom_name[4])) {
                    reference_name = "H5''";
                    flags |= any;
                    return;
                  }
                }
                break;
              }
              if (atom_name[2] == '*') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = "HO5'";
                  flags |= any | ho5prime;
                  return;
                }
                if (atom_name[3] == '1') {
                  if (rest_is_whitespace(&atom_name[4])) {
                    reference_name = " H5'";
                    flags |= any;
                    return;
                  }
                  break;
                }
                if (atom_name[3] == '2') {
                  if (rest_is_whitespace(&atom_name[4])) {
                    reference_name = "H5''";
                    flags |= any;
                    return;
                  }
                }
                break;
              }
              if (atom_name[2] == 'M' || atom_name[2] == 'm') {
                if (atom_name[3] == '1') {
                  if (rest_is_whitespace(&atom_name[4])) {
                    reference_name = " H71";
                    flags |= dt;
                    return;
                  }
                  break;
                }
                if (atom_name[3] == '2') {
                  if (rest_is_whitespace(&atom_name[4])) {
                    reference_name = " H72";
                    flags |= dt;
                    return;
                  }
                  break;
                }
                if (atom_name[3] == '3') {
                  if (rest_is_whitespace(&atom_name[4])) {
                    reference_name = " H73";
                    flags |= dt;
                    return;
                  }
                }
                break;
              }
              if (atom_name[2] == 'T' || atom_name[2] == 't') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = "HO5'";
                  flags |= any | ho5prime;
                  return;
                }
              }
              break;

            case '6':
              if (rest_is_whitespace(&atom_name[2])) {
                reference_name = " H6 ";
                flags |= c | u | dc | dt;
                return;
              }
              if (atom_name[2] == '1') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " H61";
                  flags |= a | da;
                  return;
                }
                break;
              }
              if (atom_name[2] == '2') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " H62";
                  flags |= a | da;
                  return;
                }
              }
              break;

            case '7':
              if (atom_name[2] == '1') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " H71";
                  flags |= dt;
                  return;
                }
                break;
              }
              if (atom_name[2] == '2') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " H72";
                  flags |= dt;
                  return;
                }
                break;
              }
              if (atom_name[2] == '3') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " H73";
                  flags |= dt;
                  return;
                }
              }
              break;

            case '8':
              if (!rest_is_whitespace(&atom_name[2])) break;
              reference_name = " H8 ";
              flags |= a | g | da | dg;
              return;

            case 'O':
            case 'o':
              if (atom_name[2] == '2') {
                if (atom_name[3] == '\'' || atom_name[3] == '*') {
                  if (rest_is_whitespace(&atom_name[4])) {
                    reference_name = "HO2'";
                    flags |= a | c | g | u | ho2prime;
                    return;
                  }
                }
                break;
              }
              if (atom_name[2] == '3') {
                if (atom_name[3] == '\'' || atom_name[3] == '*') {
                  if (rest_is_whitespace(&atom_name[4])) {
                    reference_name = "HO3'";
                    flags |= any | ho3prime;
                    return;
                  }
                }
                break;
              }
              if (atom_name[2] == '5') {
                if (atom_name[3] == '\'' || atom_name[3] == '*') {
                  if (rest_is_whitespace(&atom_name[4])) {
                    reference_name = "HO5'";
                    flags |= any | ho5prime;
                    return;
                  }
                }
                break;
              }
              if (atom_name[2] == 'P' || atom_name[2] == 'p') {
                if (atom_name[3] == '2') {
                  if (rest_is_whitespace(&atom_name[4])) {
                    reference_name = "HOP2";
                    flags |= any | phosphate_group;
                    return;
                  }
                  break;
                }
                if (atom_name[3] == '3') {
                  if (rest_is_whitespace(&atom_name[4])) {
                    reference_name = "HOP3";
                    flags |= any | phosphate_group | op3_or_hop3;
                    return;
                  }
                }
              }
              break;

            default:
              break;
          }
          flags = none;
          break;

        case 'N':
        case 'n':
          switch (atom_name[1])
          {
            case '1':
              if (!rest_is_whitespace(&atom_name[2])) break;
              reference_name = " N1 ";
              flags |= any;
              return;

            case '2':
              if (!rest_is_whitespace(&atom_name[2])) break;
              reference_name = " N2 ";
              flags |= g | dg;
              return;

            case '3':
              if (!rest_is_whitespace(&atom_name[2])) break;
              reference_name = " N3 ";
              flags |= any;
              return;

            case '4':
              if (!rest_is_whitespace(&atom_name[2])) break;
              reference_name = " N4 ";
              flags |= c | dc;
              return;

            case '6':
              if (!rest_is_whitespace(&atom_name[2])) break;
              reference_name = " N6 ";
              flags |= a | da;
              return;

            case '7':
              if (!rest_is_whitespace(&atom_name[2])) break;
              reference_name = " N7 ";
              flags |= a | g | da | dg;
              return;

            case '9':
              if (!rest_is_whitespace(&atom_name[2])) break;
              reference_name = " N9 ";
              flags |= a | g | da | dg;
              return;

            default:
              break;
          }
          break;

        case 'O':
        case 'o':
          switch (atom_name[1])
          {
            case '1':
              if (atom_name[2] == 'P' || atom_name[2] == 'p') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " OP1";
                  flags |= any | phosphate_group;
                  return;
                }
              }
              break;

            case '2':
              if (rest_is_whitespace(&atom_name[2])) {
                reference_name = " O2 ";
                flags |= c | u | dc | dt;
                return;
              }
              if (atom_name[2] == '\'' || atom_name[2] == '*') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " O2'";
                  flags |= a | c | g | u | o2prime;
                  return;
                }
                break;
              }
              if (atom_name[2] == 'P' || atom_name[2] == 'p') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " OP2";
                  flags |= any | phosphate_group;
                  return;
                }
              }
              break;

            case '3':
              if (atom_name[2] == '\'' || atom_name[2] == '*') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " O3'";
                  flags |= any;
                  return;
                }
                break;
              }
              if (   atom_name[2] == 'P' || atom_name[2] == 'p'
                  || atom_name[2] == 'T' || atom_name[2] == 't') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " OP3";
                  flags |= any | phosphate_group | op3_or_hop3;
                  return;
                }
              }
              break;

            case '4':
              if (rest_is_whitespace(&atom_name[2])) {
                reference_name = " O4 ";
                flags |= u | dt;
                return;
              }
              if (atom_name[2] == '\'' || atom_name[2] == '*') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " O4'";
                  flags |= any;
                  return;
                }
              }
              break;

            case '5':
              if (atom_name[2] == '\'' || atom_name[2] == '*') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " O5'";
                  flags |= any;
                  return;
                }
                break;
              }
              if (atom_name[2] == 'T' || atom_name[2] == 't') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " OP3";
                  flags |= any | phosphate_group | op3_or_hop3;
                  return;
                }
              }
              break;

            case '6':
              if (!rest_is_whitespace(&atom_name[2])) break;
              reference_name = " O6 ";
              flags |= g | dg;
              return;

            case 'P':
            case 'p':
              if (atom_name[2] == '1') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " OP1";
                  flags |= any | phosphate_group;
                  return;
                }
                break;
              }
              if (atom_name[2] == '2') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " OP2";
                  flags |= any | phosphate_group;
                  return;
                }
                break;
              }
              if (atom_name[2] == '3') {
                if (rest_is_whitespace(&atom_name[3])) {
                  reference_name = " OP3";
                  flags |= any | phosphate_group | op3_or_hop3;
                  return;
                }
              }
              break;

            default:
              break;
          }
          break;

        case 'P':
        case 'p':
          if (rest_is_whitespace(&atom_name[1])) {
            reference_name = " P  ";
            flags |= any | phosphate_group;
            return;
          }
          break;

        default:
          break;
      }
    }

    std::string
    compatible_residue_names() const
    {
      using namespace info_flags;
      std::string result;
      if (flags & any_bit) {
        result += " ANY";
      }
      else {
        if (flags & a) result += " A";
        if (flags & c) result += " C";
        if (flags & g) result += " G";
        if (flags & u) result += " U";
        if (flags & da) result += " DA";
        if (flags & dc) result += " DC";
        if (flags & dg) result += " DG";
        if (flags & dt) result += " DT";
      }
      if (result.size() == 0) return "None";
      return result.substr(1);
    }

    bool
    is_compatible_with(const char* residue_name) const
    {
      using namespace info_flags;
      char letter = residue_name[0];
      bool rna = true;
      if (letter == 'D') {
        rna = false;
        letter = residue_name[1];
      }
      switch (letter)
      {
        case 'A':
          if (rna) return (flags & a) && residue_name[1] == '\0';
          return flags & da           && residue_name[2] == '\0';
        case 'C':
          if (rna) return (flags & c) && residue_name[1] == '\0';
          return flags & dc           && residue_name[2] == '\0';
        case 'G':
          if (rna) return (flags & g) && residue_name[1] == '\0';
          return flags & dg           && residue_name[2] == '\0';
        case 'U':
          if (rna) return (flags & u) && residue_name[1] == '\0';
          break;
        case 'T':
          if (rna) break;
          return flags & dt           && residue_name[2] == '\0';
        default:
          break;
      }
      return false;
    }

    bool
    is_hydrogen() const
    {
      return flags & info_flags::hydrogen;
    }

    bool
    is_deuterium() const
    {
      return flags & info_flags::deuterium;
    }

    bool
    is_o2prime() const
    {
      return flags & info_flags::o2prime;
    }

    bool
    is_ho2prime() const
    {
      return flags & info_flags::ho2prime;
    }

    bool
    is_h2primeprime() const
    {
      return flags & info_flags::h2primeprime;
    }

    bool
    is_in_phosphate_group() const
    {
      return flags & info_flags::phosphate_group;
    }

    bool
    is_op3_or_hop3() const
    {
      return flags & info_flags::op3_or_hop3;
    }

    bool
    is_ho5prime() const
    {
      return flags & info_flags::ho5prime;
    }

    bool
    is_ho3prime() const
    {
      return flags & info_flags::ho3prime;
    }

    bool
    change_h2primeprime_to_ho2prime()
    {
      if (!is_h2primeprime()) return false;
      reference_name = "HO2'";
      using namespace info_flags;
      flags = (is_deuterium() ? deuterium : none)
            | a | c | g | u | hydrogen | ho2prime;
      return true;
    }

    bool
    change_ho5prime_to_hop3()
    {
      if (!is_ho5prime()) return false;
      reference_name = "HOP3";
      using namespace info_flags;
      flags = (is_deuterium() ? deuterium : none)
            | any | hydrogen | phosphate_group | op3_or_hop3;
      return true;
    }

    void
    change_to_unknown()
    {
      reference_name = 0;
      flags = info_flags::none;
    }

    private:

      bool
      rest_is_whitespace(const char* s) const
      {
        while (*s) {
          if (!std::isspace(*s++)) return false;
        }
        return true;
      }
  };

}}} // namespace iotbx::pdb::rna_dna_atom_names

#endif // IOTBX_PDB_RNA_DNA_ATOM_NAMES_H
