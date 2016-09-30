#include <cctbx/eltbx/neutron.h>
#include <cctbx/eltbx/basic.h>

namespace cctbx { namespace eltbx { namespace neutron {

namespace detail {
namespace {

    /*
      Neutron bound scattering lengths & cross-sections

      Data from: http://www.ncnr.nist.gov/resources/n-lengths/list.html

      All of this data was taken from the Special Feature section of
      neutron scattering lengths and cross sections of the elements and
      their isotopes in Neutron News, Vol. 3, No. 3, 1992, pp. 29-37.
     */
    const raw_record_neutron_news_1992 table[] =
    {
// BEGIN_COMPILED_IN_REFERENCE_DATA
      {"H", -3.7390, 0, 0.3326},
      {"D", 6.671, 0, 0.000519},
      {"He", 3.26, 0, 0.00747},
      {"Li",   -1.90, 0, 70.5},
      {"Li1+", -1.90, 0, 70.5}, // copy
      {"Be",   7.79, 0, 0.0076},
      {"Be2+", 7.79, 0, 0.0076}, // copy
      {"B", 5.30, -0.213, 767.},
      {"C", 6.6460, 0, 0.0035},
      {"N", 9.36, 0, 1.9},
      {"O",   5.803, 0, 0.00019},
      {"O1-", 5.803, 0, 0.00019}, // copy
      {"O2-", 5.803, 0, 0.00019}, // copy
      {"F",   5.654, 0, 0.0096},
      {"F1-", 5.654, 0, 0.0096}, // copy
      {"Ne", 4.566, 0, 0.039},
      {"Na",   3.63, 0, 0.53},
      {"Na1+", 3.63, 0, 0.53}, // copy
      {"Mg",   5.375, 0, 0.063},
      {"Mg2+", 5.375, 0, 0.063}, // copy
      {"Al",   3.449, 0, 0.231},
      {"Al3+", 3.449, 0, 0.231}, // copy
      {"Si",   4.1491, 0, 0.171},
      {"Si4+", 4.1491, 0, 0.171}, // copy
      {"P", 5.13, 0, 0.172},
      {"S", 2.847, 0, 0.53},
      {"Cl",   9.5770, 0, 33.5},
      {"Cl1-", 9.5770, 0, 33.5}, // copy
      {"Ar", 1.909, 0, 0.675},
      {"K",   3.67, 0, 2.1},
      {"K1+", 3.67, 0, 2.1}, // copy
      {"Ca",   4.70, 0, 0.43},
      {"Ca2+", 4.70, 0, 0.43}, // copy
      {"Sc",   12.29, 0, 27.5},
      {"Sc3+", 12.29, 0, 27.5}, // copy
      {"Ti",   -3.438, 0, 6.09},
      {"Ti2+", -3.438, 0, 6.09}, // copy
      {"Ti3+", -3.438, 0, 6.09}, // copy
      {"Ti4+", -3.438, 0, 6.09}, // copy
      {"V",   -0.3824, 0, 5.08},
      {"V2+", -0.3824, 0, 5.08}, // copy
      {"V3+", -0.3824, 0, 5.08}, // copy
      {"V5+", -0.3824, 0, 5.08}, // copy
      {"Cr",   3.635, 0, 3.05},
      {"Cr2+", 3.635, 0, 3.05}, // copy
      {"Cr3+", 3.635, 0, 3.05}, // copy
      {"Mn",   -3.73, 0, 13.3},
      {"Mn2+", -3.73, 0, 13.3}, // copy
      {"Mn3+", -3.73, 0, 13.3}, // copy
      {"Mn4+", -3.73, 0, 13.3}, // copy
      {"Fe",   9.45, 0, 2.56},
      {"Fe2+", 9.45, 0, 2.56}, // copy
      {"Fe3+", 9.45, 0, 2.56}, // copy
      {"Co",   2.49, 0, 37.18},
      {"Co2+", 2.49, 0, 37.18}, // copy
      {"Co3+", 2.49, 0, 37.18}, // copy
      {"Ni",   10.3, 0, 4.49},
      {"Ni2+", 10.3, 0, 4.49}, // copy
      {"Ni3+", 10.3, 0, 4.49}, // copy
      {"Cu",   7.718, 0, 3.78},
      {"Cu1+", 7.718, 0, 3.78}, // copy
      {"Cu2+", 7.718, 0, 3.78}, // copy
      {"Zn",   5.680, 0, 1.11},
      {"Zn2+", 5.680, 0, 1.11}, // copy
      {"Ga",   7.288, 0, 2.75},
      {"Ga3+", 7.288, 0, 2.75}, // copy
      {"Ge",   8.185, 0, 2.2},
      {"Ge4+", 8.185, 0, 2.2}, // copy
      {"As", 6.58, 0, 4.5},
      {"Se", 7.970, 0, 11.7},
      {"Br",   6.795, 0, 6.9},
      {"Br1-", 6.795, 0, 6.9}, // copy
      {"Kr", 7.81, 0, 25.},
      {"Rb",   7.09, 0, 0.38},
      {"Rb1+", 7.09, 0, 0.38}, // copy
      {"Sr",   7.02, 0, 1.28},
      {"Sr2+", 7.02, 0, 1.28}, // copy
      {"Y",   7.75, 0, 1.28},
      {"Y3+", 7.75, 0, 1.28}, // copy
      {"Zr",   7.16, 0, 0.185},
      {"Zr4+", 7.16, 0, 0.185}, // copy
      {"Nb",   7.054, 0, 1.15},
      {"Nb3+", 7.054, 0, 1.15}, // copy
      {"Nb5+", 7.054, 0, 1.15}, // copy
      {"Mo",   6.715, 0, 2.48},
      {"Mo3+", 6.715, 0, 2.48}, // copy
      {"Mo5+", 6.715, 0, 2.48}, // copy
      {"Mo6+", 6.715, 0, 2.48}, // copy
      {"Tc", 6.8, 0, 20.},
      {"Ru",   7.03, 0, 2.56},
      {"Ru3+", 7.03, 0, 2.56}, // copy
      {"Ru4+", 7.03, 0, 2.56}, // copy
      {"Rh",   5.88, 0, 144.8},
      {"Rh3+", 5.88, 0, 144.8}, // copy
      {"Rh4+", 5.88, 0, 144.8}, // copy
      {"Pd",   5.91, 0, 6.9},
      {"Pd2+", 5.91, 0, 6.9}, // copy
      {"Pd4+", 5.91, 0, 6.9}, // copy
      {"Ag",   5.922, 0, 63.3},
      {"Ag1+", 5.922, 0, 63.3}, // copy
      {"Ag2+", 5.922, 0, 63.3}, // copy
      {"Cd",   4.87, -0.70, 2520.},
      {"Cd2+", 4.87, -0.70, 2520.}, // copy
      {"In",   4.065, -0.0539, 193.8},
      {"In3+", 4.065, -0.0539, 193.8}, // copy
      {"Sn",   6.225, 0, 0.626},
      {"Sn2+", 6.225, 0, 0.626}, // copy
      {"Sn4+", 6.225, 0, 0.626}, // copy
      {"Sb",   5.57, 0, 4.91},
      {"Sb3+", 5.57, 0, 4.91},
      {"Sb5+", 5.57, 0, 4.91},
      {"Te", 5.80, 0, 4.7},
      {"I",   5.28, 0, 6.15},
      {"I1-", 5.28, 0, 6.15}, // copy
      {"Xe", 4.92, 0, 23.9},
      {"Cs",   5.42, 0, 29.0},
      {"Cs1+", 5.42, 0, 29.0}, // copy
      {"Ba",   5.07, 0, 1.1},
      {"Ba2+", 5.07, 0, 1.1}, // copy
      {"La",   8.24, 0, 8.97},
      {"La3+", 8.24, 0, 8.97}, // copy
      {"Ce",   4.84, 0, 0.63},
      {"Ce3+", 4.84, 0, 0.63}, // copy
      {"Ce4+", 4.84, 0, 0.63}, // copy
      {"Pr",   4.58, 0, 11.5},
      {"Pr3+", 4.58, 0, 11.5}, // copy
      {"Pr4+", 4.58, 0, 11.5}, // copy
      {"Nd",   7.69, 0, 50.5},
      {"Nd3+", 7.69, 0, 50.5}, // copy
      {"Pm",   12.6, 0, 168.4},
      {"Pm3+", 12.6, 0, 168.4}, // copy
      {"Sm",   0.80, -1.65, 5922.},
      {"Sm3+", 0.80, -1.65, 5922.}, // copy
      {"Eu",   7.22, -1.26, 4530.},
      {"Eu2+", 7.22, -1.26, 4530.}, // copy
      {"Eu3+", 7.22, -1.26, 4530.}, // copy
      {"Gd",   6.5, -13.82, 49700.},
      {"Gd3+", 6.5, -13.82, 49700.}, // copy
      {"Tb",   7.38, 0, 23.4},
      {"Tb3+", 7.38, 0, 23.4}, // copy
      {"Dy",   16.9, -0.276, 994.},
      {"Dy3+", 16.9, -0.276, 994.}, // copy
      {"Ho",   8.01, 0, 64.7},
      {"Ho3+", 8.01, 0, 64.7}, // copy
      {"Er",   7.79, 0, 159.},
      {"Er3+", 7.79, 0, 159.}, // copy
      {"Tm",   7.07, 0, 100.},
      {"Tm3+", 7.07, 0, 100.}, // copy
      {"Yb",   12.43, 0, 34.8},
      {"Yb2+", 12.43, 0, 34.8}, // copy
      {"Yb3+", 12.43, 0, 34.8}, // copy
      {"Lu",   7.21, 0, 74.},
      {"Lu3+", 7.21, 0, 74.}, // copy
      {"Hf",   7.7, 0, 104.1},
      {"Hf4+", 7.7, 0, 104.1}, // copy
      {"Ta",   6.91, 0, 20.6},
      {"Ta5+", 6.91, 0, 20.6}, // copy
      {"W",   4.86, 0, 18.3},
      {"W6+", 4.86, 0, 18.3}, // copy
      {"Re", 9.2, 0, 89.7},
      {"Os",   10.7, 0, 16},
      {"Os4+", 10.7, 0, 16}, // copy
      {"Ir",   10.6, 0, 425.},
      {"Ir3+", 10.6, 0, 425.}, // copy
      {"Ir4+", 10.6, 0, 425.}, // copy
      {"Pt",   9.60, 0, 10.3},
      {"Pt2+", 9.60, 0, 10.3}, // copy
      {"Pt4+", 9.60, 0, 10.3}, // copy
      {"Au",   7.63, 0, 98.65},
      {"Au1+", 7.63, 0, 98.65}, // copy
      {"Au3+", 7.63, 0, 98.65}, // copy
      {"Hg",   12.692, 0, 372.3},
      {"Hg1+", 12.692, 0, 372.3}, // copy
      {"Hg2+", 12.692, 0, 372.3}, // copy
      {"Tl",   8.776, 0, 3.43},
      {"Tl1+", 8.776, 0, 3.43}, // copy
      {"Tl3+", 8.776, 0, 3.43}, // copy
      {"Pb",   9.405, 0, 0.171},
      {"Pb2+", 9.405, 0, 0.171}, // copy
      {"Pb4+", 9.405, 0, 0.171}, // copy
      {"Bi",   8.532, 0, 0.0338},
      {"Bi3+", 8.532, 0, 0.0338}, // copy
      {"Bi5+", 8.532, 0, 0.0338}, // copy
      {"Th",   10.31, 0, 7.37},
      {"Th4+", 10.31, 0, 7.37}, // copy
      {"U",   8.417, 0, 7.57},
      {"U3+", 8.417, 0, 7.57}, // copy
      {"U4+", 8.417, 0, 7.57}, // copy
      {"U6+", 8.417, 0, 7.57}, // copy
      {0, 0, 0, 0}
// END_COMPILED_IN_REFERENCE_DATA
    };

  const raw_record_neutron_news_1992*
  find_record(std::string const& work_label, bool exact)
  {
    int m = 0;
    const raw_record_neutron_news_1992* matching_record = 0;
    for (const raw_record_neutron_news_1992* r=table; r->label; r++) {
      int i = basic::match_labels(work_label, r->label);
      if (i < 0) return r;
      if (i > m) {
        m = i;
        matching_record = r;
      }
    }
    if (exact || !matching_record) {
      throw std::invalid_argument("Unknown element label:" + work_label);
    }
    return matching_record;
  }

} // namespace <anonymous>
} // namespace detail

  neutron_news_1992_table::neutron_news_1992_table(std::string const& label,
                                                   bool exact)
  {
    std::string work_label = basic::strip_label(label, exact);
    record_ = detail::find_record(work_label, exact);
  }

  neutron_news_1992_table_iterator::neutron_news_1992_table_iterator()
  :
    current_("H", true)
  {}

  neutron_news_1992_table
  neutron_news_1992_table_iterator::next()
  {
    neutron_news_1992_table result = current_;
    if (current_.is_valid()) current_.record_++;
    return result;
  }

}}} // namespace cctbx::eltbx::neutron
