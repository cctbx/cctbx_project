#include <cctbx/eltbx/covalent_radii.h>
#include <cctbx/eltbx/basic.h>

namespace cctbx { namespace eltbx { namespace covalent_radii {

namespace detail {
namespace {

/*! Table of covalent radii

    Reference:
      Beatriz Cordero, Verónica Gómez, Ana E. Platero-Prats, Marc Revés,
      Jorge Echeverría, Eduard Cremades, Flavia Barragán and Santiago Alvarez.
      Covalent radii revisited. Dalton Trans., 2008, 2832-2838,
      doi:10.1039/b801115j
*/
    const raw_record table[] = {
// BEGIN_COMPILED_IN_REFERENCE_DATA
      {"H",  0.31, 0.05},
      {"D",  0.31, 0.05},
      {"He", 0.28, 0.00},
      {"Li", 1.28, 0.07},
      {"Be", 0.96, 0.03},
      {"B",  0.84, 0.03},
      {"C",  0.76, 0.01}, /* sp3 value */
      {"N",  0.71, 0.01},
      {"O",  0.66, 0.02},
      {"F",  0.57, 0.03},
      {"Ne", 0.58, 0.00},
      {"Na", 1.66, 0.09},
      {"Mg", 1.41, 0.07},
      {"Al", 1.21, 0.04},
      {"Si", 1.11, 0.02},
      {"P",  1.07, 0.03},
      {"S",  1.05, 0.03},
      {"Cl", 1.02, 0.04},
      {"Ar", 1.06, 0.10},
      {"K",  2.03, 0.12},
      {"Ca", 1.76, 0.10},
      {"Sc", 1.70, 0.07},
      {"Ti", 1.60, 0.08},
      {"V",  1.53, 0.08},
      {"Cr", 1.39, 0.05},
      {"Mn", 1.39, 0.05}, /* low-spin value */
      {"Fe", 1.32, 0.03}, /* low-spin value */
      {"Co", 1.26, 0.03}, /* low-spin value */
      {"Ni", 1.24, 0.04},
      {"Cu", 1.32, 0.04},
      {"Zn", 1.22, 0.04},
      {"Ga", 1.22, 0.03},
      {"Ge", 1.20, 0.04},
      {"As", 1.19, 0.04},
      {"Se", 1.20, 0.04},
      {"Br", 1.20, 0.03},
      {"Kr", 1.16, 0.04},
      {"Rb", 2.20, 0.09},
      {"Sr", 1.95, 0.10},
      {"Y",  1.90, 0.07},
      {"Zr", 1.75, 0.07},
      {"Nb", 1.64, 0.06},
      {"Mo", 1.54, 0.05},
      {"Tc", 1.47, 0.07},
      {"Ru", 1.46, 0.07},
      {"Rh", 1.42, 0.07},
      {"Pd", 1.39, 0.06},
      {"Ag", 1.45, 0.05},
      {"Cd", 1.44, 0.09},
      {"In", 1.42, 0.05},
      {"Sn", 1.39, 0.04},
      {"Sb", 1.39, 0.05},
      {"Te", 1.38, 0.04},
      {"I",  1.39, 0.03},
      {"Xe", 1.40, 0.09},
      {"Cs", 2.44, 0.11},
      {"Ba", 2.15, 0.11},
      {"La", 2.07, 0.08},
      {"Ce", 2.04, 0.09},
      {"Pr", 2.03, 0.07},
      {"Nd", 2.01, 0.06},
      {"Pm", 1.99, 0.00},
      {"Sm", 1.98, 0.08},
      {"Eu", 1.98, 0.06},
      {"Gd", 1.96, 0.06},
      {"Tb", 1.94, 0.05},
      {"Dy", 1.92, 0.07},
      {"Ho", 1.92, 0.07},
      {"Er", 1.89, 0.06},
      {"Tm", 1.90, 0.10},
      {"Yb", 1.87, 0.08},
      {"Lu", 1.87, 0.08},
      {"Hf", 1.75, 0.10},
      {"Ta", 1.70, 0.08},
      {"W",  1.62, 0.07},
      {"Re", 1.51, 0.07},
      {"Os", 1.44, 0.04},
      {"Ir", 1.41, 0.06},
      {"Pt", 1.36, 0.05},
      {"Au", 1.36, 0.06},
      {"Hg", 1.32, 0.05},
      {"Tl", 1.45, 0.07},
      {"Pb", 1.46, 0.05},
      {"Bi", 1.48, 0.04},
      {"Po", 1.40, 0.04},
      {"At", 1.50, 0.00},
      {"Rn", 1.50, 0.00},
      {"Fr", 2.60, 0.00},
      {"Ra", 2.21, 0.02},
      {"Ac", 2.15, 0.00},
      {"Th", 2.06, 0.06},
      {"Pa", 2.00, 0.00},
      {"U",  1.96, 0.07},
      {"Np", 1.90, 0.01},
      {"Pu", 1.87, 0.01},
      {"Am", 1.80, 0.06},
      {"Cm", 1.69, 0.03},
      {0, 0}
// END_COMPILED_IN_REFERENCE_DATA
    };

  const raw_record*
  find_record(std::string const& work_label, bool exact)
  {
    int m = 0;
    const raw_record* matching_record = 0;
    for (const raw_record* record = table; record->label; record++) {
      int i = basic::match_labels(work_label, record->label);
      if (i < 0) return record;
      if (i > m) {
        m = i;
        matching_record = record;
      }
    }
    if (exact || !matching_record) {
      throw error("Unknown atom label.");
    }
    return matching_record;
  }

} // namespace <anonymous>
} // namespace detail

  table::table(std::string const& label, bool exact)
  {
    std::string work_label = basic::strip_label(label, exact);
    record_ = detail::find_record(work_label, exact);
  }

  table_iterator::table_iterator()
  :
    current_("H", true)
  {}

  table
  table_iterator::next()
  {
    table result = current_;
    if (current_.is_valid()) current_.record_++;
    return result;
  }

}}} // namespace cctbx::eltbx::covalent_radii
