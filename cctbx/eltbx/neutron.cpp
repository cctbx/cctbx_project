// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
               Based on C code contributed by Vincent Favre-Nicolin.
 */

#include <cctbx/eltbx/basic.h>
#include <cctbx/eltbx/neutron.h>

namespace eltbx {
  namespace tables {

    /*
      Neutron bound scattering lengths & cross-sections

      Data from: http://www.ncnr.nist.gov/resources/n-lengths/list.html

      All of this data was taken from the Special Feature section of
      neutron scattering lengths and cross sections of the elements and
      their isotopes in Neutron News, Vol. 3, No. 3, 1992, pp. 29-37.
     */
    const detail::RawNeutronNews1992Record RawNeutronNews1992Records[] =
    {
// BEGIN_COMPILED_IN_REFERENCE_DATA
      {"H", -3.7390, 0, 0.3326},
      {"D", 6.671, 0, 0.000519},
      {"He", 3.26, 0, 0.00747},
      {"Li", -1.90, 0, 70.5},
      {"Be", 7.79, 0, 0.0076},
      {"B", 5.30, -0.213, 767.},
      {"C", 6.6460, 0, 0.0035},
      {"N", 9.36, 0, 1.9},
      {"O", 5.803, 0, 0.00019},
      {"F", 5.654, 0, 0.0096},
      {"Ne", 4.566, 0, 0.039},
      {"Na", 3.63, 0, 0.53},
      {"Mg", 5.375, 0, 0.063},
      {"Al", 3.449, 0, 0.231},
      {"Si", 4.1491, 0, 0.171},
      {"P", 5.13, 0, 0.172},
      {"S", 2.847, 0, 0.53},
      {"Cl", 9.5770, 0, 33.5},
      {"Ar", 1.909, 0, 0.675},
      {"K", 3.67, 0, 2.1},
      {"Ca", 4.70, 0, 0.43},
      {"Sc", 12.29, 0, 27.5},
      {"Ti", -3.438, 0, 6.09},
      {"V", -0.3824, 0, 5.08},
      {"Cr", 3.635, 0, 3.05},
      {"Mn", -3.73, 0, 13.3},
      {"Fe", 9.45, 0, 2.56},
      {"Co", 2.49, 0, 37.18},
      {"Ni", 10.3, 0, 4.49},
      {"Cu", 7.718, 0, 3.78},
      {"Zn", 5.680, 0, 1.11},
      {"Ga", 7.288, 0, 2.75},
      {"Ge", 8.185, 0, 2.2},
      {"As", 6.58, 0, 4.5},
      {"Se", 7.970, 0, 11.7},
      {"Br", 6.795, 0, 6.9},
      {"Kr", 7.81, 0, 25.},
      {"Rb", 7.09, 0, 0.38},
      {"Sr", 7.02, 0, 1.28},
      {"Y", 7.75, 0, 1.28},
      {"Zr", 7.16, 0, 0.185},
      {"Nb", 7.054, 0, 1.15},
      {"Mo", 6.715, 0, 2.48},
      {"Tc", 6.8, 0, 20.},
      {"Ru", 7.03, 0, 2.56},
      {"Rh", 5.88, 0, 144.8},
      {"Pd", 5.91, 0, 6.9},
      {"Ag", 5.922, 0, 63.3},
      {"Cd", 4.87, -0.70, 2520.},
      {"In", 4.065, -0.0539, 193.8},
      {"Sn", 6.225, 0, 0.626},
      {"Sb", 5.57, 0, 4.91},
      {"Te", 5.80, 0, 4.7},
      {"I", 5.28, 0, 6.15},
      {"Xe", 4.92, 0, 23.9},
      {"Cs", 5.42, 0, 29.0},
      {"Ba", 5.07, 0, 1.1},
      {"La", 8.24, 0, 8.97},
      {"Ce", 4.84, 0, 0.63},
      {"Pr", 4.58, 0, 11.5},
      {"Nd", 7.69, 0, 50.5},
      {"Pm", 12.6, 0, 168.4},
      {"Sm", 0.80, -1.65, 5922.},
      {"Eu", 7.22, -1.26, 4530.},
      {"Gd", 6.5, -13.82, 49700.},
      {"Tb", 7.38, 0, 23.4},
      {"Dy", 16.9, -0.276, 994.},
      {"Ho", 8.01, 0, 64.7},
      {"Er", 7.79, 0, 159.},
      {"Tm", 7.07, 0, 100.},
      {"Yb", 12.43, 0, 34.8},
      {"Lu", 7.21, 0, 74.},
      {"Hf", 7.7, 0, 104.1},
      {"Ta", 6.91, 0, 20.6},
      {"W", 4.86, 0, 18.3},
      {"Re", 9.2, 0, 89.7},
      {"Os", 10.7, 0, 16},
      {"Ir", 10.6, 0, 425.},
      {"Pt", 9.60, 0, 10.3},
      {"Au", 7.63, 0, 98.65},
      {"Hg", 12.692, 0, 372.3},
      {"Tl", 8.776, 0, 3.43},
      {"Pb", 9.405, 0, 0.171},
      {"Bi", 8.532, 0, 0.0338},
      {"Th", 10.31, 0, 7.37},
      {"U", 8.417, 0, 7.57},
      {0, 0, 0, 0}
// END_COMPILED_IN_REFERENCE_DATA
    };

  } // namespace tables
} // namespace eltbx

namespace {

  const eltbx::detail::RawNeutronNews1992Record*
  FindEntry(const std::string& WorkLabel, bool Exact)
  {
    int m = 0;
    const eltbx::detail::RawNeutronNews1992Record* mEntry = 0;
    for (const eltbx::detail::RawNeutronNews1992Record*
         Entry = eltbx::tables::RawNeutronNews1992Records;
         Entry->Symbol;
         Entry++)
    {
      int i = eltbx::MatchLabels(WorkLabel, Entry->Symbol);
      if (i < 0) return Entry;
      if (i > m) {
        m = i;
        mEntry = Entry;
      }
    }
    if (Exact || !mEntry) {
      throw eltbx::error("Unknown element symbol.");
    }
    return mEntry;
  }

} // namespace <anonymous>

namespace eltbx {

  NeutronNews1992Record::NeutronNews1992Record(const std::string& Label,
                                               bool Exact)
  {
    std::string WorkLabel = StripLabel(Label, Exact);
    m_RawEntry = FindEntry(WorkLabel, Exact);
  }

} // namespace eltbx
