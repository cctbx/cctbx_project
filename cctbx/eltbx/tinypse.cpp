// $Id$

#include <cctbx/eltbx/basic.h>
#include <cctbx/eltbx/tinypse.h>

namespace eltbx {
  namespace tables {

    /*
      CRC Handbook of Chemistry & Physics, 63rd edition, 1982-1983
      CRC Handbook of Chemistry & Physics, 70th edition, 1989-1990
     */
    const detail::TinyPSE_RawEntry TinyPSE_RawTable[] =
    {
      {   1, "H",  "hydrogen",       1.008 },
      {   1, "D",  "deuterium",      2.000 },
      {   2, "He", "helium",         4.003 },
      {   3, "Li", "lithium",        6.941 },
      {   4, "Be", "beryllium",      9.012 },
      {   5, "B",  "boron",         10.811 },
      {   6, "C",  "carbon",        12.011 },
      {   7, "N",  "nitrogen",      14.007 },
      {   8, "O",  "oxygen",        15.999 },
      {   9, "F",  "fluorine",      18.998 },
      {  10, "Ne", "neon",          20.180 },
      {  11, "Na", "sodium",        22.990 },
      {  12, "Mg", "magnesium",     24.305 },
      {  13, "Al", "aluminium",     26.982 },
      {  14, "Si", "silicon",       28.086 },
      {  15, "P",  "phosphorus",    30.974 },
      {  16, "S",  "sulphur",       32.066 },
      {  17, "Cl", "chlorine",      35.452 },
      {  18, "Ar", "argon",         39.948 },
      {  19, "K",  "potassium",     39.098 },
      {  20, "Ca", "calcium",       40.078 },
      {  21, "Sc", "scandium",      44.956 },
      {  22, "Ti", "titanium",      47.883 },
      {  23, "V",  "vanadium",      50.941 },
      {  24, "Cr", "chromium",      51.996 },
      {  25, "Mn", "manganese",     54.938 },
      {  26, "Fe", "iron",          55.847 },
      {  27, "Co", "cobalt",        58.933 },
      {  28, "Ni", "nickel",        58.691 },
      {  29, "Cu", "copper",        63.546 },
      {  30, "Zn", "zinc",          65.392 },
      {  31, "Ga", "gallium",       69.723 },
      {  32, "Ge", "germanium",     72.612 },
      {  33, "As", "arsenic",       74.922 },
      {  34, "Se", "selenium",      78.963 },
      {  35, "Br", "bromine",       79.904 },
      {  36, "Kr", "krypton",       83.801 },
      {  37, "Rb", "rubidium",      85.468 },
      {  38, "Sr", "strontium",     87.621 },
      {  39, "Y",  "yttrium",       88.906 },
      {  40, "Zr", "zirconium",     91.224 },
      {  41, "Nb", "niobium",       92.906 },
      {  42, "Mo", "molybdenum",    95.941 },
      {  43, "Tc", "technetium",    98.000 },
      {  44, "Ru", "ruthenium",    101.072 },
      {  45, "Rh", "rhodium",      102.905 },
      {  46, "Pd", "palladium",    106.421 },
      {  47, "Ag", "silver",       107.868 },
      {  48, "Cd", "cadmium",      112.411 },
      {  49, "In", "indium",       114.821 },
      {  50, "Sn", "tin",          118.710 },
      {  51, "Sb", "antimony",     121.753 },
      {  52, "Te", "tellurium",    127.603 },
      {  53, "I",  "iodine",       126.904 },
      {  54, "Xe", "xenon",        131.292 },
      {  55, "Cs", "caesium",      132.905 },
      {  56, "Ba", "barium",       137.327 },
      {  57, "La", "lanthanum",    138.906 },
      {  58, "Ce", "cerium",       140.115 },
      {  59, "Pr", "praseodymium", 140.908 },
      {  60, "Nd", "neodymium",    144.243 },
      {  61, "Pm", "promethium",   145.000 },
      {  62, "Sm", "samarium",     150.363 },
      {  63, "Eu", "europium",     151.965 },
      {  64, "Gd", "gadolinium",   157.253 },
      {  65, "Tb", "terbium",      158.925 },
      {  66, "Dy", "dysprosium",   162.503 },
      {  67, "Ho", "holmium",      164.930 },
      {  68, "Er", "erbium",       167.263 },
      {  69, "Tm", "thulium",      168.934 },
      {  70, "Yb", "ytterbium",    173.043 },
      {  71, "Lu", "lutetium",     174.967 },
      {  72, "Hf", "hafnium",      178.492 },
      {  73, "Ta", "tantalum",     180.948 },
      {  74, "W",  "tungsten",     183.853 },
      {  75, "Re", "rhenium",      186.207 },
      {  76, "Os", "osmium",       190.210 },
      {  77, "Ir", "iridium",      192.223 },
      {  78, "Pt", "platinum",     195.083 },
      {  79, "Au", "gold",         196.967 },
      {  80, "Hg", "mercury",      200.593 },
      {  81, "Tl", "thallium",     204.383 },
      {  82, "Pb", "lead",         207.210 },
      {  83, "Bi", "bismuth",      208.980 },
      {  84, "Po", "polonium",     209.000 },
      {  85, "At", "astatine",     210.000 },
      {  86, "Rn", "radon",        222.000 },
      {  87, "Fr", "francium",     223.000 },
      {  88, "Ra", "radium",       226.025 },
      {  89, "Ac", "actinium",     227.028 },
      {  90, "Th", "thorium",      232.038 },
      {  91, "Pa", "protactinium", 231.035 },
      {  92, "U",  "uranium",      238.028 },
      {  93, "Np", "neptunium",    237.048 },
      {  94, "Pu", "plutonium",    244.000 },
      {  95, "Am", "americium",    243.000 },
      {  96, "Cm", "curium",       247.000 },
      {  97, "Bk", "berkelium",    247.000 },
      {  98, "Cf", "californium",  251.000 },
      {  99, "Es", "einsteinium",  254.000 },
      { 100, "Fm", "fermium",      257.000 },
      { 101, "Md", "mendelevium",  258.000 },
      { 102, "No", "nobelium",     259.000 },
      { 103, "Lr", "lawrencium",   260.000 },
      { 0, 0, 0, 0. }
    };

  } // namespace tables
} // namespace eltbx

namespace {

  const eltbx::detail::TinyPSE_RawEntry*
  FindEntry(const std::string& WorkLabel, bool Exact)
  {
    int m = 0;
    const eltbx::detail::TinyPSE_RawEntry* mEntry = 0;
    for (const eltbx::detail::TinyPSE_RawEntry*
         Entry = eltbx::tables::TinyPSE_RawTable; Entry->Symbol; Entry++)
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

  TinyPSE::TinyPSE(const std::string& Label, bool Exact)
  {
    std::string WorkLabel = StripLabel(Label, Exact);
    m_RawEntry = FindEntry(WorkLabel, Exact);
  }

  TinyPSE::TinyPSE(int Z)
  {
    for (m_RawEntry = eltbx::tables::TinyPSE_RawTable;
         m_RawEntry->Symbol;
         m_RawEntry++) {
      if (m_RawEntry->Z == Z) return;
    }
    throw eltbx::error("Atomic number out of range.");
  }

} // namespace eltbx
