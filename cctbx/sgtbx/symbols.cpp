#include <ctype.h>
#include <string>
#include <stdio.h>
#include <cstddef>
#include <cctbx/error.h>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/reference_settings.h>
#include <cctbx/sgtbx/change_of_basis_op.h>
#include <map>

using std::string;

namespace cctbx { namespace sgtbx {

  namespace {

    static const char* not_recognized = "Space group symbol not recognized: ";

    inline
    string
    to_str(int value)
    {
      char buf[256];
      sprintf(buf, "%d", value);
      return string(buf);
    }

  }

  namespace symbols { namespace tables {

    static const char* monoclinic_sg_number_as_hm_list[][2] = {
               { 0, 0 },
               { 0, 0 },
               { 0, 0 },
      /*  3 */ { "P121",    "P112" },
      /*  4 */ { "P1211",   "P1121" },
      /*  5 */ { "C121",    "B112" },
      /*  6 */ { "P1m1",    "P11m" },
      /*  7 */ { "P1c1",    "P11b" },
      /*  8 */ { "C1m1",    "B11m" },
      /*  9 */ { "C1c1",    "B11b" },
      /* 10 */ { "P12/m1",  "P112/m" },
      /* 11 */ { "P121/m1", "P1121/m" },
      /* 12 */ { "C12/m1",  "B112/m" },
      /* 13 */ { "P12/c1",  "P112/b" },
      /* 14 */ { "P121/c1", "P1121/b" },
      /* 15 */ { "C12/c1",  "B112/b" },
    };

    static const char* schoenflies_list[] = {
// BEGIN_COMPILED_IN_REFERENCE_DATA
      0,
      "C1^1", "Ci^1", "C2^1", "C2^2", "C2^3", "Cs^1", "Cs^2", "Cs^3",
      "Cs^4", "C2h^1", "C2h^2", "C2h^3", "C2h^4", "C2h^5", "C2h^6",
      "D2^1", "D2^2", "D2^3", "D2^4", "D2^5", "D2^6", "D2^7", "D2^8",
      "D2^9", "C2v^1", "C2v^2", "C2v^3", "C2v^4", "C2v^5", "C2v^6",
      "C2v^7", "C2v^8", "C2v^9", "C2v^10", "C2v^11", "C2v^12",
      "C2v^13", "C2v^14", "C2v^15", "C2v^16", "C2v^17", "C2v^18",
      "C2v^19", "C2v^20", "C2v^21", "C2v^22", "D2h^1", "D2h^2",
      "D2h^3", "D2h^4", "D2h^5", "D2h^6", "D2h^7", "D2h^8", "D2h^9",
      "D2h^10", "D2h^11", "D2h^12", "D2h^13", "D2h^14", "D2h^15",
      "D2h^16", "D2h^17", "D2h^18", "D2h^19", "D2h^20", "D2h^21",
      "D2h^22", "D2h^23", "D2h^24", "D2h^25", "D2h^26", "D2h^27",
      "D2h^28", "C4^1", "C4^2", "C4^3", "C4^4", "C4^5", "C4^6",
      "S4^1", "S4^2", "C4h^1", "C4h^2", "C4h^3", "C4h^4", "C4h^5",
      "C4h^6", "D4^1", "D4^2", "D4^3", "D4^4", "D4^5", "D4^6",
      "D4^7", "D4^8", "D4^9", "D4^10", "C4v^1", "C4v^2", "C4v^3",
      "C4v^4", "C4v^5", "C4v^6", "C4v^7", "C4v^8", "C4v^9", "C4v^10",
      "C4v^11", "C4v^12", "D2d^1", "D2d^2", "D2d^3", "D2d^4",
      "D2d^5", "D2d^6", "D2d^7", "D2d^8", "D2d^9", "D2d^10",
      "D2d^11", "D2d^12", "D4h^1", "D4h^2", "D4h^3", "D4h^4",
      "D4h^5", "D4h^6", "D4h^7", "D4h^8", "D4h^9", "D4h^10",
      "D4h^11", "D4h^12", "D4h^13", "D4h^14", "D4h^15", "D4h^16",
      "D4h^17", "D4h^18", "D4h^19", "D4h^20", "C3^1", "C3^2", "C3^3",
      "C3^4", "C3i^1", "C3i^2", "D3^1", "D3^2", "D3^3", "D3^4",
      "D3^5", "D3^6", "D3^7", "C3v^1", "C3v^2", "C3v^3", "C3v^4",
      "C3v^5", "C3v^6", "D3d^1", "D3d^2", "D3d^3", "D3d^4", "D3d^5",
      "D3d^6", "C6^1", "C6^2", "C6^3", "C6^4", "C6^5", "C6^6",
      "C3h^1", "C6h^1", "C6h^2", "D6^1", "D6^2", "D6^3", "D6^4",
      "D6^5", "D6^6", "C6v^1", "C6v^2", "C6v^3", "C6v^4", "D3h^1",
      "D3h^2", "D3h^3", "D3h^4", "D6h^1", "D6h^2", "D6h^3", "D6h^4",
      "T^1", "T^2", "T^3", "T^4", "T^5", "Th^1", "Th^2", "Th^3",
      "Th^4", "Th^5", "Th^6", "Th^7", "O^1", "O^2", "O^3", "O^4",
      "O^5", "O^6", "O^7", "O^8", "Td^1", "Td^2", "Td^3", "Td^4",
      "Td^5", "Td^6", "Oh^1", "Oh^2", "Oh^3", "Oh^4", "Oh^5", "Oh^6",
      "Oh^7", "Oh^8", "Oh^9", "Oh^10",
// END_COMPILED_IN_REFERENCE_DATA
    };

    struct short_mono_hm_dict_entry
    {
      const char  *shrt;
      const char  *full;
    };

    static const short_mono_hm_dict_entry vol_i_short_mono_hm_dict[] = {
      { "P2",    "P112" },
      { "P21",   "P1121" },
      { "B2",    "B112" },
      { "C2",    "C121" },
      { "Pm",    "P11m" },
      { "Pb",    "P11b" },
      { "Pc",    "P1c1" },
      { "Bm",    "B11m" },
      { "Cm",    "C1m1" },
      { "Bb",    "B11b" },
      { "Cc",    "C1c1" },
      { "P2/m",  "P112/m" },
      { "P21/m", "P1121/m" },
      { "B2/m",  "B112/m" },
      { "C2/m",  "C12/m1" },
      { "P2/b",  "P112/b" },
      { "P2/c",  "P12/c1" },
      { "P21/b", "P1121/b" },
      { "P21/c", "P121/c1" },
      { "B2/b",  "B112/b" },
      { "C2/c",  "C12/c1" },
      { 0, 0 },
    };

    static const short_mono_hm_dict_entry vol_a_short_mono_hm_dict[] = {
      { "P2",    "P121" },
      { "P21",   "P1211" },
      { "C2",    "C121" },
      { "A2",    "A121" },
      { "I2",    "I121" },
      { "Pm",    "P1m1" },
      { "Pc",    "P1c1" },
      { "Pn",    "P1n1" },
      { "Pa",    "P1a1" },
      { "Cm",    "C1m1" },
      { "Am",    "A1m1" },
      { "Im",    "I1m1" },
      { "Cc",    "C1c1" },
      { "An",    "A1n1" },
      { "Ia",    "I1a1" },
      { "Aa",    "A1a1" },
      { "Cn",    "C1n1" },
      { "Ic",    "I1c1" },
      { "P2/m",  "P12/m1" },
      { "P21/m", "P121/m1" },
      { "C2/m",  "C12/m1" },
      { "A2/m",  "A12/m1" },
      { "I2/m",  "I12/m1" },
      { "P2/c",  "P12/c1" },
      { "P2/n",  "P12/n1" },
      { "P2/a",  "P12/a1" },
      { "P21/c", "P121/c1" },
      { "P21/n", "P121/n1" },
      { "P21/a", "P121/a1" },
      { "C2/c",  "C12/c1" },
      { "A2/n",  "A12/n1" },
      { "I2/a",  "I12/a1" },
      { "A2/a",  "A12/a1" },
      { "C2/n",  "C12/n1" },
      { "I2/c",  "I12/c1" },
      { 0, 0 },
    };

    // Qualifiers
    namespace qf {
      static const char* a   = "a";
      static const char* b   = "b";
      static const char* c   = "c";
      static const char* a1  = "a1";
      static const char* b1  = "b1";
      static const char* c1  = "c1";
      static const char* a2  = "a2";
      static const char* b2  = "b2";
      static const char* c2  = "c2";
      static const char* a3  = "a3";
      static const char* b3  = "b3";
      static const char* c3  = "c3";
      static const char* ma1 = "-a1";
      static const char* mb1 = "-b1";
      static const char* mc1 = "-c1";
      static const char* ma2 = "-a2";
      static const char* mb2 = "-b2";
      static const char* mc2 = "-c2";
      static const char* ma3 = "-a3";
      static const char* mb3 = "-b3";
      static const char* mc3 = "-c3";
      static const char* bamc = "ba-c";
      static const char* cab  = "cab";
      static const char* mcba = "-cba";
      static const char* bca  = "bca";
      static const char* amcb = "a-cb";
    }

    /* This table corresponds to a table in the
       International Tables for Crystallography, Volume B, 2001.
       The Hall symbols were generated with the
       STARTX module of the Xtal System of Crystallographic Programs,
       version 3.7 (http://xtal.crystal.uwa.edu.au/).
     */
    static const main_symbol_dict_entry main_symbol_dict[] = {
// BEGIN_COMPILED_IN_REFERENCE_DATA
      {   1, 0,        "P 1",         " P 1\0" },
      {   2, 0,        "P -1",        "-P 1\0" },
      {   3, qf::b,    "P 1 2 1",     " P 2y\0" },
      {   3, qf::c,    "P 1 1 2",     " P 2\0" },
      {   3, qf::a,    "P 2 1 1",     " P 2x\0" },
      {   4, qf::b,    "P 1 21 1",    " P 2yb\0" },
      {   4, qf::c,    "P 1 1 21",    " P 2c\0" },
      {   4, qf::a,    "P 21 1 1",    " P 2xa\0" },
      {   5, qf::b1,   "C 1 2 1",     " C 2y\0" },
      {   5, qf::b2,   "A 1 2 1",     " A 2y\0" },
      {   5, qf::b3,   "I 1 2 1",     " I 2y\0" },
      {   5, qf::c1,   "A 1 1 2",     " A 2\0" },
      {   5, qf::c2,   "B 1 1 2",     " B 2\0" },
      {   5, qf::c3,   "I 1 1 2",     " I 2\0" },
      {   5, qf::a1,   "B 2 1 1",     " B 2x\0" },
      {   5, qf::a2,   "C 2 1 1",     " C 2x\0" },
      {   5, qf::a3,   "I 2 1 1",     " I 2x\0" },
      {   6, qf::b,    "P 1 m 1",     " P -2y\0" },
      {   6, qf::c,    "P 1 1 m",     " P -2\0" },
      {   6, qf::a,    "P m 1 1",     " P -2x\0" },
      {   7, qf::b1,   "P 1 c 1",     " P -2yc\0" },
      {   7, qf::b2,   "P 1 n 1",     " P -2yac\0" },
      {   7, qf::b3,   "P 1 a 1",     " P -2ya\0" },
      {   7, qf::c1,   "P 1 1 a",     " P -2a\0" },
      {   7, qf::c2,   "P 1 1 n",     " P -2ab\0" },
      {   7, qf::c3,   "P 1 1 b",     " P -2b\0" },
      {   7, qf::a1,   "P b 1 1",     " P -2xb\0" },
      {   7, qf::a2,   "P n 1 1",     " P -2xbc\0" },
      {   7, qf::a3,   "P c 1 1",     " P -2xc\0" },
      {   8, qf::b1,   "C 1 m 1",     " C -2y\0" },
      {   8, qf::b2,   "A 1 m 1",     " A -2y\0" },
      {   8, qf::b3,   "I 1 m 1",     " I -2y\0" },
      {   8, qf::c1,   "A 1 1 m",     " A -2\0" },
      {   8, qf::c2,   "B 1 1 m",     " B -2\0" },
      {   8, qf::c3,   "I 1 1 m",     " I -2\0" },
      {   8, qf::a1,   "B m 1 1",     " B -2x\0" },
      {   8, qf::a2,   "C m 1 1",     " C -2x\0" },
      {   8, qf::a3,   "I m 1 1",     " I -2x\0" },
      {   9, qf::b1,   "C 1 c 1",     " C -2yc\0" },
      {   9, qf::b2,   "A 1 n 1",     " A -2yab\0" },
      {   9, qf::b3,   "I 1 a 1",     " I -2ya\0" },
      {   9, qf::mb1,  "A 1 a 1",     " A -2ya\0" },
      {   9, qf::mb2,  "C 1 n 1",     " C -2yac\0" },
      {   9, qf::mb3,  "I 1 c 1",     " I -2yc\0" },
      {   9, qf::c1,   "A 1 1 a",     " A -2a\0" },
      {   9, qf::c2,   "B 1 1 n",     " B -2ab\0" },
      {   9, qf::c3,   "I 1 1 b",     " I -2b\0" },
      {   9, qf::mc1,  "B 1 1 b",     " B -2b\0" },
      {   9, qf::mc2,  "A 1 1 n",     " A -2ab\0" },
      {   9, qf::mc3,  "I 1 1 a",     " I -2a\0" },
      {   9, qf::a1,   "B b 1 1",     " B -2xb\0" },
      {   9, qf::a2,   "C n 1 1",     " C -2xac\0" },
      {   9, qf::a3,   "I c 1 1",     " I -2xc\0" },
      {   9, qf::ma1,  "C c 1 1",     " C -2xc\0" },
      {   9, qf::ma2,  "B n 1 1",     " B -2xab\0" },
      {   9, qf::ma3,  "I b 1 1",     " I -2xb\0" },
      {  10, qf::b,    "P 1 2/m 1",   "-P 2y\0" },
      {  10, qf::c,    "P 1 1 2/m",   "-P 2\0" },
      {  10, qf::a,    "P 2/m 1 1",   "-P 2x\0" },
      {  11, qf::b,    "P 1 21/m 1",  "-P 2yb\0" },
      {  11, qf::c,    "P 1 1 21/m",  "-P 2c\0" },
      {  11, qf::a,    "P 21/m 1 1",  "-P 2xa\0" },
      {  12, qf::b1,   "C 1 2/m 1",   "-C 2y\0" },
      {  12, qf::b2,   "A 1 2/m 1",   "-A 2y\0" },
      {  12, qf::b3,   "I 1 2/m 1",   "-I 2y\0" },
      {  12, qf::c1,   "A 1 1 2/m",   "-A 2\0" },
      {  12, qf::c2,   "B 1 1 2/m",   "-B 2\0" },
      {  12, qf::c3,   "I 1 1 2/m",   "-I 2\0" },
      {  12, qf::a1,   "B 2/m 1 1",   "-B 2x\0" },
      {  12, qf::a2,   "C 2/m 1 1",   "-C 2x\0" },
      {  12, qf::a3,   "I 2/m 1 1",   "-I 2x\0" },
      {  13, qf::b1,   "P 1 2/c 1",   "-P 2yc\0" },
      {  13, qf::b2,   "P 1 2/n 1",   "-P 2yac\0" },
      {  13, qf::b3,   "P 1 2/a 1",   "-P 2ya\0" },
      {  13, qf::c1,   "P 1 1 2/a",   "-P 2a\0" },
      {  13, qf::c2,   "P 1 1 2/n",   "-P 2ab\0" },
      {  13, qf::c3,   "P 1 1 2/b",   "-P 2b\0" },
      {  13, qf::a1,   "P 2/b 1 1",   "-P 2xb\0" },
      {  13, qf::a2,   "P 2/n 1 1",   "-P 2xbc\0" },
      {  13, qf::a3,   "P 2/c 1 1",   "-P 2xc\0" },
      {  14, qf::b1,   "P 1 21/c 1",  "-P 2ybc\0" },
      {  14, qf::b2,   "P 1 21/n 1",  "-P 2yn\0" },
      {  14, qf::b3,   "P 1 21/a 1",  "-P 2yab\0" },
      {  14, qf::c1,   "P 1 1 21/a",  "-P 2ac\0" },
      {  14, qf::c2,   "P 1 1 21/n",  "-P 2n\0" },
      {  14, qf::c3,   "P 1 1 21/b",  "-P 2bc\0" },
      {  14, qf::a1,   "P 21/b 1 1",  "-P 2xab\0" },
      {  14, qf::a2,   "P 21/n 1 1",  "-P 2xn\0" },
      {  14, qf::a3,   "P 21/c 1 1",  "-P 2xac\0" },
      {  15, qf::b1,   "C 1 2/c 1",   "-C 2yc\0" },
      {  15, qf::b2,   "A 1 2/n 1",   "-A 2yab\0" },
      {  15, qf::b3,   "I 1 2/a 1",   "-I 2ya\0" },
      {  15, qf::mb1,  "A 1 2/a 1",   "-A 2ya\0" },
      {  15, qf::mb2,  "C 1 2/n 1",   "-C 2yac\0" },
      {  15, qf::mb3,  "I 1 2/c 1",   "-I 2yc\0" },
      {  15, qf::c1,   "A 1 1 2/a",   "-A 2a\0" },
      {  15, qf::c2,   "B 1 1 2/n",   "-B 2ab\0" },
      {  15, qf::c3,   "I 1 1 2/b",   "-I 2b\0" },
      {  15, qf::mc1,  "B 1 1 2/b",   "-B 2b\0" },
      {  15, qf::mc2,  "A 1 1 2/n",   "-A 2ab\0" },
      {  15, qf::mc3,  "I 1 1 2/a",   "-I 2a\0" },
      {  15, qf::a1,   "B 2/b 1 1",   "-B 2xb\0" },
      {  15, qf::a2,   "C 2/n 1 1",   "-C 2xac\0" },
      {  15, qf::a3,   "I 2/c 1 1",   "-I 2xc\0" },
      {  15, qf::ma1,  "C 2/c 1 1",   "-C 2xc\0" },
      {  15, qf::ma2,  "B 2/n 1 1",   "-B 2xab\0" },
      {  15, qf::ma3,  "I 2/b 1 1",   "-I 2xb\0" },
      {  16, 0,        "P 2 2 2",     " P 2 2\0" },
      {  17, 0,        "P 2 2 21",    " P 2c 2\0" },
      {  17, qf::cab,  "P 21 2 2",    " P 2a 2a\0" },
      {  17, qf::bca,  "P 2 21 2",    " P 2 2b\0" },
      {  18, 0,        "P 21 21 2",   " P 2 2ab\0" },
      {  18, qf::cab,  "P 2 21 21",   " P 2bc 2\0" },
      {  18, qf::bca,  "P 21 2 21",   " P 2ac 2ac\0" },
      {  19, 0,        "P 21 21 21",  " P 2ac 2ab\0" },
      {  20, 0,        "C 2 2 21",    " C 2c 2\0" },
      {  20, qf::cab,  "A 21 2 2",    " A 2a 2a\0" },
      {  20, qf::bca,  "B 2 21 2",    " B 2 2b\0" },
      {  21, 0,        "C 2 2 2",     " C 2 2\0" },
      {  21, qf::cab,  "A 2 2 2",     " A 2 2\0" },
      {  21, qf::bca,  "B 2 2 2",     " B 2 2\0" },
      {  22, 0,        "F 2 2 2",     " F 2 2\0" },
      {  23, 0,        "I 2 2 2",     " I 2 2\0" },
      {  24, 0,        "I 21 21 21",  " I 2b 2c\0" },
      {  25, 0,        "P m m 2",     " P 2 -2\0" },
      {  25, qf::cab,  "P 2 m m",     " P -2 2\0" },
      {  25, qf::bca,  "P m 2 m",     " P -2 -2\0" },
      {  26, 0,        "P m c 21",    " P 2c -2\0" },
      {  26, qf::bamc, "P c m 21",    " P 2c -2c\0" },
      {  26, qf::cab,  "P 21 m a",    " P -2a 2a\0" },
      {  26, qf::mcba, "P 21 a m",    " P -2 2a\0" },
      {  26, qf::bca,  "P b 21 m",    " P -2 -2b\0" },
      {  26, qf::amcb, "P m 21 b",    " P -2b -2\0" },
      {  27, 0,        "P c c 2",     " P 2 -2c\0" },
      {  27, qf::cab,  "P 2 a a",     " P -2a 2\0" },
      {  27, qf::bca,  "P b 2 b",     " P -2b -2b\0" },
      {  28, 0,        "P m a 2",     " P 2 -2a\0" },
      {  28, qf::bamc, "P b m 2",     " P 2 -2b\0" },
      {  28, qf::cab,  "P 2 m b",     " P -2b 2\0" },
      {  28, qf::mcba, "P 2 c m",     " P -2c 2\0" },
      {  28, qf::bca,  "P c 2 m",     " P -2c -2c\0" },
      {  28, qf::amcb, "P m 2 a",     " P -2a -2a\0" },
      {  29, 0,        "P c a 21",    " P 2c -2ac\0" },
      {  29, qf::bamc, "P b c 21",    " P 2c -2b\0" },
      {  29, qf::cab,  "P 21 a b",    " P -2b 2a\0" },
      {  29, qf::mcba, "P 21 c a",    " P -2ac 2a\0" },
      {  29, qf::bca,  "P c 21 b",    " P -2bc -2c\0" },
      {  29, qf::amcb, "P b 21 a",    " P -2a -2ab\0" },
      {  30, 0,        "P n c 2",     " P 2 -2bc\0" },
      {  30, qf::bamc, "P c n 2",     " P 2 -2ac\0" },
      {  30, qf::cab,  "P 2 n a",     " P -2ac 2\0" },
      {  30, qf::mcba, "P 2 a n",     " P -2ab 2\0" },
      {  30, qf::bca,  "P b 2 n",     " P -2ab -2ab\0" },
      {  30, qf::amcb, "P n 2 b",     " P -2bc -2bc\0" },
      {  31, 0,        "P m n 21",    " P 2ac -2\0" },
      {  31, qf::bamc, "P n m 21",    " P 2bc -2bc\0" },
      {  31, qf::cab,  "P 21 m n",    " P -2ab 2ab\0" },
      {  31, qf::mcba, "P 21 n m",    " P -2 2ac\0" },
      {  31, qf::bca,  "P n 21 m",    " P -2 -2bc\0" },
      {  31, qf::amcb, "P m 21 n",    " P -2ab -2\0" },
      {  32, 0,        "P b a 2",     " P 2 -2ab\0" },
      {  32, qf::cab,  "P 2 c b",     " P -2bc 2\0" },
      {  32, qf::bca,  "P c 2 a",     " P -2ac -2ac\0" },
      {  33, 0,        "P n a 21",    " P 2c -2n\0" },
      {  33, qf::bamc, "P b n 21",    " P 2c -2ab\0" },
      {  33, qf::cab,  "P 21 n b",    " P -2bc 2a\0" },
      {  33, qf::mcba, "P 21 c n",    " P -2n 2a\0" },
      {  33, qf::bca,  "P c 21 n",    " P -2n -2ac\0" },
      {  33, qf::amcb, "P n 21 a",    " P -2ac -2n\0" },
      {  34, 0,        "P n n 2",     " P 2 -2n\0" },
      {  34, qf::cab,  "P 2 n n",     " P -2n 2\0" },
      {  34, qf::bca,  "P n 2 n",     " P -2n -2n\0" },
      {  35, 0,        "C m m 2",     " C 2 -2\0" },
      {  35, qf::cab,  "A 2 m m",     " A -2 2\0" },
      {  35, qf::bca,  "B m 2 m",     " B -2 -2\0" },
      {  36, 0,        "C m c 21",    " C 2c -2\0" },
      {  36, qf::bamc, "C c m 21",    " C 2c -2c\0" },
      {  36, qf::cab,  "A 21 m a",    " A -2a 2a\0" },
      {  36, qf::mcba, "A 21 a m",    " A -2 2a\0" },
      {  36, qf::bca,  "B b 21 m",    " B -2 -2b\0" },
      {  36, qf::amcb, "B m 21 b",    " B -2b -2\0" },
      {  37, 0,        "C c c 2",     " C 2 -2c\0" },
      {  37, qf::cab,  "A 2 a a",     " A -2a 2\0" },
      {  37, qf::bca,  "B b 2 b",     " B -2b -2b\0" },
      {  38, 0,        "A m m 2",     " A 2 -2\0" },
      {  38, qf::bamc, "B m m 2",     " B 2 -2\0" },
      {  38, qf::cab,  "B 2 m m",     " B -2 2\0" },
      {  38, qf::mcba, "C 2 m m",     " C -2 2\0" },
      {  38, qf::bca,  "C m 2 m",     " C -2 -2\0" },
      {  38, qf::amcb, "A m 2 m",     " A -2 -2\0" },
      {  39, 0,        "A b m 2",     " A 2 -2b\0" },
      {  39, qf::bamc, "B m a 2",     " B 2 -2a\0" },
      {  39, qf::cab,  "B 2 c m",     " B -2a 2\0" },
      {  39, qf::mcba, "C 2 m b",     " C -2a 2\0" },
      {  39, qf::bca,  "C m 2 a",     " C -2a -2a\0" },
      {  39, qf::amcb, "A c 2 m",     " A -2b -2b\0" },
      {  40, 0,        "A m a 2",     " A 2 -2a\0" },
      {  40, qf::bamc, "B b m 2",     " B 2 -2b\0" },
      {  40, qf::cab,  "B 2 m b",     " B -2b 2\0" },
      {  40, qf::mcba, "C 2 c m",     " C -2c 2\0" },
      {  40, qf::bca,  "C c 2 m",     " C -2c -2c\0" },
      {  40, qf::amcb, "A m 2 a",     " A -2a -2a\0" },
      {  41, 0,        "A b a 2",     " A 2 -2ab\0" },
      {  41, qf::bamc, "B b a 2",     " B 2 -2ab\0" },
      {  41, qf::cab,  "B 2 c b",     " B -2ab 2\0" },
      {  41, qf::mcba, "C 2 c b",     " C -2ac 2\0" },
      {  41, qf::bca,  "C c 2 a",     " C -2ac -2ac\0" },
      {  41, qf::amcb, "A c 2 a",     " A -2ab -2ab\0" },
      {  42, 0,        "F m m 2",     " F 2 -2\0" },
      {  42, qf::cab,  "F 2 m m",     " F -2 2\0" },
      {  42, qf::bca,  "F m 2 m",     " F -2 -2\0" },
      {  43, 0,        "F d d 2",     " F 2 -2d\0" },
      {  43, qf::cab,  "F 2 d d",     " F -2d 2\0" },
      {  43, qf::bca,  "F d 2 d",     " F -2d -2d\0" },
      {  44, 0,        "I m m 2",     " I 2 -2\0" },
      {  44, qf::cab,  "I 2 m m",     " I -2 2\0" },
      {  44, qf::bca,  "I m 2 m",     " I -2 -2\0" },
      {  45, 0,        "I b a 2",     " I 2 -2c\0" },
      {  45, qf::cab,  "I 2 c b",     " I -2a 2\0" },
      {  45, qf::bca,  "I c 2 a",     " I -2b -2b\0" },
      {  46, 0,        "I m a 2",     " I 2 -2a\0" },
      {  46, qf::bamc, "I b m 2",     " I 2 -2b\0" },
      {  46, qf::cab,  "I 2 m b",     " I -2b 2\0" },
      {  46, qf::mcba, "I 2 c m",     " I -2c 2\0" },
      {  46, qf::bca,  "I c 2 m",     " I -2c -2c\0" },
      {  46, qf::amcb, "I m 2 a",     " I -2a -2a\0" },
      {  47, 0,        "P m m m",     "-P 2 2\0" },
      {  48, 0,        "P n n n",     " P 2 2 -1n\0-P 2ab 2bc\0" },
      {  49, 0,        "P c c m",     "-P 2 2c\0" },
      {  49, qf::cab,  "P m a a",     "-P 2a 2\0" },
      {  49, qf::bca,  "P b m b",     "-P 2b 2b\0" },
      {  50, 0,        "P b a n",     " P 2 2 -1ab\0-P 2ab 2b\0" },
      {  50, qf::cab,  "P n c b",     " P 2 2 -1bc\0-P 2b 2bc\0" },
      {  50, qf::bca,  "P c n a",     " P 2 2 -1ac\0-P 2a 2c\0" },
      {  51, 0,        "P m m a",     "-P 2a 2a\0" },
      {  51, qf::bamc, "P m m b",     "-P 2b 2\0" },
      {  51, qf::cab,  "P b m m",     "-P 2 2b\0" },
      {  51, qf::mcba, "P c m m",     "-P 2c 2c\0" },
      {  51, qf::bca,  "P m c m",     "-P 2c 2\0" },
      {  51, qf::amcb, "P m a m",     "-P 2 2a\0" },
      {  52, 0,        "P n n a",     "-P 2a 2bc\0" },
      {  52, qf::bamc, "P n n b",     "-P 2b 2n\0" },
      {  52, qf::cab,  "P b n n",     "-P 2n 2b\0" },
      {  52, qf::mcba, "P c n n",     "-P 2ab 2c\0" },
      {  52, qf::bca,  "P n c n",     "-P 2ab 2n\0" },
      {  52, qf::amcb, "P n a n",     "-P 2n 2bc\0" },
      {  53, 0,        "P m n a",     "-P 2ac 2\0" },
      {  53, qf::bamc, "P n m b",     "-P 2bc 2bc\0" },
      {  53, qf::cab,  "P b m n",     "-P 2ab 2ab\0" },
      {  53, qf::mcba, "P c n m",     "-P 2 2ac\0" },
      {  53, qf::bca,  "P n c m",     "-P 2 2bc\0" },
      {  53, qf::amcb, "P m a n",     "-P 2ab 2\0" },
      {  54, 0,        "P c c a",     "-P 2a 2ac\0" },
      {  54, qf::bamc, "P c c b",     "-P 2b 2c\0" },
      {  54, qf::cab,  "P b a a",     "-P 2a 2b\0" },
      {  54, qf::mcba, "P c a a",     "-P 2ac 2c\0" },
      {  54, qf::bca,  "P b c b",     "-P 2bc 2b\0" },
      {  54, qf::amcb, "P b a b",     "-P 2b 2ab\0" },
      {  55, 0,        "P b a m",     "-P 2 2ab\0" },
      {  55, qf::cab,  "P m c b",     "-P 2bc 2\0" },
      {  55, qf::bca,  "P c m a",     "-P 2ac 2ac\0" },
      {  56, 0,        "P c c n",     "-P 2ab 2ac\0" },
      {  56, qf::cab,  "P n a a",     "-P 2ac 2bc\0" },
      {  56, qf::bca,  "P b n b",     "-P 2bc 2ab\0" },
      {  57, 0,        "P b c m",     "-P 2c 2b\0" },
      {  57, qf::bamc, "P c a m",     "-P 2c 2ac\0" },
      {  57, qf::cab,  "P m c a",     "-P 2ac 2a\0" },
      {  57, qf::mcba, "P m a b",     "-P 2b 2a\0" },
      {  57, qf::bca,  "P b m a",     "-P 2a 2ab\0" },
      {  57, qf::amcb, "P c m b",     "-P 2bc 2c\0" },
      {  58, 0,        "P n n m",     "-P 2 2n\0" },
      {  58, qf::cab,  "P m n n",     "-P 2n 2\0" },
      {  58, qf::bca,  "P n m n",     "-P 2n 2n\0" },
      {  59, 0,        "P m m n",     " P 2 2ab -1ab\0-P 2ab 2a\0" },
      {  59, qf::cab,  "P n m m",     " P 2bc 2 -1bc\0-P 2c 2bc\0" },
      {  59, qf::bca,  "P m n m",     " P 2ac 2ac -1ac\0-P 2c 2a\0" },
      {  60, 0,        "P b c n",     "-P 2n 2ab\0" },
      {  60, qf::bamc, "P c a n",     "-P 2n 2c\0" },
      {  60, qf::cab,  "P n c a",     "-P 2a 2n\0" },
      {  60, qf::mcba, "P n a b",     "-P 2bc 2n\0" },
      {  60, qf::bca,  "P b n a",     "-P 2ac 2b\0" },
      {  60, qf::amcb, "P c n b",     "-P 2b 2ac\0" },
      {  61, 0,        "P b c a",     "-P 2ac 2ab\0" },
      {  61, qf::bamc, "P c a b",     "-P 2bc 2ac\0" },
      {  62, 0,        "P n m a",     "-P 2ac 2n\0" },
      {  62, qf::bamc, "P m n b",     "-P 2bc 2a\0" },
      {  62, qf::cab,  "P b n m",     "-P 2c 2ab\0" },
      {  62, qf::mcba, "P c m n",     "-P 2n 2ac\0" },
      {  62, qf::bca,  "P m c n",     "-P 2n 2a\0" },
      {  62, qf::amcb, "P n a m",     "-P 2c 2n\0" },
      {  63, 0,        "C m c m",     "-C 2c 2\0" },
      {  63, qf::bamc, "C c m m",     "-C 2c 2c\0" },
      {  63, qf::cab,  "A m m a",     "-A 2a 2a\0" },
      {  63, qf::mcba, "A m a m",     "-A 2 2a\0" },
      {  63, qf::bca,  "B b m m",     "-B 2 2b\0" },
      {  63, qf::amcb, "B m m b",     "-B 2b 2\0" },
      {  64, 0,        "C m c a",     "-C 2ac 2\0" },
      {  64, qf::bamc, "C c m b",     "-C 2ac 2ac\0" },
      {  64, qf::cab,  "A b m a",     "-A 2ab 2ab\0" },
      {  64, qf::mcba, "A c a m",     "-A 2 2ab\0" },
      {  64, qf::bca,  "B b c m",     "-B 2 2ab\0" },
      {  64, qf::amcb, "B m a b",     "-B 2ab 2\0" },
      {  65, 0,        "C m m m",     "-C 2 2\0" },
      {  65, qf::cab,  "A m m m",     "-A 2 2\0" },
      {  65, qf::bca,  "B m m m",     "-B 2 2\0" },
      {  66, 0,        "C c c m",     "-C 2 2c\0" },
      {  66, qf::cab,  "A m a a",     "-A 2a 2\0" },
      {  66, qf::bca,  "B b m b",     "-B 2b 2b\0" },
      {  67, 0,        "C m m a",     "-C 2a 2\0" },
      {  67, qf::bamc, "C m m b",     "-C 2a 2a\0" },
      {  67, qf::cab,  "A b m m",     "-A 2b 2b\0" },
      {  67, qf::mcba, "A c m m",     "-A 2 2b\0" },
      {  67, qf::bca,  "B m c m",     "-B 2 2a\0" },
      {  67, qf::amcb, "B m a m",     "-B 2a 2\0" },
      {  68, 0,        "C c c a",     " C 2 2 -1ac\0-C 2a 2ac\0" },
      {  68, qf::bamc, "C c c b",     " C 2 2 -1ac\0-C 2a 2c\0" },
      {  68, qf::cab,  "A b a a",     " A 2 2 -1ab\0-A 2a 2b\0" },
      {  68, qf::mcba, "A c a a",     " A 2 2 -1ab\0-A 2ab 2b\0" },
      {  68, qf::bca,  "B b c b",     " B 2 2 -1ab\0-B 2ab 2b\0" },
      {  68, qf::amcb, "B b a b",     " B 2 2 -1ab\0-B 2b 2ab\0" },
      {  69, 0,        "F m m m",     "-F 2 2\0" },
      {  70, 0,        "F d d d",     " F 2 2 -1d\0-F 2uv 2vw\0" },
      {  71, 0,        "I m m m",     "-I 2 2\0" },
      {  72, 0,        "I b a m",     "-I 2 2c\0" },
      {  72, qf::cab,  "I m c b",     "-I 2a 2\0" },
      {  72, qf::bca,  "I c m a",     "-I 2b 2b\0" },
      {  73, 0,        "I b c a",     "-I 2b 2c\0" },
      {  73, qf::bamc, "I c a b",     "-I 2a 2b\0" },
      {  74, 0,        "I m m a",     "-I 2b 2\0" },
      {  74, qf::bamc, "I m m b",     "-I 2a 2a\0" },
      {  74, qf::cab,  "I b m m",     "-I 2c 2c\0" },
      {  74, qf::mcba, "I c m m",     "-I 2 2b\0" },
      {  74, qf::bca,  "I m c m",     "-I 2 2a\0" },
      {  74, qf::amcb, "I m a m",     "-I 2c 2\0" },
      {  75, 0,        "P 4",         " P 4\0" },
      {  76, 0,        "P 41",        " P 4w\0" },
      {  77, 0,        "P 42",        " P 4c\0" },
      {  78, 0,        "P 43",        " P 4cw\0" },
      {  79, 0,        "I 4",         " I 4\0" },
      {  80, 0,        "I 41",        " I 4bw\0" },
      {  81, 0,        "P -4",        " P -4\0" },
      {  82, 0,        "I -4",        " I -4\0" },
      {  83, 0,        "P 4/m",       "-P 4\0" },
      {  84, 0,        "P 42/m",      "-P 4c\0" },
      {  85, 0,        "P 4/n",       " P 4ab -1ab\0-P 4a\0" },
      {  86, 0,        "P 42/n",      " P 4n -1n\0-P 4bc\0" },
      {  87, 0,        "I 4/m",       "-I 4\0" },
      {  88, 0,        "I 41/a",      " I 4bw -1bw\0-I 4ad\0" },
      {  89, 0,        "P 4 2 2",     " P 4 2\0" },
      {  90, 0,        "P 4 21 2",    " P 4ab 2ab\0" },
      {  91, 0,        "P 41 2 2",    " P 4w 2c\0" },
      {  92, 0,        "P 41 21 2",   " P 4abw 2nw\0" },
      {  93, 0,        "P 42 2 2",    " P 4c 2\0" },
      {  94, 0,        "P 42 21 2",   " P 4n 2n\0" },
      {  95, 0,        "P 43 2 2",    " P 4cw 2c\0" },
      {  96, 0,        "P 43 21 2",   " P 4nw 2abw\0" },
      {  97, 0,        "I 4 2 2",     " I 4 2\0" },
      {  98, 0,        "I 41 2 2",    " I 4bw 2bw\0" },
      {  99, 0,        "P 4 m m",     " P 4 -2\0" },
      { 100, 0,        "P 4 b m",     " P 4 -2ab\0" },
      { 101, 0,        "P 42 c m",    " P 4c -2c\0" },
      { 102, 0,        "P 42 n m",    " P 4n -2n\0" },
      { 103, 0,        "P 4 c c",     " P 4 -2c\0" },
      { 104, 0,        "P 4 n c",     " P 4 -2n\0" },
      { 105, 0,        "P 42 m c",    " P 4c -2\0" },
      { 106, 0,        "P 42 b c",    " P 4c -2ab\0" },
      { 107, 0,        "I 4 m m",     " I 4 -2\0" },
      { 108, 0,        "I 4 c m",     " I 4 -2c\0" },
      { 109, 0,        "I 41 m d",    " I 4bw -2\0" },
      { 110, 0,        "I 41 c d",    " I 4bw -2c\0" },
      { 111, 0,        "P -4 2 m",    " P -4 2\0" },
      { 112, 0,        "P -4 2 c",    " P -4 2c\0" },
      { 113, 0,        "P -4 21 m",   " P -4 2ab\0" },
      { 114, 0,        "P -4 21 c",   " P -4 2n\0" },
      { 115, 0,        "P -4 m 2",    " P -4 -2\0" },
      { 116, 0,        "P -4 c 2",    " P -4 -2c\0" },
      { 117, 0,        "P -4 b 2",    " P -4 -2ab\0" },
      { 118, 0,        "P -4 n 2",    " P -4 -2n\0" },
      { 119, 0,        "I -4 m 2",    " I -4 -2\0" },
      { 120, 0,        "I -4 c 2",    " I -4 -2c\0" },
      { 121, 0,        "I -4 2 m",    " I -4 2\0" },
      { 122, 0,        "I -4 2 d",    " I -4 2bw\0" },
      { 123, 0,        "P 4/m m m",   "-P 4 2\0" },
      { 124, 0,        "P 4/m c c",   "-P 4 2c\0" },
      { 125, 0,        "P 4/n b m",   " P 4 2 -1ab\0-P 4a 2b\0" },
      { 126, 0,        "P 4/n n c",   " P 4 2 -1n\0-P 4a 2bc\0" },
      { 127, 0,        "P 4/m b m",   "-P 4 2ab\0" },
      { 128, 0,        "P 4/m n c",   "-P 4 2n\0" },
      { 129, 0,        "P 4/n m m",   " P 4ab 2ab -1ab\0-P 4a 2a\0" },
      { 130, 0,        "P 4/n c c",   " P 4ab 2n -1ab\0-P 4a 2ac\0" },
      { 131, 0,        "P 42/m m c",  "-P 4c 2\0" },
      { 132, 0,        "P 42/m c m",  "-P 4c 2c\0" },
      { 133, 0,        "P 42/n b c",  " P 4n 2c -1n\0-P 4ac 2b\0" },
      { 134, 0,        "P 42/n n m",  " P 4n 2 -1n\0-P 4ac 2bc\0" },
      { 135, 0,        "P 42/m b c",  "-P 4c 2ab\0" },
      { 136, 0,        "P 42/m n m",  "-P 4n 2n\0" },
      { 137, 0,        "P 42/n m c",  " P 4n 2n -1n\0-P 4ac 2a\0" },
      { 138, 0,        "P 42/n c m",  " P 4n 2ab -1n\0-P 4ac 2ac\0" },
      { 139, 0,        "I 4/m m m",   "-I 4 2\0" },
      { 140, 0,        "I 4/m c m",   "-I 4 2c\0" },
      { 141, 0,        "I 41/a m d",  " I 4bw 2bw -1bw\0-I 4bd 2\0" },
      { 142, 0,        "I 41/a c d",  " I 4bw 2aw -1bw\0-I 4bd 2c\0" },
      { 143, 0,        "P 3",         " P 3\0" },
      { 144, 0,        "P 31",        " P 31\0" },
      { 145, 0,        "P 32",        " P 32\0" },
      { 146, 0,        "R 3",         " R 3\0 P 3*\0" },
      { 147, 0,        "P -3",        "-P 3\0" },
      { 148, 0,        "R -3",        "-R 3\0-P 3*\0" },
      { 149, 0,        "P 3 1 2",     " P 3 2\0" },
      { 150, 0,        "P 3 2 1",     " P 3 2\"\0" },
      { 151, 0,        "P 31 1 2",    " P 31 2 (0 0 4)\0" },
      { 152, 0,        "P 31 2 1",    " P 31 2\"\0" },
      { 153, 0,        "P 32 1 2",    " P 32 2 (0 0 2)\0" },
      { 154, 0,        "P 32 2 1",    " P 32 2\"\0" },
      { 155, 0,        "R 3 2",       " R 3 2\"\0 P 3* 2\0" },
      { 156, 0,        "P 3 m 1",     " P 3 -2\"\0" },
      { 157, 0,        "P 3 1 m",     " P 3 -2\0" },
      { 158, 0,        "P 3 c 1",     " P 3 -2\"c\0" },
      { 159, 0,        "P 3 1 c",     " P 3 -2c\0" },
      { 160, 0,        "R 3 m",       " R 3 -2\"\0 P 3* -2\0" },
      { 161, 0,        "R 3 c",       " R 3 -2\"c\0 P 3* -2n\0" },
      { 162, 0,        "P -3 1 m",    "-P 3 2\0" },
      { 163, 0,        "P -3 1 c",    "-P 3 2c\0" },
      { 164, 0,        "P -3 m 1",    "-P 3 2\"\0" },
      { 165, 0,        "P -3 c 1",    "-P 3 2\"c\0" },
      { 166, 0,        "R -3 m",      "-R 3 2\"\0-P 3* 2\0" },
      { 167, 0,        "R -3 c",      "-R 3 2\"c\0-P 3* 2n\0" },
      { 168, 0,        "P 6",         " P 6\0" },
      { 169, 0,        "P 61",        " P 61\0" },
      { 170, 0,        "P 65",        " P 65\0" },
      { 171, 0,        "P 62",        " P 62\0" },
      { 172, 0,        "P 64",        " P 64\0" },
      { 173, 0,        "P 63",        " P 6c\0" },
      { 174, 0,        "P -6",        " P -6\0" },
      { 175, 0,        "P 6/m",       "-P 6\0" },
      { 176, 0,        "P 63/m",      "-P 6c\0" },
      { 177, 0,        "P 6 2 2",     " P 6 2\0" },
      { 178, 0,        "P 61 2 2",    " P 61 2 (0 0 5)\0" },
      { 179, 0,        "P 65 2 2",    " P 65 2 (0 0 1)\0" },
      { 180, 0,        "P 62 2 2",    " P 62 2 (0 0 4)\0" },
      { 181, 0,        "P 64 2 2",    " P 64 2 (0 0 2)\0" },
      { 182, 0,        "P 63 2 2",    " P 6c 2c\0" },
      { 183, 0,        "P 6 m m",     " P 6 -2\0" },
      { 184, 0,        "P 6 c c",     " P 6 -2c\0" },
      { 185, 0,        "P 63 c m",    " P 6c -2\0" },
      { 186, 0,        "P 63 m c",    " P 6c -2c\0" },
      { 187, 0,        "P -6 m 2",    " P -6 2\0" },
      { 188, 0,        "P -6 c 2",    " P -6c 2\0" },
      { 189, 0,        "P -6 2 m",    " P -6 -2\0" },
      { 190, 0,        "P -6 2 c",    " P -6c -2c\0" },
      { 191, 0,        "P 6/m m m",   "-P 6 2\0" },
      { 192, 0,        "P 6/m c c",   "-P 6 2c\0" },
      { 193, 0,        "P 63/m c m",  "-P 6c 2\0" },
      { 194, 0,        "P 63/m m c",  "-P 6c 2c\0" },
      { 195, 0,        "P 2 3",       " P 2 2 3\0" },
      { 196, 0,        "F 2 3",       " F 2 2 3\0" },
      { 197, 0,        "I 2 3",       " I 2 2 3\0" },
      { 198, 0,        "P 21 3",      " P 2ac 2ab 3\0" },
      { 199, 0,        "I 21 3",      " I 2b 2c 3\0" },
      { 200, 0,        "P m -3",      "-P 2 2 3\0" },
      { 201, 0,        "P n -3",      " P 2 2 3 -1n\0-P 2ab 2bc 3\0" },
      { 202, 0,        "F m -3",      "-F 2 2 3\0" },
      { 203, 0,        "F d -3",      " F 2 2 3 -1d\0-F 2uv 2vw 3\0" },
      { 204, 0,        "I m -3",      "-I 2 2 3\0" },
      { 205, 0,        "P a -3",      "-P 2ac 2ab 3\0" },
      { 206, 0,        "I a -3",      "-I 2b 2c 3\0" },
      { 207, 0,        "P 4 3 2",     " P 4 2 3\0" },
      { 208, 0,        "P 42 3 2",    " P 4n 2 3\0" },
      { 209, 0,        "F 4 3 2",     " F 4 2 3\0" },
      { 210, 0,        "F 41 3 2",    " F 4d 2 3\0" },
      { 211, 0,        "I 4 3 2",     " I 4 2 3\0" },
      { 212, 0,        "P 43 3 2",    " P 4acd 2ab 3\0" },
      { 213, 0,        "P 41 3 2",    " P 4bd 2ab 3\0" },
      { 214, 0,        "I 41 3 2",    " I 4bd 2c 3\0" },
      { 215, 0,        "P -4 3 m",    " P -4 2 3\0" },
      { 216, 0,        "F -4 3 m",    " F -4 2 3\0" },
      { 217, 0,        "I -4 3 m",    " I -4 2 3\0" },
      { 218, 0,        "P -4 3 n",    " P -4n 2 3\0" },
      { 219, 0,        "F -4 3 c",    " F -4a 2 3\0" },
      { 220, 0,        "I -4 3 d",    " I -4bd 2c 3\0" },
      { 221, 0,        "P m -3 m",    "-P 4 2 3\0" },
      { 222, 0,        "P n -3 n",    " P 4 2 3 -1n\0-P 4a 2bc 3\0" },
      { 223, 0,        "P m -3 n",    "-P 4n 2 3\0" },
      { 224, 0,        "P n -3 m",    " P 4n 2 3 -1n\0-P 4bc 2bc 3\0" },
      { 225, 0,        "F m -3 m",    "-F 4 2 3\0" },
      { 226, 0,        "F m -3 c",    "-F 4a 2 3\0" },
      { 227, 0,        "F d -3 m",    " F 4d 2 3 -1d\0-F 4vw 2vw 3\0" },
      { 228, 0,        "F d -3 c",    " F 4d 2 3 -1ad\0-F 4ud 2vw 3\0" },
      { 229, 0,        "I m -3 m",    "-I 4 2 3\0" },
      { 230, 0,        "I a -3 d",    "-I 4bd 2c 3\0" },
      { 0, 0, 0, 0 },
// END_COMPILED_IN_REFERENCE_DATA
    };

  } // namespace tables

  // remove whitespace and underscores, map to lower case
  string pre_process_symbol(string const& raw_symbol)
  {
    string result;
    for(std::size_t i=0;i<raw_symbol.size();i++) {
      const char r = raw_symbol[i];
      if (!isspace(r) && r != '_') result += tolower(r);
    }
    return result;
  }

  char strip_extension(string& work_symbol)
  {
    char ext = '\0';
    string::size_type stop = work_symbol.find(':');
    if (stop != string::npos) {
      string e = work_symbol.substr(stop + 1);
      if (e.size() == 1) {
        ext = e[0];
      }
      else if (e == "o1" || e == "o2") {
        ext = e[1];
      }
    }
    else {
      // check if last character is S, Z, R or H
      if (work_symbol.size() > 0) {
        string::size_type i = work_symbol.size() - 1;
        switch (work_symbol[i]) {
          case 's':
          case 'z':
          case 'r':
          case 'h':
            ext = work_symbol[i];
            stop = i;
            break;
        }
        // check if last two characters are O1 or O2
        if (ext == '\0' && work_symbol.size() > 1) {
          string::size_type i = work_symbol.size() - 2;
          string e = work_symbol.substr(i);
          if (e == "o1" || e == "o2") {
            ext = e[1];
            stop = i;
          }
        }
      }
    }

    if (ext != '\0') {
        switch (ext) {
        case '1': break;
        case 's': ext = '1'; break;
        case '2': break;
        case 'z': ext = '2'; break;
        case 'r': ext = 'R'; break;
        case 'h': ext = 'H'; break;
        default:
          ext = '\0';
      }
    }

    if (ext != '\0') work_symbol.erase(stop);

    return ext;
  }

  // remove parentheses, e.g. "P2(1)2(1)2(1)" -> "P212121"
  void remove_screw_component_parentheses(string& work_symbol)
  {
    static const int  rot_numbers[] = { 2, 3, 4, 6 };
    static const char rot_symbols[] = "2346";
    static const char scr_symbols[] = "012345";
    string pat = "r(s)";
    int ir;
    for(ir=0;ir<4;ir++) {
      pat[0] = rot_symbols[ir];
      int is;
      for (is = 1; is < rot_numbers[ir]; is++) {
        pat[2] = scr_symbols[is];
        for (;;) {
          string::size_type i = work_symbol.find(pat);
          if (i == string::npos) break;
          work_symbol.erase(i + 3, 1);
          work_symbol.erase(i + 1, 1);
        }
      }
    }
  }

  string
  split_off_cb_symbol(string& work_symbol)
  {
    string result;
    string::size_type sz = work_symbol.size();
    if (sz != 0 && work_symbol[sz-1] == ')') {
      string::size_type i = work_symbol.rfind("(");
      if (i != string::npos && i >= 2) {
        result = work_symbol.substr(i+1, sz-i-2);
        work_symbol.resize(i);
      }
    }
    return result;
  }

  string remove_spaces(string const& inp)
  {
    string result;
    for(std::size_t i=0;i<inp.size();i++) {
      if (inp[i] != ' ') result += inp[i];
    }
    return result;
  }

  int cmp_schoenflies_symbols(string const& from_table,
                              string const& work_symbol)
  {
    if (from_table.size() != work_symbol.size()) return -1;
    for(std::size_t i=0;i<from_table.size();i++) {
        if (    from_table[i] != work_symbol[i]
            && (from_table[i] != '^' || isalpha(work_symbol[i])
                                     || isdigit(work_symbol[i])))
          return -1;
    }
    return 0;
  }

  int schoenflies_as_sg_number(string const& work_symbol)
  {
    for (int sg_number = 1; sg_number <= 230; sg_number++) {
      if (cmp_schoenflies_symbols(tables::schoenflies_list[sg_number],
                                  work_symbol) == 0)
        return sg_number;
    }
    return 0;
  }

  string get_std_table_id(string const& table_id)
  {
    string std_table_id = remove_spaces(table_id);
    if (std_table_id.size() > 0) {
      std_table_id[0] = toupper(std_table_id[0]);
      if (   std_table_id == "I"
          || std_table_id == "I1952"
          || std_table_id == "1")
        std_table_id = "I1952";
      else if (std_table_id == "A" || std_table_id == "A1983")
        std_table_id = "A1983";
      else
        throw error("table_id not recognized: " + table_id);
    }
    return std_table_id;
  }

  bool
  short_mono_hm_as_full_mono_hm(
    string const& std_table_id,
    string& work_symbol)
  {
    const tables::short_mono_hm_dict_entry* entry;
    if (std_table_id == "I1952") entry = tables::vol_i_short_mono_hm_dict;
    else                         entry = tables::vol_a_short_mono_hm_dict;
    for (int i = 0; entry[i].shrt; i++) {
      if (work_symbol == entry[i].shrt) {
        work_symbol = entry[i].full;
        return true;
      }
    }
    return false;
  }

  void
  ad_hoc_1992_symbol_as_a1983_symbol(
    string& work_symbol)
  {
    // see documentation of space_group_type::lookup_symbol()
    // simlar table in space_group_type.cpp, space_group_type::lookup_symbol()
    static char const* adh_a38_pairs[] = {
      // No. 39
      "Aem2", "Abm2",
      "Bme2", "Bma2",
      "B2em", "B2cm",
      "C2me", "C2mb",
      "Cm2e", "Cm2a",
      "Ae2m", "Ac2m",
      // No. 41
      "Aea2", "Aba2",
      "Bbe2", "Bba2",
      "B2eb", "B2cb",
      "C2ce", "C2cb",
      "Cc2e", "Cc2a",
      "Ae2a", "Ac2a",
      // No. 64
      "Cmce", "Cmca",
      "Ccme", "Ccmb",
      "Aema", "Abma",
      "Aeam", "Acam",
      "Bbem", "Bbcm",
      "Bmeb", "Bmab",
      // No. 67
      "Cmme", "Cmma", // ambiguous: Cmmb
      "Aemm", "Abmm", // ambiguous: Acmm
      "Bmem", "Bmcm", // ambiguous: Bmam
      // No. 68
      "Ccce", "Ccca", // ambiguous: Cccb
      "Aeaa", "Abaa", // ambiguous: Acaa
      "Bbeb", "Bbcb"  // ambiguous: Bbab
    };
    typedef std::map<std::string, const char*> map_t;
    static map_t adh_a38_map;
    if (adh_a38_map.size() == 0) {
      std::size_t n = sizeof(adh_a38_pairs) / sizeof(const char*);
      for(int i=0;i<n;i+=2) {
        adh_a38_map[adh_a38_pairs[i]] = adh_a38_pairs[i+1];
      }
      CCTBX_ASSERT(adh_a38_map.size()*2 == n);
    }
    map_t::const_iterator match = adh_a38_map.find(work_symbol);
    if (match != adh_a38_map.end()) work_symbol = match->second;
  }

  const tables::main_symbol_dict_entry*
  find_main_symbol_dict_entry(string const& work_symbol)
  {
    const tables::main_symbol_dict_entry* entry;
    for (entry = tables::main_symbol_dict; entry->sg_number; entry++) {
      string hm = remove_spaces(entry->hermann_mauguin);
      if (hm == work_symbol) return entry;
      string::size_type bar = hm.find('-');
      if (bar != string::npos) {
        // reverse "-N" to "N-", e.g. "P -1" -> "P 1-"
        hm[bar] = hm[bar + 1];
        hm[bar + 1] = '-';
        if (hm == work_symbol) return entry;
        if (   (entry->sg_number >= 200 && entry->sg_number <= 206)
            || (entry->sg_number >= 221 && entry->sg_number <= 230)) {
          hm.erase(bar + 1, 1); // remove '-', e.g. "P m -3 m" -> "P m 3 m"
          if (hm == work_symbol) return entry;
        }
      }
    }
    return 0;
  }

  const tables::main_symbol_dict_entry*
  sg_number_to_main_symbol_dict_entry(int sg_number,
                                      string const& std_table_id)
  {
    if (sg_number < 1 || sg_number > 230) {
      throw error("Space group number out of range: " + to_str(sg_number));
    }
    if (sg_number < 3 || sg_number > 15) {
      const tables::main_symbol_dict_entry* entry;
      for (entry = tables::main_symbol_dict; entry->sg_number; entry++)
        if (entry->sg_number == sg_number)
          return entry;
      throw CCTBX_INTERNAL_ERROR();
    }
    int i = 0;
    if (std_table_id == "I1952") i = 1;
    string mono_hm = tables::monoclinic_sg_number_as_hm_list[sg_number][i];
    const tables::main_symbol_dict_entry*
      result = find_main_symbol_dict_entry(mono_hm);
    CCTBX_ASSERT(result != 0);
    return result;
  }

  const char* select_hall(const tables::main_symbol_dict_entry* entry,
                          char& work_extension,
                          string const& std_table_id)
  {
    const char* hall2 = &entry->hall[string(entry->hall).size() + 1];
    if (string(hall2).size() == 0) {
      if (work_extension == '\0') return entry->hall;
    }
    else {
      if (entry->hermann_mauguin[0] == 'R') {
        if (work_extension == '\0') {
          if (std_table_id == "I1952") work_extension = 'R';
          else                         work_extension = 'H';
        }
        if      (work_extension == 'H') return entry->hall;
        else if (work_extension == 'R') return hall2;
      }
      else {
        if (work_extension == '\0') {
          if (std_table_id.size() == 0) work_extension = '2';
          else                          work_extension = '1';
        }
        if      (work_extension == '1') return entry->hall;
        else if (work_extension == '2') return hall2;
      }
    }
    return 0;
  }

  } // namespace symbols

  using namespace symbols;

  bool
  space_group_symbols::set_all(
    const symbols::tables::main_symbol_dict_entry* entry,
    char work_extension,
    string const& std_table_id)
  {
    const char* table_hall = select_hall(entry, work_extension, std_table_id);
    if (table_hall == 0) return false;
    CCTBX_ASSERT(   work_extension == '\0'
                 || work_extension == '1' || work_extension == '2'
                 || work_extension == 'H' || work_extension == 'R');
    number_ = entry->sg_number;
    schoenflies_ = symbols::tables::schoenflies_list[entry->sg_number];
    qualifier_ = string(entry->qualifier ? entry->qualifier : "");
    hermann_mauguin_ = entry->hermann_mauguin;
    extension_ = work_extension;
    change_of_basis_symbol_ = "";
    universal_hermann_mauguin_ = hermann_mauguin_;
    if (extension_ != '\0') {
      universal_hermann_mauguin_ += " :";
      universal_hermann_mauguin_ += extension_;
    }
    hall_ = table_hall;
    return true;
  }

  int space_group_symbols::hall_pass_through(string const& symbol)
  {
    string::const_iterator s;
    for (s = symbol.begin(); s != symbol.end(); s++) {
      if (!isspace(*s)) break;
    }
    for (int i = 0; i < 4; i++, s++) {
      if (s == symbol.end()) return 0;
      if (tolower(*s) != "hall"[i]) return 0;
    }
    if (*s != ':' && !isspace(*s)) return 0;
    for (; s != symbol.end(); s++) if (!isspace(*s)) break;
    if (s == symbol.end()) return 0;
    if (*s == ':') {
      for (s++; s != symbol.end(); s++) if (!isspace(*s)) break;
      if (s == symbol.end()) return 0;
    }
    hall_ = string(s, symbol.end());
    return 1;
  }

  void space_group_symbols::clear()
  {
    number_ = 0;
    schoenflies_ = "";
    qualifier_ = "";
    hermann_mauguin_ = "";
    extension_ = '\0';
    change_of_basis_symbol_ = "";
    universal_hermann_mauguin_ = "";
    hall_ = "";
  }

  namespace {

    bool is_hall(string const& table_id)
    {
      string::const_iterator s;
      for (s = table_id.begin(); s != table_id.end(); s++) {
        if (!isspace(*s)) break;
      }
      for (std::size_t i = 0; i < 4; i++, s++) {
        if (s == table_id.end()) return false;
        if (tolower(*s) != "hall"[i]) return false;
      }
      for (;; s++) {
        if (s == table_id.end()) return true;
        if (!isspace(*s)) break;
      }
      return false;
    }

  }

  space_group_symbols::space_group_symbols(string const& symbol,
                                           string const& table_id)
  {
    clear();
    if (is_hall(table_id)) {
      hall_ = symbol;
      return;
    }
    string std_table_id = get_std_table_id(table_id);
    if (std_table_id.size() == 0) {
      if (hall_pass_through(symbol) != 0) return;
    }
    string work_symbol = pre_process_symbol(symbol);
    work_symbol[0] = toupper(work_symbol[0]);
    if (std_table_id.size() == 0 && work_symbol[0] == 'H') {
      if (work_symbol == "H3") {
        work_symbol = "R3";
      }
      else if (work_symbol == "H32") {
        work_symbol = "R32";
      }
    }
    remove_screw_component_parentheses(work_symbol);
    change_of_basis_op cb_op(0, 0);
    std::string cb_mx_symbol = split_off_cb_symbol(work_symbol);
    if (cb_mx_symbol.size() != 0) {
      cb_op = change_of_basis_op(cb_mx_symbol);
    }
    char work_extension = strip_extension(work_symbol);
    const symbols::tables::main_symbol_dict_entry* entry = 0;
    int sg_no;
    char xtrac;
    int n = sscanf(work_symbol.c_str(), "%d%c", &sg_no, &xtrac);
    if (n == 1) {
      entry = sg_number_to_main_symbol_dict_entry(sg_no, std_table_id);
    }
    else {
      sg_no = schoenflies_as_sg_number(work_symbol);
      if (sg_no != 0) {
        try {
          entry = sg_number_to_main_symbol_dict_entry(sg_no, std_table_id);
        }
        catch (error const&) {
          throw CCTBX_INTERNAL_ERROR();
        }
      }
    }
    if (entry == 0) {
      if (!short_mono_hm_as_full_mono_hm(std_table_id, work_symbol)) {
        ad_hoc_1992_symbol_as_a1983_symbol(work_symbol);
      }
      entry = find_main_symbol_dict_entry(work_symbol);
      if (entry == 0) {
        throw error(not_recognized + symbol);
      }
    }
    if (!set_all(entry, work_extension, std_table_id)) {
      throw error(not_recognized + symbol);
    }
    if (cb_op.is_valid()) {
      change_of_basis_symbol_ = cb_op.as_xyz();
      universal_hermann_mauguin_ += " (" + change_of_basis_symbol_ + ")";
      cb_mx_symbol = split_off_cb_symbol(hall_);
      if (cb_mx_symbol.size() != 0) {
        hall_.resize(hall_.size()-1); // strip space
        // All tabulated Hall symbol change-of-basis operators are of
        // the form (0 0 i), with i = {1,2,4,5}
        CCTBX_ASSERT(cb_mx_symbol.size() == 5);
        CCTBX_ASSERT(cb_mx_symbol.substr(0, 4) == "0 0 ");
        int t12;
        switch (cb_mx_symbol[4]) {
          case '1': t12 = 1; break;
          case '2': t12 = 2; break;
          case '4': t12 = 4; break;
          case '5': t12 = 5; break;
          default:
            throw CCTBX_INTERNAL_ERROR();
        }
        cb_op = change_of_basis_op(rt_mx(
          cb_op.c().r(),
          cb_op.c().t().plus(cb_op.c().r().multiply(tr_vec(0,0,t12,12)))));
      }
      cb_op.mod_short_in_place();
      if (!cb_op.is_identity_op()) {
        hall_ += " (" + cb_op.as_xyz() + ")";
      }
    }
  }

  space_group_symbols::space_group_symbols(
    int space_group_number,
    string const& extension,
    string const& table_id)
  {
    clear();
    string std_table_id = get_std_table_id(table_id);
    string work_symbol = pre_process_symbol(extension);
    if (work_symbol.size() != 0 && work_symbol[0] != ':') {
      work_symbol.insert(0, ":");
    }
    string inp_work_symbol = work_symbol;
    char work_extension = strip_extension(work_symbol);
    if (work_symbol.size() != 0) {
      throw error(
        not_recognized + to_str(space_group_number) + inp_work_symbol);
    }
    const symbols::tables::main_symbol_dict_entry* entry = 0;
    entry = sg_number_to_main_symbol_dict_entry(
      space_group_number, std_table_id);
    if (!set_all(entry, work_extension, std_table_id)) {
      throw error(
        not_recognized + to_str(space_group_number) + inp_work_symbol);
    }
  }

  space_group_symbols::space_group_symbols(
    const symbols::tables::main_symbol_dict_entry* entry,
    char extension)
  {
    clear();
    if (entry->sg_number) {
      CCTBX_ASSERT(set_all(entry, extension, ""));
    }
  }

  matrix_group::code
  space_group_symbols::point_group_type() const
  {
    std::size_t space_group_number = number();
    CCTBX_ASSERT(space_group_number >= 1);
    CCTBX_ASSERT(space_group_number <= 230);
    return reference_settings::matrix_group_code_table(space_group_number);
  }

  space_group_symbol_iterator::space_group_symbol_iterator()
  : entry_(symbols::tables::main_symbol_dict), n_hall_(1), i_hall_(0)
  {}

  space_group_symbols space_group_symbol_iterator::next()
  {
    const symbols::tables::main_symbol_dict_entry* current_entry = entry_;
    char extension = '\0';
    if (entry_->sg_number) {
      if (n_hall_ == 2) {
        if (143 <= entry_->sg_number && entry_->sg_number < 168) {
          extension = "HR"[i_hall_];
        }
        else {
          extension = "12"[i_hall_];
        }
      }
      i_hall_++;
      if (i_hall_ == n_hall_) {
        entry_++;
        n_hall_ = 0;
        i_hall_ = 0;
        if (entry_->sg_number) {
          n_hall_ = 1;
          const char* hall2 = &entry_->hall[string(entry_->hall).size() + 1];
          if (string(hall2).size() != 0) n_hall_ = 2;
        }
      }
    }
    return space_group_symbols(current_entry, extension);
  }

}} // namespace cctbx::sgtbx
