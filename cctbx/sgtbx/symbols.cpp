// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 02: six Hall symbols corrected (R.W. Grosse-Kunstleve)
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <ctype.h> // cannot use cctype b/o non-conforming compilers
#include <string>
#include <stdio.h>
#include <cctbx/sgtbx/basic.h>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/basic/define_range.h>

using std::string;

namespace {
  static const sgtbx::error
  symbol_not_recognized("Space group symbol not recognized.");
}

namespace sgtbx {
  namespace symbols {
    namespace tables {

      static const char* Monoclinic_SgNumber_as_HM_List[][2] = {
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

      static const char* Schoenflies_List[] = {
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

      struct Short_Mono_HM_Dict_Entry {
        const char  *shrt;
        const char  *full;
      };

      static const Short_Mono_HM_Dict_Entry VolI_Short_Mono_HM_Dict[] = {
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

      static const Short_Mono_HM_Dict_Entry VolA_Short_Mono_HM_Dict[] = {
        { "P2",    "P121" },
        { "P21",   "P1211" },
        { "C2",    "C121" },
        { "Pm",    "P1m1" },
        { "Pc",    "P1c1" },
        { "Cm",    "C1m1" },
        { "Cc",    "C1c1" },
        { "P2/m",  "P12/m1" },
        { "P21/m", "P121/m1" },
        { "C2/m",  "C12/m1" },
        { "P2/c",  "P12/c1" },
        { "P21/c", "P121/c1" },
        { "C2/c",  "C12/c1" },
        { 0, 0 },
      };

      // Qualifiers
      namespace Qf {
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
        static const char* MC = "MC"; // Multiple cell C
        static const char* MF = "MF"; // Multiple cell F
        static const char* TH = "TH"; // Triple cell H
      }

      /* This tables corresponds to a table in the
         International Tables for Crystallography, Volume B, 2001.
         The Hall symbols were generated with the
         STARTX module of the Xtal System of Crystallographic Programs,
         version 3.7 (http://xtal.crystal.uwa.edu.au/).
         10-Aug-2001, 02-Oct-2001:
         Additional 126 settings from the "Multiple cell C or F" and
         the "Triple cell H" columns of Table 4.3.1 in Int. Tab. A1983.
         The multiple cell settings were generated by applying
         the following change-of-basis matrices to the standard
         settings:
           Multiple cell C or F: 1/2*x+1/2*y,-1/2*x+1/2*y,z
           Triple cell H: 2/3*x-1/3*y,1/3*x+1/3*y,z
         The Hall symbols for "Multiple cell C or F" were computed
         with Xtal 3.7. The Hall symbols for "Triple cell H" were
         derived manually.
         Reference for "Multiple cell C or F" and "Triple cell H":
           Internationale Tabellen zur Bestimmung von Kristallstrukturen,
           Erster Band, Gruppentheoretische Tafeln, 1935.
           The change-of-basis matrices are defined on page 30.
         The settings below are identical to that used
         by ALTWYK (http://ylp.icpet.nrc.ca/altwyk/).
       */
      static const Main_Symbol_Dict_Entry Main_Symbol_Dict[] = {
// BEGIN_COMPILED_IN_REFERENCE_DATA
        {   1, 0,        "P 1",         " P 1\0" },
        {   2, 0,        "P -1",        "-P 1\0" },
        {   3, Qf::b,    "P 1 2 1",     " P 2y\0" },
        {   3, Qf::c,    "P 1 1 2",     " P 2\0" },
        {   3, Qf::a,    "P 2 1 1",     " P 2x\0" },
        {   4, Qf::b,    "P 1 21 1",    " P 2yb\0" },
        {   4, Qf::c,    "P 1 1 21",    " P 2c\0" },
        {   4, Qf::a,    "P 21 1 1",    " P 2xa\0" },
        {   5, Qf::b1,   "C 1 2 1",     " C 2y\0" },
        {   5, Qf::b2,   "A 1 2 1",     " A 2y\0" },
        {   5, Qf::b3,   "I 1 2 1",     " I 2y\0" },
        {   5, Qf::c1,   "A 1 1 2",     " A 2\0" },
        {   5, Qf::c2,   "B 1 1 2",     " B 2\0" },
        {   5, Qf::c3,   "I 1 1 2",     " I 2\0" },
        {   5, Qf::a1,   "B 2 1 1",     " B 2x\0" },
        {   5, Qf::a2,   "C 2 1 1",     " C 2x\0" },
        {   5, Qf::a3,   "I 2 1 1",     " I 2x\0" },
        {   6, Qf::b,    "P 1 m 1",     " P -2y\0" },
        {   6, Qf::c,    "P 1 1 m",     " P -2\0" },
        {   6, Qf::a,    "P m 1 1",     " P -2x\0" },
        {   7, Qf::b1,   "P 1 c 1",     " P -2yc\0" },
        {   7, Qf::b2,   "P 1 n 1",     " P -2yac\0" },
        {   7, Qf::b3,   "P 1 a 1",     " P -2ya\0" },
        {   7, Qf::c1,   "P 1 1 a",     " P -2a\0" },
        {   7, Qf::c2,   "P 1 1 n",     " P -2ab\0" },
        {   7, Qf::c3,   "P 1 1 b",     " P -2b\0" },
        {   7, Qf::a1,   "P b 1 1",     " P -2xb\0" },
        {   7, Qf::a2,   "P n 1 1",     " P -2xbc\0" },
        {   7, Qf::a3,   "P c 1 1",     " P -2xc\0" },
        {   8, Qf::b1,   "C 1 m 1",     " C -2y\0" },
        {   8, Qf::b2,   "A 1 m 1",     " A -2y\0" },
        {   8, Qf::b3,   "I 1 m 1",     " I -2y\0" },
        {   8, Qf::c1,   "A 1 1 m",     " A -2\0" },
        {   8, Qf::c2,   "B 1 1 m",     " B -2\0" },
        {   8, Qf::c3,   "I 1 1 m",     " I -2\0" },
        {   8, Qf::a1,   "B m 1 1",     " B -2x\0" },
        {   8, Qf::a2,   "C m 1 1",     " C -2x\0" },
        {   8, Qf::a3,   "I m 1 1",     " I -2x\0" },
        {   9, Qf::b1,   "C 1 c 1",     " C -2yc\0" },
        {   9, Qf::b2,   "A 1 n 1",     " A -2yab\0" },
        {   9, Qf::b3,   "I 1 a 1",     " I -2ya\0" },
        {   9, Qf::mb1,  "A 1 a 1",     " A -2ya\0" },
        {   9, Qf::mb2,  "C 1 n 1",     " C -2yac\0" },
        {   9, Qf::mb3,  "I 1 c 1",     " I -2yc\0" },
        {   9, Qf::c1,   "A 1 1 a",     " A -2a\0" },
        {   9, Qf::c2,   "B 1 1 n",     " B -2ab\0" },
        {   9, Qf::c3,   "I 1 1 b",     " I -2b\0" },
        {   9, Qf::mc1,  "B 1 1 b",     " B -2b\0" },
        {   9, Qf::mc2,  "A 1 1 n",     " A -2ab\0" },
        {   9, Qf::mc3,  "I 1 1 a",     " I -2a\0" },
        {   9, Qf::a1,   "B b 1 1",     " B -2xb\0" },
        {   9, Qf::a2,   "C n 1 1",     " C -2xac\0" },
        {   9, Qf::a3,   "I c 1 1",     " I -2xc\0" },
        {   9, Qf::ma1,  "C c 1 1",     " C -2xc\0" },
        {   9, Qf::ma2,  "B n 1 1",     " B -2xab\0" },
        {   9, Qf::ma3,  "I b 1 1",     " I -2xb\0" },
        {  10, Qf::b,    "P 1 2/m 1",   "-P 2y\0" },
        {  10, Qf::c,    "P 1 1 2/m",   "-P 2\0" },
        {  10, Qf::a,    "P 2/m 1 1",   "-P 2x\0" },
        {  11, Qf::b,    "P 1 21/m 1",  "-P 2yb\0" },
        {  11, Qf::c,    "P 1 1 21/m",  "-P 2c\0" },
        {  11, Qf::a,    "P 21/m 1 1",  "-P 2xa\0" },
        {  12, Qf::b1,   "C 1 2/m 1",   "-C 2y\0" },
        {  12, Qf::b2,   "A 1 2/m 1",   "-A 2y\0" },
        {  12, Qf::b3,   "I 1 2/m 1",   "-I 2y\0" },
        {  12, Qf::c1,   "A 1 1 2/m",   "-A 2\0" },
        {  12, Qf::c2,   "B 1 1 2/m",   "-B 2\0" },
        {  12, Qf::c3,   "I 1 1 2/m",   "-I 2\0" },
        {  12, Qf::a1,   "B 2/m 1 1",   "-B 2x\0" },
        {  12, Qf::a2,   "C 2/m 1 1",   "-C 2x\0" },
        {  12, Qf::a3,   "I 2/m 1 1",   "-I 2x\0" },
        {  13, Qf::b1,   "P 1 2/c 1",   "-P 2yc\0" },
        {  13, Qf::b2,   "P 1 2/n 1",   "-P 2yac\0" },
        {  13, Qf::b3,   "P 1 2/a 1",   "-P 2ya\0" },
        {  13, Qf::c1,   "P 1 1 2/a",   "-P 2a\0" },
        {  13, Qf::c2,   "P 1 1 2/n",   "-P 2ab\0" },
        {  13, Qf::c3,   "P 1 1 2/b",   "-P 2b\0" },
        {  13, Qf::a1,   "P 2/b 1 1",   "-P 2xb\0" },
        {  13, Qf::a2,   "P 2/n 1 1",   "-P 2xbc\0" },
        {  13, Qf::a3,   "P 2/c 1 1",   "-P 2xc\0" },
        {  14, Qf::b1,   "P 1 21/c 1",  "-P 2ybc\0" },
        {  14, Qf::b2,   "P 1 21/n 1",  "-P 2yn\0" },
        {  14, Qf::b3,   "P 1 21/a 1",  "-P 2yab\0" },
        {  14, Qf::c1,   "P 1 1 21/a",  "-P 2ac\0" },
        {  14, Qf::c2,   "P 1 1 21/n",  "-P 2n\0" },
        {  14, Qf::c3,   "P 1 1 21/b",  "-P 2bc\0" },
        {  14, Qf::a1,   "P 21/b 1 1",  "-P 2xab\0" },
        {  14, Qf::a2,   "P 21/n 1 1",  "-P 2xn\0" },
        {  14, Qf::a3,   "P 21/c 1 1",  "-P 2xac\0" },
        {  15, Qf::b1,   "C 1 2/c 1",   "-C 2yc\0" },
        {  15, Qf::b2,   "A 1 2/n 1",   "-A 2yab\0" },
        {  15, Qf::b3,   "I 1 2/a 1",   "-I 2ya\0" },
        {  15, Qf::mb1,  "A 1 2/a 1",   "-A 2ya\0" },
        {  15, Qf::mb2,  "C 1 2/n 1",   "-C 2yac\0" },
        {  15, Qf::mb3,  "I 1 2/c 1",   "-I 2yc\0" },
        {  15, Qf::c1,   "A 1 1 2/a",   "-A 2a\0" },
        {  15, Qf::c2,   "B 1 1 2/n",   "-B 2ab\0" },
        {  15, Qf::c3,   "I 1 1 2/b",   "-I 2b\0" },
        {  15, Qf::mc1,  "B 1 1 2/b",   "-B 2b\0" },
        {  15, Qf::mc2,  "A 1 1 2/n",   "-A 2ab\0" },
        {  15, Qf::mc3,  "I 1 1 2/a",   "-I 2a\0" },
        {  15, Qf::a1,   "B 2/b 1 1",   "-B 2xb\0" },
        {  15, Qf::a2,   "C 2/n 1 1",   "-C 2xac\0" },
        {  15, Qf::a3,   "I 2/c 1 1",   "-I 2xc\0" },
        {  15, Qf::ma1,  "C 2/c 1 1",   "-C 2xc\0" },
        {  15, Qf::ma2,  "B 2/n 1 1",   "-B 2xab\0" },
        {  15, Qf::ma3,  "I 2/b 1 1",   "-I 2xb\0" },
        {  16, 0,        "P 2 2 2",     " P 2 2\0" },
        {  17, 0,        "P 2 2 21",    " P 2c 2\0" },
        {  17, Qf::cab,  "P 21 2 2",    " P 2a 2a\0" },
        {  17, Qf::bca,  "P 2 21 2",    " P 2 2b\0" },
        {  18, 0,        "P 21 21 2",   " P 2 2ab\0" },
        {  18, Qf::cab,  "P 2 21 21",   " P 2bc 2\0" },
        {  18, Qf::bca,  "P 21 2 21",   " P 2ac 2ac\0" },
        {  19, 0,        "P 21 21 21",  " P 2ac 2ab\0" },
        {  20, 0,        "C 2 2 21",    " C 2c 2\0" },
        {  20, Qf::cab,  "A 21 2 2",    " A 2a 2a\0" },
        {  20, Qf::bca,  "B 2 21 2",    " B 2 2b\0" },
        {  21, 0,        "C 2 2 2",     " C 2 2\0" },
        {  21, Qf::cab,  "A 2 2 2",     " A 2 2\0" },
        {  21, Qf::bca,  "B 2 2 2",     " B 2 2\0" },
        {  22, 0,        "F 2 2 2",     " F 2 2\0" },
        {  23, 0,        "I 2 2 2",     " I 2 2\0" },
        {  24, 0,        "I 21 21 21",  " I 2b 2c\0" },
        {  25, 0,        "P m m 2",     " P 2 -2\0" },
        {  25, Qf::cab,  "P 2 m m",     " P -2 2\0" },
        {  25, Qf::bca,  "P m 2 m",     " P -2 -2\0" },
        {  26, 0,        "P m c 21",    " P 2c -2\0" },
        {  26, Qf::bamc, "P c m 21",    " P 2c -2c\0" },
        {  26, Qf::cab,  "P 21 m a",    " P -2a 2a\0" },
        {  26, Qf::mcba, "P 21 a m",    " P -2 2a\0" },
        {  26, Qf::bca,  "P b 21 m",    " P -2 -2b\0" },
        {  26, Qf::amcb, "P m 21 b",    " P -2b -2\0" },
        {  27, 0,        "P c c 2",     " P 2 -2c\0" },
        {  27, Qf::cab,  "P 2 a a",     " P -2a 2\0" },
        {  27, Qf::bca,  "P b 2 b",     " P -2b -2b\0" },
        {  28, 0,        "P m a 2",     " P 2 -2a\0" },
        {  28, Qf::bamc, "P b m 2",     " P 2 -2b\0" },
        {  28, Qf::cab,  "P 2 m b",     " P -2b 2\0" },
        {  28, Qf::mcba, "P 2 c m",     " P -2c 2\0" },
        {  28, Qf::bca,  "P c 2 m",     " P -2c -2c\0" },
        {  28, Qf::amcb, "P m 2 a",     " P -2a -2a\0" },
        {  29, 0,        "P c a 21",    " P 2c -2ac\0" },
        {  29, Qf::bamc, "P b c 21",    " P 2c -2b\0" },
        {  29, Qf::cab,  "P 21 a b",    " P -2b 2a\0" },
        {  29, Qf::mcba, "P 21 c a",    " P -2ac 2a\0" },
        {  29, Qf::bca,  "P c 21 b",    " P -2bc -2c\0" },
        {  29, Qf::amcb, "P b 21 a",    " P -2a -2ab\0" },
        {  30, 0,        "P n c 2",     " P 2 -2bc\0" },
        {  30, Qf::bamc, "P c n 2",     " P 2 -2ac\0" },
        {  30, Qf::cab,  "P 2 n a",     " P -2ac 2\0" },
        {  30, Qf::mcba, "P 2 a n",     " P -2ab 2\0" },
        {  30, Qf::bca,  "P b 2 n",     " P -2ab -2ab\0" },
        {  30, Qf::amcb, "P n 2 b",     " P -2bc -2bc\0" },
        {  31, 0,        "P m n 21",    " P 2ac -2\0" },
        {  31, Qf::bamc, "P n m 21",    " P 2bc -2bc\0" },
        {  31, Qf::cab,  "P 21 m n",    " P -2ab 2ab\0" },
        {  31, Qf::mcba, "P 21 n m",    " P -2 2ac\0" },
        {  31, Qf::bca,  "P n 21 m",    " P -2 -2bc\0" },
        {  31, Qf::amcb, "P m 21 n",    " P -2ab -2\0" },
        {  32, 0,        "P b a 2",     " P 2 -2ab\0" },
        {  32, Qf::cab,  "P 2 c b",     " P -2bc 2\0" },
        {  32, Qf::bca,  "P c 2 a",     " P -2ac -2ac\0" },
        {  33, 0,        "P n a 21",    " P 2c -2n\0" },
        {  33, Qf::bamc, "P b n 21",    " P 2c -2ab\0" },
        {  33, Qf::cab,  "P 21 n b",    " P -2bc 2a\0" },
        {  33, Qf::mcba, "P 21 c n",    " P -2n 2a\0" },
        {  33, Qf::bca,  "P c 21 n",    " P -2n -2ac\0" },
        {  33, Qf::amcb, "P n 21 a",    " P -2ac -2n\0" },
        {  34, 0,        "P n n 2",     " P 2 -2n\0" },
        {  34, Qf::cab,  "P 2 n n",     " P -2n 2\0" },
        {  34, Qf::bca,  "P n 2 n",     " P -2n -2n\0" },
        {  35, 0,        "C m m 2",     " C 2 -2\0" },
        {  35, Qf::cab,  "A 2 m m",     " A -2 2\0" },
        {  35, Qf::bca,  "B m 2 m",     " B -2 -2\0" },
        {  36, 0,        "C m c 21",    " C 2c -2\0" },
        {  36, Qf::bamc, "C c m 21",    " C 2c -2c\0" },
        {  36, Qf::cab,  "A 21 m a",    " A -2a 2a\0" },
        {  36, Qf::mcba, "A 21 a m",    " A -2 2a\0" },
        {  36, Qf::bca,  "B b 21 m",    " B -2 -2b\0" },
        {  36, Qf::amcb, "B m 21 b",    " B -2b -2\0" },
        {  37, 0,        "C c c 2",     " C 2 -2c\0" },
        {  37, Qf::cab,  "A 2 a a",     " A -2a 2\0" },
        {  37, Qf::bca,  "B b 2 b",     " B -2b -2b\0" },
        {  38, 0,        "A m m 2",     " A 2 -2\0" },
        {  38, Qf::bamc, "B m m 2",     " B 2 -2\0" },
        {  38, Qf::cab,  "B 2 m m",     " B -2 2\0" },
        {  38, Qf::mcba, "C 2 m m",     " C -2 2\0" },
        {  38, Qf::bca,  "C m 2 m",     " C -2 -2\0" },
        {  38, Qf::amcb, "A m 2 m",     " A -2 -2\0" },
        {  39, 0,        "A b m 2",     " A 2 -2b\0" },
        {  39, Qf::bamc, "B m a 2",     " B 2 -2a\0" },
        {  39, Qf::cab,  "B 2 c m",     " B -2a 2\0" },
        {  39, Qf::mcba, "C 2 m b",     " C -2a 2\0" },
        {  39, Qf::bca,  "C m 2 a",     " C -2a -2a\0" },
        {  39, Qf::amcb, "A c 2 m",     " A -2b -2b\0" },
        {  40, 0,        "A m a 2",     " A 2 -2a\0" },
        {  40, Qf::bamc, "B b m 2",     " B 2 -2b\0" },
        {  40, Qf::cab,  "B 2 m b",     " B -2b 2\0" },
        {  40, Qf::mcba, "C 2 c m",     " C -2c 2\0" },
        {  40, Qf::bca,  "C c 2 m",     " C -2c -2c\0" },
        {  40, Qf::amcb, "A m 2 a",     " A -2a -2a\0" },
        {  41, 0,        "A b a 2",     " A 2 -2ab\0" },
        {  41, Qf::bamc, "B b a 2",     " B 2 -2ab\0" },
        {  41, Qf::cab,  "B 2 c b",     " B -2ab 2\0" },
        {  41, Qf::mcba, "C 2 c b",     " C -2ac 2\0" },
        {  41, Qf::bca,  "C c 2 a",     " C -2ac -2ac\0" },
        {  41, Qf::amcb, "A c 2 a",     " A -2ab -2ab\0" },
        {  42, 0,        "F m m 2",     " F 2 -2\0" },
        {  42, Qf::cab,  "F 2 m m",     " F -2 2\0" },
        {  42, Qf::bca,  "F m 2 m",     " F -2 -2\0" },
        {  43, 0,        "F d d 2",     " F 2 -2d\0" },
        {  43, Qf::cab,  "F 2 d d",     " F -2d 2\0" },
        {  43, Qf::bca,  "F d 2 d",     " F -2d -2d\0" },
        {  44, 0,        "I m m 2",     " I 2 -2\0" },
        {  44, Qf::cab,  "I 2 m m",     " I -2 2\0" },
        {  44, Qf::bca,  "I m 2 m",     " I -2 -2\0" },
        {  45, 0,        "I b a 2",     " I 2 -2c\0" },
        {  45, Qf::cab,  "I 2 c b",     " I -2a 2\0" },
        {  45, Qf::bca,  "I c 2 a",     " I -2b -2b\0" },
        {  46, 0,        "I m a 2",     " I 2 -2a\0" },
        {  46, Qf::bamc, "I b m 2",     " I 2 -2b\0" },
        {  46, Qf::cab,  "I 2 m b",     " I -2b 2\0" },
        {  46, Qf::mcba, "I 2 c m",     " I -2c 2\0" },
        {  46, Qf::bca,  "I c 2 m",     " I -2c -2c\0" },
        {  46, Qf::amcb, "I m 2 a",     " I -2a -2a\0" },
        {  47, 0,        "P m m m",     "-P 2 2\0" },
        {  48, 0,        "P n n n",     " P 2 2 -1n\0-P 2ab 2bc\0" },
        {  49, 0,        "P c c m",     "-P 2 2c\0" },
        {  49, Qf::cab,  "P m a a",     "-P 2a 2\0" },
        {  49, Qf::bca,  "P b m b",     "-P 2b 2b\0" },
        {  50, 0,        "P b a n",     " P 2 2 -1ab\0-P 2ab 2b\0" },
        {  50, Qf::cab,  "P n c b",     " P 2 2 -1bc\0-P 2b 2bc\0" },
        {  50, Qf::bca,  "P c n a",     " P 2 2 -1ac\0-P 2a 2c\0" },
        {  51, 0,        "P m m a",     "-P 2a 2a\0" },
        {  51, Qf::bamc, "P m m b",     "-P 2b 2\0" },
        {  51, Qf::cab,  "P b m m",     "-P 2 2b\0" },
        {  51, Qf::mcba, "P c m m",     "-P 2c 2c\0" },
        {  51, Qf::bca,  "P m c m",     "-P 2c 2\0" },
        {  51, Qf::amcb, "P m a m",     "-P 2 2a\0" },
        {  52, 0,        "P n n a",     "-P 2a 2bc\0" },
        {  52, Qf::bamc, "P n n b",     "-P 2b 2n\0" },
        {  52, Qf::cab,  "P b n n",     "-P 2n 2b\0" },
        {  52, Qf::mcba, "P c n n",     "-P 2ab 2c\0" },
        {  52, Qf::bca,  "P n c n",     "-P 2ab 2n\0" },
        {  52, Qf::amcb, "P n a n",     "-P 2n 2bc\0" },
        {  53, 0,        "P m n a",     "-P 2ac 2\0" },
        {  53, Qf::bamc, "P n m b",     "-P 2bc 2bc\0" },
        {  53, Qf::cab,  "P b m n",     "-P 2ab 2ab\0" },
        {  53, Qf::mcba, "P c n m",     "-P 2 2ac\0" },
        {  53, Qf::bca,  "P n c m",     "-P 2 2bc\0" },
        {  53, Qf::amcb, "P m a n",     "-P 2ab 2\0" },
        {  54, 0,        "P c c a",     "-P 2a 2ac\0" },
        {  54, Qf::bamc, "P c c b",     "-P 2b 2c\0" },
        {  54, Qf::cab,  "P b a a",     "-P 2a 2b\0" },
        {  54, Qf::mcba, "P c a a",     "-P 2ac 2c\0" },
        {  54, Qf::bca,  "P b c b",     "-P 2bc 2b\0" },
        {  54, Qf::amcb, "P b a b",     "-P 2b 2ab\0" },
        {  55, 0,        "P b a m",     "-P 2 2ab\0" },
        {  55, Qf::cab,  "P m c b",     "-P 2bc 2\0" },
        {  55, Qf::bca,  "P c m a",     "-P 2ac 2ac\0" },
        {  56, 0,        "P c c n",     "-P 2ab 2ac\0" },
        {  56, Qf::cab,  "P n a a",     "-P 2ac 2bc\0" },
        {  56, Qf::bca,  "P b n b",     "-P 2bc 2ab\0" },
        {  57, 0,        "P b c m",     "-P 2c 2b\0" },
        {  57, Qf::bamc, "P c a m",     "-P 2c 2ac\0" },
        {  57, Qf::cab,  "P m c a",     "-P 2ac 2a\0" },
        {  57, Qf::mcba, "P m a b",     "-P 2b 2a\0" },
        {  57, Qf::bca,  "P b m a",     "-P 2a 2ab\0" },
        {  57, Qf::amcb, "P c m b",     "-P 2bc 2c\0" },
        {  58, 0,        "P n n m",     "-P 2 2n\0" },
        {  58, Qf::cab,  "P m n n",     "-P 2n 2\0" },
        {  58, Qf::bca,  "P n m n",     "-P 2n 2n\0" },
        {  59, 0,        "P m m n",     " P 2 2ab -1ab\0-P 2ab 2a\0" },
        {  59, Qf::cab,  "P n m m",     " P 2bc 2 -1bc\0-P 2c 2bc\0" },
        {  59, Qf::bca,  "P m n m",     " P 2ac 2ac -1ac\0-P 2c 2a\0" },
        {  60, 0,        "P b c n",     "-P 2n 2ab\0" },
        {  60, Qf::bamc, "P c a n",     "-P 2n 2c\0" },
        {  60, Qf::cab,  "P n c a",     "-P 2a 2n\0" },
        {  60, Qf::mcba, "P n a b",     "-P 2bc 2n\0" },
        {  60, Qf::bca,  "P b n a",     "-P 2ac 2b\0" },
        {  60, Qf::amcb, "P c n b",     "-P 2b 2ac\0" },
        {  61, 0,        "P b c a",     "-P 2ac 2ab\0" },
        {  61, Qf::bamc, "P c a b",     "-P 2bc 2ac\0" },
        {  62, 0,        "P n m a",     "-P 2ac 2n\0" },
        {  62, Qf::bamc, "P m n b",     "-P 2bc 2a\0" },
        {  62, Qf::cab,  "P b n m",     "-P 2c 2ab\0" },
        {  62, Qf::mcba, "P c m n",     "-P 2n 2ac\0" },
        {  62, Qf::bca,  "P m c n",     "-P 2n 2a\0" },
        {  62, Qf::amcb, "P n a m",     "-P 2c 2n\0" },
        {  63, 0,        "C m c m",     "-C 2c 2\0" },
        {  63, Qf::bamc, "C c m m",     "-C 2c 2c\0" },
        {  63, Qf::cab,  "A m m a",     "-A 2a 2a\0" },
        {  63, Qf::mcba, "A m a m",     "-A 2 2a\0" },
        {  63, Qf::bca,  "B b m m",     "-B 2 2b\0" },
        {  63, Qf::amcb, "B m m b",     "-B 2b 2\0" },
        {  64, 0,        "C m c a",     "-C 2ac 2\0" },
        {  64, Qf::bamc, "C c m b",     "-C 2ac 2ac\0" },
        {  64, Qf::cab,  "A b m a",     "-A 2ab 2ab\0" },
        {  64, Qf::mcba, "A c a m",     "-A 2 2ab\0" },
        {  64, Qf::bca,  "B b c m",     "-B 2 2ab\0" },
        {  64, Qf::amcb, "B m a b",     "-B 2ab 2\0" },
        {  65, 0,        "C m m m",     "-C 2 2\0" },
        {  65, Qf::cab,  "A m m m",     "-A 2 2\0" },
        {  65, Qf::bca,  "B m m m",     "-B 2 2\0" },
        {  66, 0,        "C c c m",     "-C 2 2c\0" },
        {  66, Qf::cab,  "A m a a",     "-A 2a 2\0" },
        {  66, Qf::bca,  "B b m b",     "-B 2b 2b\0" },
        {  67, 0,        "C m m a",     "-C 2a 2\0" },
        {  67, Qf::bamc, "C m m b",     "-C 2a 2a\0" },
        {  67, Qf::cab,  "A b m m",     "-A 2b 2b\0" },
        {  67, Qf::mcba, "A c m m",     "-A 2 2b\0" },
        {  67, Qf::bca,  "B m c m",     "-B 2 2a\0" },
        {  67, Qf::amcb, "B m a m",     "-B 2a 2\0" },
        {  68, 0,        "C c c a",     " C 2 2 -1ac\0-C 2a 2ac\0" },
        {  68, Qf::bamc, "C c c b",     " C 2 2 -1ac\0-C 2a 2c\0" },
        {  68, Qf::cab,  "A b a a",     " A 2 2 -1ab\0-A 2a 2b\0" },
        {  68, Qf::mcba, "A c a a",     " A 2 2 -1ab\0-A 2ab 2b\0" },
        {  68, Qf::bca,  "B b c b",     " B 2 2 -1ab\0-B 2ab 2b\0" },
        {  68, Qf::amcb, "B b a b",     " B 2 2 -1ab\0-B 2b 2ab\0" },
        {  69, 0,        "F m m m",     "-F 2 2\0" },
        {  70, 0,        "F d d d",     " F 2 2 -1d\0-F 2uv 2vw\0" },
        {  71, 0,        "I m m m",     "-I 2 2\0" },
        {  72, 0,        "I b a m",     "-I 2 2c\0" },
        {  72, Qf::cab,  "I m c b",     "-I 2a 2\0" },
        {  72, Qf::bca,  "I c m a",     "-I 2b 2b\0" },
        {  73, 0,        "I b c a",     "-I 2b 2c\0" },
        {  73, Qf::bamc, "I c a b",     "-I 2a 2b\0" },
        {  74, 0,        "I m m a",     "-I 2b 2\0" },
        {  74, Qf::bamc, "I m m b",     "-I 2a 2a\0" },
        {  74, Qf::cab,  "I b m m",     "-I 2c 2c\0" },
        {  74, Qf::mcba, "I c m m",     "-I 2 2b\0" },
        {  74, Qf::bca,  "I m c m",     "-I 2 2a\0" },
        {  74, Qf::amcb, "I m a m",     "-I 2c 2\0" },
        {  75, 0,        "P 4",         " P 4\0" },
        {  75, Qf::MC,   "C 4",         " C 4\0" },
        {  76, 0,        "P 41",        " P 4w\0" },
        {  76, Qf::MC,   "C 41",        " C 4w\0" },
        {  77, 0,        "P 42",        " P 4c\0" },
        {  77, Qf::MC,   "C 42",        " C 4c\0" },
        {  78, 0,        "P 43",        " P 4cw\0" },
        {  78, Qf::MC,   "C 43",        " C 4cw\0" },
        {  79, 0,        "I 4",         " I 4\0" },
        {  79, Qf::MF,   "F 4",         " F 4\0" },
        {  80, 0,        "I 41",        " I 4bw\0" },
        {  80, Qf::MF,   "F 41",        " F 4d\0" },
        {  81, 0,        "P -4",        " P -4\0" },
        {  81, Qf::MC,   "C -4",        " C -4\0" },
        {  82, 0,        "I -4",        " I -4\0" },
        {  82, Qf::MF,   "F -4",        " F -4\0" },
        {  83, 0,        "P 4/m",       "-P 4\0" },
        {  83, Qf::MC,   "C 4/m",       "-C 4\0" },
        {  84, 0,        "P 42/m",      "-P 4c\0" },
        {  84, Qf::MC,   "C 42/m",      "-C 4c\0" },
        {  85, 0,        "P 4/n",       " P 4ab -1ab\0-P 4a\0" },
        {  85, Qf::MC,   "C 4/a",       " C 4a -1a\0-C 4auv\0" },
        {  86, 0,        "P 42/n",      " P 4n -1n\0-P 4bc\0" },
        {  86, Qf::MC,   "C 42/a",      " C 4ac -1ac\0-C 4wd\0" },
        {  87, 0,        "I 4/m",       "-I 4\0" },
        {  87, Qf::MF,   "F 4/m",       "-F 4\0" },
        {  88, 0,        "I 41/a",      " I 4bw -1bw\0-I 4ad\0" },
        {  88, Qf::MF,   "F 41/d",      " F 4d -1d\0-F 4vw\0" },
        {  89, 0,        "P 4 2 2",     " P 4 2\0" },
        {  89, Qf::MC,   "C 4 2 2",     " C 4 2\0" },
        {  90, 0,        "P 4 21 2",    " P 4ab 2ab\0" },
        {  90, Qf::MC,   "C 4 2 21",    " C 4a 2\0" },
        {  91, 0,        "P 41 2 2",    " P 4w 2c\0" },
        {  91, Qf::MC,   "C 41 2 2",    " C 4w 2cw\0" },
        {  92, 0,        "P 41 21 2",   " P 4abw 2nw\0" },
        {  92, Qf::MC,   "C 41 2 21",   " C 4aw 2\0" },
        {  93, 0,        "P 42 2 2",    " P 4c 2\0" },
        {  93, Qf::MC,   "C 42 2 2",    " C 4c 2c\0" },
        {  94, 0,        "P 42 21 2",   " P 4n 2n\0" },
        {  94, Qf::MC,   "C 42 2 21",   " C 4ac 2\0" },
        {  95, 0,        "P 43 2 2",    " P 4cw 2c\0" },
        {  95, Qf::MC,   "C 43 2 2",    " C 4cw 2w\0" },
        {  96, 0,        "P 43 21 2",   " P 4nw 2abw\0" },
        {  96, Qf::MC,   "C 43 2 21",   " C 4acw 2\0" },
        {  97, 0,        "I 4 2 2",     " I 4 2\0" },
        {  97, Qf::MF,   "F 4 2 2",     " F 4 2\0" },
        {  98, 0,        "I 41 2 2",    " I 4bw 2bw\0" },
        {  98, Qf::MF,   "F 41 2 2",    " F 4d 2\0" },
        {  99, 0,        "P 4 m m",     " P 4 -2\0" },
        {  99, Qf::MC,   "C 4 m m",     " C 4 -2\0" },
        { 100, 0,        "P 4 b m",     " P 4 -2ab\0" },
        { 100, Qf::MC,   "C 4 m g1",    " C 4 -2a\0" },
        { 101, 0,        "P 42 c m",    " P 4c -2c\0" },
        { 101, Qf::MC,   "C 42 m c",    " C 4c -2\0" },
        { 102, 0,        "P 42 n m",    " P 4n -2n\0" },
        { 102, Qf::MC,   "C 42 m g2",   " C 4ac -2\0" },
        { 103, 0,        "P 4 c c",     " P 4 -2c\0" },
        { 103, Qf::MC,   "C 4 c c",     " C 4 -2c\0" },
        { 104, 0,        "P 4 n c",     " P 4 -2n\0" },
        { 104, Qf::MC,   "C 4 c g2",    " C 4 -2ac\0" },
        { 105, 0,        "P 42 m c",    " P 4c -2\0" },
        { 105, Qf::MC,   "C 42 c m",    " C 4c -2c\0" },
        { 106, 0,        "P 42 b c",    " P 4c -2ab\0" },
        { 106, Qf::MC,   "C 42 c g1",   " C 4c -2ac\0" },
        { 107, 0,        "I 4 m m",     " I 4 -2\0" },
        { 107, Qf::MF,   "F 4 m m",     " F 4 -2\0" },
        { 108, 0,        "I 4 c m",     " I 4 -2c\0" },
        { 108, Qf::MF,   "F 4 m c",     " F 4 -2a\0" },
        { 109, 0,        "I 41 m d",    " I 4bw -2\0" },
        { 109, Qf::MF,   "F 41 d m",    " F 4d -2d\0" },
        { 110, 0,        "I 41 c d",    " I 4bw -2c\0" },
        { 110, Qf::MF,   "F 41 d c",    " F 4d -2ad\0" },
        { 111, 0,        "P -4 2 m",    " P -4 2\0" },
        { 111, Qf::MC,   "C -4 m 2",    " C -4 -2\0" },
        { 112, 0,        "P -4 2 c",    " P -4 2c\0" },
        { 112, Qf::MC,   "C -4 c 2",    " C -4 -2c\0" },
        { 113, 0,        "P -4 21 m",   " P -4 2ab\0" },
        { 113, Qf::MC,   "C -4 m 21",   " C -4 -2a\0" },
        { 114, 0,        "P -4 21 c",   " P -4 2n\0" },
        { 114, Qf::MC,   "C -4 c 21",   " C -4 -2ac\0" },
        { 115, 0,        "P -4 m 2",    " P -4 -2\0" },
        { 115, Qf::MC,   "C -4 2 m",    " C -4 2\0" },
        { 116, 0,        "P -4 c 2",    " P -4 -2c\0" },
        { 116, Qf::MC,   "C -4 2 c",    " C -4 2c\0" },
        { 117, 0,        "P -4 b 2",    " P -4 -2ab\0" },
        { 117, Qf::MC,   "C -4 2 g1",   " C -4 2a\0" },
        { 118, 0,        "P -4 n 2",    " P -4 -2n\0" },
        { 118, Qf::MC,   "C -4 2 g2",   " C -4 2ac\0" },
        { 119, 0,        "I -4 m 2",    " I -4 -2\0" },
        { 119, Qf::MF,   "F -4 2 m",    " F -4 2\0" },
        { 120, 0,        "I -4 c 2",    " I -4 -2c\0" },
        { 120, Qf::MF,   "F -4 2 c",    " F -4 2a\0" },
        { 121, 0,        "I -4 2 m",    " I -4 2\0" },
        { 121, Qf::MF,   "F -4 m 2",    " F -4 -2\0" },
        { 122, 0,        "I -4 2 d",    " I -4 2bw\0" },
        { 122, Qf::MF,   "F -4 d 2",    " F -4 -2d\0" },
        { 123, 0,        "P 4/m m m",   "-P 4 2\0" },
        { 123, Qf::MC,   "C 4/m m m",   "-C 4 2\0" },
        { 124, 0,        "P 4/m c c",   "-P 4 2c\0" },
        { 124, Qf::MC,   "C 4/m c c",   "-C 4 2c\0" },
        { 125, 0,        "P 4/n b m",   " P 4 2 -1ab\0-P 4a 2b\0" },
        { 125, Qf::MC,   "C 4/a m g1",  " C 4 2 -1a\0-C 4auv 2\0" },
        { 126, 0,        "P 4/n n c",   " P 4 2 -1n\0-P 4a 2bc\0" },
        { 126, Qf::MC,   "C 4/a c g2",  " C 4 2 -1ac\0-C 4auv 2c\0" },
        { 127, 0,        "P 4/m b m",   "-P 4 2ab\0" },
        { 127, Qf::MC,   "C 4/m m g1",  "-C 4 2a\0" },
        { 128, 0,        "P 4/m n c",   "-P 4 2n\0" },
        { 128, Qf::MC,   "C 4/m c g2",  "-C 4 2ac\0" },
        { 129, 0,        "P 4/n m m",   " P 4ab 2ab -1ab\0-P 4a 2a\0" },
        { 129, Qf::MC,   "C 4/a m m",   " C 4a 2 -1a\0-C 4auv 2a\0" },
        { 130, 0,        "P 4/n c c",   " P 4ab 2n -1ab\0-P 4a 2ac\0" },
        { 130, Qf::MC,   "C 4/a c c",   " C 4a 2c -1a\0-C 4auv 2ac\0" },
        { 131, 0,        "P 42/m m c",  "-P 4c 2\0" },
        { 131, Qf::MC,   "C 42/m c m",  "-C 4c 2c\0" },
        { 132, 0,        "P 42/m c m",  "-P 4c 2c\0" },
        { 132, Qf::MC,   "C 42/m m c",  "-C 4c 2\0" },
        { 133, 0,        "P 42/n b c",  " P 4n 2c -1n\0-P 4ac 2b\0" },
        { 133, Qf::MC,   "C 42/a c g1", " C 4ac 2a -1ac\0-C 4acuv 2c\0" },
        { 134, 0,        "P 42/n n m",  " P 4n 2 -1n\0-P 4ac 2bc\0" },
        { 134, Qf::MC,   "C 42/a m g2", " C 4ac 2ac -1ac\0-C 4acuv 2\0" },
        { 135, 0,        "P 42/m b c",  "-P 4c 2ab\0" },
        { 135, Qf::MC,   "C 42/m c g1", "-C 4c 2ac\0" },
        { 136, 0,        "P 42/m n m",  "-P 4n 2n\0" },
        { 136, Qf::MC,   "C 42/m m g2", "-C 4ac 2\0" },
        { 137, 0,        "P 42/n m c",  " P 4n 2n -1n\0-P 4ac 2a\0" },
        { 137, Qf::MC,   "C 42/a c m",  " C 4ac 2 -1ac\0-C 4acuv 2ac\0" },
        { 138, 0,        "P 42/n c m",  " P 4n 2ab -1n\0-P 4ac 2ac\0" },
        { 138, Qf::MC,   "C 42/a m c",  " C 4ac 2c -1ac\0-C 4acuv 2a\0" },
        { 139, 0,        "I 4/m m m",   "-I 4 2\0" },
        { 139, Qf::MF,   "F 4/m m m",   "-F 4 2\0" },
        { 140, 0,        "I 4/m c m",   "-I 4 2c\0" },
        { 140, Qf::MF,   "F 4/m m c",   "-F 4 2a\0" },
        { 141, 0,        "I 41/a m d",  " I 4bw 2bw -1bw\0-I 4bd 2\0" },
        { 141, Qf::MF,   "F 41/d d m",  " F 4d 2 -1d\0-F 4ud 2ud\0" },
        { 142, 0,        "I 41/a c d",  " I 4bw 2aw -1bw\0-I 4bd 2c\0" },
        { 142, Qf::MF,   "F 41/d d c",  " F 4d 2a -1d\0-F 4ud 2vw\0" },
        { 143, 0,        "P 3",         " P 3\0" },
        { 143, Qf::TH,   "H 3",         " H 3\0" },
        { 144, 0,        "P 31",        " P 31\0" },
        { 144, Qf::TH,   "H 31",        " H 31\0" },
        { 145, 0,        "P 32",        " P 32\0" },
        { 145, Qf::TH,   "H 32",        " H 32\0" },
        { 146, 0,        "R 3",         " R 3\0 P 3*\0" },
        { 147, 0,        "P -3",        "-P 3\0" },
        { 147, Qf::TH,   "H -3",        "-H 3\0" },
        { 148, 0,        "R -3",        "-R 3\0-P 3*\0" },
        { 149, 0,        "P 3 1 2",     " P 3 2\0" },
        { 149, Qf::TH,   "H 3 2 1",     " H 3 2\"\0" },
        { 150, 0,        "P 3 2 1",     " P 3 2\"\0" },
        { 150, Qf::TH,   "H 3 1 2",     " H 3 2\0" },
        { 151, 0,        "P 31 1 2",    " P 31 2 (0 0 4)\0" },
        { 151, Qf::TH,   "H 31 2 1",    " H 31 2\"\0" },
        { 152, 0,        "P 31 2 1",    " P 31 2\"\0" },
        { 152, Qf::TH,   "H 31 1 2",    " H 31 2 (0 0 2)\0" },
        { 153, 0,        "P 32 1 2",    " P 32 2 (0 0 2)\0" },
        { 153, Qf::TH,   "H 32 2 1",    " H 32 2\"\0" },
        { 154, 0,        "P 32 2 1",    " P 32 2\"\0" },
        { 154, Qf::TH,   "H 32 1 2",    " H 32 2 (0 0 4)\0" },
        { 155, 0,        "R 3 2",       " R 3 2\"\0 P 3* 2\0" },
        { 156, 0,        "P 3 m 1",     " P 3 -2\"\0" },
        { 156, Qf::TH,   "H 3 1 m",     " H 3 -2\0" },
        { 157, 0,        "P 3 1 m",     " P 3 -2\0" },
        { 157, Qf::TH,   "H 3 m 1",     " H 3 -2\"\0" },
        { 158, 0,        "P 3 c 1",     " P 3 -2\"c\0" },
        { 158, Qf::TH,   "H 3 1 c",     " H 3 -2c\0" },
        { 159, 0,        "P 3 1 c",     " P 3 -2c\0" },
        { 159, Qf::TH,   "H 3 c 1",     " H 3 -2\"c\0" },
        { 160, 0,        "R 3 m",       " R 3 -2\"\0 P 3* -2\0" },
        { 161, 0,        "R 3 c",       " R 3 -2\"c\0 P 3* -2n\0" },
        { 162, 0,        "P -3 1 m",    "-P 3 2\0" },
        { 162, Qf::TH,   "H -3 m 1",    "-H 3 2\"\0" },
        { 163, 0,        "P -3 1 c",    "-P 3 2c\0" },
        { 163, Qf::TH,   "H -3 c 1",    "-H 3 2\"c\0" },
        { 164, 0,        "P -3 m 1",    "-P 3 2\"\0" },
        { 164, Qf::TH,   "H -3 1 m",    "-H 3 2\0" },
        { 165, 0,        "P -3 c 1",    "-P 3 2\"c\0" },
        { 165, Qf::TH,   "H -3 1 c",    "-H 3 2c\0" },
        { 166, 0,        "R -3 m",      "-R 3 2\"\0-P 3* 2\0" },
        { 167, 0,        "R -3 c",      "-R 3 2\"c\0-P 3* 2n\0" },
        { 168, 0,        "P 6",         " P 6\0" },
        { 168, Qf::TH,   "H 6",         " H 6\0" },
        { 169, 0,        "P 61",        " P 61\0" },
        { 169, Qf::TH,   "H 61",        " H 61\0" },
        { 170, 0,        "P 65",        " P 65\0" },
        { 170, Qf::TH,   "H 65",        " H 65\0" },
        { 171, 0,        "P 62",        " P 62\0" },
        { 171, Qf::TH,   "H 62",        " H 62\0" },
        { 172, 0,        "P 64",        " P 64\0" },
        { 172, Qf::TH,   "H 64",        " H 64\0" },
        { 173, 0,        "P 63",        " P 6c\0" },
        { 173, Qf::TH,   "H 63",        " H 6c\0" },
        { 174, 0,        "P -6",        " P -6\0" },
        { 174, Qf::TH,   "H -6",        " H -6\0" },
        { 175, 0,        "P 6/m",       "-P 6\0" },
        { 175, Qf::TH,   "H 6/m",       "-H 6\0" },
        { 176, 0,        "P 63/m",      "-P 6c\0" },
        { 176, Qf::TH,   "H 63/m",      "-H 6c\0" },
        { 177, 0,        "P 6 2 2",     " P 6 2\0" },
        { 177, Qf::TH,   "H 6 2 2",     " H 6 2\0" },
        { 178, 0,        "P 61 2 2",    " P 61 2 (0 0 5)\0" },
        { 178, Qf::TH,   "H 61 2 2",    " H 61 2 (0 0 4)\0" },
        { 179, 0,        "P 65 2 2",    " P 65 2 (0 0 1)\0" },
        { 179, Qf::TH,   "H 65 2 2",    " H 65 2 (0 0 2)\0" },
        { 180, 0,        "P 62 2 2",    " P 62 2 (0 0 4)\0" },
        { 180, Qf::TH,   "H 62 2 2",    " H 62 2 (0 0 2)\0" },
        { 181, 0,        "P 64 2 2",    " P 64 2 (0 0 2)\0" },
        { 181, Qf::TH,   "H 64 2 2",    " H 64 2 (0 0 4)\0" },
        { 182, 0,        "P 63 2 2",    " P 6c 2c\0" },
        { 182, Qf::TH,   "H 63 2 2",    " H 6c 2\0" },
        { 183, 0,        "P 6 m m",     " P 6 -2\0" },
        { 183, Qf::TH,   "H 6 m m",     " H 6 -2\0" },
        { 184, 0,        "P 6 c c",     " P 6 -2c\0" },
        { 184, Qf::TH,   "H 6 c c",     " H 6 -2c\0" },
        { 185, 0,        "P 63 c m",    " P 6c -2\0" },
        { 185, Qf::TH,   "H 63 m c",    " H 6c -2c\0" },
        { 186, 0,        "P 63 m c",    " P 6c -2c\0" },
        { 186, Qf::TH,   "H 63 c m",    " H 6c -2\0" },
        { 187, 0,        "P -6 m 2",    " P -6 2\0" },
        { 187, Qf::TH,   "H -6 2 m",    " H -6 -2\0" },
        { 188, 0,        "P -6 c 2",    " P -6c 2\0" },
        { 188, Qf::TH,   "H -6 2 c",    " H -6c -2c\0" },
        { 189, 0,        "P -6 2 m",    " P -6 -2\0" },
        { 189, Qf::TH,   "H -6 m 2",    " H -6 2\0" },
        { 190, 0,        "P -6 2 c",    " P -6c -2c\0" },
        { 190, Qf::TH,   "H -6 c 2",    " H -6c 2\0" },
        { 191, 0,        "P 6/m m m",   "-P 6 2\0" },
        { 191, Qf::TH,   "H 6/m m m",   "-H 6 2\0" },
        { 192, 0,        "P 6/m c c",   "-P 6 2c\0" },
        { 192, Qf::TH,   "H 6/m c c",   "-H 6 2c\0" },
        { 193, 0,        "P 63/m c m",  "-P 6c 2\0" },
        { 193, Qf::TH,   "H 63/m m c",  "-H 6c 2c\0" },
        { 194, 0,        "P 63/m m c",  "-P 6c 2c\0" },
        { 194, Qf::TH,   "H 63/m c m",  "-H 6c 2\0" },
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
    string PreProcessSymbol(const string& RawSymbol)
    {
      string result;
      rangei(RawSymbol.size()) {
        const char r = RawSymbol[i];
        if (!isspace(r) && r != '_') result += tolower(r);
      }
      return result;
    }

    char StripExtension(string& WorkSymbol)
    {
      char Ext = '\0';
      string::size_type stop = WorkSymbol.find(':');
      if (stop != string::npos) {
        string e = WorkSymbol.substr(stop + 1);
        if (e.size() == 1) {
          Ext = e[0];
        }
        else if (e == "o1" || e == "o2") {
          Ext = e[1];
        }
      }
      else {
        // check if last character is S, Z, R or H
        if (WorkSymbol.size() > 0) {
          string::size_type i = WorkSymbol.size() - 1;
          switch (WorkSymbol[i]) {
            case 's':
            case 'z':
            case 'r':
            case 'h':
              Ext = WorkSymbol[i];
              stop = i;
              break;
          }
          // check if last two characters are O1 or O2
          if (Ext == '\0' && WorkSymbol.size() > 1) {
            string::size_type i = WorkSymbol.size() - 2;
            string e = WorkSymbol.substr(i);
            if (e == "o1" || e == "o2") {
              Ext = e[1];
              stop = i;
            }
          }
        }
      }

      if (Ext != '\0') {
          switch (Ext) {
          case '1': break;
          case 's': Ext = '1'; break;
          case '2': break;
          case 'z': Ext = '2'; break;
          case 'r': Ext = 'R'; break;
          case 'h': Ext = 'H'; break;
          default:
            Ext = '\0';
        }
      }

      if (Ext != '\0') WorkSymbol.erase(stop);

      return Ext;
    }

    // remove parentheses, e.g. "P2(1)2(1)2(1)" -> "P212121"
    void RemoveParentheses(string& WorkSymbol)
    {
      static const int  RotNumbers[] = { 2, 3, 4, 6 };
      static const char RotSymbols[] = "2346";
      static const char ScrSymbols[] = "012345";
      string pat = "r(s)";
      int ir;
      for(ir=0;ir<4;ir++) {
        pat[0] = RotSymbols[ir];
        int is;
        for (is = 1; is < RotNumbers[ir]; is++) {
          pat[2] = ScrSymbols[is];
          for (;;) {
            string::size_type i = WorkSymbol.find(pat);
            if (i == string::npos) break;
            WorkSymbol.erase(i + 3, 1);
            WorkSymbol.erase(i + 1, 1);
          }
        }
      }
    }

    string RemoveSpaces(const string& inp)
    {
      string result;
      rangei(inp.size()) if (inp[i] != ' ') result += inp[i];
      return result;
    }

    int CmpSchoenfliesSymbols(const string& FromTable,
                              const string& WorkSymbol)
    {
      if (FromTable.size() != WorkSymbol.size()) return -1;
      rangei(FromTable.size()) {
          if (    FromTable[i] != WorkSymbol[i]
              && (FromTable[i] != '^' || isalpha(WorkSymbol[i])
                                      || isdigit(WorkSymbol[i])))
            return -1;
      }
      return 0;
    }

    int Schoenflies_as_SgNumber(const string& WorkSymbol)
    {
      for (int SgNumber = 1; SgNumber <= 230; SgNumber++) {
        if (CmpSchoenfliesSymbols(tables::Schoenflies_List[SgNumber],
                                  WorkSymbol) == 0)
          return SgNumber;
      }
      return 0;
    }

    string GetStdTableId(const string& TableId)
    {
      string StdTableId = RemoveSpaces(TableId);
      if (StdTableId.size() > 0) {
        StdTableId[0] = toupper(StdTableId[0]);
        if (StdTableId == "I" || StdTableId == "I1952" || StdTableId == "1")
          StdTableId = "I1952";
        else if (StdTableId == "A" || StdTableId == "A1983")
          StdTableId = "A1983";
        else
          throw error("TableId not recognized.");
      }
      return StdTableId;
    }

    void ShortMonoHM_as_FullMonoHM(const string& StdTableId,
                                   string& WorkSymbol)
    {
      const tables::Short_Mono_HM_Dict_Entry* Entry;
      if (StdTableId == "I1952") Entry = tables::VolI_Short_Mono_HM_Dict;
      else                       Entry = tables::VolA_Short_Mono_HM_Dict;
      for (int i = 0; Entry[i].shrt; i++) {
        if (WorkSymbol == Entry[i].shrt) {
          WorkSymbol = Entry[i].full;
          return;
        }
      }
    }

    const tables::Main_Symbol_Dict_Entry*
    find_Main_Symbol_Dict_Entry(const string& WorkSymbol)
    {
      const tables::Main_Symbol_Dict_Entry* Entry;
      for (Entry = tables::Main_Symbol_Dict; Entry->SgNumber; Entry++) {
        string HM = RemoveSpaces(Entry->Hermann_Mauguin);
        if (HM == WorkSymbol) return Entry;
        string::size_type bar = HM.find('-');
        if (bar != string::npos) {
          // reverse "-N" to "N-", e.g. "P -1" -> "P 1-"
          HM[bar] = HM[bar + 1];
          HM[bar + 1] = '-';
          if (HM == WorkSymbol) return Entry;
          if (   (Entry->SgNumber >= 200 && Entry->SgNumber <= 206)
              || (Entry->SgNumber >= 221 && Entry->SgNumber <= 230)) {
            HM.erase(bar + 1, 1); // remove '-', e.g. "P m -3 m" -> "P m 3 m"
            if (HM == WorkSymbol) return Entry;
          }
        }
      }
      throw symbol_not_recognized;
    }

    const tables::Main_Symbol_Dict_Entry*
    SgNumber_to_Main_Symbol_Dict_Entry(int SgNumber, const string& StdTableId)
    {
      if (SgNumber < 1 || SgNumber > 230) {
        throw error("Space group number out of range.");
      }
      if (SgNumber < 3 || SgNumber > 15) {
        const tables::Main_Symbol_Dict_Entry* Entry;
        for (Entry = tables::Main_Symbol_Dict; Entry->SgNumber; Entry++)
          if (Entry->SgNumber == SgNumber)
            return Entry;
        throw cctbx_internal_error();
      }
      int i = 0;
      if (StdTableId == "I1952") i = 1;
      string MonoHM = tables::Monoclinic_SgNumber_as_HM_List[SgNumber][i];
      try {
        return find_Main_Symbol_Dict_Entry(MonoHM);
      }
      catch (const error&) {
        throw cctbx_internal_error();
      }
    }

    const char* SelectHall(const tables::Main_Symbol_Dict_Entry* Entry,
                           char& WorkExtension,
                           const string& StdTableId)
    {
      const char* Hall2 = &Entry->Hall[string(Entry->Hall).size() + 1];
      if (string(Hall2).size() == 0) {
        if (WorkExtension == '\0') return Entry->Hall;
      }
      else {
        if (Entry->Hermann_Mauguin[0] == 'R') {
          if (WorkExtension == '\0') {
            if (StdTableId == "I1952") WorkExtension = 'R';
            else                       WorkExtension = 'H';
          }
          if      (WorkExtension == 'H') return Entry->Hall;
          else if (WorkExtension == 'R') return Hall2;
        }
        else {
          if (WorkExtension == '\0') {
            if (StdTableId.size() == 0) WorkExtension = '2';
            else                        WorkExtension = '1';
          }
          if      (WorkExtension == '1') return Entry->Hall;
          else if (WorkExtension == '2') return Hall2;
        }
      }
      throw symbol_not_recognized;
    }

  } // namespace symbols

  using namespace symbols;

  void SpaceGroupSymbols::SetAll(const tables::Main_Symbol_Dict_Entry* Entry,
                                 char WorkExtension,
                                 const string& StdTableId)
  {
    const char* TableHall = SelectHall(Entry, WorkExtension, StdTableId);
    cctbx_assert(   WorkExtension == '\0'
                 || WorkExtension == '1' || WorkExtension == '2'
                 || WorkExtension == 'H' || WorkExtension == 'R');
    m_SgNumber        = Entry->SgNumber;
    m_Schoenflies     = tables::Schoenflies_List[Entry->SgNumber];
    m_Qualifier       = string(Entry->Qualifier ? Entry->Qualifier : "");
    m_Hermann_Mauguin = Entry->Hermann_Mauguin;
    m_Extension       = WorkExtension;
    if (m_Extension == '\0') m_ExtendedHermann_Mauguin = m_Hermann_Mauguin;
    else m_ExtendedHermann_Mauguin = m_Hermann_Mauguin + " :" + m_Extension;
    m_Hall            = TableHall;
  }

  int SpaceGroupSymbols::HallPassThrough(const string& Symbol)
  {
    string::const_iterator s;
    for (s = Symbol.begin(); s != Symbol.end(); s++) {
      if (!isspace(*s)) break;
    }
    for (int i = 0; i < 4; i++, s++) {
      if (s == Symbol.end()) return 0;
      if (tolower(*s) != "hall"[i]) return 0;
    }
    if (*s == ':') s++;
    else if (!isspace(*s)) return 0;
    for (; s != Symbol.end(); s++) if (!isspace(*s)) break;
    if (s == Symbol.end()) return 0;
    m_Hall = string(s, Symbol.end());
    return 1;
  }

  void SpaceGroupSymbols::Clear()
  {
    m_SgNumber = 0;
    m_Schoenflies = "";
    m_Qualifier = "";
    m_Hermann_Mauguin = "";
    m_Extension = '\0';
    m_Hall = "";
  }

  namespace detail {

    bool isHall(const string& TableId)
    {
      string::const_iterator s;
      for (s = TableId.begin(); s != TableId.end(); s++) {
        if (!isspace(*s)) break;
      }
      for (std::size_t i = 0; i < 4; i++, s++) {
        if (s == TableId.end()) return false;
        if (tolower(*s) != "hall"[i]) return false;
      }
      for (;; s++) {
        if (s == TableId.end()) return true;
        if (!isspace(*s)) break;
      }
      return false;
    }

  }

  SpaceGroupSymbols::SpaceGroupSymbols(const string& Symbol,
                                       const string& TableId)
  {
    Clear();
    if (detail::isHall(TableId)) {
      m_Hall = Symbol;
      return;
    }
    string StdTableId = GetStdTableId(TableId);
    if (StdTableId.size() == 0) {
      if (HallPassThrough(Symbol) != 0) return;
    }
    string WorkSymbol = PreProcessSymbol(Symbol);
    char WorkExtension = StripExtension(WorkSymbol);
    WorkSymbol[0] = toupper(WorkSymbol[0]);
    RemoveParentheses(WorkSymbol);

    const tables::Main_Symbol_Dict_Entry* Entry = 0;
    int SgNo;
    char xtrac;
    int n = sscanf(WorkSymbol.c_str(), "%d%c", &SgNo, &xtrac);
    if (n == 1) {
      Entry = SgNumber_to_Main_Symbol_Dict_Entry(SgNo, StdTableId);
    }
    else {
      SgNo = Schoenflies_as_SgNumber(WorkSymbol);
      if (SgNo != 0) {
        try {
          Entry = SgNumber_to_Main_Symbol_Dict_Entry(SgNo, StdTableId);
        }
        catch (const error&) {
          throw cctbx_internal_error();
        }
      }
    }
    if (Entry == 0) {
      ShortMonoHM_as_FullMonoHM(StdTableId, WorkSymbol);
      Entry = find_Main_Symbol_Dict_Entry(WorkSymbol);
    }
    SetAll(Entry, WorkExtension, StdTableId);
  }

  SpaceGroupSymbols::SpaceGroupSymbols(int SgNumber, const string& Extension,
                                       const string& TableId)
  {
    Clear();
    string StdTableId = GetStdTableId(TableId);
    string WorkSymbol = PreProcessSymbol(Extension);
    if (WorkSymbol.size() != 0 && WorkSymbol[0] != ':')
      WorkSymbol.insert(0, ":");
    char WorkExtension = StripExtension(WorkSymbol);
    if (WorkSymbol.size() != 0) throw symbol_not_recognized;

    const tables::Main_Symbol_Dict_Entry* Entry = 0;
    Entry = SgNumber_to_Main_Symbol_Dict_Entry(SgNumber, StdTableId);
    SetAll(Entry, WorkExtension, StdTableId);
  }

  SpaceGroupSymbols::SpaceGroupSymbols(
    const symbols::tables::Main_Symbol_Dict_Entry* Entry,
    char Extension)
  {
    Clear();
    if (Entry->SgNumber) SetAll(Entry, Extension, "");
  }

  SpaceGroupSymbolIterator::SpaceGroupSymbolIterator()
    : Entry(symbols::tables::Main_Symbol_Dict), nHall(1), iHall(0) {}

  SpaceGroupSymbols SpaceGroupSymbolIterator::next()
  {
    const symbols::tables::Main_Symbol_Dict_Entry* CurrentEntry = Entry;
    char Extension = '\0';
    if (Entry->SgNumber) {
      if (nHall == 2) {
        if (143 <= Entry->SgNumber && Entry->SgNumber < 168) {
          Extension = "HR"[iHall];
        }
        else {
          Extension = "12"[iHall];
        }
      }
      iHall++;
      if (iHall == nHall) {
        Entry++;
        nHall = 0;
        iHall = 0;
        if (Entry->SgNumber) {
          nHall = 1;
          const char* Hall2 = &Entry->Hall[string(Entry->Hall).size() + 1];
          if (string(Hall2).size() != 0) nHall = 2;
        }
      }
    }
    return SpaceGroupSymbols(CurrentEntry, Extension);
  }

} // namespace sgtbx

/*
http://www.iucr.org/iucr-top/comm/cnom/symsym/node7.html
Space group No. 39 41 64 67 68
Symbol in ITA83:  Abm2 Aba2 Cmca Cmma Ccca
New symbol: Aem2 Aea2 Cmce Cmme Ccce.
*/
