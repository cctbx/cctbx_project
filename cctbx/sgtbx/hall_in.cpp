// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctype>
#include <stdio.h>
#include <string.h>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/sgtbx/tables.h>
#include <cctbx/basic/define_range.h>

namespace sgtbx {
  namespace hall {

    bool IsHSymSpace(char c)
    {
      if (c == '\0') return 0;
      if (c == '_') return 1;
      return isspace(c);
    }

    bool IsHSymChar(char c)
    {
      if (c == '\0') return 0;
      return ! IsHSymSpace(c);
    }

    int GetAbsOrder(char c)
    {
      if (c == '1') return 1;
      if (c == '2') return 2;
      if (c == '3') return 3;
      if (c == '4') return 4;
      if (c == '6') return 6;
      return 0;
    }

    int GetScrew(char c)
    {
      if (c == '1') return 1;
      if (c == '2') return 2;
      if (c == '3') return 3;
      if (c == '4') return 4;
      if (c == '5') return 5;
      return 0;
    }

    char GetRefAxis(char c)
    {
      c = tolower(c);
      if (   c == 'x'
          || c == 'y'
          || c == 'z') return c;
      return '\0';
    }

    char GetDirCode(char c)
    {
      if (   c == '\''
          || c ==  '"'
          || c ==  '*') return c;
      if (   c == ','
          || c == '.') return '\'';
      if (   c == ';'
          || c == ':') return  '"';
      return '\0';
    }

    const TrVec& GetTr(char Sym)
    {
      struct TrMap {
        char  Sym;
        TrVec Vec;
        inline TrMap(const char S, const TrVec& V) : Sym(S), Vec(V) {}
      };
      static const TrMap HallTrs[] =
      {
        TrMap('a', TrVec12(6, 0, 0)),
        TrMap('b', TrVec12(0, 6, 0)),
        TrMap('c', TrVec12(0, 0, 6)),
        TrMap('n', TrVec12(6, 6, 6)),
        TrMap('u', TrVec12(3, 0, 0)),
        TrMap('v', TrVec12(0, 3, 0)),
        TrMap('w', TrVec12(0, 0, 3)),
        TrMap('d', TrVec12(3, 3, 3)),
      };
      const int nHallTrs = sizeof HallTrs / sizeof (*HallTrs);
      Sym = tolower(Sym);
      for (int i = 0; i < nHallTrs; i++)
        if (HallTrs[i].Sym == Sym) return HallTrs[i].Vec;
      static TrVec null(0);
      return null;
    }

    const RotMx GetRMx(bool Improper, int AbsOrder,
                       char RefAxis, char DirCode)
    {
      struct TabRMxEntry
      {
        int          Order;
        char         DirCode;
        const RotMx* RMx;
      };
      using namespace tables::RotationMatrices;
      const TabRMxEntry TabRMx[] = {
        { 1, '\0', &R_1_000 },
        { 2, '\0', &R_2_001 },
        { 2, '\'', &R_2_1b0 },
        { 2,  '"', &R_2_110 },
        { 3, '\0', &R_3_001 },
        { 3,  '*', &R_3_111 },
        { 4, '\0', &R_4_001 },
        { 6, '\0', &R_6_001 },
      };
      const int n = (sizeof TabRMx / sizeof (*TabRMx));
      rangei(n) {
        if (   TabRMx[i].Order == AbsOrder
            && TabRMx[i].DirCode == DirCode) {
          RotMx R;
          if (! Improper) R =  (*TabRMx[i].RMx);
          else            R = -(*TabRMx[i].RMx);
          if (RefAxis == 'x') return R_3_111 * R * R_3i111;
          if (RefAxis == 'y') return R_3i111 * R * R_3_111;
          return R;
        }
      }
      throw cctbx_internal_error();
    }

    TrVec ParseShortCBO(parse_string& HSym, const char* StopChars, int TBF)
    {
      cctbx_assert(TBF % 12 == 0);
      TrVec result(TBF);
      for (int Row = 0; Row < 3; Row++) {
        while (IsHSymSpace(HSym())) HSym.skip();
        if (Row && HSym() == ',') {
          HSym.skip();
          while (IsHSymSpace(HSym())) HSym.skip();
        }
        if (strchr(StopChars, HSym()))
          return TrVec(0);
        int i = 1;
        int n = sscanf(HSym.peek(), "%d%n", &result[Row], &i);
        HSym.skip(i - 1);
        if (n != 1) return TrVec(0);
        HSym.skip();
        result[Row] *= (TBF / 12);
      }
      return result;
    }

  } // namespace hall

  int SgOps::ParseHallSymbolCBOp(parse_string& HSym, ChOfBasisOp& CBOp,
                                 bool Pedantic, bool NoCType)
  {
    using namespace hall;

    int nAddedMx = 0;

    // Interpret the lattice type code.
    if (! NoCType) {
      while (IsHSymSpace(HSym())) HSym.skip();
      if (HSym() == '-') {
        expandInv(TrVec(STBF));
        HSym.skip();
        nAddedMx++;
      }
      if (HSym() == '\0') throw error("Lattice type not specified.");
      nAddedMx += expandConventionalCentringType(HSym());
      HSym.skip();
    }

    const char char_after_lattice_type_symbol = HSym();
    while (IsHSymSpace(HSym())) HSym.skip();
    if (HSym() == '\0' || HSym() == '(') {
      if (Pedantic) throw error("Matrix symbol expected.");
      if (HSym() == '\0') return nAddedMx;
    }
    if (   ! NoCType
        && Pedantic
        && ! IsHSymSpace(char_after_lattice_type_symbol))
      throw error("Space expected after lattice type symbol.");

    // Loop over the matrix symbols.
    int  iMxSym = 0;
    int  FirstAbsOrder = 0;
    char FirstRefAxis  = '\0';
    while (HSym() != '\0' && HSym() != '(')
    {
      bool Improper = false;
      int  AbsOrder =  0;
      int  Screw = 0;
      char RefAxis = '\0';
      char DirCode = '\0';
      RotMx SMx_R;
      TrVec SMx_T;

      if (HSym() == '-') {
        Improper = true;
        HSym.skip();
        if (! IsHSymChar(HSym())) {
          throw error("Incomplete matrix symbol.");
        }
      }
            AbsOrder = GetAbsOrder(HSym());
      if (! AbsOrder) {
        throw error("Improper symbol for rotational order.");
      }
      HSym.skip();

          Screw = GetScrew(HSym());
      if (Screw) {
        if (Screw >= AbsOrder) {
          throw error("Improper screw translation.");
        }
        HSym.skip();
      }

      while (IsHSymChar(HSym()))
      {
        if (  RefAxis == '\0') {
              RefAxis = GetRefAxis(HSym());
          if (RefAxis != '\0') {
            if (    AbsOrder == 1
                || (AbsOrder == 3 && DirCode == '*')) {
              throw error("Inconsistent matrix symbol.");
            }
            HSym.skip();
            continue;
          }
        }
        else if (GetRefAxis(HSym()) != '\0') {
          throw error("Multiple axis symbols.");
        }

        if (  DirCode == '\0') {
              DirCode = GetDirCode(HSym());
          if (DirCode != '\0') {
            if (   ! (AbsOrder == 2 && (   DirCode ==  '"'
                                        || DirCode == '\''))
                && ! (AbsOrder == 3 && DirCode == '*')) {
              throw error("Inconsistent matrix symbol.");
            }
            if (Screw) {
              throw error("Screw translation for non-principal direction.");
            }
            HSym.skip();
            continue;
          }
        }
        else if (GetDirCode(HSym()) != '\0') {
          throw error("Multiple axis symbols.");
        }

        const TrVec& HTr = GetTr(HSym());
        if (HTr.isValid()) {
          SMx_T = (SMx_T + HTr).ModPositive();
          HSym.skip();
          continue;
        }

        if (HSym() == '(') {
          if (Pedantic) {
            throw error("Space expected before change-of-basis operator.");
          }
          break;
        }

        throw error("Malformed matrix symbol.");
      }

      if (RefAxis == '\0') {
        if      (iMxSym == 0) {
          if (      AbsOrder != 1
              && ! (AbsOrder == 3 && DirCode == '*'))
            RefAxis = 'z';
        }
        else if (iMxSym == 1) {
          if      (AbsOrder == 2) {
            if      (FirstAbsOrder == 2 || FirstAbsOrder == 4) {
              if (DirCode == '\0') {
                RefAxis = 'x';
              }
            }
            else if (FirstAbsOrder == 3 || FirstAbsOrder == 6) {
              if (DirCode == '\0') {
                DirCode = '\'';
              }
              RefAxis = FirstRefAxis;
            }
          }
          else if (   AbsOrder == 3
                   && (FirstAbsOrder == 2 || FirstAbsOrder == 4)
                   && DirCode == '\0') {
            DirCode = '*';
          }
        }
        else if (iMxSym == 2) {
          if (AbsOrder == 3 && DirCode == '\0') {
            DirCode = '*';
          }
        }
      }

      if (RefAxis == '\0' && (   DirCode ==  '"'
                              || DirCode == '\'')) {
        RefAxis = 'z';
      }

      if (RefAxis == '\0' && AbsOrder != 1 && DirCode != '*') {
        throw error("Need explicit axis symbol.");
      }

      SMx_R = GetRMx(Improper, AbsOrder, RefAxis, DirCode);

      if (Screw) {
        int i;
        switch (RefAxis) {
          case 'x': i = 0; break;
          case 'y': i = 1; break;
          default:  i = 2; break;
        }
        cctbx_assert(SMx_T.BF() * Screw % AbsOrder == 0);
        SMx_T[i] += SMx_T.BF() * Screw / AbsOrder;
      }

      expandSMx(RTMx(SMx_R, SMx_T));

      if (iMxSym == 0) {
        FirstAbsOrder = AbsOrder;
        FirstRefAxis  = RefAxis;
      }
      iMxSym++;

      if (Improper || AbsOrder != 1)
        nAddedMx++;

      while (IsHSymSpace(HSym())) HSym.skip();
    }

    // Interpret the change-of-basis operator.
    if (HSym() == '(') {
      HSym.skip();
      HSym.set_mark();
      RTMx V = RTMx(ParseShortCBO(HSym, ")", CTBF), CRBF);
      if (! V.isValid()) {
        HSym.go_to_mark();
        V = RTMx(HSym, ")", CRBF, CTBF);
        if (! V.isValid()) {
          throw error("Malformed change-of-basis operator.");
        }
      }
      while (IsHSymSpace(HSym())) HSym.skip();
      if (HSym() != ')') {
        throw error(
          "Closing parenthesis expected after change-of-basis operator");
      }
      try {
        CBOp = ChOfBasisOp(V);
      }
      catch (error) {
        throw error("Change-of-basis operator is not invertible.");
      }
      HSym.skip();
    }

    while (IsHSymSpace(HSym())) HSym.skip();
    if (HSym() != '\0') {
      throw error("Unexpected extra character.");
    }

    return nAddedMx;
  }

  int SgOps::ParseHallSymbol(parse_string& HSym, bool Pedantic, bool NoCType)
  {
    ChOfBasisOp CBOp(0, 0);
    int nAddedMx = ParseHallSymbolCBOp(HSym, CBOp, Pedantic, NoCType);
    if (CBOp.isValid()) {
      SgOps NewSgOps = ChangeBasis(CBOp);
      *this = NewSgOps;
    }
    return nAddedMx;
  }

  SgOps::SgOps(parse_string& HSym, bool Pedantic, bool NoCType,
               bool NoExpand)
    : m_NoExpand(NoExpand)
  {
    reset();
    ParseHallSymbol(HSym, Pedantic, NoCType);
  }

  SgOps::SgOps(const std::string& HSym, bool Pedantic, bool NoCType,
               bool NoExpand)
    : m_NoExpand(NoExpand)
  {
    reset();
    parse_string HSymPS(HSym);
    ParseHallSymbol(HSymPS, Pedantic, NoCType);
  }

  SgOps::SgOps(const char* HSym, bool Pedantic, bool NoCType,
               bool NoExpand)
    : m_NoExpand(NoExpand)
  {
    reset();
    parse_string HSymPS(HSym);
    ParseHallSymbol(HSymPS, Pedantic, NoCType);
  }

} // namespace sgtbx
