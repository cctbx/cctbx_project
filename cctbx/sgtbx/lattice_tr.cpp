// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctype>
#include <cctbx/sgtbx/groups.h>

namespace sgtbx {
  namespace lattice {

    const CentringTypeMap* getConventionalCentringType(char Symbol)
    {
      Symbol = toupper(Symbol);
      if (Symbol == 'Q') return 0;
      for (const CentringTypeMap* map = tables::ConventionalCentringTypeMap;
           map->Symbol != '\0';
           map++) {
        if (map->Symbol == Symbol) return map;
      }
      return 0;
    }

    const RotMx& getConventionalZ2PMatrix(int Symbol)
    {
      using namespace tables::ConventionalZ2PMatrices;
      switch (Symbol) {
        case 'P': return PP;
        case 'A': return AP;
        case 'B': return BP;
        case 'C': return CP;
        case 'I': return IP;
        case 'R': return RP;
        case 'H': return HP;
        case 'F': return FP;
        default: break;
      }
      static RotMx null(0);
      return null;
    }

  } // namespace lattice

  char TrOps::getConventionalCentringTypeSymbol() const
  {
    for (const lattice::CentringTypeMap*
         map = lattice::tables::ConventionalCentringTypeMap;
         map->Symbol != '\0';
         map++)
    {
      if (map->nTrs == nVects()) {
        int nMatch = 0;
        int Match[4];
        cctbx_assert(map->nTrs <= sizeof Match / sizeof (*Match));
        int i, j;
        for(i=0;i<map->nTrs;i++) Match[i] = 0;
        for(i=0;i<map->nTrs;i++) {
          for(j=0;j<map->nTrs;j++) {
            if (Match[j] == 0 && map->Trs[i] == m_Vects[j]) {
              Match[j] = 1;
              nMatch++;
              break;
            }
          }
        }
        if (nMatch == map->nTrs) return map->Symbol;
      }
    }
    return '\0';
  }

  ChOfBasisOp TrOps::getConventionalZ2POp(int RBF, int TBF) const
  {
    char ZSymbol = getConventionalCentringTypeSymbol();
    const RotMx& Z2PMatrix = lattice::getConventionalZ2PMatrix(ZSymbol);
    if (! Z2PMatrix.isValid()) return ChOfBasisOp(0, 0);
    return ChOfBasisOp(RTMx(Z2PMatrix.newBaseFactor(RBF), TrVec(TBF)));
  }

  ChOfBasisOp SgOps::ConstructZ2POp(int RBF, int TBF) const
  {
    throw cctbx_not_implemented();
  }

  ChOfBasisOp SgOps::getZ2POp(int RBF, int TBF) const
  {
    ChOfBasisOp CBOp = m_LTr.getConventionalZ2POp(RBF, TBF);
    if (CBOp.isValid()) return CBOp;
    return ConstructZ2POp(RBF, TBF);
  }

} // namespace sgtbx
