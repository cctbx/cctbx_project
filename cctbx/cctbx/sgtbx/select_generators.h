// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Sep 10: created from fragments in type.cpp, seminvariant.cpp (rwgk)
 */

#ifndef CCTBX_SGTBX_SELECT_GENERATORS_H
#define CCTBX_SGTBX_SELECT_GENERATORS_H

#include <cctbx/sgtbx/groups.h>

namespace sgtbx {
  namespace detail {

    struct AnyGenerators {
      AnyGenerators() : nGen(0) {}
      AnyGenerators(const SgOps& sgo);
      inline int nAll() const {
        if (ZInvT.isValid()) return nGen + 1;
        return nGen;
      }
      void setPrimitive();
      ChOfBasisOp Z2POp;
      TrVec ZInvT;
      TrVec PInvT;
      int nGen;
      RTMx ZGen[2];
      RTMx PGen[2];
    };

    struct StdGenerators: AnyGenerators {
      StdGenerators(const SgOps& WorkSgOps,
                    const tables::MatrixGroup::Code& PG_MGC);
    };

  } // namespace detail
} // namespace sgtbx

#endif // CCTBX_SGTBX_SELECT_GENERATORS_H
