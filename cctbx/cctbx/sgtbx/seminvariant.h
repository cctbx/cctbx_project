// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Sep 02: start port of sglite/sgss.c (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_SEMINVARIANT_H
#define CCTBX_SGTBX_SEMINVARIANT_H

#include <cctbx/static_vector.h>
#include <cctbx/sgtbx/groups.h>

namespace sgtbx {

  struct ssVM {
    Vec3 V;
    int M;
  };

  class StructureSeminvariant {
    public:
      StructureSeminvariant() {}
      StructureSeminvariant(const SgOps& sgo);
      inline std::size_t size() const { return m_VM.size(); }
      inline const ssVM& VM(std::size_t i) const {
        if (i >= size()) throw error_index();
        return m_VM[i];
      }
      inline const Vec3& V(std::size_t i) const { return VM(i).V; }
      inline int M(std::size_t i) const { return VM(i).M; }
      bool is_ss(const Miller::Index& H) const;
      Vec3 get_uvw(const Miller::Index& H) const;
    private:
      cctbx::static_vector<ssVM, 3> m_VM;
  };

} // namespace sgtbx

#endif // CCTBX_SGTBX_SEMINVARIANT_H
