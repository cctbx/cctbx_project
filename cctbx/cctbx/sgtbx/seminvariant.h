// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Sep 02: start port of sglite/sgss.c (R.W. Grosse-Kunstleve)
 */

#include <cctbx/static_vector.h>
#include <cctbx/sgtbx/groups.h>

namespace sgtbx {

  class ssVM {
    public:
      Vec3 V;
      int M;
      inline void zero_out() {
        V.assign(0);
        M = -1;
      }
  };

  namespace detail {
    class AnyGenerators;
    class DiscrList;
    const int mDiscrList = 8;
  }

  class StructureSeminvariant {
    public:
      StructureSeminvariant() {}
      StructureSeminvariant(const SgOps& sgo);
      inline std::size_t size() const { return m_size; }
      inline const ssVM& VM(std::size_t i) const {
        if (i >= size()) throw error_index();
        return m_VM[i];
      }
    private:
      boost::array<ssVM, 3> m_VM;
      std::size_t m_size;
      void GetContNullSpace(const detail::AnyGenerators& Gen);
      void BestVectors(
        const SgOps& sgo,
        cctbx::static_vector<detail::DiscrList,
                             detail::mDiscrList>& DiscrLst) const;
  };

} // namespace sgtbx
