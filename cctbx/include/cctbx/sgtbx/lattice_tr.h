/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_LATTICE_TR_H
#define CCTBX_SGTBX_LATTICE_TR_H

#include <cctbx/sgtbx/rot_mx.h>

namespace cctbx { namespace sgtbx { namespace lattice_tr {

  namespace conventional_centring_types {

    struct table_entry
    {
      char symbol;
      std::size_t n_translations;
      const tr_vec* t;
    };

    const table_entry* table();

    const table_entry* get(char symbol);
  }

  namespace conventional_z2p_matrices {

    rot_mx const& get(char symbol);

  }

}}} // namespace cctbx::sgtbx::lattice_tr

#endif // CCTBX_SGTBX_LATTICE_TR_H
