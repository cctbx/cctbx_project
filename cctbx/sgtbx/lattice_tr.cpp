/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2001 Jul: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 Apr: Implementation of ConstructZ2POp() (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/sgtbx/lattice_tr.h>
#include <ctype.h>

namespace cctbx { namespace sgtbx { namespace lattice_tr {

  namespace conventional_centring_types {

    namespace {

      static const tr_vec p[] = { tr_vec_12(0, 0, 0),
                                };
      static const tr_vec a[] = { tr_vec_12(0, 0, 0),
                                  tr_vec_12(0, 6, 6),
                                };
      static const tr_vec b[] = { tr_vec_12(0, 0, 0),
                                  tr_vec_12(6, 0, 6),
                                };
      static const tr_vec c[] = { tr_vec_12(0, 0, 0),
                                  tr_vec_12(6, 6, 0),
                                };
      static const tr_vec i[] = { tr_vec_12(0, 0, 0),
                                  tr_vec_12(6, 6, 6),
                                };
      static const tr_vec r[] = { tr_vec_12(0, 0, 0),
                                  tr_vec_12(8, 4, 4),
                                  tr_vec_12(4, 8, 8),
                                };
      static const tr_vec q[] = { tr_vec_12(0, 0, 0), // reverse setting
                                  tr_vec_12(4, 8, 4), // for internal use
                                  tr_vec_12(8, 4, 8), // only
                                };
      static const tr_vec h[] = { tr_vec_12(0, 0, 0),
                                  tr_vec_12(8, 4, 0),
                                  tr_vec_12(4, 8, 0),
                                };
      static const tr_vec f[] = { tr_vec_12(0, 0, 0),
                                  tr_vec_12(0, 6, 6),
                                  tr_vec_12(6, 0, 6),
                                  tr_vec_12(6, 6, 0),
                                };

    } // namespace <anonymous>

    const table_entry* table()
    {
      static const table_entry table_[] = {
        { 'P', 1, p },
        { 'A', 2, a },
        { 'B', 2, b },
        { 'C', 2, c },
        { 'I', 2, i },
        { 'R', 3, r },
        { 'Q', 3, q }, // reverse setting
        { 'H', 3, h },
        { 'F', 4, f },
        { '\0', 0, 0 },
      };
      return table_;
    }

    const table_entry* get(char symbol)
    {
      symbol = toupper(symbol);
      if (symbol == 'Q') return 0;
      for (const table_entry* e = table(); e->symbol != '\0'; e++) {
        if (e->symbol == symbol) return e;
      }
      return 0;
    }

  } // namespace conventional_centring_types

  namespace conventional_z2p_matrices {

    namespace {

      // change-of-basis matrices (coordinate transformations)

      static const rot_mx pp(
         1,  0,  0,
         0,  1,  0,
         0,  0,  1 );
      static const rot_mx ap(
        -1,  0,  0,
         0, -1,  1,
         0,  1,  1 );
      static const rot_mx bp(
        -1,  0,  1,
         0, -1,  0,
         1,  0,  1 );
      static const rot_mx cp(
         1,  1,  0,
         1, -1,  0,
         0,  0, -1 );
      static const rot_mx ip(
         0,  1,  1,
         1,  0,  1,
         1,  1,  0 );
      static const rot_mx rp(
         1,  0,  1,
        -1,  1,  1,
         0, -1,  1 );
      static const rot_mx hp(
         1,  1,  0,
        -1,  2,  0,
         0,  0,  1 );
      static const rot_mx fp(
        -1,  1,  1,
         1, -1,  1,
         1,  1, -1 );

    } // namespace <anonymous>

    rot_mx const& get(char symbol)
    {
      switch (symbol) {
        case 'P': return pp;
        case 'A': return ap;
        case 'B': return bp;
        case 'C': return cp;
        case 'I': return ip;
        case 'R': return rp;
        case 'H': return hp;
        case 'F': return fp;
        default: break;
      }
      static rot_mx null(0);
      return null;
    }

  } // namespace conventional_z2p_matrices

}}} // namespace cctbx::sgtbx::lattice_tr
