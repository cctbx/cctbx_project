/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2001 May: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
               Based on C code contributed by Vincent Favre-Nicolin.
 */

#include <cctbx/eltbx/fp_fdp.h>
#include <cstddef>

namespace cctbx { namespace eltbx {

  namespace anomalous {

    fp_fdp
    interpolate(const label_z_e_fp_fdp* label_z_e_fp_fdp_, double energy)
    {
      float fp = fp_fdp_undefined;
      float fdp = fp_fdp_undefined;
      const e_fp_fdp* data = label_z_e_fp_fdp_->data;
      float energy1 = data[0].e;
      float energy2 = data[1].e;
      std::size_t i;
      for(i = 2; energy2 > 0. && energy2 < energy; i++) {
        energy1 = energy2;
        energy2 = data[i].e;
      }
      if (energy >= energy1 && energy2 > 0.) {
        float f = (energy - energy1) / (energy2 - energy1);
        if (   data[i-2].fp != fp_fdp_undefined
            && data[i-1].fp != fp_fdp_undefined) {
          fp = data[i-2].fp + f * (data[i-1].fp - data[i-2].fp);
        }
        fdp = data[i-2].fdp + f * (data[i-1].fdp - data[i-2].fdp);
      }
      return fp_fdp(fp, fdp);
    }

  } // namespace anomalous

}} // namespace cctbx::eltbx
