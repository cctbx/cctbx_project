/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Refactored parts of cctbx/sgtbx/matrix.h (rwgk)
     2001 Jul: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_ROT_MX_INFO_H
#define CCTBX_SGTBX_ROT_MX_INFO_H

#include <cctbx/sgtbx/rot_mx.h>

namespace cctbx { namespace sgtbx {

  //! Class for grouping information about rotation matrices.
  /*! Determing the sense of rotation requires the rotation type
      and Eigenvector. Therefore it is most efficient to group
      all these properties so that the intermediate results
      can also be used.
   */
  class rot_mx_info
  {
    public:
      //! Default constructor. Some data members are not initialized!
      rot_mx_info() {}

      //! Determination of all properties.
      rot_mx_info(rot_mx const& r);

      //! Rotation-part type (1, 2, 3, 4, 6, -1, -2=m, -3, -4, -6)
      int type() const { return type_; }

      //! Axis direction (Eigenvector) of the corresponding proper rotation.
      /*! Only defined if abs(type()) > 1.<br>
          For type() > 0, the proper rotation is defined as r.<br>
          For type() < 0, the proper rotation is defined as -r.
       */
      sg_vec3 const& ev() const { return ev_; }

      //! Sense of rotation with respect to the axis direction.
      /*! Only defined if abs(type()) > 1.
       */
      int sense() const { return sense_; }

    private:
      int type_;
      sg_vec3 ev_;
      int sense_;
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_ROT_MX_INFO_H
