// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created from fragments of cctbx/sgtbx/miller_asu.h (rwgk)
 */

#ifndef CCTBX_MILLER_ASU_H
#define CCTBX_MILLER_ASU_H

#include <cctbx/sgtbx/miller_asu.h>

namespace cctbx { namespace miller {

  template <typename DataType>
  class map_to_asu
  {
    public:
      map_to_asu() {}

      map_to_asu(
        const sgtbx::SpaceGroupInfo& sginfo,
        bool friedel_flag,
        af::shared<miller::Index> miller_indices,
        af::shared<DataType> data_array,
        bool in_place = false)
        : friedel_flag_(friedel_flag),
          asu_(sginfo)
      {
        cctbx_assert(miller_indices.size() == data_array.size());
        const sgtbx::SpaceGroup& sgops = sginfo.SgOps();
        for(std::size_t i=0;i<miller_indices.size();i++) {
          AsymIndex ai(sgops, asu_, miller_indices[i]);
          IndexTableLayoutAdaptor ila = ai.one_column(friedel_flag);
          if (in_place) {
            miller_indices[i] = ila.H();
            data_array[i] = ila.complex_eq(data_array[i]);
          }
          else {
            asym_miller_indices_.push_back(ila.H());
            asym_data_array_.push_back(ila.complex_eq(data_array[i]));
          }
        }
        if (in_place) {
          asym_miller_indices_ = miller_indices;
          asym_data_array_ = data_array;
        }
      }

      bool friedel_flag() const { return friedel_flag_; }

      sgtbx::ReciprocalSpaceASU const& asu() const { return asu_; }

      af::shared<miller::Index> asym_miller_indices() const {
        return asym_miller_indices_;
      }

      af::shared<DataType> asym_data_array() const {
        return asym_data_array_;
      }
    protected:
      bool friedel_flag_;
      sgtbx::ReciprocalSpaceASU asu_;
      af::shared<miller::Index> asym_miller_indices_;
      af::shared<DataType> asym_data_array_;
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_ASU_H
