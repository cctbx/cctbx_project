// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created from fragments of cctbx/miller.h (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MILLER_SPAN_H
#define CCTBX_MILLER_SPAN_H

#include <cctbx/miller.h>
#include <cctbx/math/utils.h>
#include <cctbx/array_family/misc_functions.h>
#include <cctbx/array_family/shared.h>

namespace cctbx { namespace miller {

  /*! \brief Determine min(H[i]) and max(H[i])+1, i=1..3,
       for an array of Miller indices.
   */
  struct index_span : af::tiny<af::tiny<int, 2>, 3>
  {
    typedef af::tiny<int, 2> min_end;
    typedef af::tiny<min_end, 3> base_class;

    index_span() {}

    index_span(af::shared<Index> index_array)
    {
      this->fill(min_end(0, 0));
      if (index_array.size()) {
        for(std::size_t j=0;j<3;j++) {
          (*this)[j] = min_end().fill(index_array[0][j]);
        }
      }
      for(std::size_t i=1;i<index_array.size();i++) {
        for(std::size_t j=0;j<3;j++) {
          math::update_min((*this)[j][0], index_array[i][j]);
          math::update_max((*this)[j][1], index_array[i][j]);
        }
      }
      for(std::size_t j=0;j<3;j++) (*this)[j][1]++;
    }

    af::int3 min() const
    {
      af::int3 result;
      for(std::size_t j=0;j<3;j++) result[j] = (*this)[j][0];
      return result;
    }

    af::int3 max() const
    {
      af::int3 result;
      for(std::size_t j=0;j<3;j++) result[j] = (*this)[j][1] - 1;
      return result;
    }

    af::int3 abs_range() const
    {
      af::int3 result;
      std::size_t j;
      for(j=0;j<3;j++) {
        result[j] = fn::absolute((*this)[j][0]);
        math::update_max(result[j], fn::absolute((*this)[j][1]-1));
      }
      for(j=0;j<3;j++) result[j] += 1;
      return result;
    }

    af::int3 map_grid() const
    {
      af::int3 result = abs_range();
      for(std::size_t j=0;j<3;j++) {
        result[j] = (result[j] - 1) * 2 + 1;
      }
      return result;
    }

    bool is_in_domain(Index const& h) const
    {
      for(std::size_t j=0;j<3;j++) {
        if (!((*this)[j][0] <= h[j] && h[j] < (*this)[j][1])) return false;
      }
      return true;
    }

    std::size_t pack(Index const& h) const
    {
      return ((h[0] - (*this)[0][0]) * range_((*this)[1])
            + (h[1] - (*this)[1][0])) * range_((*this)[2])
            + (h[2] - (*this)[2][0]);
    }

    private:
      static int range_(min_end const& span)
      {
        return span[1] - span[0];
      }
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_SPAN_H
