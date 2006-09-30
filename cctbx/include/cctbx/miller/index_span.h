#ifndef CCTBX_MILLER_INDEX_SPAN_H
#define CCTBX_MILLER_INDEX_SPAN_H

#include <cctbx/miller.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/tiny_types.h>

namespace cctbx { namespace miller {

  /*! \brief Determines min(indices[i]) and max(indices[i])+1, i=1..3,
       for an array of Miller indices.
   */
  class index_span : af::tiny<af::tiny<int, 2>, 3>
  {
    public:
      typedef af::tiny<int, 2> min_end;
      typedef af::tiny<min_end, 3> base_class;

      //! Default constructor. Some data members are not initialized!
      index_span() {}

      //! Determines the min and max elements in the array of indices.
      index_span(af::const_ref<index<> > const& indices);

      //! Constructor with specified min and max values
      index_span(af::tiny< af::tiny<int, 3>, 2 > const& min_max_miller)
      {
        *this = min_max_miller;
        (*this)[0][1]++;
        (*this)[1][1]++;
        (*this)[2][1]++;
      }

      //! The min elements found in the array passed to the constructor.
      af::int3
      min() const;

      //! The max elements found in the array passed to the constructor.
      af::int3
      max() const;

      //! Maximum of abs(min()) and abs(max()) + 1.
      af::int3
      abs_range() const;

      /*! \brief Dimensions of 3-dimensional array for storing the indices
          found in the 1-dimensional array passed to the constructor.
       */
      /*! Formula used: (abs_range()-1) * 2 + 1
       */
      af::int3
      map_grid() const;

      //! Tests if min() <= h <= max().
      bool
      is_in_domain(index<> const& h) const;

      //! Computes a 1-dimensional index for h.
      /*! Useful for generating fast lookup maps.
          <p>
          Not available in Python.
       */
      std::size_t
      pack(index<> const& h) const
      {
        return ((h[0] - (*this)[0][0]) * range_((*this)[1])
              + (h[1] - (*this)[1][0])) * range_((*this)[2])
              + (h[2] - (*this)[2][0]);
      }

      //! Computes 1-dimensional indices for given Miller indices.
      af::shared<std::size_t>
      pack(af::const_ref<index<> > const& miller_indices) const
      {
        af::shared<std::size_t> result((af::reserve(miller_indices.size())));
        for(std::size_t i=0;i<miller_indices.size();i++) {
          result.push_back(pack(miller_indices[i]));
        }
        return result;
      }

    private:
      static int
      range_(min_end const& span) { return span[1] - span[0]; }
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_INDEX_SPAN_H
