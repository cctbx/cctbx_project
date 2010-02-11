#ifndef CCTBX_MASKS_FLOOD_FILL_H
#define CCTBX_MASKS_FLOOD_FILL_H

#include <cctbx/import_scitbx_af.h>
#include <scitbx/array_family/accessors/c_grid_periodic.h>
#include <scitbx/array_family/tiny_algebra.h>

#include <stack>

namespace cctbx { namespace masks {

  namespace af = scitbx::af;

  //! Separates a mask into individual voids.
  /*!
     Takes a mask of 0's and 1's where 0 indicates outside a void and 1 is
     inside the void.  Individual voids will be labelled as 2, 3, etc...

     This is not the most efficient implementation, but the time spent in
     this function is small relative to that taken computing the mask in the
     first place.
   */
  template <typename DataType, typename FloatType>
  class flood_fill
  {
  public:
    flood_fill(af::ref<DataType, af::c_grid_periodic<3> > const & data)
      :
    gridding_n_real(data.accessor())
    {
      std::stack<index_t> stack;
      DataType target = 1;
      DataType replacement = 2;
      for (std::size_t i=0; i<gridding_n_real[0]; i++) {
        for (std::size_t j=0; j<gridding_n_real[1]; j++) {
          for (std::size_t k=0; k<gridding_n_real[2]; k++) {
            if (data(i,j,k) == target) {
              stack.push(index_t(i,j,k));
              data(i,j,k) = replacement;
              accumulated_indices.push_back(index_t(0,0,0));
              grid_points_per_void.push_back(0);
              while (!stack.empty()) {
                index_t index = stack.top();
                stack.pop();
                accumulated_indices[accumulated_indices.size() - 1] += index;
                grid_points_per_void[grid_points_per_void.size() - 1]++;
                for (std::size_t i=0; i<3; i++) {
                  index_t index_1 = index;
                  index_1[i]++;
                  if (data(index_1) == target) {
                    data(index_1) = replacement;
                    stack.push(index_1);
                  }
                  index_1[i] = index[i]-1;
                  if (data(index_1) == target) {
                    data(index_1) = replacement;
                    stack.push(index_1);
                  }
                }
              }
              replacement++;
            }
          }
        }
      }
    }

    unsigned n_voids()
    {
      return grid_points_per_void.size();
    }

    //! Provides the average indices for each void.
    /*! Move coordinates (if necessary) so the fractional coordinates will be
        between -1 and 1. Where there is a channel, an index along the
        direction of the channel may drift many unit cells along. In the case
        of a channel this index (along the channel) is not really applicable
        anyway.
     */
    af::shared<scitbx::vec3<FloatType> > averaged_indices()
    {
      af::shared<scitbx::vec3<FloatType> > result(
        af::reserve(accumulated_indices.size()));
      for (std::size_t i=0; i<n_voids(); i++) {
        result.push_back(scitbx::vec3<FloatType>(accumulated_indices[i]) /
          ((FloatType)grid_points_per_void[i]));
        for (std::size_t j=0; j<3; j++) {
          while (result[i][j] > gridding_n_real[j]) {
            result[i][j] -= gridding_n_real[j];
          }
          while (result[i][j] < -gridding_n_real[j]) {
            result[i][j] += gridding_n_real[j];
          }
        }
      }
      return result;
    }

    //! Provides the average fractional coordinates for each void.
    af::shared<scitbx::vec3<FloatType> > averaged_frac_coords()
    {
      af::shared<scitbx::vec3<FloatType> > result = averaged_indices();
      for (std::size_t i=0; i<result.size(); i++) {
        for (std::size_t j=0; j<3 ; j++) {
          result[i][j] /= (FloatType)gridding_n_real[j];
        }
      }
      return result;
    }

    af::shared<int> grid_points_per_void;

  private:
    typedef af::c_grid_periodic<3>::index_type index_t;
    index_t const gridding_n_real;
    af::shared<index_t> accumulated_indices;
  };

}} // namespace cctbx::masks

#endif // CCTBX_MASKS_FLOOD_FILL_H
