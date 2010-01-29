#ifndef CCTBX_MASKS_FLOOD_FILL_H
#define CCTBX_MASKS_FLOOD_FILL_H

#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/accessors/c_grid.h>

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
  template <typename DataType>
  void flood_fill(af::ref<DataType, af::c_grid<3> > const & data)
  {
    typedef af::c_grid<3>::index_type index_t;
    index_t const& gridding_n_real = data.accessor();
    DataType replacement = 2;
    DataType target = 1;
    std::stack<index_t> stack;
    for (std::size_t i=0; i<gridding_n_real[0]; i++) {
      for (std::size_t j=0; j<gridding_n_real[1]; j++) {
        for (std::size_t k=0; k<gridding_n_real[2]; k++) {
          if (data(i,j,k) == target) {
            stack.push(index_t(i,j,k));
            while (!stack.empty()) {
              index_t index = stack.top();
              stack.pop();
              data(index) = replacement;
              for (std::size_t i=0; i<3; i++) {
                index_t index_1 = index;
                index_1[i]++;
                if (index_1[i]==gridding_n_real[i]) index_1[i]=0;
                if (data(index_1) == target) {
                  stack.push(index_1);
                }
                if (index[i]==0) index_1[i] = gridding_n_real[i] - 1;
                else index_1[i] = index[i]-1;
                if (data(index_1) == target) {
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

}} // namespace cctbx::masks

#endif // CCTBX_MASKS_FLOOD_FILL_H
