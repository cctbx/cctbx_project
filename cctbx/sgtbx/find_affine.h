#ifndef CCTBX_SGTBX_FIND_AFFINE_H
#define CCTBX_SGTBX_FIND_AFFINE_H

#include <cctbx/sgtbx/space_group.h>

namespace cctbx { namespace sgtbx {

  //! Find affine normalizers of triclinic and monoclinic space groups.
  /*! Do not use for higher symmetry space groups. For these it is
      much more efficient to use the tabulated affine normalizers.
   */
  class find_affine
  {
    public:
      //! Default constructor. Some data members are not initialized!
      find_affine() {}

      //! Exhaustive search.
      /*! Generates all matrices c with elements in the interval
          [-range, +range] and determinant one. For each of these
          it is tested if the space group is invariant under
          c*s*c.inverse(). If so the matrix c is stored in cb_mx_.
          If use_p1_algorithm == false (default) and the space group
          is different from P1 a much faster algorithm is used
          to generate the change-of-basis matrices.
          Set use_p1_algorithm=true for testing purposes only.
       */
      find_affine(
        space_group const& group,
        int range=2,
        bool use_p1_algorithm=false);

      //! Change-of-basis matrices that leave the space group invariant.
      af::shared<rt_mx>
      cb_mx() const { return cb_mx_; }

    protected:
      af::shared<rt_mx> cb_mx_;

      void
      p1_algorithm(space_group const& group, int range);

      void
      sg_algorithm(space_group const& group, int range);
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_FIND_AFFINE_H
