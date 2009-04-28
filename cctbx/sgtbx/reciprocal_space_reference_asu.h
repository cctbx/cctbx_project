#ifndef CCTBX_SGTBX_RECIPROCAL_SPACE_REFERENCE_ASU_H
#define CCTBX_SGTBX_RECIPROCAL_SPACE_REFERENCE_ASU_H

#include <cctbx/error.h>
#include <cctbx/miller.h>
#include <cctbx/sgtbx/group_codes.h>
#include <scitbx/array_family/tiny_types.h>

namespace cctbx { namespace sgtbx { namespace reciprocal_space {

  /*! \brief Contiguous reciprocal space asymmetric units for
      the 230 reference settings.
   */
  /*! 12 contiguous reciprocal space asymmetric units (11 Laue
      classes, two settings for Laue class -3m) are
      tabulated. The tabulated asymmetric units are
      compatible with the asymmetric units of the CCP4
      package.
      <p>
      This implementation is based on work by
      <a href="http://www.ysbl.york.ac.uk/~cowtan/"
      >Kevin Cowtan</a>.
      <p>
      See also: class asu, class miller_index_generator
   */
  class reference_asu
  {
    public:
      //! Returns one of exactly 12 Laue group codes.
      /*! The labels of the possible return codes are:<br>
            -1, 2/m, mmm, 4/m, 4/mmm, -3, -31m, -3m1, 6/m, 6/mmm, m-3, m-3m
          <p>
          For Laue class -3m there are two possible orientations of the
          mirror plane with respect to the periodic lattice.
       */
      virtual matrix_group::code
      laue_group() const
      {
        throw CCTBX_INTERNAL_ERROR();
      }

      //! Test if given Miller index is in the tabulated asymmetric unit.
      virtual bool
      is_inside(miller::index<> const& /*h*/) const
      {
        throw CCTBX_INTERNAL_ERROR();
      }

      //! Tests Friedel pair.
      /*! Returns 1 if is_inside(h), -1 if is_inside(minus_h), 0 otherwise.
       */
      int
      which(miller::index<> const& h, miller::index<> const& minus_h) const
      {
        if      (is_inside(      h)) return  1;
        else if (is_inside(minus_h)) return -1;
        return 0;
      }

      /*! Returns 1 if is_inside(h), -1 if is_inside(-h), 0 otherwise.
       */
      int
      which(miller::index<> const& h) const
      {
        if (is_inside(h)) return 1;
        miller::index<> minus_h = -h; // work around Visual C++ 7.1 bug
        if (is_inside(minus_h)) return -1;
        return 0;
      }

      //! String representation of the tabluated asymmetric unit.
      /*! Example: "h>=k and k>=0 and (k>0 or l>=0)"
       */
      virtual const char*
      as_string() const
      {
        throw CCTBX_INTERNAL_ERROR();
      }

      //! "Cut parameters" for building Miller indices.
      /*! When building (or generating) a large list of Miller indices,
          it is useful to restrict the loop over all possible indices
          to 1/2, 1/4, or 1/8 of reciprocal space, if possible.
          <p>
          The cut parameters are used in the next() method
          of the class miller_index_generator. In general it
          should be much more convenient to use that higher-level
          class rather than hand-crafting a loop for building
          Miller indices.
          <p>
          Each element of cut_parameters() is either -1 or 0. A value
          of 0 indicates that the corresponding negative half-space
          can be omitted in the loop over possible indices.
          <p>
          Friedel symmetry is implied. If the Friedel mates
          are needed explicitly, they have to be added in a
          separate step. Note that the Friedel mate appears
          explicitly only for acentric reflections (use e.g.
          !space_group::is_centric() to determine which reflections
          are acentric).
       */
      virtual af::int3
      const& cut_parameters() const
      {
        throw CCTBX_INTERNAL_ERROR();
      }

      virtual ~reference_asu() {}
  };

  bool
  is_in_reference_asu_1b(miller::index<> const& h);

  const reference_asu*
  lookup_reference_asu(matrix_group::code const& group_code);

}}} // namespace cctbx::sgtbx::reciprocal_space

#endif // CCTBX_SGTBX_RECIPROCAL_SPACE_REFERENCE_ASU_H
