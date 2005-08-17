#ifndef CCTBX_SGTBX_RECIPROCAL_SPACE_ASU_H
#define CCTBX_SGTBX_RECIPROCAL_SPACE_ASU_H

#include <cctbx/sgtbx/space_group_type.h>
#include <cctbx/sgtbx/reciprocal_space_reference_asu.h>
#include <cctbx/sgtbx/miller_ops.h>

namespace cctbx { namespace sgtbx { namespace reciprocal_space {

  //! Access to general contiguous reciprocal space asymmetric units.
  /*! class reference_reciprocal_space_asu implements 12 contiguous
      reciprocal space asymmetric units that cover the 230
      reference settings. The algorithm for the determination of the
      space group type (class space_group_type) is used to derive a
      change-of-basis matrix for the transformation of the tabulated
      asymmetric units. In this way a contiguous asymmetric unit is
      available for any arbitrary setting.
      <p>
      Python binding: reciprocal_space_asu
   */
  class asu
  {
    public:
      //! Default constructor.
      /*! Default-constructed instances will throw exceptions if
          some of the member functions are used.
       */
      asu()
      :
        is_reference_(true)
      {}

      //! Initialization.
      /*! Based on the space group number (space_group_type::number()),
          the Laue group class is derived which is in turn used
          to select the corresponding tabulated
          reference_asu.
       */
      asu(sgtbx::space_group_type const& space_group_type);

      //! Access to the selected tabulated reference_asu.
      /*! Not available in Python.
       */
      const reference_asu*
      reference() const { return reference_; }

      /*! \brief Access to the string representation of the selected
          tabulated reference_asu.
       */
      /*! Equivalent to: reference().as_string()
       */
      const char*
      reference_as_string() const { return reference_->as_string(); }

      //! Access to the change-of-basis operator.
      /*! This operator is a copy of space_group_type.cb_op() as passed to
          the constructor.
       */
      change_of_basis_op const&
      cb_op() const { return cb_op_; }

      //! Tests if the given asymmetric unit is one of the tabulated units.
      /*! This test is equivalent to the test cb_op().is_identity_op().
          If this is the case, some optimizations are used.
       */
      bool
      is_reference() const { return is_reference_; }

      //! Tests if the given Miller index is in the asymmetric unit.
      /*! The change-of-basis matrix is used to transform
          the Miller index (miller_index * cb_op().c_inv()). It is then
          tested if the result is in the tabulated reference
          asymmetric unit.
       */
      bool
      is_inside(miller::index<> const& miller_index) const
      {
        if (is_reference_) return reference_->is_inside(miller_index);
        return reference_->is_inside(miller_index * cb_op_.c_inv().r());
      }

      //! Tests Friedel pair.
      /*! Returns 1 if is_inside(miller_index),
          -1 if is_inside(minus_miller_index),
          0 otherwise.
          <p>
          Not available in Python.
       */
      int
      which(
        miller::index<> const& miller_index,
        miller::index<> const& minus_miller_index) const
      {
        if (is_reference_) {
          return reference_->which(miller_index, minus_miller_index);
        }
        return reference_->which(miller_index * cb_op_.c_inv().r());
      }

      //! Tests Friedel pair.
      /*! Returns 1 if is_inside(miller_index),
          -1 if is_inside(-miller_index), 0 otherwise.
       */
      int
      which(miller::index<> const& miller_index) const
      {
        if (is_reference_) return reference_->which(miller_index);
        return reference_->which(miller_index * cb_op_.c_inv().r());
      }

    private:
      change_of_basis_op cb_op_;
      bool is_reference_;
      const reference_asu* reference_;
  };

}}} // namespace cctbx::sgtbx::reciprocal_space

#endif // CCTBX_SGTBX_RECIPROCAL_SPACE_ASU_H
