#ifndef CCTBX_SGTBX_CHANGE_OF_BASIS_OP_H
#define CCTBX_SGTBX_CHANGE_OF_BASIS_OP_H

#include <cctbx/sgtbx/rt_mx.h>
#include <cctbx/sgtbx/utils.h>
#include <cctbx/uctbx.h>

namespace cctbx { namespace sgtbx {

  //! Change-of-basis (transformation) operator.
  /*! For ease of use, a change-of-basis matrix c() and its inverse c_inv()
      are grouped by this class.
   */
  class change_of_basis_op
  {
    public:
      //! Initializes the change-of-basis operator with c and c_inv.
      /*! The input matrices are NOT checked for consistency.
       */
      change_of_basis_op(rt_mx const& c, rt_mx const& c_inv)
      :
        c_(c), c_inv_(c_inv)
      {}

      //! Initializes the change-of-basis operator with c.
      /*! The inverse matrix c_inv() is computed by inverting c.
          An exception is thrown if c is not invertible.
       */
      explicit
      change_of_basis_op(rt_mx const& c)
      :
        c_(c), c_inv_(c.inverse())
      {}

      /*! \brief Initializes the change-of-basis operator with
          matrix given as xyz, hkl, or abc string.
       */
      /*! An exception is thrown if the given matrix is not invertible.
          <p>
          See also: constructor of class rt_mx
          <p>
          Not available in Python.
       */
      change_of_basis_op(
        parse_string& symbol,
        const char* stop_chars="",
        int r_den=cb_r_den,
        int t_den=cb_t_den);

      /*! \brief Initializes the change-of-basis operator with
          matrix given as xyz, hkl, or abc string.
       */
      /*! An exception is thrown if the given matrix is not invertible.
          <p>
          See also: constructor of class rt_mx
       */
      change_of_basis_op(
        std::string const& symbol,
        const char* stop_chars="",
        int r_den=cb_r_den,
        int t_den=cb_t_den);

      //! Initializes the change-of-basis operator with unit matrices.
      /*! The unit matrices are initialized with the rotation part
          denominator r_den and the translation part denominator t_den.
       */
      explicit
      change_of_basis_op(int r_den=cb_r_den, int t_den=cb_t_den)
      :
        c_(r_den, t_den), c_inv_(r_den, t_den)
      {}

      //! Tests if the change-of-basis operator is valid.
      /*! A change_of_basis_op is valid only if the rotation part
          denominator and the translation part denominator of both
          c() and c_inv() are not zero.
       */
      bool
      is_valid() const
      {
        return c_.is_valid() && c_inv_.is_valid();
      }

      //! Returns a new change-of-basis operator with unit matrices.
      /*! The new matrices inherit the rotation and translation part
          denominators.
       */
      change_of_basis_op
      identity_op() const
      {
        return change_of_basis_op(c_.unit_mx(), c_inv_.unit_mx());
      }

      //! Tests if the change-of-basis operator is the identity.
      bool
      is_identity_op() const
      {
        return c_.is_unit_mx() && c_inv_.is_unit_mx();
      }

      //! Returns a new copy with the denominators r_den and t_den.
      /*! An exception is thrown if the elements cannot be scaled to
          the new denominators.<br>
          r_den or t_den == 0 indicates that the corresponding old
          denominator is retained.
       */
      change_of_basis_op
      new_denominators(int r_den, int t_den) const
      {
        return change_of_basis_op(c_.new_denominators(r_den, t_den),
                                  c_inv_.new_denominators(r_den, t_den));
      }

      /*! Returns a new copy of the operator, but with the denominators
          of other.
       */
      /*! An exception is thrown if the elements cannot be scaled to
          the new denominators.
       */
      change_of_basis_op
      new_denominators(change_of_basis_op const& other) const
      {
        return change_of_basis_op(c_.new_denominators(other.c()),
                                  c_inv_.new_denominators(other.c_inv()));
      }

      //! Returns the change-of-basis matrix.
      rt_mx const&
      c() const { return c_; }

      //! Returns the inverse of the change-of-basis matrix.
      rt_mx const&
      c_inv() const { return c_inv_; }

      //! Returns c() for inv == false, and c_inv() for inv == true.
      rt_mx const&
      select(bool inv) const
      {
        if (inv) return c_inv_;
        return c_;
      }

      //! Returns a new copy with c() and c_inv() swapped.
      change_of_basis_op
      inverse() const
      {
        return change_of_basis_op(c_inv_, c_);
      }

      //! Applies modulus operation such that 0 <= x < t_den.
      /*! The operation is applied to the elements of the
          translation vectors of c() and c_inv(). The vectors are
          modified in place.
       */
      void
      mod_positive_in_place()
      {
        c_.mod_positive_in_place();
        c_inv_.mod_positive_in_place();
      }

      //! Applies modulus operation such that -t_den/2+1 < x <= t_den/2.
      /*! The operation is applied to the elements of the
          translation vectors of c() and c_inv(). The vectors are
          modified in place.
       */
      void
      mod_short_in_place()
      {
        c_.mod_short_in_place();
        c_inv_.mod_short_in_place();
      }

      //! Applies modulus operation such that -t_den/2+1 < x <= t_den/2.
      /*! The operation is applied to the elements of the
          translation vectors of c() and c_inv().
       */
      change_of_basis_op
      mod_short() const
      {
        return change_of_basis_op(c_.mod_short(), c_inv_.mod_short());
      }

      //! c().r() * r * c_inv().r(), for r with rotation part denominator 1.
      /*! The rotation part denominator of the result is 1.
          <p>
          Not available in Python.
       */
      rot_mx
      operator()(rot_mx const& r) const;

      //! c() * s * c_inv(), for s with rotation part denominator 1.
      /*! Similar to apply(s), but faster. The translation denominator
          of the result is equal to the translation denominator of s.
          <p>
          Not available in Python.
       */
      rt_mx
      operator()(rt_mx const& s) const;

      //! c() * s * c_inv(), for s with any rotation part denominator.
      /*! Similar to opertor()(). s may have any rotation part denominator
          or translation part denominator. The denominators of the result
          are made as small as possible.
          <p>
          See also: rt_mx::multiply(), rt_mx::cancel()
       */
      rt_mx
      apply(rt_mx const& s) const;

      //! Deprecated. Use cctbx::uctbx::unit_cell::change_basis(cb_op) instead.
      uctbx::unit_cell
      apply(uctbx::unit_cell const& unit_cell) const
      {
        return unit_cell.change_basis(c_inv().r());
      }

      //! Transforms a Miller index.
      /*! result = miller_index * c_inv()
       */
      miller::index<>
      apply(miller::index<> const& miller_index) const;

      //! Transforms an array of Miller indices.
      /*! result = miller_indices * c_inv()
       */
      af::shared<miller::index<> >
      apply(af::const_ref<miller::index<> > const& miller_indices) const;

      //! Flags Miller indices for which apply() leads to non-integral indices.
      af::shared<std::size_t>
      apply_results_in_non_integral_indices(
        af::const_ref<miller::index<> > const& miller_indices) const;

      //! c() * (rt_mx(rot_mx(sign_identity), t)) * c_inv()
      /*! Not available in Python.
       */
      tr_vec
      operator()(tr_vec const& t, int sign_identity) const;

      //! Transform fractional coordinates: c() * site_frac
      template <class FloatType>
      fractional<FloatType>
      operator()(fractional<FloatType> const& site_frac) const
      {
        return c_ * site_frac;
      }

      //! c() = other.c() * c(); c_inv() = c_inv() * other.c_inv();
      void
      update(change_of_basis_op const& other)
      {
        c_ = (other.c() * c_).new_denominators(other.c());
        c_inv_ = (c_inv_ * other.c_inv()).new_denominators(other.c_inv());
      }

      //! c() = (I|shift) * c(); c_inv() = c_inv() * (I|-shift);
      /*! Not available in Python.
       */
      void
      update(tr_vec const& shift)
      {
        // (I|S)*(R|T) = (R|T+S)
        c_ = rt_mx(c_.r(), c_.t() + shift);
        // (R|T)*(I|-S) = (R|T-R*S)
        c_inv_ = rt_mx(
          c_inv_.r(),
          c_inv_.t() - (c_inv_.r() * shift).new_denominator(c_inv_.t().den()));
      }

      //! Multiplication of change-of-basis operators.
      change_of_basis_op
      operator*(change_of_basis_op const& rhs)
      {
        return change_of_basis_op(
          (c() * rhs.c()).new_denominators(c()),
          (rhs.c_inv() * c_inv()).new_denominators(c_inv()));
      }

      //! Returns c() in xyz format.
      /*! See also: rt_mx::as_xyz()
       */
      std::string
      as_xyz(bool decimal=false,
             bool t_first=false,
             const char* symbol_letters="xyz",
             const char* separator=",") const
      {
        return c_.as_xyz(decimal, t_first, symbol_letters, separator);
      }

      //! Returns transpose of c_inv() in hkl format.
      /*! c_inv.t() must be zero. An exception is thrown otherwise.
          <p>
          See also: rt_mx::as_xyz(), rot_mx::as_hkl()
       */
      std::string
      as_hkl(
        bool decimal=false,
        const char* letters_hkl="hkl",
        const char* separator=",") const
      {
        CCTBX_ASSERT(c_inv_.t().is_zero());
        return c_inv_.r().as_hkl(decimal, letters_hkl, separator);
      }

      //! Returns transpose of c_inv() in abc format.
      /*! See also: rt_mx::as_xyz()
       */
      std::string
      as_abc(bool decimal=false,
             bool t_first=false,
             const char* letters_abc="abc",
             const char* separator=",") const
      {
        return rt_mx(c_inv_.r().transpose(), c_inv_.t())
          .as_xyz(decimal, t_first, letters_abc, separator);
      }

      //! Returns abc format, or xyz if it is shorter than abc.
      /*! See also: as_abc(), as_xyz()
       */
      std::string
      symbol() const
      {
        std::string result = as_abc();
        std::string alt = as_xyz();
        if (result.size() > alt.size()) return alt;
        return result;
      }

    private:
      rt_mx c_;
      rt_mx c_inv_;
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_CHANGE_OF_BASIS_OP_H
