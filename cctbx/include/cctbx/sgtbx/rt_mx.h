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

#ifndef CCTBX_SGTBX_RT_MX_H
#define CCTBX_SGTBX_RT_MX_H

#include <cctbx/sgtbx/rot_mx.h>
#include <cctbx/sgtbx/parse_string.h>
#include <scitbx/array_family/tiny_types.h>
#include <cctbx/import_scitbx_af.h>
#include <boost/rational.hpp>

namespace cctbx { namespace sgtbx {

  //! Rotation-Translation matrix.
  /*! r() is the (3x3) rotation part or rotation matrix and t()
      the (3x1) translation part or translation vector of the symmetry
      operation.  For efficiency, the elements of the matrix and the
      vector are stored as integer numbers and denominators. The actual
      value of an element is obtained by dividing the integer number by
      the corresponding denominator:<br>
      r()[i] / r().den(), t()[i] / t().den().
   */
  class rt_mx
  {
    public:
      //! Identity rotation matrix and zero translation vector.
      /*! The denominators used are r_den and t_den for the rotation part
          and the translation part, respectively.
       */
      explicit
      rt_mx(int r_den=1, int t_den=sg_t_den)
      : r_(r_den), t_(t_den)
      {}

      //! Initialization with the given rotation matrix and translation vector.
      rt_mx(rot_mx const& r, tr_vec const& t)
      : r_(r), t_(t)
      {}

      //! Initialization with the given rotation matrix.
      /*! The translation part is initialized as the zero translation
          with the supplied denominator t_den.
       */
      explicit
      rt_mx(rot_mx const& r, int t_den=sg_t_den)
      : r_(r), t_(t_den)
      {}

      //! Initialization with the given translation vector.
      /*! The rotation part is initialized as the identity matrix
          with the supplied denominator r_den.
       */
      explicit
      rt_mx(tr_vec const& t, int r_den=1)
      : r_(r_den), t_(t)
      {}

      //! Initialization with a symbolic expression, e.g. "x+1/2,y,z".
      /*! Parsing will stop at any of the characters listed in stop_chars.<br>
          The denominators of the new rt_mx are r_den and t_den.<br>
          If an error occurs, an exception is thrown. parse_string can
          be investigated to locate the input character that triggered
          the error.
       */
      rt_mx(parse_string& str_xyz,
            const char* stop_chars="",
            int r_den=1, int t_den=sg_t_den);

      //! Initialize by parsing a symbolic expression, e.g. "x+1/2,y,z".
      /*! Identical to the constructor that takes a parse_string
          as the first argument. However, if an exception is thrown
          there is no way to locate the input character that triggered
          the error.
       */
      rt_mx(std::string const& str_xyz,
            const char* stop_chars="",
            int r_den=1, int t_den=sg_t_den);

      //! Initialize with floating point rotation and translation parts.
      rt_mx(scitbx::mat3<double> const& r,
            scitbx::vec3<double> const& t,
            int r_den=1, int t_den=sg_t_den);

      //! Access to rotation part.
      rot_mx const& r() const { return r_; };

      //! Access to rotation part.
      rot_mx&       r()       { return r_; };

      //! Access to translation part.
      tr_vec const& t() const { return t_; };

      //! Access to translation part.
      tr_vec&       t()       { return t_; };

      //! Tests equality.
      bool operator==(rt_mx const& rhs) const
      {
        return (this->r_ == rhs.r_) && (this->t_ == rhs.t_);
      }

      //! Tests inequality.
      bool operator!=(rt_mx const& rhs) const
      {
        return !((*this) == rhs);
      }

      //! Test if the matrix is valid.
      /*! A rt_mx is valid only if both the rotation denominator and the
          translation denominator are not zero.
       */
      bool is_valid() const
      {
        return r_.is_valid() && t_.is_valid();
      }

      //! New unit matrix.
      /*! The new matrix inherits r().den() and t().den().
       */
      rt_mx unit_mx() const
      {
        return rt_mx(r().den(), t().den());
      }

      //! Test if the matrix is the unit matrix.
      bool is_unit_mx() const
      {
        return r_.is_unit_mx() && t_.is_zero();
      }

      //! Conversion to a symbolic expression, e.g. "x+1/2,y,z".
      /*! Constants can be formatted as fractional or decimal numbers.<br>
          E.g. "x+1/2,y,z" or "x+0.5,y,z".<br>
          The translation component can appear last or first.<br>
          E.g. "x+1/2,y,z" or "1/2+x,y,z".<br>
          letters_xyz must contain three characters that are used to
          represent x, y, and z, respectively. Typically one only
          uses either letters_xyz = "xyz" or letters_xyz = "XYZ".
          Other letters could be used, but the resulting symbolic
          expression cannot be translated back with the constructors
          above.<br>
          separator is inserted after the first and the second term
          of the symbolic expression. Typical strings used are
          separator = "," and separator = ", ".
       */
      std::string
      as_xyz(bool decimal=false,
             bool t_first=false,
             const char* letters_xyz="xyz",
             const char* separator=",") const;

      //! Copy matrix elements to a plain array of int.
      /*! The 9 elements of the rotation part are copied to
          elements 0..8 of the plain array.
          The 3 elements of the translation part are copied to
          elements 9..11 of the plain array.
          Use r().den() and t().den() to access the denominators
          for conversion to floating point representations.
          <p>
          See also: as_double_array(), as_float_array()
       */
      af::tiny<int, 12>
      as_int_array() const;

      //! Copy matrix elements to a floating point array.
      /*! The 9 elements of the rotation part are divided by r().den()
          and copied to elements 0..8 of the array.
          The 3 elements of the translation part are divided by t().den()
          and copied to elements 9..11 of the array.
          <p>
          See also: as_int_array(), as_double_array(), as_float_array()
       */
      template <class FloatType>
      af::tiny<FloatType, 12>
      as_array(scitbx::type_holder<FloatType>) const;

      //! Copy matrix elements to a plain array of float.
      /*! See also: as_array(), as_int_array(), as_double_array()
       */
      af::tiny<float, 12>
      as_float_array() const
      {
        return as_array(scitbx::type_holder<float>());
      }

      //! Copy matrix elements to a plain array of double.
      /*! See also: as_array(), as_int_array(), as_float_array()
       */
      af::tiny<double, 12>
      as_double_array() const
      {
        return as_array(scitbx::type_holder<double>());
      }

      //! Copy with the denominators r_den and t_den.
      /*! An exception is thrown if the elements cannot be scaled to
          the new denominators.<br>
          r_den == 0 or t_den == 0 indicates that the corresponding old
          denominator is retained.
       */
      rt_mx new_denominators(int r_den, int t_den=0) const;

      //! Copy of this matrix with the denominators of other.
      /*! An exception is thrown if the elements cannot be scaled to
          the new denominators.
       */
      rt_mx new_denominators(rt_mx const& other) const
      {
        return new_denominators(other.r().den(), other.t().den());
      }

      //! Multiplies the elements and the denominators by the given factors.
      /*! if (factor_r == 0) factor_t = factor_r;
       */
      rt_mx scale(int factor_r, int factor_t=0) const;

      //! Applies modulus operation such that 0 <= x < t().den().
      /*! The operation is applied to the elements of the
          translation vector. The vector is modified in place.
       */
      void mod_positive_in_place() { t_ = t_.mod_positive(); }

      //! Applies modulus operation such that 0 <= x < t().den().
      /*! The operation is applied to the elements of the
          translation vector. A new instance of rt_mx is created.
       */
      rt_mx mod_positive() const { return rt_mx(r_, t_.mod_positive()); }

      /*! \brief Applies modulus operation such that
          -t().den()/2+1 < x <= t().den()/2.
       */
      /*! The operation is applied to the elements of the
          translation vector. The vector is modified in place.
       */
      void mod_short_in_place() { t_ = t_.mod_short(); }

      /*! \brief Applies modulus operation such that
          -t().den()/2+1 < x <= t().den()/2.
       */
      /*! The operation is applied to the elements of the
          translation vector. A new instance of rt_mx is created.
       */
      rt_mx mod_short() const { return rt_mx(r_, t_.mod_short()); }

      /*! \brief Tests if the vector v is perpendicular to the axis
          direction of the rotation part.
       */
      bool is_perpendicular(sg_vec3 const& v) const;

      //! Computes the intrinsic (screw or glide) part of the translation part.
      /*! Let N be the rotation-part type of r().
          Following the procedure given by Fischer & Koch
          (International Tables for Crystallography, Volume A, 1983,
          Chapter 11), the first step in the analysis of the
          translation part t() is the decomposition into the
          intrinsic part wi and the location part wl. For this,
          (r()|t())^n = (I|t) has to be computed, where the
          rotational order n = abs(N), except for N = -1
          and N = -3. For those two cases, n = -2*N. The
          intrinsic part is obtained as wi = (1/n)*t.<br>
          See also: translation_part_info
       */
      tr_vec t_intrinsic_part() const;

      //! Computes the location part given the intrinsic part wi.
      /*! wi is the result of intrinsic_part(). The location part
          is simply the difference wl = t() - wi.<br>
          See also: translation_part_info
       */
      tr_vec t_location_part(tr_vec const& wi) const;

      //! Computes the origin shift given the location part wl.
      /*! wl is the result of location_part(). The origin shift
          is obtained by solving the equation (r()|wl)*x = x
          for x (see
          <A HREF="http://journals.iucr.org/a/issues/1999/02/02/au0146/"
          ><I>Acta Cryst.</I> 1999, <B>A55</B>:383-395</A>).
          For rotation-part types N > 1 and N = -2, the combination
          with the axis direction of r() is a convenient
          description of all fixed points. For the other rotation-part
          types, the fixed point is unique.<br>
          See also: translation_part_info
       */
      tr_vec t_origin_shift(tr_vec const& wl) const;

      //! Efficient computation of (-I|inv_t) * (r|t).
      /*! I is the identidy matrix. inv_t is the translation part
          of a centre of inversion.
       */
      rt_mx pre_multiply_inv_t(tr_vec const& inv_t)
      {
        return rt_mx(-r_, -t_ + inv_t);
      }

      //! Computes the inverse matrix.
      /*! An exception is thrown if the matrix cannot be inverted
          or if the result cannot be scaled to the rotation part
          denominator or the translation part denominator.
       */
      rt_mx inverse() const;

      /*! \brief Similar to /=, but multiplies denominators instead
          of dividing elements.
       */
      void pseudo_divide(int rhs)
      {
        r_.den() *= rhs;
        t_.den() *= rhs;
      }

      //! Unary minus.
      rt_mx operator-() const { return rt_mx(-r_, -t_); }

      //! Addition operator.
      friend rt_mx operator+(rt_mx const& lhs, rt_mx const& rhs);

      //! += operator.
      rt_mx& operator+=(rt_mx const& rhs)
      {
        r_ += rhs.r_;
        t_ += rhs.t_;
        return *this;
      }

      //! Multiplication of homogeneous rt_mx with r().den()=1.
      /*! The rotation denominator of both matrices must be equal to 1.
          The translation denominators of the two matrices must be equal.
          <p>
          The denominators of the result are equal to the denominators
          of the operands.
          <p>
          operator*() is faster than multiply().
       */
      friend rt_mx operator*(rt_mx const& lhs, rt_mx const& rhs);

      //! Addition of translation vector to translation part.
      friend rt_mx operator+(rt_mx const& lhs, tr_vec const& rhs);

      /*! \brief Refines gridding such that each grid point is
          mapped onto another grid point by the symmetry operation.
       */
      template <typename GridTupleType>
      GridTupleType
      refine_gridding(GridTupleType const& grid) const;

      //! Reduces denominators by cancellation.
      /*! The rotation part denominator and the elements of the rotation
          part are divided by their greatest common denominator.
          The same procedure is applied to the translation part.
       */
      rt_mx cancel() const;

      //! Computes the inverse matrix.
      /*! An exception is thrown if the matrix cannot be inverted.
       */
      rt_mx inverse_cancel() const;

      //! Multiplication with cancellation for general rt_mx.
      /*! Similar to opertor*(). However, the operands may have any
          rotation denominator or translation denominator.
          <p>
          The denominators of the result are made as small as possible.
          <p>
          See also: cancel()
       */
      rt_mx multiply(rt_mx const& rhs) const;

    private:
      rot_mx r_;
      tr_vec t_;
  };

  template <class FloatType>
  af::tiny<FloatType, 12>
  rt_mx::as_array(scitbx::type_holder<FloatType>) const
  {
    af::tiny<FloatType, 12> result;
    for(std::size_t i=0;i<9;i++) {
      result[i    ] = FloatType(r_[i]) / FloatType(r().den());
    }
    for(std::size_t i=0;i<3;i++) {
      result[i + 9] = FloatType(t_[i]) / FloatType(t().den());
    }
    return result;
  }

  template <typename IntType>
  IntType
  norm_denominator(IntType numerator, IntType denominator)
  {
    return denominator / boost::gcd(numerator, denominator);
  }

  template <typename GridTupleType>
  GridTupleType
  rt_mx::refine_gridding(GridTupleType const& grid) const
  {
    GridTupleType result;
    for(std::size_t ir=0;ir<3;ir++) {
      result[ir] = boost::lcm(
        grid[ir], norm_denominator(t_[ir], t_.den()));
      for(std::size_t ic=0;ic<3;ic++) {
        result[ir] = boost::lcm(
          result[ir], norm_denominator(r_(ir, ic), grid[ic]));
      }
    }
    return result;
  }

  //! Multiplication of rt_mx with a vector of double.
  /*! Python: __mul__
   */
  inline
  scitbx::vec3<double>
  operator*(rt_mx const& lhs, af::tiny_plain<double, 3> const& rhs)
  {
    rot_mx const& r = lhs.r();
    tr_vec const& t = lhs.t();
    double rd = r.den();
    double td = t.den();
    scitbx::vec3<double> result;
    for(std::size_t i=0;i<3;i++) {
      result[i] = (  r(i,0) * rhs[0]
                   + r(i,1) * rhs[1]
                   + r(i,2) * rhs[2]) / rd + t[i] / td;
    }
    return result;
  }

  //! Analysis of the translation part of a rotation-translation matrix.
  /*! Grouping of the results of rt_mx::t_intrinsic_part(),
      rt_mx::t_location_part() and rt_mx::t_origin_shift().
   */
  class translation_part_info
  {
    public:
      //! Default constructor. Some data members are not initialized!
      translation_part_info()
      {}

      /*! \brief Determination of intrinsic_part(), location_part()
          and origin_shift().
       */
      translation_part_info(rt_mx const& s);

      //! Intrinsic (srew or glide) part.
      /*! See rt_mx::t_intrinsic_part()
       */
      tr_vec const& intrinsic_part() const { return ip_; }

      //! Location part.
      /*! See rt_mx::t_location_part()
       */
      tr_vec const& location_part() const { return lp_; }

      //! Origin shift.
      /*! See rt_mx::t_origin_shift()
       */
      tr_vec const& origin_shift() const { return os_; }

    private:
      tr_vec ip_;
      tr_vec lp_;
      tr_vec os_;
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_RT_MX_H
