#ifndef CCTBX_SGTBX_SITE_SYMMETRY_H
#define CCTBX_SGTBX_SITE_SYMMETRY_H

#include <cctbx/sgtbx/site_constraints.h>
#include <cctbx/sgtbx/tensor_rank_2.h>
#include <cctbx/uctbx.h>

namespace cctbx { namespace sgtbx {

  class rt_point_group
  {
    public:
      rt_point_group() : is_valid_(true) {}

      rt_point_group(sgtbx::space_group const& sg, rt_mx const& projection);

      bool
      is_valid() const { return is_valid_; }

      af::shared<rt_mx> const&
      matrices() const { return matrices_; }

      void
      reset(rt_mx const& s);

      void
      add(rt_mx const& s);

      void
      expand(rt_mx const& s);

      bool
      try_expand(rt_mx const& s);

      rt_mx
      accumulate() const;

      sgtbx::space_group
      space_group() const;

      matrix_group::code
      type() const;

    protected:
      bool is_valid_;
      af::shared<rt_mx> matrices_;
  };


  //! Base class for site_symmetry, holding the essential results.
  class site_symmetry_ops
  {
    public:
      //! Default constructor. Some data members are not initialized!
      site_symmetry_ops()
      :
        have_site_constraints_(false),
        have_adp_constraints_(false),
        have_cartesian_adp_constraints_(false)
      {}

      /*! \brief Number of distinct symmetrically equivalent positions
          of site_symmetry::exact_site().
       */
      int
      multiplicity() const { return multiplicity_; }

      //! Special position operation.
      /*! This operation is used to compute site_symmetry::exact_site() from
          site_symmetry::original_site(). It satisfies the following two
          conditions:
          <ul>
          <li>special_op() * original_site() = exact_site()
          <li>special_op() * exact_site() = exact_site()
          </ul>
       */
      rt_mx const&
      special_op() const { return special_op_; }

      //! Symmetry matrices of site symmetry point group.
      af::shared<rt_mx> const&
      matrices() const { return matrices_; }

      //! Shorthand for: matrices().size()
      std::size_t
      n_matrices() const { return matrices_.size(); }

      //! Tests if the site symmetry is point group 1.
      bool
      is_point_group_1() const { return n_matrices() == 1; }

      //! True if special_op() and all matrices() are identical.
      bool
      operator==(site_symmetry_ops const& other) const
      {
        if (other.special_op_ != special_op_) return false;
        return other.matrices_.const_ref().all_eq(matrices_.const_ref());
      }

      /*! \brief Tests if the given anisotropic displacement parameters
          u_star are compatible with the site symmetry.
       */
      /*! The condition
          <p>
          r * u_star * r.transposed() = u_star
          <p>
          is evaluated for all rotation parts r of the site symmetry.
       */
      template <class FloatType>
      bool
      is_compatible_u_star(
        scitbx::sym_mat3<FloatType> const& u_star,
        FloatType tolerance=1e-6) const;

      /*! \brief Averages symmetrically equivalent u_star tensors to
          obtain a tensor that satisfies the symmetry constraints.
       */
      /*! The averaged tensor is equivalent to beta_inv
          of Giacovazzo, Fundamentals of Crystallography 1992,
          p. 189.
       */
      template <class FloatType>
      scitbx::sym_mat3<FloatType>
      average_u_star(scitbx::sym_mat3<FloatType> const& u_star) const
      {
        return average_tensor(matrices_.const_ref(), u_star, true);
      }

      /*! \brief Construct a new site_symmetry_ops instance with the
          denominators (special_op().den(), matrices()[0].den()) of this.
       */
      /*! Not available in Python.
       */
      site_symmetry_ops
      make_point_group_1() const
      {
        site_symmetry_ops result(multiplicity_*matrices_.size(), 1, 1);
        result.matrices_.push_back(matrices_[0]);
        return result;
      }

      //! Apply change-of-basis operator.
      site_symmetry_ops
      change_basis(change_of_basis_op const& cb_op) const
      {
        site_symmetry_ops result;
        rat new_multiplicity = cb_op.c_inv().r().determinant() * multiplicity_;
        CCTBX_ASSERT(new_multiplicity.denominator() == 1);
        result.multiplicity_ = scitbx::fn::absolute(
          new_multiplicity.numerator());
        result.special_op_ = cb_op.apply(special_op_);
        af::const_ref<rt_mx> m = matrices_.const_ref();
        result.matrices_.reserve(m.size());
        result.matrices_.push_back(m[0]); // identity matrix
        for(std::size_t i=1;i<m.size();i++) {
          result.matrices_.push_back(cb_op(m[i]));
        }
        return result;
      }

      //! Return reference to cached site_constraints instance.
      sgtbx::site_constraints<> const&
      site_constraints() const
      {
        if (!have_site_constraints_) {
          site_constraints_ = sgtbx::site_constraints<>(matrices_.const_ref());
          have_site_constraints_ = true;
        }
        return site_constraints_;
      }

      //! Return reference to cached tensor_rank_2::constraints instance.
      tensor_rank_2::constraints<> const&
      adp_constraints() const
      {
        if (!have_adp_constraints_) {
          adp_constraints_ = tensor_rank_2::constraints<>(
            matrices_.const_ref(), 1, true);
          have_adp_constraints_ = true;
        }
        return adp_constraints_;
      }

      //! Return reference to cached constraints on u_cart
      /*! The unit_cell is only used the first time this member function is
          called. On subsequent calls, the cached value is used and unit_cell is
          ignored unless recompute is true, in which case it runs like the first
          time.
      */
      tensor_rank_2::cartesian_constraints<> const&
      cartesian_adp_constraints(uctbx::unit_cell const& unit_cell,
                                bool recompute=false) const
      {
        if (recompute || !have_cartesian_adp_constraints_) {
          cartesian_adp_constraints_ = tensor_rank_2::cartesian_constraints<>(
            unit_cell,
            matrices_.const_ref());
          have_cartesian_adp_constraints_ = true;
        }
        return cartesian_adp_constraints_;
      }

      //! Support for Python's pickle facility. Do not use for other purposes.
      site_symmetry_ops(
        int multiplicity,
        rt_mx const& special_op,
        af::shared<rt_mx> const& matrices)
      :
        multiplicity_(multiplicity),
        special_op_(special_op),
        matrices_(matrices),
        have_site_constraints_(false),
        have_adp_constraints_(false),
        have_cartesian_adp_constraints_(false)
      {}

    protected:
      int multiplicity_;
      rt_mx special_op_;
      af::shared<rt_mx> matrices_;
      mutable bool have_site_constraints_;
      mutable sgtbx::site_constraints<> site_constraints_;
      mutable bool have_adp_constraints_;
      mutable tensor_rank_2::constraints<> adp_constraints_;
      mutable bool have_cartesian_adp_constraints_;
      mutable tensor_rank_2::cartesian_constraints<> cartesian_adp_constraints_;

      site_symmetry_ops(int multiplicity, int r_den, int t_den)
      :
        multiplicity_(multiplicity),
        special_op_(r_den, t_den),
        have_site_constraints_(false),
        have_adp_constraints_(false),
        have_cartesian_adp_constraints_(false)
      {}
  };

  template <class FloatType>
  bool
  site_symmetry_ops::
  is_compatible_u_star(
    scitbx::sym_mat3<FloatType> const& u_star,
    FloatType tolerance) const
  {
    FloatType scaled_tolerance = af::max_absolute(u_star) * tolerance;
    for (std::size_t i=0;i<matrices_.size();i++) {
      scitbx::mat3<FloatType>
        r = matrices_[i].r()
              .as_floating_point(scitbx::type_holder<FloatType>());
      if (!u_star.const_ref().all_approx_equal(
           u_star.tensor_transform(r).const_ref(), scaled_tolerance)) {
        return false;
      }
    }
    return true;
  }

  //! Numerically robust algorithm for the determination of site-symmetries.
  class site_symmetry : public site_symmetry_ops
  {
    public:
      //! Default constructor. Some data members are not initialized!
      site_symmetry() {}

      //! Determines the site symmetry of original_postion.
      /*! The first step in the determination of the site-symmetry
          is the computation of the oriented and translated
          point_group() for the original_postion using the given parameters.
          If the site symmetry is higher than 1, the exact_site()
          of the special position is determined based on point_group().
          If the site symmetry is 1, exact_site() is equivalent
          to original_site().

          If the distance between symmetrically equivalent sites is
          less than or equal to min_distance_sym_equiv,
          original_site() is considered to be at a special position.

          If, after moving the site to exact_site(), the distance
          between symmetrically equivalent sites is still less than or
          equal to min_distance_sym_equiv, an exception is raised if
          assert_min_distance_sym_equiv == true. This condition is
          usually the consequence of an input error: the unit cell is
          too small relative to min_distance_sym_equiv. The exception
          is not raised if assert_min_distance_sym_equiv == false.
          check_min_distance_sym_equiv() can be used to query the
          status.

          See also: class wyckoff::table
       */
      site_symmetry(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        fractional<> const& original_site,
        double min_distance_sym_equiv=0.5,
        bool assert_min_distance_sym_equiv=true);

      //! Unit cell used in the computation of the site symmetry.
      uctbx::unit_cell const&
      unit_cell() const { return unit_cell_; }

      //! Space group used in the computation of the site symmetry.
      sgtbx::space_group const&
      space_group() const { return space_group_; }

      //! Retrieves the original coordinates.
      fractional<> const&
      original_site() const { return original_site_; }

      //! Retrieves the parameter passed to the constructor.
      double
      min_distance_sym_equiv_sq() const
      {
        return min_distance_sym_equiv_sq_;
      }

      //! Retrieves the parameter passed to the constructor.
      double
      min_distance_sym_equiv() const
      {
        return std::sqrt(min_distance_sym_equiv_sq_);
      }

      //! Exact location of the special position.
      fractional<> const&
      exact_site() const { return exact_site_; }

      //! Distance^2 between original_site() and exact_site().
      /*! Not available in Python.
       */
      double
      distance_moved_sq() const
      {
        return unit_cell_.distance_sq(exact_site_, original_site_);
      }

      //! Distance between original_site() and exact_site().
      double
      distance_moved() const
      {
        return std::sqrt(distance_moved_sq());
      }

      /*! \brief Shortest distance^2 between symmetrically equivalent
          positions of exact_site().
       */
      /*! Not available in Python.
       */
      double
      shortest_distance_sq() const { return shortest_distance_sq_; }

      /*! \brief Shortest distance between symmetrically equivalent
          positions of exact_site().
       */
      double
      shortest_distance() const { return std::sqrt(shortest_distance_sq_); }

      //! Tests if shortest_distance() > min_distance_sym_equiv().
      /*! As a special case, always true if min_distance_sym_equiv() is zero.
       */
      bool
      check_min_distance_sym_equiv() const
      {
        if (min_distance_sym_equiv_sq_ == 0) return true;
        return shortest_distance_sq_ > min_distance_sym_equiv_sq_;
      }

      //! Oriented and translated site-symmetry point group.
      /*! Not available in Python.
       */
      rt_point_group const&
      point_group() const { return point_group_; }

      //! Point group type of the site-symmetry.
      matrix_group::code
      point_group_type() const { return point_group_.type(); }

      //! Computes the unique symmetry operations for a special position.
      /*! special_op() is multiplied with all symmetry operations
          of space_group(). The unique results are returned.
          <p>
          See also: sgtbx::space_group::unique()
       */
      af::shared<rt_mx>
      unique_ops();

    protected:
      uctbx::unit_cell unit_cell_;
      sgtbx::space_group space_group_;
      fractional<> original_site_;
      double min_distance_sym_equiv_sq_;

      double shortest_distance_sq_;
      rt_point_group point_group_;
      fractional<> exact_site_;

      void build_special_op();
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_SITE_SYMMETRY_H
