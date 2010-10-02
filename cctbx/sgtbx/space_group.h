#ifndef CCTBX_SGTBX_SPACE_GROUP_H
#define CCTBX_SGTBX_SPACE_GROUP_H

#include <cctbx/sgtbx/tr_group.h>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/group_codes.h>
#include <cctbx/sgtbx/phase_info.h>
#include <map>

namespace cctbx { namespace sgtbx {

  //! Specific exception to indicate failure of group multiplication.
  class error_non_crystallographic_rotation_matrix_encountered : public error
  {
    public:
      error_non_crystallographic_rotation_matrix_encountered()
      :
        error("Non-crystallographic rotation matrix encountered.")
      {}
  };

  class space_group_type; // forward declaration

  //! Space group operations.
  /*! This class is composed of the following main components:<br>
      <ol>
      <li>A list of n_ltr() centring vectors.
      <li>An integer f_inv() and the translation part inv_t() of an inversion
          operation, if present. f_inv() = 1 if no inversion operation
          exists, f_inv() = 2 otherwise.
      <li>A list of n_smx() = order_z()/(n_ltr()*f_inv()) representative
          symmetry matrices.
      </ol>
   */
  class space_group
  {
    public:
      //! Type of internal array of Seitz matrices.
      typedef af::small<rt_mx, n_max_repr_rot_mx> smx_array_type;

      //! Default constructor. Symmetry is set to P1.
      /*! With no_expand = true, group multiplication will not be
          carried out. This option is for internal use only.
          <p>
          Not available in Python.
       */
      explicit
      space_group( bool no_expand=false, int t_den=sg_t_den)
      :
        no_expand_(no_expand)
      {
        reset(t_den);
      }

      //! Initialization with symmetry encoded by a Hall symbol.
      /*! The Hall symbol parser can be instructed to be more
          or less pedantic.<br>
          With no_centring_type_symbol = true, the parser
          can also be used to interpret matrix symbols only. This
          option is for internal use only.<br>
          An exception is thrown if an error occurs. parse_string can
          be investigated to locate the input character that triggered
          the error.
       */
      explicit
      space_group(
        parse_string& hall_symbol,
        bool pedantic=false,
        bool no_centring_type_symbol=false,
        bool no_expand=false,
        int t_den=sg_t_den);

      //! Initialization with symmetry encoded by a Hall symbol.
      /*! Identical to the constructor that takes a parse_string
          as the first argument. However, if an exception is thrown
          it is not possible to locate the input character that
          triggered the error.
       */
      explicit
      space_group(
        std::string const& hall_symbol,
        bool pedantic=false,
        bool no_centring_type_symbol=false,
        bool no_expand=false,
        int t_den=sg_t_den);

      /*! Identical to the constructor that takes a parse_string
          as the first argument. However, if an exception is thrown
          it is not possible to locate the input character that
          triggered the error.
       */
      explicit
      space_group(
        const char* hall_symbol,
        bool pedantic=false,
        bool no_centring_type_symbol=false,
        bool no_expand=false,
        int t_den=sg_t_den);

      /*! \brief Initialization with the Hall symbol of the
          space_group_symbols object.
       */
      explicit
      space_group(
        sgtbx::space_group_symbols const& space_group_symbols,
        int t_den=sg_t_den);

      //! Resets the symmetry to P1.
      void
      reset(int t_den=sg_t_den);

      //! Adds a lattice translation vector to the space group.
      /*! Group multiplication is automatically performed
          (unless no_expand = true).
       */
      space_group&
      expand_ltr(tr_vec const& new_t);

      //! Adds a centre of inversion to the space group.
      /*! The matrix added is (-I, new_inv_t). This is, the rotation part
          -I is implied, and the translation part is new_inv_t.<br>
          Group multiplication is automatically performed
          (unless no_expand = true).
       */
      space_group&
      expand_inv(tr_vec const& new_inv_t);

      //! Adds a Seitz matrix to the space group.
      /*! Group multiplication is automatically performed
          (unless no_expand = true).
       */
      space_group&
      expand_smx(rt_mx const& new_smx);

      //! Adds a vector of Seitz matrices to the space group.
      /*! See also: overload for individual new_smx.
          <p>
          Not available in Python.
       */
      space_group&
      expand_smx(af::const_ref<rt_mx> const& new_smx)
      {
        for(std::size_t i=0;i<new_smx.size();i++) {
          expand_smx(new_smx[i]);
        }
        return *this;
      }

      //! Shorthand for: expand_smx(rt_mx(smx_symbol))
      space_group&
      expand_smx(std::string const& smx_symbol)
      {
        expand_smx(rt_mx(smx_symbol));
        return *this;
      }

      //! Adds lattice translation vectors to the space group.
      /*! The lattice translation vectors corresponding to the
          conventional centring type symbol are determined from
          a lookup table and added to the space group.<br>
          Group multiplication is automatically performed
          (unless no_expand = true).
       */
      std::size_t
      expand_conventional_centring_type(char symbol);

      //! Parses a Hall Symbol and adds the encoded symmetry matrices.
      /*! Similar to the constructors that take a Hall symbol as the
          first argument. However, the new symmetry matrices are added
          to the existing ones.
       */
      std::size_t
      parse_hall_symbol(
        parse_string& hall_symbol,
        bool pedantic=false,
        bool no_centring_type_symbol=false);

      //! Parses a Hall Symbol and adds the encoded symmetry matrices.
      /*! The change-of-basis operator is parsed but not applied.
          For internal use only.
          <p>
          Not available in Python.
       */
      std::size_t
      parse_hall_symbol_cb_op(
        parse_string& hall_symbol,
        change_of_basis_op& cb_op,
        bool pedantic=false,
        bool no_centring_type_symbol=false);

      //! Applies a change-of-basis operator to the symmetry operations.
      /*! The transformed space group is returned as a new object.<br>
          An exception is thrown if the change-of-basis operator is
          invalid or if the new symmetry matrices can not be
          represented as integer matrices with denominators in use.
       */
      space_group
      change_basis(change_of_basis_op const& cb_op) const;

      //! The change of origin which moves the inversion centre to the origin
      /*! The basis axes are therefore not altered */
      change_of_basis_op
      change_of_origin_realising_origin_centricity() const;

      //! Rotation part denominator of the Seitz matrices.
      int r_den() const { return smx_[0].r().den(); }
      //! Translation part denominator of the Seitz matrices.
      int
      t_den() const { return smx_[0].t().den(); }

      //! Order of the point group = f_inv() * n_smx().
      std::size_t
      order_p() const { return f_inv_ * n_smx(); }

      //! Order of the space group = n_ltr() * f_inv() * n_smx().
      /*! Python: __len__
       */
      std::size_t
      order_z() const { return ltr_.size() * f_inv_ * n_smx(); }

      //! Alias for order_z().
      std::size_t
      n_equivalent_positions() const { return order_z(); }

      //! Number of lattice translations.
      std::size_t
      n_ltr() const { return ltr_.size(); }

      //! Access to the group of lattice translation vectors.
      /*! Not available in Python.
       */
      tr_group
      const& ltr() const { return ltr_; }

      //! Access to elements of the group of lattice translation vectors.
      /*! Not available in Python.
       */
      tr_vec
      const& ltr(std::size_t i) const { return ltr_[i]; }

      //! Tests for a centre of inversion.
      bool
      is_centric() const { return f_inv_ == 2; }

      //! Tests for a centre of inversion at the origin.
      bool
      is_origin_centric() const
      {
        return is_centric() && inv_t(true).is_zero();
      }

      //! Flag for centre of inversion.
      /*! f_inv() = 1 if no inversion operation exists,<br>
          f_inv() = 2 otherwise.
       */
      std::size_t
      f_inv() const { return f_inv_; }

      //! Translation part for centre of inversion.
      /*! Not available in Python.
       */
      tr_vec
      inv_t(bool tidy=false) const
      {
        if (tidy == false) return inv_t_;
        if (!inv_t_.is_valid()) return inv_t_;
        return ltr_.tidy(inv_t_);
      }

      //! Number of representative Seitz matrices.
      std::size_t
      n_smx() const { return smx_.size(); }

      //! Access to the list of representative Seitz matrices.
      /*! Not available in Python.
       */
      smx_array_type
      const& smx() const { return smx_; }

      //! Access to elements of the list of representative Seitz matrices.
      /*! Not available in Python.
       */
      rt_mx
      const& smx(std::size_t i) const { return smx_[i]; }

      //! Returns a symmetry operation.
      /*! Usage:<pre>
          space_group sg(...);
          for (std::size_t i_ltr = 0; i_ltr < sg.n_ltr(); i_ltr++)
            for (std::size_t i_inv = 0; i_inv < sg.f_inv(); i_inv++)
              for (std::size_t i_smx = 0; i_smx < sg.n_smx(); i_smx++)
                rt_mx M = sg(i_ltr, i_inv, i_smx);</pre>
          The symmetry operations are generated from the representative
          list of Seitz matrices. The lattice translation with the index
          i_ltr is added to the translation part. If i_inv = 1, the
          representative Seitz matrix is pre-multiplied by (-I, inv_t).<br>
          Note that the translation part of the returned Seitz matrix
          is NOT modified by application of the modulus operation.
       */
      rt_mx
      operator()(
        std::size_t i_ltr,
        std::size_t i_inv,
        std::size_t i_smx) const;

      //! Returns a symmetry operation.
      /*! Usage:<pre>
          space_group sg(...);
          for (std::size_t i_op = 0; i_op < sg.order_z(); i_op++)
            rt_mx s = sg(i_op);</pre>
          The symmetry operations are generated from the representative
          list of Seitz matrices. Internally i_op is interpreted as:<br>
          i_op = ((i_ltr * f_inv()) + i_inv) * n_smx() + i_smx<br>
          The comments for
          operator()(std::size_t i_ltr, std::size_t i_inv, std::size_t i_smx)
          apply.
          <p>
          Python: __getitem__
       */
      rt_mx
      operator()(std::size_t i_op) const;

      //! Tidies the lists of representative symmetry operations in place.
      /*! The list of lattice translations is sorted in a certain order.
          If there is a centre of inversion, a certain normalized
          translation part is selected for inv_t(), and a proper
          rotation part is selected for all Seitz matrices in the
          representative list. The list of representative Seitz
          matrices is sorted, and the translation parts are
          normalized.
          <p>
          After application of make_tidy(), a given space group
          representation will always result in exactly the same
          internal representation.
       */
      space_group&
      make_tidy();

      //! True if make_tidy() was called before.
      bool
      is_tidy() const { return is_tidy_; }

      //! Tests if smx modulo 1 is an operation of the group.
      bool
      contains(rt_mx const& smx) const;

      //! Tests for equality.
      /*! Internally, make_tidy() is used, followed by essentially a
          byte-wise comparison of the objects.<br>
          Each space_group object maintains an internal flag
          indicating whether or not make_tidy() was applied already.
          If an space_group object is repeatedly used in a test for equality,
          the test will therefore be significantly faster if make_tidy()
          is applied outside the loop.
       */
      bool
      operator==(space_group const& rhs) const;

      //! Negation of test for equality.
      bool
      operator!=(space_group const& rhs) const
      {
        return !((*this) == rhs);
      }

      //! Determines the conventional centring type symbol (P,A,B,C,I,R,H,F).
      /*! If the lattice translation vectors do not correspond to
          any of the conventional centring types, the null character
          is returned.
       */
      char
      conventional_centring_type_symbol() const
      {
        return ltr_.conventional_centring_type_symbol();
      }

      //! Determines an operator for centred-to-primitive basis transformation.
      /*! Determines a change-of-basis operator that transforms the
          symmetry from a centered to a primitive setting.<br>
          For the conventional lattice centring types (P, A, B, C, I, R, H, F)
          a standard change-of-basis operator is determined by a lookup in a
          table. For unconventional lattice centring types, a change-of-basis
          operator is constructed with the algorithm published in
          <A HREF="http://journals.iucr.org/a/issues/1999/02/02/au0146/"
          ><I>Acta Cryst.</I> 1999, <B>A55</B>:383-395</A>.<br>
          The change-of-basis operator can be used as the argument to
          the change_basis() member function:<pre>
          space_group centred(...);
          space_group primitive = centred.change_basis(centred.z2p_op());</pre>
       */
      change_of_basis_op
      z2p_op(int r_den=cb_r_den, int t_den=cb_t_den) const;

      //! Constructs an operator for centred->primitive basis transformation.
      /*! Similar to z2p_op(), but does not use standard change-of-basis
          operators for conventional lattice centring types.<br>
          For internal use only.
       */
      change_of_basis_op
      construct_z2p_op(int r_den=cb_r_den, int t_den=cb_t_den) const;

      //! Tests for chirality.
      /*! A space group is chiral if all its symmetry operations
          have a positive rotation-part type (1, 2, 3, 4, 6).<br>
          If there are symmetry operations with negative rotation-part
          types (-1, -2=m, -3, -4, -6) the space group is not chiral.<br>
          There are exactly 65 chiral space groups.<br>
          Note that proteins always crystallize in a chiral space
          group.
       */
      bool
      is_chiral() const;

      /*! \brief Tests if a reflection with given Miller index is
          systematically absent.
       */
      bool
      is_sys_absent(miller::index<> const& miller_index) const
      {
        return phase_info(*this, miller_index, false).is_sys_absent();
      }

      /*! \brief Tests if a reflection with given Miller index is
          systematically absent.
       */
      /*! Overload for arrays.
       */
      af::shared<bool>
      is_sys_absent(
        af::const_ref<miller::index<> > const& miller_indices) const;

      //! Tests if a reflection with given Miller index is centric.
      /*! A reflection with the Miller index h is "centric" if
          there is a symmetry operation with rotation part r such
          that h*r = -h.<br>
          The phase of a centric reflection is restricted to two phase
          angels (modulo pi).<br>
          If the phase restriction is also needed later in the
          calculation it is more efficient to use phase_restriction().
       */
      bool
      is_centric(miller::index<> const& miller_index) const;

      //! Tests if a reflection with given Miller index is centric.
      /*! Overload for arrays.
       */
      af::shared<bool>
      is_centric(
        af::const_ref<miller::index<> > const& miller_indices) const;

      //! Determines the phase restriction for the given Miller index.
      /*! See class phase_info. The conditions for systematically
          absent reflections are NOT tested.
       */
      phase_info
      phase_restriction(miller::index<> const& miller_index) const
      {
        return phase_info(*this, miller_index, true);
      }

      //! See class phase_info::is_valid_phase().
      bool
      is_valid_phase(
        miller::index<> const& miller_index,
        double phi,
        bool deg=false,
        double tolerance=1e-5) const
      {
        return phase_restriction(miller_index)
          .is_valid_phase(phi, deg, tolerance);
      }

      //! Computes nearest phases compatible with the phase restrictions.
      /*! See also: class phase_info::nearest_valid_phase()
       */
      af::shared<double>
      nearest_valid_phases(
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<double> const& phases,
        bool deg=false) const;

      //! Determines the reflection multiplicity for the given Miller index.
      /*! The multiplicity is defined as the number of symmetrically
          equivalent but distinct reflections.<br>
          If anomalous_flag == false, a centre of inversion is added to the
          list of symmetry matrices considered in the determination of
          the multiplicity.
       */
      int
      multiplicity(
        miller::index<> const& miller_index,
        bool anomalous_flag) const;

      //! Determines the reflection multiplicity for the given Miller index.
      /*! Overload for arrays.
       */
      af::shared<int>
      multiplicity(
        af::const_ref<miller::index<> > const& miller_indices,
        bool anomalous_flag) const;

      //! Determines "epsilon" for the given Miller index.
      /*! The factor epsilon counts the number of times a Miller
          index h is mapped onto itself by symmetry. This factor
          is used for "statistical averaging" and in direct methods
          formulae.<br>
          See also: miller::sym_equiv_indices::epsilon()
       */
      int
      epsilon(miller::index<> const& miller_index) const;

      //! Determines "epsilon" for the given Miller index.
      /*! Overload for arrays.
       */
      af::shared<int>
      epsilon(af::const_ref<miller::index<> > const& miller_indices) const;

      //! Determines multiplicity of site given rational coordinates.
      int
      multiplicity(
        vec3_rat const& site) const;

      //! Computes a metrical matrix compatible with the space group symmetry.
      /*! A metrical matrix g is compatible with a given space group
          representation if the following relation holds for all
          rotation matrices r of the space group:
          <p>
          r.transposed() * g * r == g
          <p>
          This function returns the average of r.transposed() * g * r.
          The result is guaranteed to be compatible with the space
          group symmetry.
          <p>
          See also: average_unit_cell(),
                    uctbx::unit_cell::metrical_matrix()
          <p>
          Not available in Python.
       */
      uc_sym_mat3
      average_metrical_matrix(uc_sym_mat3 const& g) const
      {
        return average_tensor(smx_.const_ref(), g, false);
      }

      //! Computes a unit cell compatible with the space group symmetry.
      /*! Shorthand for:
          uctbx::unit_cell(average_metrical_matrix(ucell.metrical_matrix()));
       */
      uctbx::unit_cell
      average_unit_cell(uctbx::unit_cell const& unit_cell) const
      {
        return uctbx::unit_cell(
          average_metrical_matrix(unit_cell.metrical_matrix()));
      }

      //! Tests if a unit cell is compatible with the symmetry operations.
      bool
      is_compatible_unit_cell(
        uctbx::unit_cell const& unit_cell,
        double relative_length_tolerance=0.01,
        double absolute_angle_tolerance=1.) const
      {
        return unit_cell.is_similar_to(
          average_unit_cell(unit_cell),
          relative_length_tolerance,
          absolute_angle_tolerance);
      }

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
        return average_tensor(smx_.const_ref(), u_star, true);
      }

      //! The translation parts of the symmetry operations are set to 0.
      /*! If discard_z = false, the lattice translation vectors are not
          modified. If add_inv = true, a centre of inversion is added
          at the origin.
          <p>
          See also:
            build_derived_reflection_intensity_group(),
            build_derived_patterson_group(),
            build_derived_point_group(),
            build_derived_laue_group()
          <p>
          Not available in Python.
       */
      space_group
      build_derived_group(bool discard_z, bool add_inv) const;

      //! New instance without a centre of inversion.
      space_group
      build_derived_acentric_group() const;

      //! Builds the symmetry for reflection intensities.
      /*! The translation parts of the symmetry operations are set to 0.
          However, the lattice translation vectors are not modified.
          If anomalous_flag = false a centre of inversion is added at
          the origin.
       */
      space_group
      build_derived_reflection_intensity_group(bool anomalous_flag) const
      {
        return build_derived_group(false, !anomalous_flag);
      }

      //! Builds the corresponding Patterson space group.
      /*! The translation parts of the symmetry operations are set to 0.
          However, the lattice translation vectors are not modified.
          A centre of inversion is added at the origin.
       */
      space_group
      build_derived_patterson_group() const
      {
        return build_derived_group(false, true);
      }

      //! Builds the corresponding point group.
      /*! The translation parts of the symmetry operations are set to 0,
          and the lattice translation vectors are discarded.
       */
      space_group
      build_derived_point_group() const
      {
        return build_derived_group(true, false);
      }

      //! Builds the corresponding Laue group.
      /*! The translation parts of the symmetry operations are set to 0,
          the lattice translation vectors are discarded,
          and a centre of inversion is added at the origin.
       */
      space_group
      build_derived_laue_group() const
      {
        return build_derived_group(true, true);
      }

      //! Histogram of rotation-part types.
      /*! The histogram is accumulated in a std::map. The key
          is the rotation-part type (-6,-4,-3,-2,-1,1,2,3,4,6),
          and the value is the number of times this rotation-part
          type occurs in the list of representative symmetry
          matrices.<br>
          This function is used by the algorithm for the determination
          of the space group type and its usefulness is probably
          limited to this context.
          <p>
          Not available in Python.
       */
      std::map<int, int>
      count_rotation_part_types() const;

      //! Determines the point group type.
      /*! The code returned is a matrix group code. There are
          exactly 32 possible return values, corresponding to
          the 32 crystallographic point group types.
          <p>
          Python: returns a string representing the point group type.
       */
      matrix_group::code
      point_group_type() const;

      //! Determines the Laue group type.
      /*! The code returned is a matrix group code. There are
          exactly 11 possible return values, corresponding to
          the 11 Laue group types.
          <p>
          Python: returns a string representing the Laue group type.
       */
      matrix_group::code
      laue_group_type() const
      {
        return point_group_type().laue_group_type();
      }

      //! Determines the crystal system.
      /*! There are exactly 7 possible return values.
          <p>
          Python: returns a string representing the crystal system.
       */
      crystal_system::code
      crystal_system() const
      {
        return point_group_type().crystal_system();
      }

      //! Matches given symmetry operations with the 530 tabulated settings.
      /*! This is a light-weight alternative to using space_group_type,
          but is applicable only to the 530 tabulated settings.
          <p>
          This algorithm is not particularly optimized. Therefore
          the runtime for match_tabulated_settings() and the runtime for
          instantiating space_group_type are in general comparable.
          The main purpose of this algorithm is therefore the
          retrieval of conventionally used Hermann-Mauguin symbols.
       */
      space_group_symbols
      match_tabulated_settings() const;

      //! Refines gridding starting with grid (1,1,1).
      /*! See also: refine_gridding()
       */
      sg_vec3
      gridding() const
      {
        return refine_gridding(sg_vec3(1,1,1));
      }

      /*! \brief Refines gridding such that each grid point is
          mapped onto another grid point by all symmetry operations.
       */
      template <typename GridTupleType>
      GridTupleType
      refine_gridding(GridTupleType const& grid) const;

      //! Returns a straight list of all symmetry operations.
      /*! The operations are in the same order as returned by
          operator(std::size_t).
          <p>
          If mod > 0, rt_mx::mod_positive() is applied to all operations.<br>
          If mod < 0, rt_mx::mod_short() is applied to all operations.<br>
          If mod == 0, the translation parts are not treated in a special way.
          <p>
          If cancel == true, the translation denominators of the
          returned operations are made as small as possible by
          cancellation of factors.
       */
      af::shared<rt_mx>
      all_ops(int mod=0, bool cancel=false) const;

      /*! \brief Computes the unique symmetry operations given a
          special position operator.
       */
      /*! special_op is multiplied with all symmetry operations
          of the space group. The unique results are returned.
          <p>
          The rotation and translation denominators of the returned
          operations are made as small as possible by cancellation
          of factors.
          <p>
          See also: site_symmetry::unique_ops(),
                    wyckoff::position::unique_ops()
       */
      af::shared<rt_mx>
      unique(rt_mx const& special_op) const;

      //! Convenience method for instantiating class space_group_type.
      space_group_type type() const;

    private:
      bool no_expand_;
      std::size_t n_lsl_;
      std::size_t n_ssl_;
      std::size_t f_inv_;
      tr_group ltr_;
      tr_vec inv_t_;
      bool is_tidy_;
      smx_array_type smx_;

      void add_inv(tr_vec const& new_inv_t);

      void add_smx(rt_mx const& new_smx);
  };

  template <typename GridTupleType>
  GridTupleType
  space_group::refine_gridding(GridTupleType const& grid) const
  {
    GridTupleType prev_grid = grid;
    GridTupleType refined_grid = grid;
    for (;;) {
      for(std::size_t i=0;i<order_z();i++) {
        refined_grid = operator()(i).refine_gridding(refined_grid);
      }
      if (        af::make_const_ref(prev_grid)
          .all_eq(af::make_const_ref(refined_grid))) break;
      prev_grid = refined_grid;
    }
    return refined_grid;
  }

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_SPACE_GROUP_H
