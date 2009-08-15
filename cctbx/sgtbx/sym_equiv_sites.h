#ifndef CCTBX_SGTBX_SYM_EQUIV_SITES_H
#define CCTBX_SGTBX_SYM_EQUIV_SITES_H

#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/wyckoff.h>

namespace cctbx { namespace sgtbx {

  //! Container for symmetrically equivalent coordinates.
  template <typename FloatType=double>
  class sym_equiv_sites
  {
    public:
      //! Convenience typedef.
      typedef fractional<FloatType> frac_t;
      //! Convenience typedef.
      typedef scitbx::vec3<FloatType> coor_t;

      //! Default constructor. Some data members are not initialized!
      sym_equiv_sites() {}

      /*! \brief Computes symmetrically equivalent coordinates given
          a site symmetry.
       */
      explicit
      sym_equiv_sites(site_symmetry const& site_sym)
      :
        unit_cell_(site_sym.unit_cell()),
        space_group_(site_sym.space_group()),
        original_site_(site_sym.original_site()),
        special_op_(site_sym.special_op()),
        use_special_op_(true),
        max_accepted_tolerance_(-1)
      {
        initialize_with_special_op(site_sym.multiplicity());
        CCTBX_ASSERT(coordinates_.size() == site_sym.multiplicity());
      }

      /*! \brief Computes symmetrically equivalent coordinates given
          site_symmetry_ops.
       */
      /*! See also: site_symmetry_table
       */
      sym_equiv_sites(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        frac_t const& original_site,
        sgtbx::site_symmetry_ops const& site_symmetry_ops)
      :
        unit_cell_(unit_cell),
        space_group_(space_group),
        original_site_(original_site),
        special_op_(site_symmetry_ops.special_op()),
        use_special_op_(true),
        max_accepted_tolerance_(-1)
      {
        std::size_t multiplicity = site_symmetry_ops.multiplicity();
        initialize_with_special_op(multiplicity);
        CCTBX_ASSERT(coordinates_.size() == multiplicity);
      }

      /*! \brief Computes symmetrically equivalent coordinates given
          a Wyckoff mapping.
       */
      explicit
      sym_equiv_sites(wyckoff::mapping const& wyckoff_mapping)
      :
        unit_cell_(wyckoff_mapping.unit_cell()),
        space_group_(
          wyckoff_mapping.position().table().space_group_type().group()),
        original_site_(wyckoff_mapping.original_site()),
        special_op_(wyckoff_mapping.special_op()),
        use_special_op_(true),
        max_accepted_tolerance_(-1)
      {
        initialize_with_special_op(wyckoff_mapping.position().multiplicity());
        CCTBX_ASSERT(
          coordinates_.size() == wyckoff_mapping.position().multiplicity());
      }

      /*! \brief Computes symmetrically equivalent coordinates without a
          treatment of special positions.
       */
      /*! The symmetry operations of space_group are applied to orignal_site.
          Duplicates on special positions are not removed. This algorithm is
          therefore very fast and is suitable for structures with most
          atoms in general positions (e.g. protein structures).
          <p>
          unit_cell is not used in the computation of the equivalent
          coordinates and may therefore be omitted. However, to
          faciliate subsequent use of the sym_equiv_sites instance it
          is good practice to supply the unit cell parameters if
          available.
          <p>
          See also: site_symmetry, wyckoff::mapping
       */
      sym_equiv_sites(
        sgtbx::space_group const& space_group,
        frac_t const& original_site,
        uctbx::unit_cell const& unit_cell=uctbx::unit_cell())
      :
        unit_cell_(unit_cell),
        space_group_(space_group),
        original_site_(original_site),
        special_op_(0, 0),
        use_special_op_(false),
        max_accepted_tolerance_(-1)
      {
        initialize_trivial();
        CCTBX_ASSERT(coordinates_.size() == space_group_.order_z());
      }

      /*! \brief Computes symmetrically equivalent coordinates given
          a special position operator.
       */
      /*! In most situations it will be more convenient to use the
          constructor with the site_symmetry parameter, or the
          constructor with the wyckoff::mapping parameter.
          <p>
          See also:
            site_symmetry,
            wyckoff::mapping
       */
      sym_equiv_sites(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        frac_t const& original_site,
        rt_mx const& special_op,
        std::size_t multiplicity=0)
      :
        unit_cell_(unit_cell),
        space_group_(space_group),
        original_site_(original_site),
        special_op_(special_op),
        use_special_op_(true),
        max_accepted_tolerance_(-1)
      {
        initialize_with_special_op(multiplicity);
        CCTBX_ASSERT(space_group_.order_z() % coordinates_.size() == 0);
      }

      /*! \brief Computes symmetrically equivalent coordinates using simple
          distance calculations.
       */
      /*! This function should only be used with preprocessed coordinates
          that are known to be well defined. In general it is best to use
          class site_symmetry instead. Unless there is a large number of
          special positions in a structure, the computation via class
          site_symmetry is only marginally slower.

          If the distance between symmetrically equivalent sites is less
          than or equal to tolerance, the site x is considered to be on a
          special position. However, the equivalent sites are not moved to
          the exact location of the special position. Therefore the value
          for tolerance should in general be very small.

          The member function max_accepted_tolerance() returns the
          largest tolerance accepted in the determination of special
          positions.

          If the distance between symmetry mates is less than
          minimum_distance, but not less than or equal to tolerance, an
          exception is raised.

          To avoid numerical instabilities and potentially inaccurate
          results, minimum_distance should be strictly greater than
          tolerance. Only under this condition are the results
          guaranteed to be reliable.

          As a safeguard it is asserted that the number of symmetrically
          equivalent sites is a factor of the space group multiplicity.
       */
      sym_equiv_sites(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        frac_t const& original_site,
        FloatType const& minimum_distance=0.5,
        FloatType const& tolerance=0.01)
      :
        unit_cell_(unit_cell),
        space_group_(space_group),
        original_site_(original_site),
        special_op_(0, 0),
        use_special_op_(false),
        max_accepted_tolerance_(0)
      {
        initialize_with_simple_distance_calculations(
          minimum_distance * minimum_distance,
          tolerance * tolerance);
      }

      /*! \brief Unit cell used in the computation of the
          symmetrically equivalent sites.
       */
      uctbx::unit_cell const&
      unit_cell() const
      {
        CCTBX_ASSERT(unit_cell_.volume() > 0);
        return unit_cell_;
      }

      /*! \brief Space group used in the computation of the
          symmetrically equivalent sites.
       */
      sgtbx::space_group const&
      space_group() const { return space_group_; }

      /*! \brief Original site used in the computation of the
          symmetrically equivalent sites.
       */
      frac_t const&
      original_site() const { return original_site_; }

      /*! \brief Special position operator used in the computation of the
          symmetrically equivalent sites.
       */
      /*! Use special_op().is_valid() to test if a special position
          operator was used or not.
          <p>
          See also:
            site_symmetry,
            wyckoff::mapping
       */
      rt_mx const&
      special_op() const { return special_op_; }

      /*! \brief Larget tolerance accepted in the determination of special
          positions via simple distance calculations.
       */
      /*! This value is >= 0 only if the simple distance calculations
          are used.
       */
      FloatType
      max_accepted_tolerance() const { return max_accepted_tolerance_; }

      //! Array of computed coordinates.
      af::shared<coor_t>
      coordinates() const { return coordinates_; }

      /*! \brief Symmetry operation that gave rise to the i_coor'th element
          of coordinates().
       */
      rt_mx
      sym_op(std::size_t i_coor) const
      {
        CCTBX_ASSERT(i_coor < sym_op_indices_.size());
        return space_group_(sym_op_indices_[i_coor]);
      }

      //! Array of symmetry-operation indices corresponding to coordinates().
      /*! See also: sym_op()
       */
      af::shared<std::size_t>
      sym_op_indices() const { return sym_op_indices_; }

      //! True if coordinates.size() < space_group().order_z().
      bool
      is_special_position() const
      {
        return coordinates_.size() < space_group_.order_z();
      }

    protected:
      uctbx::unit_cell unit_cell_;
      sgtbx::space_group space_group_;
      frac_t original_site_;
      rt_mx special_op_;
      bool use_special_op_;
      FloatType max_accepted_tolerance_;
      af::shared<std::size_t> sym_op_indices_;
      af::shared<coor_t> coordinates_;

      void
      reserve(std::size_t sz)
      {
        sym_op_indices_.reserve(sz);
        coordinates_.reserve(sz);
      }

      void
      push_back(std::size_t i, coor_t const& x)
      {
        sym_op_indices_.push_back(i);
        coordinates_.push_back(x);
      }

      void initialize_trivial();

      void initialize_with_special_op(std::size_t multiplicity);

      void initialize_with_simple_distance_calculations(
        FloatType const& minimum_distance_sq,
        FloatType const& tolerance_sq);
  };

  template <typename FloatType>
  void
  sym_equiv_sites<FloatType>
  ::initialize_trivial()
  {
    reserve(space_group_.order_z());
    push_back(0, original_site_);
    std::size_t n_smx = space_group_.n_smx();
    std::size_t i = 1;
    for(std::size_t j=1;j<n_smx;j++,i++) {
      push_back(i, space_group_.smx(i) * original_site_);
    }
    if (space_group_.is_centric()) {
      coor_t inv_t = space_group_.inv_t()
                       .as_floating_point(scitbx::type_holder<FloatType>());
      for(std::size_t j=0;j<n_smx;j++,i++) {
        push_back(i, -coordinates_[j] + inv_t);
      }
    }
    std::size_t n_primitive = i;
    for(std::size_t l=1;l<space_group_.n_ltr();l++) {
      coor_t ltr = space_group_.ltr(l)
                     .as_floating_point(scitbx::type_holder<FloatType>());
      for(std::size_t j=0;j<n_primitive;j++,i++) {
        push_back(i, coordinates_[j] + ltr);
      }
    }
  }

  template <typename FloatType>
  void
  sym_equiv_sites<FloatType>
  ::initialize_with_special_op(std::size_t multiplicity)
  {
    CCTBX_ASSERT(special_op_.is_valid());
    use_special_op_ = !special_op_.is_unit_mx();
    if (!use_special_op_) {
      initialize_trivial();
    }
    else {
      std::vector<rt_mx> unique_ops;
      if (multiplicity) {
        reserve(multiplicity);
        unique_ops.reserve(multiplicity);
      }
      for(std::size_t i=0;i<space_group_.order_z();i++) {
        rt_mx s = space_group_(i).multiply(special_op_);
        rt_mx sm = s.mod_positive();
        if (   std::find(unique_ops.begin(), unique_ops.end(), sm)
            == unique_ops.end()) {
          unique_ops.push_back(sm);
          push_back(i, s * original_site_);
        }
      }
    }
  }

  template <typename FloatType>
  void
  sym_equiv_sites<FloatType>
  ::initialize_with_simple_distance_calculations(
    FloatType const& minimum_distance_sq,
    FloatType const& tolerance_sq)
  {
    CCTBX_ASSERT((tolerance_sq == 0) == (minimum_distance_sq == 0));
    CCTBX_ASSERT(tolerance_sq == 0 || tolerance_sq < minimum_distance_sq);
    push_back(0, original_site_);
    for(std::size_t i=1;i<space_group_.order_z();i++) {
      frac_t sx = space_group_(i) * original_site_;
      FloatType
        dist_sq = unit_cell_.min_mod_short_distance_sq(
          coordinates_.const_ref(), sx);
      if (tolerance_sq == 0) {
        push_back(i, sx);
      }
      else if (dist_sq >= tolerance_sq) {
        if (dist_sq < minimum_distance_sq) {
          throw error(
            "Special position not well defined: use class site_symmetry.");
        }
        push_back(i, sx);
      }
      else {
        scitbx::math::update_max(max_accepted_tolerance_, dist_sq);
      }
    }
    if (space_group_.order_z() % coordinates_.size() != 0) {
      throw error("Numerical instability: use class site_symmetry.");
    }
    max_accepted_tolerance_ = std::sqrt(max_accepted_tolerance_);
  }

  //! Default value for min_sym_equiv_distance_info constructor.
  static const af::tiny<bool, 3>
    no_continuous_allowed_shifts = af::tiny<bool, 3>(false,false,false);

  /*! \brief Algorithm for finding the shortest distance
      between two sites under application of space group
      symmetry.
      <p>
      This class may be used to move groups of sites, such as
      residues in proteins, or entire molecules. For example, with
      reference to the constructor of this class, let the center of
      mass of one residue be the "reference site", and the center of
      mass of a second residue "other." Use the member function apply()
      to move the center of mass of the second residue to the
      symmetrically equivalent position with the shortest distance to
      the reference site.
   */
  /*! See also:
        sym_equiv_sites,
        uctbx::unit_cell::min_mod_short_distance()
   */
  template <typename FloatType=double>
  class min_sym_equiv_distance_info
  {
    public:
      //! Convenience typedef.
      typedef fractional<FloatType> frac_t;
      //! Convenience typedef.
      typedef scitbx::vec3<FloatType> coor_t;

      //! Default constructor. Some data members are not initialized!
      min_sym_equiv_distance_info() {}

      /*! \brief Computes the shortest distance between other and
          reference_sites.coordinates()[0].
       */
      /*! If principal_continuous_allowed_origin_shift_flags is
          supplied, the three flags signal a continuous allowed
          origin shift along a,b,c, respectively. These shifts
          are considered in the distance calculations. The
          member function continuous_shifts() returns the
          shifts that were applied in the calculation of the
          shortest distance.
       */
      min_sym_equiv_distance_info(
        sym_equiv_sites<FloatType> const& reference_sites,
        frac_t const& other,
        af::tiny<bool, 3> const&
          principal_continuous_allowed_origin_shift_flags
            = no_continuous_allowed_shifts)
      {
        init(
          reference_sites,
          af::const_ref<coor_t>(&other, 1),
          principal_continuous_allowed_origin_shift_flags);
      }

      /*! \brief Computes the shortest distance between others and
          reference_sites.coordinates()[0].
       */
      min_sym_equiv_distance_info(
        sym_equiv_sites<FloatType> const& reference_sites,
        af::const_ref<coor_t> const& others,
        af::tiny<bool, 3> const&
          principal_continuous_allowed_origin_shift_flags
            = no_continuous_allowed_shifts)
      {
        CCTBX_ASSERT(others.size() > 0);
        init(
          reference_sites,
          others,
          principal_continuous_allowed_origin_shift_flags);
      }

      //! Index to coordinates in others array as passed to the constructor.
      /*! The index is zero if the non-array constructor was used.
       */
      std::size_t
      i_other() const { return i_other_; }

      //! Symmetry operation that gives rise to the shortest distance.
      /*! The symmetry operation is defined by the equation:
          <p>
            reference_sites.coordinates()[0]
              = sym_op() * other + continuous_shifts() + diff()
          <p>
          whith reference_sites and other as passed to the constructor.
          <p>
          See also: apply()
       */
      rt_mx const&
      sym_op() const { return sym_op_; }

      //! Continuous allowed origin shift applied in the distance calculation.
      frac_t const&
      continuous_shifts() const { return continuous_shifts_; }

      //! Difference vector according to the shortest distance.
      /*! The difference vector is defined by the equation:
          <p>
            reference_sites.coordinates()[0]
              = sym_op() * other + continuous_shifts() + diff()
          <p>
          whith reference_sites and other as passed to the constructor.
       */
      frac_t const&
      diff() const { return diff_; }

      //! Shortest distance^2.
      /*! Not available in Python.
       */
      FloatType
      dist_sq() const { return dist_sq_; }

      //! Shortest distance.
      FloatType
      dist() const { return std::sqrt(dist_sq_); }

      //! Function for moving groups of sites.
      /*! Applies sym_op() and continuous_shifts() to the coordinates
          in the array x.
          <p>
          Formula applied:
          <p>
          result[i] = sym_op() * x[i] + continuous_shifts()
          <p>
          The calculation is optimized if the array of
          continuous_shifts().is_zero() .
       */
      // Visual C++ 7 Internal Compiler Error if moved out of class body.
      af::shared<coor_t>
      apply(af::const_ref<coor_t> const& sites_frac) const
      {
        af::shared<coor_t> result((af::reserve(sites_frac.size())));
        if (continuous_shifts().is_zero()) {
          for(std::size_t i=0;i<sites_frac.size();i++) {
            result.push_back(sym_op_ * sites_frac[i]);
          }
        }
        else {
          for(std::size_t i=0;i<sites_frac.size();i++) {
            result.push_back(sym_op_ * sites_frac[i] + continuous_shifts_);
          }
        }
        return result;
      }

    protected:
      std::size_t i_other_;
      rt_mx sym_op_;
      frac_t continuous_shifts_;
      frac_t diff_;
      FloatType dist_sq_;

      void
      init(
        sym_equiv_sites<FloatType> const& reference_sites,
        af::const_ref<coor_t> const& others,
        af::tiny<bool, 3> const&
          principal_continuous_allowed_origin_shift_flags);
  };

  template <typename FloatType>
  void
  min_sym_equiv_distance_info<FloatType>::
  init(
    sym_equiv_sites<FloatType> const& reference_sites,
    af::const_ref<coor_t> const& others,
    af::tiny<bool, 3> const& continuous_shift_flags)
  {
    uctbx::unit_cell const& unit_cell = reference_sites.unit_cell();
    af::const_ref<coor_t>
      reference_coordinates = reference_sites.coordinates().const_ref();
    bool no_continuous_shifts = continuous_shift_flags.all_eq(false);
    std::size_t min_i_coor = 0;
    frac_t min_unit_shifts;
    FloatType min_dist_sq = -1;
    for(std::size_t i_corr=0;i_corr<reference_coordinates.size();i_corr++) {
      for(std::size_t i_other=0;i_other<others.size();i_other++) {
        frac_t other = others[i_other];
        frac_t diff_raw = other - reference_coordinates[i_corr];
        frac_t diff_mod = diff_raw.mod_short();
        FloatType dist_sq;
        if (no_continuous_shifts) {
          dist_sq = unit_cell.length_sq(diff_mod);
        }
        else {
          frac_t diff_proj;
          for (std::size_t i=0;i<3;i++) {
            diff_proj[i] = (continuous_shift_flags[i] ? 0 : diff_mod[i]);
          }
          dist_sq = unit_cell.length_sq(diff_proj);
        }
        if (min_dist_sq > dist_sq || min_dist_sq == -1) {
          i_other_ = i_other;
          min_i_coor = i_corr;
          min_unit_shifts = diff_raw - diff_mod;
          min_dist_sq = dist_sq;
        }
      }
    }
    CCTBX_ASSERT(min_dist_sq != -1);
    rt_mx s = reference_sites.sym_op(min_i_coor);
    sym_op_ = (s + tr_vec(min_unit_shifts.unit_shifts(), 1)
                   .scale(s.t().den()))
              .inverse();
    diff_ = reference_coordinates[0] - sym_op_ * others[i_other_];
    if (no_continuous_shifts) {
      continuous_shifts_ = frac_t(0,0,0);
    }
    else {
      frac_t diff_proj;
      for (std::size_t j=0;j<3;j++) {
        diff_proj[j] = (continuous_shift_flags[j] ? 0 : diff_[j]);
      }
      continuous_shifts_ = diff_ - diff_proj;
      diff_ = diff_proj;
    }
    dist_sq_ = unit_cell.length_sq(diff_);
    CCTBX_ASSERT(
      dist_sq_ <= min_dist_sq + unit_cell.longest_vector_sq() * 1.e-6);
  }

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_SYM_EQUIV_SITES_H
