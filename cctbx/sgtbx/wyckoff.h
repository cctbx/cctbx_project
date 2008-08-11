#ifndef CCTBX_SGTBX_WYCKOFF_H
#define CCTBX_SGTBX_WYCKOFF_H

#include <cctbx/sgtbx/site_symmetry.h>
#include <cctbx/sgtbx/space_group_type.h>

namespace cctbx { namespace sgtbx {

//! Wyckoff tables.
namespace wyckoff {

  class table; // forward declaration

  //! Information for Wyckoff positions.
  /*! See also: wyckoff::table, wyckoff::mapping, sym_equiv_coordinates
      <p>
      Python binding: wyckoff_position
   */
  class position
  {
    public:
      //! Default constructor. Some data members are not initialized!
      position() {}

      //! Back-reference to the Wyckoff table.
      /*! For efficiency, the table is not copied. This is, the
          wyckoff::table instance must exist as long as this
          reference is used.
          <p>
          Not available in Python.
       */
      wyckoff::table const&
      table() const { return *table_; }

      //! The multiplicity of the Wyckoff position.
      int
      multiplicity() const { return multiplicity_; }

      //! The Wyckoff letter according to the International Tables.
      /*! The Wyckoff letter "alpha" is represented by the character "@".
       */
      char
      letter() const { return letter_; }

      //! A representative special position operation.
      rt_mx const&
      special_op() const { return special_op_; }

      //! Point group type of the site-symmetry of the Wyckoff position.
      matrix_group::code
      point_group_type() const;

      //! Computes the unique symmetry operations for a Wyckoff position.
      /*! special_op() is multiplied with all symmetry operations
          of space_group(). The unique results are returned.
          <p>
          See also: space_group::unique()
       */
      af::shared<rt_mx>
      unique_ops(const sgtbx::space_group& space_group);

    private:
      wyckoff::table const* table_;
      int multiplicity_;
      char letter_;
      rt_mx special_op_;

      friend class wyckoff::table;

      position(wyckoff::table const* table,
               int multiplicity,
               char letter,
               rt_mx const& special_op);
  };

  //! Pair of wyckoff::position and a symmetry operation.
  /*! The positions are defined by a representative special
      position operation. These are tabulated according to the
      conventions of the International Tables for Crystallography,
      Volume A, 1983. A representative special position operation
      only applies to one particular point, axis, or plane.
      To determine the Wyckoff position of an aritrary point,
      a symmetry operation (sym_op()) is constructed that maps the
      original_site() to the point, axis, or plane for which
      the representative operation is defined.
      <p>
      Both the rotation part and the translation part of sym_op() are
      unique for a given position. In particular, the modulus operation
      that can be applied to regular symmetry operations must not be
      applied to sym_op().
      <p>
      See also: wyckoff::table, wyckoff::position, sym_equiv_coordinates
      <p>
      Python binding: wyckoff_mapping
   */
  class mapping
  {
    public:
      //! Default constructor. Calling certain methods may cause crashes!
      mapping()
      :
        position_(0)
      {}

      //! Unit cell used in the computation of the Wyckoff mapping.
      uctbx::unit_cell const&
      unit_cell() const { return unit_cell_; }

      //! Original site used to compute the Wyckoff position().
      fractional<>
      original_site() const { return original_site_; }

      //! Reference to entry in wyckoff::table.
      /*! For efficiency, the entry is not copied. This is, the
          wyckoff::table instance must exist as long as this
          reference is used.
       */
      wyckoff::position const&
      position() const { return *position_; }

      /*! Symmetry operation that maps the original_site()
          to the representative_site().
       */
      /*! See class details and wyckoff::table::mapping().
       */
      rt_mx const&
      sym_op() const { return sym_op_; }

      //! Exact location of the representative Wyckoff position.
      /*! Formula used: position().special_op() * sym_op() * original_site()
          <p>
          The result can be used directly with the representative
          Wyckoff position to compute sym_equiv_coordinates. This
          is useful for repeated, efficient computation of
          symmetry mates.
       */
      fractional<>
      representative_site() const
      {
        return position_->special_op() * (sym_op_ * original_site_);
      }

      //! Exact location of the special position.
      /*! This is equivalent to site_symmetry::exact_site().
          <p>
          Formula used: sym_op().inverse_cancel() * representative_site()
       */
      fractional<>
      exact_site() const
      {
        return sym_op_.inverse_cancel() * representative_site();
      }

      //! Distance^2 between original_site() and exact_site().
      /*! Not available in Python.
       */
      double
      distance_moved_sq() const
      {
        return unit_cell_.distance_sq(exact_site(), original_site_);
      }

      //! Distance between original_site() and exact_site().
      double
      distance_moved() const
      {
        return std::sqrt(distance_moved_sq());
      }

      //! Special position operator that maps original_site() to exact_site().
      /*! This is equivalent to site_symmetry::special_op().
       */
      rt_mx
      special_op() const
      {
        return sym_op_.inverse_cancel()
               .multiply(position_->special_op())
               .multiply(sym_op_);
      }

    private:
      uctbx::unit_cell unit_cell_;
      fractional<> original_site_;
      const wyckoff::position* position_;
      rt_mx sym_op_;

      friend class table;

      mapping(uctbx::unit_cell const& unit_cell,
              fractional<> const& original_site,
              wyckoff::position const& position,
              rt_mx const& sym_op)
      : unit_cell_(unit_cell),
        original_site_(original_site),
        position_(&position),
        sym_op_(sym_op)
      {}
  };

  //! Table of Wyckoff positions.
  /*! This class represents the information for all Wyckoff positions
      for a given space group.
      <p>
      See also: wyckoff::position, wyckoff::mapping, site_symmetry,
                sym_equiv_coordinates
      <p>
      Python binding: wyckoff_table
   */
  class table
  {
    public:
      //! Initializes a table of size() zero.
      /*! Not available in Python.
       */
      table() {}

      //! Constructor.
      /*! The Wyckoff positions for the 230 reference settings are
          tabulated. sgtbx::space_group_type::cb_op() is used to
          transform the tabulated settings to the given setting.
       */
      explicit
      table(sgtbx::space_group_type const& space_group_type);

      //! The space group information as passed to the constructor.
      sgtbx::space_group_type const&
      space_group_type() const { return space_group_type_; }

      //! The number of Wyckoff positions for the given space group.
      /*! This number varies between 1 and 27.
          <p>
          Shorthand for positions().size()
       */
      std::size_t
      size() const { return positions_.size(); }

      //! Returns an entry from the list of Wyckoff positions().
      wyckoff::position const&
      position(std::size_t i) const
      {
        CCTBX_ASSERT(i < positions_.size());
        return positions_[i];
      }

      //! Retrieves the Wyckoff position with the given letter.
      /*! An exception is thrown if the letter is not valid.
          <p>
          The Wyckoff letter "alpha" must be given as the character "@".
          <p>
          Retrieving a Wyckoff position by letter is the slowest form
          of access. For repeated access, use lookup_index() to map
          the Wyckoff letter to an index.
          For a given space group type, there is a setting independent
          one-to-one correspondence between letters and indices.
       */
      wyckoff::position const&
      position(char letter) const
      {
        return positions_[lookup_index(letter)];
      }

      //! Returns the list of Wyckoff positions.
      /*! Not available in Python.
       */
      af::shared<wyckoff::position>
      positions() const { return positions_; }

      //! Determines the Wyckoff position index for a given Wyckoff letter.
      /*! An exception is thrown if the letter is not valid.
          <p>
          The Wyckoff letter "alpha" must be given as the character "@".
          <p>
          For a given space group type, there is a setting independent
          one-to-one correspondence between letters and indices.
       */
      std::size_t
      lookup_index(char letter) const;

      //! Determines the Wyckoff position using site_symmetry::special_op().
      /*! Use this overload for maximum reliability and flexibility.
          In general, determining site_symmetry::special_op() followed
          by using this overload is also faster than the alternative.
          Usage:<pre>
          uctbx::unit_cell ucell = ...;
          sgtbx::space_group sg = ...;
          fractional<> x = ...;
          sgtbx::site_symmetry ss(ucell, sg, x);
          sgtbx::wyckoff::mapping wm = mapping(ss);</pre>
       */
      wyckoff::mapping
      mapping(sgtbx::site_symmetry const& site_symmetry) const;

      //! Determines the Wyckoff position of x.
      /*! This overload is included mainly for debugging purposes.
          It is highly recommended to use the alternative algorithm.
       */
      wyckoff::mapping
      mapping(
        uctbx::unit_cell const& unit_cell,
        fractional<> const& original_site,
        double special_position_radius=0.5) const;

    private:
      sgtbx::space_group_type space_group_type_;
      af::shared<wyckoff::position> positions_;
  };

}}} // namespace cctbx::sgtbx::wyckoff

#endif // CCTBX_SGTBX_WYCKOFF_H
