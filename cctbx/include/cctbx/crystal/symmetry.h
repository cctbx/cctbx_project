#ifndef CCTBX_CRYSTAL_SYMMETRY_H
#define CCTBX_CRYSTAL_SYMMETRY_H

#include <cctbx/sgtbx/space_group_type.h>

namespace cctbx { namespace crystal {

  //! Grouping of unit cell and space group.
  /*! This class is *not* directly related to cctbx.crystal.symmetry
      in Python.
   */
  class symmetry
  {
    public:
      //! Default constructor. Some data members are not initialized!
      symmetry() {}

      //! Initialization with unit cell and space group.
      symmetry(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group)
      :
        unit_cell_(unit_cell),
        space_group_(space_group)
      {
        space_group_.make_tidy();
      }

      //! Unit cell as passed to the constructor.
      uctbx::unit_cell const&
      unit_cell() const { return unit_cell_; }

      //! Space group as passed to the constructor.
      sgtbx::space_group const&
      space_group() const { return space_group_; }

      //! New symmetry object in a new basis system.
      symmetry
      change_basis(sgtbx::change_of_basis_op const& cb_op) const
      {
        return symmetry(
          unit_cell_.change_basis(cb_op),
          space_group_.change_basis(cb_op));
      }

    protected:
      uctbx::unit_cell unit_cell_;
      sgtbx::space_group space_group_;
  };

}} // namespace cctbx::crystal

#endif // CCTBX_CRYSTAL_SYMMETRY_H
