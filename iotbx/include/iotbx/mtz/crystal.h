#ifndef IOTBX_MTZ_CRYSTAL_H
#define IOTBX_MTZ_CRYSTAL_H

#include <iotbx/mtz/object.h>

namespace iotbx { namespace mtz {

  //! Safe access to CMtz::MTZXTAL* owned by an iotbx::mtz::object .
  /*! A crystal object contains (owns) an iotbx::mtz::object and
      an integer index into the list of crystals owned by the
      object.
   */
  class crystal
  {
    public:
      //! Not available from Python.
      crystal() {}

      /*! \brief Initialization given a mtz_object and an integer
          index into the list of crystals owned by the mtz_object.
       */
      /*! An exception is thrown if i_crystal is out of range.
       */
      crystal(object const& mtz_object, int i_crystal);

      //! The contained iotbx::mtz::object instance.
      object
      mtz_object() const { return mtz_object_; }

      //! Integer index into the list of crystals owned by mtz_object() .
      int
      i_crystal() const { return i_crystal_; }

      //! Raw C pointer. Not available from Python.
      /*! An exception is thrown if i_crystal() is out of range.
       */
      CMtz::MTZXTAL*
      ptr() const;

      //! Read-only access.
      int
      id() const { return ptr()->xtalid; }

      //! Write access.
      crystal&
      set_id(int id);

      //! Read-only access.
      const char*
      name() const { return ptr()->xname; }

      //! Read-only access.
      const char*
      project_name() const { return ptr()->pname; }

      //! Read-only access.
      af::small<double, 6>
      unit_cell_parameters() const;

      //! Read-only access.
      cctbx::uctbx::unit_cell
      unit_cell() const
      {
        return cctbx::uctbx::unit_cell(unit_cell_parameters());
      }

      //! Read-only access.
      int
      n_datasets() const { return CMtz::MtzNsetsInXtal(ptr()); }

      //! Read-only access.
      af::shared<dataset>
      datasets() const;

      //! Adds a new dataset the this crystal.
      dataset
      add_dataset(
        const char *name,
        double wavelength);

    protected:
      object mtz_object_;
      int i_crystal_;
  };

}} // namespace iotbx::mtz

#endif // IOTBX_MTZ_CRYSTAL_H
