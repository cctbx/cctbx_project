#ifndef IOTBX_MTZ_CRYSTAL_H
#define IOTBX_MTZ_CRYSTAL_H

#include <iotbx/mtz/object.h>

namespace iotbx { namespace mtz {

  class crystal
  {
    public:
      crystal() {}

      crystal(object const& mtz_object, int i_crystal);

      object
      mtz_object() const { return mtz_object_; }

      int
      i_crystal() const { return i_crystal_; }

      CMtz::MTZXTAL*
      ptr() const;

      int
      id() const { return ptr()->xtalid; }

      crystal&
      set_id(int id);

      const char*
      name() const { return ptr()->xname; }

      const char*
      project_name() const { return ptr()->pname; }

      af::small<double, 6>
      unit_cell_parameters() const;

      cctbx::uctbx::unit_cell
      unit_cell() const
      {
        return cctbx::uctbx::unit_cell(unit_cell_parameters());
      }

      int
      n_datasets() const { return CMtz::MtzNsetsInXtal(ptr()); }

      af::shared<dataset>
      datasets() const;

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
