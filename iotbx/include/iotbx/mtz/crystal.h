#ifndef IOTBX_MTZ_CRYSTAL_H
#define IOTBX_MTZ_CRYSTAL_H

#include <iotbx/mtz/object.h>

namespace iotbx { namespace mtz {

  class crystal
  {
    public:
      crystal() {}

      crystal(object const& mtz_object, int i_crystal)
      :
        mtz_object_(mtz_object),
        i_crystal_(i_crystal)
      {
        CCTBX_ASSERT(i_crystal >= 0);
        CCTBX_ASSERT(i_crystal < mtz_object.n_crystals());
      }

      object
      mtz_object() const { return mtz_object_; }

      int
      i_crystal() const { return i_crystal_; }

      CMtz::MTZXTAL*
      ptr() const
      {
        CCTBX_ASSERT(mtz_object_.n_crystals() > i_crystal_);
        return CMtz::MtzIxtal(mtz_object_.ptr(), i_crystal_);
      }

      int
      id() const { return ptr()->xtalid; }

      const char*
      name() const { return ptr()->xname; }

      const char*
      project_name() const { return ptr()->pname; }

      af::small<double, 6>
      unit_cell_parameters() const
      {
        af::small<double, 6> result;
        float* cell = ptr()->cell;
        for(std::size_t i=0;i<6;i++) result.push_back(cell[i]);
        return result;
      }

      cctbx::uctbx::unit_cell
      unit_cell() const
      {
        return cctbx::uctbx::unit_cell(unit_cell_parameters());
      }

      int
      n_datasets() const { return CMtz::MtzNsetsInXtal(ptr()); }

      af::shared<dataset>
      datasets() const;

    protected:
      object mtz_object_;
      int i_crystal_;
  };

  inline
  af::shared<crystal>
  object::crystals() const
  {
    af::shared<crystal> result((af::reserve(n_crystals())));
    for(int i_crystal=0;i_crystal<n_crystals();i_crystal++) {
      result.push_back(crystal(*this, i_crystal));
    }
    return result;
  }

}} // namespace iotbx::mtz

#endif // IOTBX_MTZ_CRYSTAL_H
