#ifndef IOTBX_MTZ_CRYSTAL_H
#define IOTBX_MTZ_CRYSTAL_H

#include <iotbx/mtz/object.h>

namespace iotbx { namespace mtz {

  class crystal
  {
    public:
      crystal() {}

      crystal(mtz::object const& object, int i_crystal)
      :
        object_(object),
        i_crystal_(i_crystal)
      {
        CCTBX_ASSERT(i_crystal >= 0);
        CCTBX_ASSERT(i_crystal < object.n_crystals());
      }

      mtz::object
      object() const { return object_; }

      int
      i_crystal() const { return i_crystal_; }

    protected:
      mtz::object object_;
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
