#ifndef IOTBX_MTZ_DATASET_H
#define IOTBX_MTZ_DATASET_H

#include <iotbx/mtz/crystal.h>

namespace iotbx { namespace mtz {

  class dataset
  {
    public:
      dataset() {}

      dataset(crystal const& mtz_crystal, int i_dataset);

      crystal
      mtz_crystal() const { return mtz_crystal_; }

      int
      i_dataset() const { return i_dataset_; }

      object
      mtz_object() const { return mtz_crystal_.mtz_object(); }

      CMtz::MTZSET*
      ptr() const;

      int
      id() const { return ptr()->setid; }

      dataset&
      set_id(int id);

      const char*
      name() const { return ptr()->dname; }

      float
      wavelength() const { return ptr()->wavelength; }

      int
      n_batches() const
      {
        return CMtz::MtzNbatchesInSet(mtz_object().ptr(), ptr());
      }

      af::shared<batch>
      batches() const;

      batch
      add_batch();

      int
      n_columns() const { return CMtz::MtzNcolsInSet(ptr()); }

      af::shared<column>
      columns() const;

      column
      add_column(
        const char *label,
        const char *type);

    protected:
      crystal mtz_crystal_;
      int i_dataset_;
  };

}} // namespace iotbx::mtz

#endif // IOTBX_MTZ_DATASET_H
