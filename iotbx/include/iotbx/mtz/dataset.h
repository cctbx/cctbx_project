#ifndef IOTBX_MTZ_DATASET_H
#define IOTBX_MTZ_DATASET_H

#include <iotbx/mtz/crystal.h>

namespace iotbx { namespace mtz {

  class dataset
  {
    public:
      dataset() {}

      dataset(crystal const& mtz_crystal, int i_dataset)
      :
        mtz_crystal_(mtz_crystal),
        i_dataset_(i_dataset)
      {
        CCTBX_ASSERT(i_dataset >= 0);
        CCTBX_ASSERT(i_dataset < mtz_crystal.n_datasets());
      }

      crystal
      mtz_crystal() const { return mtz_crystal_; }

      int
      i_dataset() const { return i_dataset_; }

      object
      mtz_object() const { return mtz_crystal_.mtz_object(); }

      CMtz::MTZSET*
      ptr() const
      {
        CCTBX_ASSERT(mtz_crystal_.n_datasets() > i_dataset_);
        return CMtz::MtzIsetInXtal(mtz_crystal_.ptr(), i_dataset_);
      }

      int
      id() const { return ptr()->setid; }

      dataset&
      set_id(int id)
      {
        if (ptr()->setid != id) {
          CMtz::MTZ* p = mtz_crystal().mtz_object().ptr();
          CCTBX_ASSERT(p->refs_in_memory);
          for(int i=0;i<p->nxtal;i++) {
            for(int j=0;j<p->xtal[i]->nset;j++) {
              CCTBX_ASSERT(p->xtal[i]->set[j]->setid != id);
            }
          }
          ptr()->setid = id;
        }
        return *this;
      }

      const char*
      name() const { return ptr()->dname; }

      float
      wavelength() const { return ptr()->wavelength; }

      int
      n_columns() const { return CMtz::MtzNcolsInSet(ptr()); }

      af::shared<column>
      columns() const;

      inline
      column
      add_column(
        const char *label,
        const char *type);

    protected:
      crystal mtz_crystal_;
      int i_dataset_;
  };

  inline
  af::shared<dataset>
  crystal::datasets() const
  {
    af::shared<dataset> result((af::reserve(n_datasets())));
    for(int i_dataset=0;i_dataset<n_datasets();i_dataset++) {
      result.push_back(dataset(*this, i_dataset));
    }
    return result;
  }

  inline
  dataset
  crystal::add_dataset(
    const char *name,
    double wavelength)
  {
    int i_dataset = n_datasets();
    CMtz::MTZSET* dataset_ptr = CMtz::MtzAddDataset(
      mtz_object().ptr(), ptr(), name, static_cast<float>(wavelength));
    CCTBX_ASSERT(dataset_ptr != 0);
    CCTBX_ASSERT(n_datasets() == i_dataset+1);
    dataset result(*this, i_dataset);
    CCTBX_ASSERT(result.ptr() == dataset_ptr);
    return result;
  }

}} // namespace iotbx::mtz

#endif // IOTBX_MTZ_DATASET_H
