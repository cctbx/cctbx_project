#include <iotbx/mtz/batch.h>
#include <iotbx/mtz/column.h>

namespace iotbx { namespace mtz {

  dataset::dataset(crystal const& mtz_crystal, int i_dataset)
  :
    mtz_crystal_(mtz_crystal),
    i_dataset_(i_dataset)
  {
    CCTBX_ASSERT(i_dataset >= 0);
    CCTBX_ASSERT(i_dataset < mtz_crystal.n_datasets());
  }

  CMtz::MTZSET*
  dataset::ptr() const
  {
    CCTBX_ASSERT(mtz_crystal_.n_datasets() > i_dataset_);
    return CMtz::MtzIsetInXtal(mtz_crystal_.ptr(), i_dataset_);
  }

  dataset&
  dataset::set_id(int id)
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

  af::shared<batch>
  dataset::batches() const
  {
    af::shared<batch> result;
    int n_orig_bat = mtz_object().ptr()->n_orig_bat;
    CMtz::MTZBAT* batch_ptr = mtz_object().ptr()->batch;
    for(int i_batch=0;
            batch_ptr!=0;
            i_batch++, batch_ptr=batch_ptr->next) {
      if (   batch_ptr->nbsetid == id()
          && i_batch >= n_orig_bat) {
        result.push_back(batch(mtz_object(), i_batch));
      }
    }
    return result;
  }

  batch
  dataset::add_batch()
  {
    batch result = mtz_object().add_batch();
    result.set_nbsetid(id());
    return result;
  }

  af::shared<column>
  dataset::columns() const
  {
    af::shared<column> result((af::reserve(n_columns())));
    for(int i_column=0;i_column<n_columns();i_column++) {
      result.push_back(column(*this, i_column));
    }
    return result;
  }

  column
  dataset::add_column(
    const char *label,
    const char *type)
  {
    CMtz::MTZCOL* column_ptr = 0;
    CCTBX_ASSERT(label != 0);
    CCTBX_ASSERT(strlen(label) < sizeof(column_ptr->label));
    CCTBX_ASSERT(type != 0);
    CCTBX_ASSERT(strlen(type) < sizeof(column_ptr->type));
    CCTBX_ASSERT(!mtz_object().has_column(label));
    int i_column = n_columns();
    column_ptr = CMtz::MtzAddColumn(mtz_object().ptr(), ptr(), label, type);
    CCTBX_ASSERT(column_ptr != 0);
    CCTBX_ASSERT(n_columns() == i_column+1);
    column result(*this, i_column);
    CCTBX_ASSERT(result.ptr() == column_ptr);
    return result;
  }

}} // namespace iotbx::mtz
