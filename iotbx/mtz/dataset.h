#ifndef IOTBX_MTZ_DATASET_H
#define IOTBX_MTZ_DATASET_H

#include <iotbx/mtz/crystal.h>

namespace iotbx { namespace mtz {

  //! Safe access to CMtz::MTZSET* owned by an iotbx::mtz::crystal .
  /*! A dataset object contains (owns) an iotbx::mtz::crystal and
      an integer index into the list of datasets owned by the
      crystal.
   */
  class dataset
  {
    public:
      //! Not available from Python.
      dataset() {}

      /*! \brief Initialization given a mtz_crystal and an integer
          index into the list of datasets owned by the mtz_crystal.
       */
      /*! An exception is thrown if i_dataset is out of range.
       */
      dataset(crystal const& mtz_crystal, int i_dataset);

      //! The contained iotbx::mtz::crystal instance.
      crystal
      mtz_crystal() const { return mtz_crystal_; }

      //! Integer index into the list of datasets owned by mtz_crystal() .
      int
      i_dataset() const { return i_dataset_; }

      //! Shorthand for mtz_crystal().mtz_object() .
      object
      mtz_object() const { return mtz_crystal_.mtz_object(); }

      //! Raw C pointer. Not available from Python.
      /*! An exception is thrown if i_dataset() is out of range.
       */
      CMtz::MTZSET*
      ptr() const;

      //! Read-only access.
      int
      id() const { return ptr()->setid; }

      //! Write access.
      dataset&
      set_id(int id);

      //! Read-only access.
      const char*
      name() const { return ptr()->dname; }

      //! Write access.
      /*! An exception is thrown if the new_name is used already by
          another dataset in the mtz_crystal().
       */
      dataset&
      set_name(const char* new_name);

      //! Read-only access.
      float
      wavelength() const { return ptr()->wavelength; }

      //! Write access.
      dataset&
      set_wavelength(float new_wavelength)
      {
        ptr()->wavelength = new_wavelength;
        return *this;
      }

      //! Read-only access.
      int
      n_batches() const
      {
        return CMtz::MtzNbatchesInSet(mtz_object().ptr(), ptr());
      }

      //! Read-only access.
      af::shared<batch>
      batches() const;

      //! Adds a batch to this dataset.
      batch
      add_batch();

      //! Read-only access.
      int
      n_columns() const { return CMtz::MtzNcolsInSet(ptr()); }

      //! Read-only access.
      af::shared<column>
      columns() const;

      //! Adds a column to this dataset.
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
