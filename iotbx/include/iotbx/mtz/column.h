#ifndef IOTBX_MTZ_COLUMN_H
#define IOTBX_MTZ_COLUMN_H

#include <iotbx/mtz/dataset.h>

namespace iotbx { namespace mtz {

  class column
  {
    public:
      column() {}

      column(dataset const& mtz_dataset, int i_column);

      dataset
      mtz_dataset() const { return mtz_dataset_; }

      int
      i_column() const { return i_column_; }

      crystal
      mtz_crystal() const { return mtz_dataset_.mtz_crystal(); }

      object
      mtz_object() const { return mtz_crystal().mtz_object(); }

      CMtz::MTZCOL*
      ptr() const;

      const char*
      label() const { return ptr()->label; }

      const char*
      type() const { return ptr()->type; }

      bool
      is_active() const { return (ptr()->active != 0); }

      int
      array_size() const { return column_array_size(ptr()); }

      int
      array_capacity() const { return column_array_capacity(ptr()); }

      std::string
      path() const;

      column
      get_other(const char* label) const
      {
        return mtz_object().get_column(label);
      }

      float const&
      float_datum(int i) const { return ptr()->ref[i]; }

      bool
      is_ccp4_nan(int i) const { return mtz::is_ccp4_nan(float_datum(i)); }

      int
      int_datum(int i) const { return static_cast<int>(float_datum(i)); }

      int
      n_valid_values() const;

      af::shared<int>
      set_reals(
        af::const_ref<cctbx::miller::index<> > const& miller_indices,
        af::const_ref<double> const& data);

      void
      set_reals(
        af::const_ref<int> const& mtz_reflection_indices,
        af::const_ref<double> const& data);

    protected:
      dataset mtz_dataset_;
      int i_column_;
  };

  class hkl_columns : public af::tiny<column, 3>
  {
    public:
      hkl_columns() {}

      hkl_columns(column const& h, column const& k, column const& l)
      :
        af::tiny<column, 3>(h, k, l)
      {}

      cctbx::miller::index<>
      get_miller_index(int iref) const
      {
        return cctbx::miller::index<>(
          this->elems[0].int_datum(iref),
          this->elems[1].int_datum(iref),
          this->elems[2].int_datum(iref));
      }
  };

}} // namespace iotbx::mtz

#endif // IOTBX_MTZ_COLUMN_H
