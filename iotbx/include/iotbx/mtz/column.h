#ifndef IOTBX_MTZ_COLUMN_H
#define IOTBX_MTZ_COLUMN_H

#include <iotbx/mtz/dataset.h>

namespace iotbx { namespace mtz {

  class column
  {
    public:
      column() {}

      column(dataset const& mtz_dataset, int i_column)
      :
        mtz_dataset_(mtz_dataset),
        i_column_(i_column)
      {
        CCTBX_ASSERT(i_column >= 0);
        CCTBX_ASSERT(i_column < mtz_dataset.n_columns());
      }

      dataset
      mtz_dataset() const { return mtz_dataset_; }

      int
      i_column() const { return i_column_; }

      crystal
      mtz_crystal() const { return mtz_dataset().mtz_crystal(); }

      object
      mtz_object() const { return mtz_crystal().mtz_object(); }

      CMtz::MTZCOL*
      ptr() const
      {
        CCTBX_ASSERT(mtz_dataset_.n_columns() > i_column_);
        return CMtz::MtzIcolInSet(mtz_dataset_.ptr(), i_column_);
      }

      const char*
      label() const { return ptr()->label; }

      const char*
      type() const { return ptr()->type; }

      bool
      is_active() const { return (ptr()->active != 0); }

      std::string
      path() const
      {
        boost::shared_ptr<char>
          p(CMtz::MtzColPath(mtz_object().ptr(), ptr()), free);
        return std::string(p.get());
      }

      column
      lookup_other(const char* label) const
      {
        return mtz_object().lookup_column(label);
      }

      float const&
      float_datum(int i) const { return ptr()->ref[i]; }

      bool
      is_ccp4_nan(int i) const { return mtz::is_ccp4_nan(float_datum(i)); }

      int
      int_datum(int i) const { return static_cast<int>(float_datum(i)); }

      af::shared<cctbx::miller::index<> >
      valid_indices() const
      {
        int n_reflections = mtz_object().n_reflections();
        af::shared<cctbx::miller::index<> >
          result((af::reserve(n_reflections)));
        column h(lookup_other("H"));
        column k(lookup_other("K"));
        column l(lookup_other("L"));
        for(int i=0;i<n_reflections;i++) {
          if (!is_ccp4_nan(i)) {
            result.push_back(cctbx::miller::index<>(
              h.int_datum(i), k.int_datum(i), l.int_datum(i)));
          }
        }
        return result;
      }

      af::shared<double>
      valid_values()
      {
        int n_reflections = mtz_object().n_reflections();
        af::shared<double> result((af::reserve(n_reflections)));
        for(int i=0;i<n_reflections;i++) {
          if (!is_ccp4_nan(i)) {
            result.push_back(float_datum(i));
          }
        }
        return result;
      }

      af::shared<int>
      valid_integers()
      {
        if (   strcmp(type(), "H") != 0
            && strcmp(type(), "B") != 0
            && strcmp(type(), "Y") != 0
            && strcmp(type(), "I") != 0) {
          throw cctbx::error(
            std::string("Not a MTZ integer column: ") + label());
        }
        int n_reflections = mtz_object().n_reflections();
        af::shared<int> result((af::reserve(n_reflections)));
        for(int i=0;i<n_reflections;i++) {
          if (!is_ccp4_nan(i)) {
            result.push_back(int_datum(i));
          }
        }
        return result;
      }

    protected:
      dataset mtz_dataset_;
      int i_column_;
  };

  inline
  af::shared<column>
  dataset::columns() const
  {
    af::shared<column> result((af::reserve(n_columns())));
    for(int i_column=0;i_column<n_columns();i_column++) {
      result.push_back(column(*this, i_column));
    }
    return result;
  }

  inline
  column
  object::lookup_column(const char* label) const
  {
    for(int i_crystal=0;i_crystal<n_crystals();i_crystal++) {
      crystal x(*this, i_crystal);
      for(int i_dataset=0;i_dataset<x.n_datasets();i_dataset++) {
        dataset s(x, i_dataset);
        for(int i_column=0;i_column<s.n_columns();i_column++) {
          column c(s, i_column);
          if (CMtz::MtzPathMatch(label, c.label())) {
            return c;
          }
        }
      }
    }
    throw cctbx::error(std::string("Unknown MTZ column label: ") + label);
  }

}} // namespace iotbx::mtz

#endif // IOTBX_MTZ_COLUMN_H
