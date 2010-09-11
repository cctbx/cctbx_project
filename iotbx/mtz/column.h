#ifndef IOTBX_MTZ_COLUMN_H
#define IOTBX_MTZ_COLUMN_H

#include <iotbx/mtz/dataset.h>

namespace iotbx { namespace mtz {

  //! Safe access to CMtz::MTZCOL* owned by an iotbx::mtz::dataset .
  /*! A column object contains (owns) an iotbx::mtz::dataset and
      an integer index into the list of columns owned by the
      dataset.
   */
  class column
  {
    public:
      //! Not available from Python.
      column() {}

      /*! \brief Initialization given a mtz_dataset and an integer
          index into the list of columns owned by the mtz_dataset.
       */
      /*! An exception is thrown if i_column is out of range.
       */
      column(dataset const& mtz_dataset, int i_column);

      //! The contained iotbx::mtz::dataset instance.
      dataset
      mtz_dataset() const { return mtz_dataset_; }

      //! Integer index into the list of columns owned by mtz_dataset() .
      int
      i_column() const { return i_column_; }

      //! Shorthand for mtz_dataset().mtz_crystal() .
      crystal
      mtz_crystal() const { return mtz_dataset_.mtz_crystal(); }

      //! Shorthand for mtz_crystal().mtz_object() .
      object
      mtz_object() const { return mtz_crystal().mtz_object(); }

      //! Raw C pointer. Not available from Python.
      /*! An exception is thrown if i_column() is out of range.
       */
      CMtz::MTZCOL*
      ptr() const;

      //! Read-only access.
      const char*
      label() const { return ptr()->label; }

      //! Write access.
      /*! An exception is thrown if the new_label is used already by
          any other column in the mtz_object().
       */
      column&
      set_label(const char* new_label);

      //! Read-only access.
      const char*
      type() const { return ptr()->type; }

      //! Write access.
      column&
      set_type(const char* new_type);

      //! Test.
      bool
      is_active() const { return (ptr()->active != 0); }

      //! Read-only access.
      const char*
      source() const
      {
#if defined(CCP4_MTZDATA) && CCP4_MTZDATA >= 20100630
        return ptr()->colsource;
#else
        return 0;
#endif
      }

      //! Write access.
      column&
      set_source(const char* new_source);

      //! Read-only access.
      const char*
      group_name() const
      {
#if defined(CCP4_MTZDATA) && CCP4_MTZDATA >= 20100630
        return ptr()->grpname;
#else
        return 0;
#endif
      }

      //! Write access.
      column&
      set_group_name(const char* new_group_name);

      //! Read-only access.
      const char*
      group_type() const
      {
#if defined(CCP4_MTZDATA) && CCP4_MTZDATA >= 20100630
        return ptr()->grptype;
#else
        return 0;
#endif
      }

      //! Write access.
      column&
      set_group_type(const char* new_group_type);

      //! Read-only access.
      int
      group_position() const
      {
#if defined(CCP4_MTZDATA) && CCP4_MTZDATA >= 20100630
        return ptr()->grpposn;
#else
        return -1;
#endif
      }

      //! Write access.
      column&
      set_group_position(int new_group_position);

      //! Read-only access.
      int
      array_size() const { return column_array_size(ptr()); }

      //! Read-only access.
      int
      array_capacity() const { return column_array_capacity(ptr()); }

      //! Read-only access.
      std::string
      path() const;

      //! Shorthand for mtz_object().get_column(label) .
      column
      get_other(const char* label) const
      {
        return mtz_object().get_column(label);
      }

      //! Read-only access. Not available from Python.
      float const&
      float_datum(int i) const { return ptr()->ref[i]; }

      //! Write access. Not available from Python.
      float&
      float_datum(int i)       { return ptr()->ref[i]; }

      //! Test. Not available from Python.
      bool
      is_ccp4_nan(int i) const { return mtz::is_ccp4_nan(float_datum(i)); }

      //! Read-only access. Not available from Python.
      int
      int_datum(int i) const { return static_cast<int>(float_datum(i)); }

      //! Count of all values in this column which are not "not-a-number."
      int
      n_valid_values() const;

      //! All values in this column which are not "not-a-number."
      af::shared<float>
      extract_valid_values() const;

      //! Bool selection of elements which are not "not-a-number."
      af::shared<bool>
      selection_valid() const;

      //! Copy of values in this column.
      /*! "not-a-number" elements are replaced by not_a_number_substitute.
       */
      af::shared<float>
      extract_values(float not_a_number_substitute=0) const;

      //! Copy values into this column.
      /*! selection_valid[i] == false indicates "not-a-number".
       */
      void
      set_values(
        af::const_ref<float> const& values,
        af::const_ref<bool> const& selection_valid) const;

      //! Copy values into this column.
      void
      set_values(
        af::const_ref<float> const& values) const
      {
        set_values(values, af::const_ref<bool>(0,0));
      }

      //! Write access given Miller indices and corresponding data.
      af::shared<int>
      set_reals(
        af::const_ref<cctbx::miller::index<> > const& miller_indices,
        af::const_ref<double> const& data);

      /*! \brief Write access given Miller indices and
          mtz_reflection_indices returned from a previous call
          of the other set_reals() overload.
       */
      void
      set_reals(
        af::const_ref<int> const& mtz_reflection_indices,
        af::const_ref<double> const& data);

    protected:
      dataset mtz_dataset_;
      int i_column_;
  };

  //! Result type for object::get_hkl_columns() .
  /*! Not available from Python.

      The family of object::extract_* functions provides safe
      access to Miller indices and corresponding data.
   */
  class hkl_columns : public af::tiny<column, 3>
  {
    public:
      //! Not available from Python.
      hkl_columns() {}

      //! Not available from Python.
      hkl_columns(column const& h, column const& k, column const& l)
      :
        af::tiny<column, 3>(h, k, l)
      {}

      //! Not available from Python.
      cctbx::miller::index<>
      get_miller_index(int iref) const
      {
        return cctbx::miller::index<>(
          this->elems[0].int_datum(iref),
          this->elems[1].int_datum(iref),
          this->elems[2].int_datum(iref));
      }

      //! Not available from Python.
      void
      replace_miller_index(
        int iref,
        cctbx::miller::index<> const& miller_index)
      {
        for(int i=0;i<3;i++) {
          this->elems[i].float_datum(iref)=static_cast<float>(miller_index[i]);
        }
      }
  };

}} // namespace iotbx::mtz

#endif // IOTBX_MTZ_COLUMN_H
