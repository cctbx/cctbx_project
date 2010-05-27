#ifndef IOTBX_MTZ_OBJECT_H
#define IOTBX_MTZ_OBJECT_H

#if defined(__sgi) && defined(__host_mips)
# include <sys/stat.h>
# include <math.h>
#endif

#include <cmtzlib.h>
#include <ccp4_array.h>
#include <ccp4_errno.h>
#include <cctbx/sgtbx/space_group.h>
#include <cctbx/hendrickson_lattman.h>
#include <boost/shared_ptr.hpp>

namespace iotbx {

//! C++ wrapper for CCP4 CMtz library.
/*! Please start with the documentation for iotbx::mtz::object .

    A simple Python example:
<pre>
      from iotbx import mtz
      import sys
      mtz_object = mtz.object(file_name=sys.argv[1])
      print mtz_object.title()
      for crystal in mtz_object.crystals():
        for dataset in crystal.datasets():
          for column in dataset.columns():
            print column.label()
</pre>
 */
namespace mtz {

  namespace af = scitbx::af;

  //! Returns a list of sizes of the C structs defined in mtzdata.h.
  /*! Used for consistency checking in the Python layer.
      See also:
        expected_cmtz_struct_sizes in
          $IOTBX_DIST/iotbx/mtz/__init__.py
   */
  af::shared<std::size_t>
  cmtz_struct_sizes();

  //! Thin wrapper, mainly to facilitate Python bindings.
  inline
  int
  ccp4_liberr_verbosity(int level)
  {
    return CCP4::ccp4_liberr_verbosity(level);
  }

  //! Wrapper for CCP4::ccp4_utils_isnan() .
  inline
  bool
  is_ccp4_nan(float const& datum)
  {
    return CCP4::ccp4_utils_isnan(
      reinterpret_cast<const union float_uint_uchar*>(&datum));
  }

  //! Complex sequence of type casts.
  inline
  ccp4array_base*
  ccp4array_base_ptr(CMtz::MTZCOL* ptr)
  {
    ccp4_ptr* p = reinterpret_cast<ccp4_ptr*>(&ptr->ref);
    return reinterpret_cast<ccp4array_base*>(
      reinterpret_cast<ccp4_byteptr>(*p) - sizeof(ccp4array_base));
  }

  //! Access to ccp4array_base->size .
  inline
  int
  column_array_size(CMtz::MTZCOL* ptr)
  {
    return ccp4array_base_ptr(ptr)->size;
  }

  //! Access to ccp4array_base->capacity .
  inline
  int
  column_array_capacity(CMtz::MTZCOL* ptr)
  {
    return ccp4array_base_ptr(ptr)->capacity;
  }

  //! Result type for family of object::extract_* functions.
  template <typename DataType>
  struct data_group
  {
    //! Not available from Python.
    data_group() {}

    //! Not available from Python.
    data_group(bool anomalous_flag_, std::size_t size)
    :
      anomalous_flag(anomalous_flag_)
    {
      mtz_reflection_indices.reserve(size);
      indices.reserve(size);
      data.reserve(size);
    }

    bool anomalous_flag;
    af::shared<int> mtz_reflection_indices;
    af::shared<cctbx::miller::index<> > indices;
    af::shared<DataType> data;
  };

  //! Result type for family of object::extract_* functions.
  typedef data_group<int> integer_group;
  //! Result type for family of object::extract_* functions.
  typedef data_group<double> real_group;
  //! Result type for family of object::extract_* functions.
  typedef data_group<cctbx::hendrickson_lattman<> > hl_group;

  //! Result type for family of object::extract_* functions.
  struct observations_group : real_group
  {
    //! Not available from Python.
    observations_group() {}

    //! Not available from Python.
    observations_group(bool anomalous_flag, std::size_t size)
    :
      real_group(anomalous_flag, size)
    {
      sigmas.reserve(size);
    }

    af::shared<double> sigmas;
  };

  //! Result type for family of object::extract_* functions.
  struct complex_group
  {
    //! Not available from Python.
    complex_group() {}

    //! Not available from Python.
    complex_group(bool anomalous_flag_, std::size_t size)
    :
      anomalous_flag(anomalous_flag_)
    {
      mtz_reflection_indices.reserve(size);
      indices.reserve(size);
      data.reserve(size);
    }

    bool anomalous_flag;
    af::shared<int> mtz_reflection_indices;
    af::shared<cctbx::miller::index<> > indices;
    af::shared<std::complex<double> > data;
  };

  class column;
  class hkl_columns;
  class dataset;
  class crystal;
  class batch;

  //! Wrapper for CMtz::MTZ* .
  /*! The life-time of the wrapped CMtz::MTZ "object" is
      managed by a boost::shared_ptr<CMtz::MTZ>. I.e.
      CMtz::MtzFree() is called automatically if object instances
      go out of scope. Access to the raw C ptr() is provided
      for completeness. However, all important functionality
      is available through safer interfaces which should
      be preferred.

      All safe interfaces are availabe in Python.
      Interfaces marked as "Not available from Python"
      should also not normally be used from C++
      (with the exception of default constructors).

      Basic Python regression tests:
        $IOTBX_DIST/include/iotbx/mtz/tst_ext.py

      Advanced Python regression tests:
        $IOTBX_DIST/iotbx/mtz/tst.py

      The show_summary() method injected in
        $IOTBX_DIST/iotbx/mtz/__init__.py
      demonstrates how to traverse the data hierarchy of
      object, crystal, dataset, column.

      Copying an object is inexpensive because the only
      data member is a boost::shared_ptr<CMtz::MTZ>.

      This wrapper was tested only with refs_in_memory = true.

      See also:
        http://www.ccp4.ac.uk/dist/html/C_library/cmtz_page.html
   */
  class object
  {
    public:
      //! Not available from Python.
      object();

      //! Wrapper for CMtz::MtzGet() .
      object(const char* file_name);

      //! Access to raw C pointer. Not available in Python.
      CMtz::MTZ*
      ptr() const { return ptr_.get(); }

      //! Read-only access.
      std::string
      title() const;

      //! Write access.
      object&
      set_title(const char* title, bool append=false);

      //! Read-only access.
      af::shared<std::string>
      history() const;

      //! Write access. Adds multiple lines.
      object&
      add_history(af::const_ref<std::string> const& lines);

      //! Write access. Adds one line.
      object&
      add_history(const char* line);

      //! Read-only access.
      std::string
      space_group_name() const { return ptr()->mtzsymm.spcgrpname; }

      //! Write access.
      object&
      set_space_group_name(const char* name);

      //! Read-only access.
      int
      space_group_number() const { return ptr()->mtzsymm.spcgrp; }

      //! Write access.
      object&
      set_space_group_number(int number)
      {
        ptr()->mtzsymm.spcgrp = number;
        return *this;
      }

      //! Read-only access.
      std::string
      point_group_name() const { return ptr()->mtzsymm.pgname; }

      //! Write access.
      object&
      set_point_group_name(const char* name);

      //! Read-only access.
      char
      lattice_centring_type() const
      {
        return ptr()->mtzsymm.symtyp;
      }

      //! Write access.
      object&
      set_lattice_centring_type(char symbol)
      {
        ptr()->mtzsymm.symtyp = symbol;
        return *this;
      }

      //! Read-only access.
      int
      n_symmetry_matrices() const { return ptr()->mtzsymm.nsym; }

      //! Read-only access.
      char
      space_group_confidence() const
      {
#if defined(CCP4_MTZDATA) && CCP4_MTZDATA > 20100418
        return ptr()->mtzsymm.spg_confidence;
#else
        return '\0';
#endif
      }

      //! Read-only access.
      cctbx::sgtbx::space_group
      space_group() const;

      //! Write access.
      object&
      set_space_group(cctbx::sgtbx::space_group const& space_group);

      //! Pre-allocates memory for reflection arrays.
      void
      reserve(int capacity);

      /*! \brief Allocates memory for new_nref and fill all new data
          slots with "not-a-number" bit patterns.
       */
      /*! Not available in Python.
       */
      void
      adjust_column_array_sizes(int new_nref);

      //! Read-only access.
      int
      n_batches() const { return CMtz::MtzNbat(ptr()); }

      //! Read-only access.
      af::shared<batch>
      batches() const;

      //! Adds a new batch to this object (with default values).
      batch
      add_batch();

      //! Wrapper for CMtz::sort_batches() .
      void
      sort_batches()
      {
        CMtz::MTZ* p = ptr();
        p->batch = CMtz::sort_batches(p->batch, n_batches());
      }

      //! Read-only access.
      int
      n_reflections() const { return ptr()->nref; }

      //! Maximum and minimum resolution, typically in Angstroms.
      af::tiny<double, 2>
      max_min_resolution() const;

      //! Read-only access.
      int
      n_crystals() const { return CMtz::MtzNxtal(ptr()); }

      //! Read-only access.
      int
      n_active_crystals() const { return CMtz::MtzNumActiveXtal(ptr()); }

      //! List of all crystals owned by this object.
      af::shared<crystal>
      crystals() const;

      //! Adds a crystal to this object.
      crystal
      add_crystal(
        const char* name,
        const char* project_name,
        af::double6 const& unit_cell_parameters);

      //! Adds a crystal to this object.
      crystal
      add_crystal(
        const char* name,
        const char* project_name,
        cctbx::uctbx::unit_cell const& unit_cell);

      //! Test.
      bool
      has_crystal(const char* name) const;

      //! Test.
      bool
      has_column(const char* label) const;

      //! Retrieves a column owned by this object.
      /*! An exception is thrown if the column label is unknown.
       */
      column
      get_column(const char* label) const;

      //! Retrieves columns H, K, L.
      /*! An exception is thrown if any of the columns does not exist.

          Not available in Python.
       */
      hkl_columns
      get_hkl_columns() const;

      //! Copies Miller indices from columns H, K, L.
      af::shared<cctbx::miller::index<> >
      extract_miller_indices() const;

      //! Overwrites Miller indices in columns H, K, L.
      /*! The miller_indices.size() must be equal to n_reflections().
       */
      void
      replace_miller_indices(
        af::const_ref<cctbx::miller::index<> > const& miller_indices);

      //! Read-only access.
      integer_group
      extract_integers(
        const char* column_label) const;

      //! Read-only access.
      af::shared<int>
      extract_integers(
        af::const_ref<int> const& mtz_reflection_indices,
        const char* column_label) const;

      //! Read-only access.
      integer_group
      extract_integers_anomalous(
        const char* column_label_plus,
        const char* column_label_minus) const;

      //! Read-only access.
      real_group
      extract_reals(
        const char* column_label) const;

      //! Read-only access.
      af::shared<double>
      extract_reals(
        af::const_ref<int> const& mtz_reflection_indices,
        const char* column_label) const;

      //! Read-only access.
      real_group
      extract_reals_anomalous(
        const char* column_label_plus,
        const char* column_label_minus) const;

      //! Read-only access.
      hl_group
      extract_hendrickson_lattman(
        const char* column_label_a,
        const char* column_label_b,
        const char* column_label_c,
        const char* column_label_d) const;

      //! Read-only access.
      hl_group
      extract_hendrickson_lattman_ab_only(
        const char* column_label_a,
        const char* column_label_b) const;

      //! Read-only access.
      hl_group
      extract_hendrickson_lattman_anomalous(
        const char* column_label_a_plus,
        const char* column_label_b_plus,
        const char* column_label_c_plus,
        const char* column_label_d_plus,
        const char* column_label_a_minus,
        const char* column_label_b_minus,
        const char* column_label_c_minus,
        const char* column_label_d_minus) const;

      //! Read-only access.
      hl_group
      extract_hendrickson_lattman_anomalous_ab_only(
        const char* column_label_a_plus,
        const char* column_label_b_plus,
        const char* column_label_a_minus,
        const char* column_label_b_minus) const;

      //! Read-only access.
      observations_group
      extract_observations(
        const char* column_label_data,
        const char* column_label_sigmas) const;

      //! Read-only access.
      observations_group
      extract_observations_anomalous(
        const char* column_label_data_plus,
        const char* column_label_sigmas_plus,
        const char* column_label_data_minus,
        const char* column_label_sigmas_minus) const;

      //! Read-only access.
      /*! http://www.ccp4.ac.uk/dist/html/mtzMADmod.html
            F(+) = F + 0.5*D
            F(-) = F - 0.5*D
            SIGF(+) = sqrt( SIGF**2 + 0.25*SIGD**2 )
            SIGF(-) = SIGF(+)
       */
      observations_group
      extract_delta_anomalous(
        const char* column_label_f_data,
        const char* column_label_f_sigmas,
        const char* column_label_d_data,
        const char* column_label_d_sigmas) const;

      //! Read-only access.
      complex_group
      extract_complex(
        const char* column_label_ampl,
        const char* column_label_phi) const;

      //! Read-only access.
      complex_group
      extract_complex_anomalous(
        const char* column_label_ampl_plus,
        const char* column_label_phi_plus,
        const char* column_label_ampl_minus,
        const char* column_label_phi_minus) const;

      //! Wrapper for CMtz::MtzPut() .
      void
      write(const char* file_name);

      //! Read-only access to "not-a-number" value.
      /*! Not available in Python.
       */
      const float&
      not_a_number_value() { return not_a_number_value_.f; }

    protected:
      boost::shared_ptr<CMtz::MTZ> ptr_;

      static void
      ptr_deleter(CMtz::MTZ* ptr);

      union float_uint_uchar not_a_number_value_;

      void
      init_not_a_number_value();
  };

}} // namespace iotbx::mtz

#endif // IOTBX_MTZ_OBJECT_H
