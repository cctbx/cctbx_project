#ifndef IOTBX_MTZ_OBJECT_H
#define IOTBX_MTZ_OBJECT_H

#include <cmtzlib.h>
#include <ccp4_array.h>
#include <ccp4_errno.h>
#include <cctbx/sgtbx/space_group.h>
#include <cctbx/hendrickson_lattman.h>
#include <boost/shared_ptr.hpp>

namespace iotbx { namespace mtz {

  namespace af = scitbx::af;

  af::shared<std::size_t>
  cmtz_struct_sizes();

  inline
  int
  ccp4_liberr_verbosity(int level)
  {
    return CCP4::ccp4_liberr_verbosity(level);
  }

  inline
  bool
  is_ccp4_nan(float const& datum)
  {
    return CCP4::ccp4_utils_isnan(
      reinterpret_cast<const union float_uint_uchar*>(&datum));
  }

  inline
  ccp4array_base*
  ccp4array_base_ptr(CMtz::MTZCOL* ptr)
  {
    ccp4_ptr* p = reinterpret_cast<ccp4_ptr*>(&ptr->ref);
    return reinterpret_cast<ccp4array_base*>(
      reinterpret_cast<ccp4_byteptr>(*p) - sizeof(ccp4array_base));
  }

  inline
  int
  column_array_size(CMtz::MTZCOL* ptr)
  {
    return ccp4array_base_ptr(ptr)->size;
  }

  inline
  int
  column_array_capacity(CMtz::MTZCOL* ptr)
  {
    return ccp4array_base_ptr(ptr)->capacity;
  }

  template <typename DataType>
  struct data_group
  {
    data_group() {}

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

  typedef data_group<int> integer_group;
  typedef data_group<double> real_group;
  typedef data_group<cctbx::hendrickson_lattman<> > hl_group;

  struct observations_group : real_group
  {
    observations_group() {}

    observations_group(bool anomalous_flag, std::size_t size)
    :
      real_group(anomalous_flag, size)
    {
      sigmas.reserve(size);
    }

    af::shared<double> sigmas;
  };

  struct complex_group
  {
    complex_group() {}

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

  class object
  {
    public:
      object();

      object(af::const_ref<int> const& n_datasets_for_each_crystal);

      object(const char* file_name);

      CMtz::MTZ*
      ptr() const { return ptr_.get(); }

      std::string
      title() const;

      object&
      set_title(const char* title, bool append=false);

      af::shared<std::string>
      history() const;

      object&
      add_history(af::const_ref<std::string> const& lines);

      object&
      add_history(const char* line);

      std::string
      space_group_name() const { return ptr()->mtzsymm.spcgrpname; }

      object&
      set_space_group_name(const char* name);

      int
      space_group_number() const { return ptr()->mtzsymm.spcgrp; }

      object&
      set_space_group_number(int number)
      {
        ptr()->mtzsymm.spcgrp = number;
        return *this;
      }

      std::string
      point_group_name() const { return ptr()->mtzsymm.pgname; }

      object&
      set_point_group_name(const char* name);

      char
      lattice_centring_type() const
      {
        return ptr()->mtzsymm.symtyp;
      }

      object&
      set_lattice_centring_type(char symbol)
      {
        ptr()->mtzsymm.symtyp = symbol;
        return *this;
      }

      int
      n_symmetry_matrices() const { return ptr()->mtzsymm.nsym; }

      cctbx::sgtbx::space_group
      space_group() const;

      object&
      set_space_group(cctbx::sgtbx::space_group const& space_group);

      void
      reserve(int capacity);

      void
      adjust_column_array_sizes(int new_nref);

      int
      n_batches() const { return CMtz::MtzNbat(ptr()); }

      af::shared<batch>
      batches() const;

      batch
      add_batch();

      void
      sort_batches()
      {
        CMtz::MTZ* p = ptr();
        p->batch = CMtz::sort_batches(p->batch, n_batches());
      }

      int
      n_reflections() const { return ptr()->nref; }

      af::tiny<double, 2>
      max_min_resolution() const;

      int
      n_crystals() const { return CMtz::MtzNxtal(ptr()); }

      int
      n_active_crystals() const { return CMtz::MtzNumActiveXtal(ptr()); }

      af::shared<crystal>
      crystals() const;

      crystal
      add_crystal(
        const char* name,
        const char* project_name,
        af::double6 const& unit_cell_parameters);

      crystal
      add_crystal(
        const char* name,
        const char* project_name,
        cctbx::uctbx::unit_cell const& unit_cell);

      bool
      has_column(const char* label) const;

      column
      get_column(const char* label) const;

      hkl_columns
      get_hkl_columns() const;

      integer_group
      extract_integers(
        const char* column_label);

      af::shared<int>
      extract_integers(
        af::const_ref<int> const& mtz_reflection_indices,
        const char* column_label);

      integer_group
      extract_integers_anomalous(
        const char* column_label_plus,
        const char* column_label_minus);

      real_group
      extract_reals(
        const char* column_label);

      af::shared<double>
      extract_reals(
        af::const_ref<int> const& mtz_reflection_indices,
        const char* column_label);

      real_group
      extract_reals_anomalous(
        const char* column_label_plus,
        const char* column_label_minus);

      hl_group
      extract_hendrickson_lattman(
        const char* column_label_a,
        const char* column_label_b,
        const char* column_label_c,
        const char* column_label_d);

      hl_group
      extract_hendrickson_lattman_anomalous(
        const char* column_label_a_plus,
        const char* column_label_b_plus,
        const char* column_label_c_plus,
        const char* column_label_d_plus,
        const char* column_label_a_minus,
        const char* column_label_b_minus,
        const char* column_label_c_minus,
        const char* column_label_d_minus);

      observations_group
      extract_observations(
        const char* column_label_data,
        const char* column_label_sigmas);

      observations_group
      extract_observations_anomalous(
        const char* column_label_data_plus,
        const char* column_label_sigmas_plus,
        const char* column_label_data_minus,
        const char* column_label_sigmas_minus);

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
        const char* column_label_d_sigmas);

      complex_group
      extract_complex(
        const char* column_label_ampl,
        const char* column_label_phi);

      complex_group
      extract_complex_anomalous(
        const char* column_label_ampl_plus,
        const char* column_label_phi_plus,
        const char* column_label_ampl_minus,
        const char* column_label_phi_minus);

      void
      write(const char* file_name);

    protected:
      boost::shared_ptr<CMtz::MTZ> ptr_;

      static void
      ptr_deleter(CMtz::MTZ* ptr);
  };

}} // namespace iotbx::mtz

#endif // IOTBX_MTZ_OBJECT_H
