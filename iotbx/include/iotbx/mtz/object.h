#ifndef IOTBX_MTZ_OBJECT_H
#define IOTBX_MTZ_OBJECT_H

#include <cmtzlib.h>
#include <ccp4_array.h>
#include <ccp4_errno.h>
#include <cctbx/sgtbx/space_group.h>

namespace iotbx { namespace mtz {

  namespace af = scitbx::af;

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
      indices.reserve(size);
      data.reserve(size);
    }

    bool anomalous_flag;
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
      indices.reserve(size);
      data.reserve(size);
    }

    bool anomalous_flag;
    af::shared<cctbx::miller::index<> > indices;
    af::shared<std::complex<double> > data;
  };

  class column;
  class hkl_columns;
  class dataset;
  class crystal;

  class object
  {
    public:
      object()
      :
        ptr_(CMtz::MtzMalloc(0, 0), ptr_deleter)
      {
        if (ptr_.get() == 0) throw cctbx::error("MtzMalloc failed.");
        ptr_->refs_in_memory = true;
      }

      object(af::const_ref<int> const& n_datasets_for_each_crystal)
      :
        ptr_(CMtz::MtzMalloc(
            n_datasets_for_each_crystal.size(),
            const_cast<int*>(&*n_datasets_for_each_crystal.begin())),
          ptr_deleter)
      {
        if (ptr_.get() == 0) throw cctbx::error("MtzMalloc failed.");
        ptr_->refs_in_memory = true;
      }

      object(const char* file_name)
      :
        ptr_(CMtz::MtzGet(file_name, true), ptr_deleter)
      {
        if (ptr_.get() == 0) {
          throw cctbx::error(std::string("MTZ file read error: ") + file_name);
        }
      }

      CMtz::MTZ*
      ptr() const { return ptr_.get(); }

      std::string
      title() const
      {
        char result[sizeof(ptr()->title)];
        int title_length = CMtz::ccp4_lrtitl(ptr(), result);
        return std::string(result, title_length);
      }

      object&
      set_title(const char* title, bool append=false)
      {
        int set_title_success = CMtz::ccp4_lwtitl(ptr(), title, append);
        CCTBX_ASSERT(set_title_success);
        ptr()->title[sizeof(ptr()->title)-1] = '\0';
        return *this;
      }

      af::shared<std::string>
      history() const
      {
        CMtz::MTZ* p = ptr();
        af::shared<std::string> result((af::reserve(p->histlines)));
        for(int i=0;i<p->histlines;i++) {
          const char* line = p->hist+MTZRECORDLENGTH*i;
          int j = 0;
          for(;j<MTZRECORDLENGTH;j++) {
            if (line[j] == '\0') break;
          }
          result.push_back(std::string(line, j));
        }
        return result;
      }

      object&
      add_history(af::const_ref<std::string> const& lines)
      {
        boost::shared_ptr<char> buffer(
          CMtz::MtzCallocHist(lines.size()), CMtz::MtzFreeHist);
        for(std::size_t i=0;i<lines.size();i++) {
          strncpy(
            buffer.get()+i*MTZRECORDLENGTH,
            lines[i].c_str(),
            std::min(
              static_cast<std::size_t>(MTZRECORDLENGTH),
              lines[i].size()));
        }
        int add_history_success = CMtz::MtzAddHistory(
          ptr(),
          reinterpret_cast<char (*)[MTZRECORDLENGTH]>(buffer.get()),
          lines.size());
        CCTBX_ASSERT(add_history_success);
        return *this;
      }

      object&
      add_history(const char* line)
      {
        return add_history(af::tiny<std::string, 1>(line).const_ref());
      }

      std::string
      space_group_name() const { return ptr()->mtzsymm.spcgrpname; }

      object&
      set_space_group_name(const char* name)
      {
        char* target = ptr()->mtzsymm.spcgrpname;
        const unsigned target_size = sizeof(ptr()->mtzsymm.spcgrpname);
        strncpy(target, name, target_size-1);
        target[target_size-1] = '\0';
        return *this;
      }

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
      set_point_group_name(const char* name)
      {
        char* target = ptr()->mtzsymm.pgname;
        const unsigned target_size = sizeof(ptr()->mtzsymm.pgname);
        strncpy(target, name, target_size-1);
        target[target_size-1] = '\0';
        return *this;
      }

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

      cctbx::sgtbx::space_group
      space_group() const
      {
        CMtz::MTZ* p = ptr();
        cctbx::sgtbx::space_group result;
        scitbx::mat3<double> rm;
        scitbx::vec3<double> tv;
        for(int im=0;im<p->mtzsymm.nsym;im++) {
          for (int ir=0;ir<3;ir++) {
            for (int ic=0;ic<3;ic++) {
              rm(ir,ic) = p->mtzsymm.sym[im][ir][ic];
            }
            tv[ir] = p->mtzsymm.sym[im][ir][3];
          }
          result.expand_smx(cctbx::sgtbx::rt_mx(rm, tv));
        }
        return result;
      }

      object&
      set_space_group(cctbx::sgtbx::space_group const& space_group)
      {
        CMtz::MTZ* p = ptr();
        CCTBX_ASSERT(sizeof(p->mtzsymm.sym) / sizeof(*p->mtzsymm.sym)
                  >= space_group.order_z());
        p->mtzsymm.nsymp = static_cast<int>(space_group.order_p());
        p->mtzsymm.nsym = static_cast<int>(space_group.order_z());
        for (int im=0;im<p->mtzsymm.nsym;im++) {
          cctbx::sgtbx::rt_mx sm = space_group(im).mod_positive();
          cctbx::sgtbx::rot_mx rm = sm.r();
          cctbx::sgtbx::tr_vec tv = sm.t();
          scitbx::mat3<int> rm_num = rm.num();
          float rm_den = rm.den();
          scitbx::vec3<int> tv_num = tv.num();
          float tv_den = tv.den();
          for (int ir=0;ir<3;ir++) {
            for (int ic=0;ic<3;ic++) {
              p->mtzsymm.sym[im][ir][ic] = rm_num(ir,ic)/rm_den;
            }
            p->mtzsymm.sym[im][ir][3] = tv_num[ir]/tv_den;
          }
          for (int ic=0;ic<3;ic++) {
            p->mtzsymm.sym[im][3][ic] = 0.;
          }
          p->mtzsymm.sym[im][3][3] = 1.;
        }
        return *this;
      }

      void
      reserve(int capacity)
      {
        CMtz::MTZ* p = ptr();
        if (!p->refs_in_memory) return;
        for(int i=0;i<p->nxtal;i++) {
          for(int j=0;j<p->xtal[i]->nset;j++) {
            for(int k=0;k<p->xtal[i]->set[j]->ncol;k++) {
              capacity = std::max(
                capacity,
                column_array_capacity(p->xtal[i]->set[j]->col[k]));
            }
          }
        }
        for(int i=0;i<p->nxtal;i++) {
          for(int j=0;j<p->xtal[i]->nset;j++) {
            for(int k=0;k<p->xtal[i]->set[j]->ncol;k++) {
              ccp4array_reserve(p->xtal[i]->set[j]->col[k]->ref, capacity);
            }
          }
        }
      }

      void
      adjust_column_array_sizes(int new_nref)
      {
        CMtz::MTZ* p = ptr();
        if (!p->refs_in_memory) return;
        if (new_nref > p->nref) {
          reserve(new_nref);
          for(int i=0;i<p->nxtal;i++) {
            for(int j=0;j<p->xtal[i]->nset;j++) {
              for(int k=0;k<p->xtal[i]->set[j]->ncol;k++) {
                CMtz::MTZCOL* col_k = p->xtal[i]->set[j]->col[k];
                int old_size = column_array_size(col_k);
                if (new_nref > old_size) {
                  ccp4array_resize(col_k->ref, new_nref);
                  for(int iref=old_size;iref<new_nref;iref++) {
                    *(reinterpret_cast<union float_uint_uchar*>(
                      &col_k->ref[iref])) = CCP4::ccp4_nan();
                  }
                }
              }
            }
          }
        }
      }

      int
      n_batches() const { return CMtz::MtzNbat(ptr()); }

      int
      n_reflections() const { return ptr()->nref; }

      inline
      af::tiny<double, 2>
      max_min_resolution() const;

      int
      n_crystals() const { return CMtz::MtzNxtal(ptr()); }

      int
      n_active_crystals() const { return CMtz::MtzNumActiveXtal(ptr()); }

      inline
      af::shared<crystal>
      crystals() const;

      inline
      crystal
      add_crystal(
        const char* name,
        const char* project_name,
        af::double6 const& unit_cell_parameters);

      inline
      crystal
      add_crystal(
        const char* name,
        const char* project_name,
        cctbx::uctbx::unit_cell const& unit_cell);

      inline
      bool
      has_column(const char* label) const;

      inline
      column
      get_column(const char* label) const;

      inline
      hkl_columns
      get_hkl_columns() const;

      inline
      integer_group
      extract_integers(
        const char* column_label);

      inline
      integer_group
      extract_integers_anomalous(
        const char* column_label_plus,
        const char* column_label_minus);

      inline
      real_group
      extract_reals(
        const char* column_label);

      inline
      real_group
      extract_reals_anomalous(
        const char* column_label_plus,
        const char* column_label_minus);

      inline
      hl_group
      extract_hls(
        const char* column_label_a,
        const char* column_label_b,
        const char* column_label_c,
        const char* column_label_d);

      inline
      hl_group
      extract_hls_anomalous(
        const char* column_label_a_plus,
        const char* column_label_b_plus,
        const char* column_label_c_plus,
        const char* column_label_d_plus,
        const char* column_label_a_minus,
        const char* column_label_b_minus,
        const char* column_label_c_minus,
        const char* column_label_d_minus);

      inline
      observations_group
      extract_observations(
        const char* column_label_data,
        const char* column_label_sigmas);

      inline
      observations_group
      extract_observations_anomalous(
        const char* column_label_data_plus,
        const char* column_label_sigmas_plus,
        const char* column_label_data_minus,
        const char* column_label_sigmas_minus);

      inline
      observations_group
      extract_delta_anomalous(
        const char* column_label_f_data,
        const char* column_label_f_sigmas,
        const char* column_label_d_data,
        const char* column_label_d_sigmas);

      inline
      complex_group
      extract_complex(
        const char* column_label_ampl,
        const char* column_label_phi);

      inline
      complex_group
      extract_complex_anomalous(
        const char* column_label_ampl_plus,
        const char* column_label_phi_plus,
        const char* column_label_ampl_minus,
        const char* column_label_phi_minus);

      void
      write(const char* file_name)
      {
        if (!CMtz::MtzPut(ptr(), file_name)) {
          throw cctbx::error("MTZ write failed.");
        }
      }

    protected:
      boost::shared_ptr<CMtz::MTZ> ptr_;

      static void
      ptr_deleter(CMtz::MTZ* ptr)
      {
        if (ptr != 0) CCTBX_ASSERT(CMtz::MtzFree(ptr));
      }
  };

}} // namespace iotbx::mtz

#endif // IOTBX_MTZ_OBJECT_H
