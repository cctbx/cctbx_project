#ifndef IOTBX_MTZ_OBJECT_H
#define IOTBX_MTZ_OBJECT_H

#include <cmtzlib.h>
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
    return CCP4::ccp4_utils_isnan((union float_uint_uchar *) &datum);
  }

  template <typename DataType>
  struct data_group
  {
    data_group() {}

    data_group(std::size_t size)
    {
      indices.reserve(size);
      data.reserve(size);
    }

    af::shared<cctbx::miller::index<> > indices;
    af::shared<DataType> data;
  };

  typedef data_group<int> integer_group;
  typedef data_group<double> real_group;
  typedef data_group<cctbx::hendrickson_lattman<> > hl_group;

  struct observations_group : real_group
  {
    observations_group() {}

    observations_group(std::size_t size)
    :
      real_group(size)
    {
      sigmas.reserve(size);
    }

    af::shared<double> sigmas;
  };

  struct complex_group
  {
    complex_group() {}

    complex_group(std::size_t size)
    {
      indices.reserve(size);
      data.reserve(size);
    }

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
        ptr_(CMtz::MtzMalloc(0, 0), ptr_deleter())
      {
        if (ptr_.get() == 0) throw cctbx::error("MtzMalloc failed.");
        ptr_->refs_in_memory = true;
      }

      object(af::const_ref<int> const& n_datasets_for_each_crystal)
      :
        ptr_(CMtz::MtzMalloc(
            n_datasets_for_each_crystal.size(),
            const_cast<int*>(&*n_datasets_for_each_crystal.begin())),
          ptr_deleter())
      {
        if (ptr_.get() == 0) throw cctbx::error("MtzMalloc failed.");
        ptr_->refs_in_memory = true;
      }

      object(const char* file_name)
      :
        ptr_(CMtz::MtzGet(file_name, true), ptr_deleter())
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

      int
      n_batches() const { return CMtz::MtzNbat(ptr()); }

      int
      n_reflections() const { return CMtz::MtzNref(ptr()); }

      af::tiny<double, 2>
      max_min_resolution() const
      {
        if (ptr()->nxtal == 0) return af::tiny<double, 2>(-1., -1.);
        float max_resolution;
        float min_resolution;
        CCTBX_ASSERT(
          CMtz::MtzResLimits(ptr(), &min_resolution, &max_resolution));
        return af::tiny<double, 2>(max_resolution, min_resolution);
      }

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
      column
      lookup_column(const char* label) const;

      inline
      hkl_columns
      lookup_hkl_columns() const;

      inline
      integer_group
      extract_integers(
        const char* column_label);

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

    protected:
      boost::shared_ptr<CMtz::MTZ> ptr_;

      struct ptr_deleter
      {
        void
        operator()(CMtz::MTZ* ptr) const
        {
          CCTBX_ASSERT(CMtz::MtzFree(ptr));
        }
      };
  };

}} // namespace iotbx::mtz

#endif // IOTBX_MTZ_OBJECT_H
