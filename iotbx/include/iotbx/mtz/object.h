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

  struct array_group
  {
    array_group() {}

    array_group(std::size_t size)
    {
      indices.reserve(size);
      data.reserve(size);
      sigmas.reserve(size);
    }

    af::shared<cctbx::miller::index<> > indices;
    af::shared<double> data;
    af::shared<double> sigmas;
  };

  class column;
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
        std::vector<char> result(strlen(ptr()->title)+1);
        int len = CMtz::ccp4_lrtitl(ptr(), &*result.begin());
        CCTBX_ASSERT(len < result.size());
        return std::string(&*result.begin(), len);
      }

      af::shared<std::string>
      history() const
      {
        CMtz::MTZ* p = ptr();
        af::shared<std::string> result((af::reserve(p->histlines)));
        for(int i=0;i<p->histlines;i++) {
          result.push_back(
            std::string(p->hist+MTZRECORDLENGTH*i, MTZRECORDLENGTH));
        }
        return result;
      }

      std::string
      space_group_name() const { return ptr()->mtzsymm.spcgrpname; }

      std::string
      point_group_name() const { return ptr()->mtzsymm.pgname; }

      cctbx::sgtbx::space_group
      space_group() const
      {
        cctbx::sgtbx::space_group result;
        scitbx::mat3<double> r;
        scitbx::vec3<double> t;
        for(int i=0;i<ptr()->mtzsymm.nsym;i++) {
          for (int p=0;p<3;p++) {
            for (int q=0;q<3;q++) {
              r(p,q) = ptr()->mtzsymm.sym[i][p][q];
            }
            t[p] = ptr()->mtzsymm.sym[i][p][3];
          }
          result.expand_smx(cctbx::sgtbx::rt_mx(r, t));
        }
        return result;
      }

      int
      n_batches() const { return CMtz::MtzNbat(ptr()); }

      int
      n_reflections() const { return CMtz::MtzNref(ptr()); }

      int
      space_group_number() const { return CMtz::MtzSpacegroupNumber(ptr()); }

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

      af::shared<crystal>
      crystals() const;

      inline
      column
      lookup_column(const char* label) const;

      inline
      af::shared<cctbx::miller::index<> >
      valid_indices(const char* column_label) const;

      inline
      af::shared<double>
      valid_values(const char* column_label) const;

      inline
      af::shared<int>
      valid_integers(const char* column_label) const;

      inline
      af::shared<cctbx::miller::index<> >
      valid_indices_anomalous(
        const char* column_label_plus,
        const char* column_label_minus) const;

      inline
      af::shared<double>
      valid_values_anomalous(
        const char* column_label_plus,
        const char* column_label_minus) const;

      af::shared<std::complex<double> >
      valid_complex(
        const char* column_label_ampl,
        const char* column_label_phi);

      af::shared<std::complex<double> >
      valid_complex_anomalous(
        const char* column_label_ampl_plus,
        const char* column_label_phi_plus,
        const char* column_label_ampl_minus,
        const char* column_label_phi_minus);

      array_group
      valid_delta_anomalous(
        const char* column_label_f_data,
        const char* column_label_f_sigmas,
        const char* column_label_d_data,
        const char* column_label_d_sigmas);

      af::shared<cctbx::hendrickson_lattman<> >
      valid_hl(
        const char* column_label_a,
        const char* column_label_b,
        const char* column_label_c,
        const char* column_label_d);

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
