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
      mtz_crystal() const { return mtz_dataset_.mtz_crystal(); }

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

  inline
  af::shared<cctbx::miller::index<> >
  object::valid_indices(const char* column_label) const
  {
    return lookup_column(column_label).valid_indices();
  }

  inline
  af::shared<double>
  object::valid_values(const char* column_label) const
  {
    return lookup_column(column_label).valid_values();
  }

  inline
  af::shared<int>
  object::valid_integers(const char* column_label) const
  {
    return lookup_column(column_label).valid_integers();
  }

  inline
  af::shared<cctbx::miller::index<> >
  object::valid_indices_anomalous(
    const char* column_label_plus,
    const char* column_label_minus) const
  {
    column h(lookup_column("H"));
    column k(lookup_column("K"));
    column l(lookup_column("L"));
    column values_plus(lookup_column(column_label_plus));
    column values_minus(lookup_column(column_label_minus));
    af::shared<cctbx::miller::index<> >
      result((af::reserve(2*n_reflections())));
    for(int i=0;i<n_reflections();i++) {
      if (!values_plus.is_ccp4_nan(i)) {
        result.push_back(cctbx::miller::index<>(
          h.int_datum(i), k.int_datum(i), l.int_datum(i)));
      }
      if (!values_minus.is_ccp4_nan(i)) {
        result.push_back(cctbx::miller::index<>(
          -h.int_datum(i), -k.int_datum(i), -l.int_datum(i)));
      }
    }
    return result;
  }

  inline
  af::shared<double>
  object::valid_values_anomalous(
    const char* column_label_plus,
    const char* column_label_minus) const
  {
    column values_plus(lookup_column(column_label_plus));
    column values_minus(lookup_column(column_label_minus));
    af::shared<double> result((af::reserve(2*n_reflections())));
    for(int i=0;i<n_reflections();i++) {
      if (!values_plus.is_ccp4_nan(i)) {
        result.push_back(values_plus.float_datum(i));
      }
      if (!values_minus.is_ccp4_nan(i)) {
        result.push_back(values_minus.float_datum(i));
      }
    }
    return result;
  }

  namespace detail {
    inline
    std::complex<double>
    polar_deg(double ampl, double phi)
    {
      return std::polar(ampl, phi*scitbx::constants::pi_180);
    }
  }

  inline
  af::shared<std::complex<double> >
  object::valid_complex(
    const char* column_label_ampl,
    const char* column_label_phi)
  {
    column values_ampl(lookup_column(column_label_ampl));
    column values_phi(lookup_column(column_label_phi));
    af::shared<std::complex<double> > result((af::reserve(n_reflections())));
    for(int i=0;i<n_reflections();i++) {
      if (values_ampl.is_ccp4_nan(i) != values_phi.is_ccp4_nan(i)) {
        throw cctbx::error(std::string(
          "Unexpected NAN while extracting complex array from columns: ")
          + column_label_ampl + ", " + column_label_phi);
      }
      if (!values_ampl.is_ccp4_nan(i)) {
        result.push_back(detail::polar_deg(
          values_ampl.float_datum(i), values_phi.float_datum(i)));
      }
    }
    return result;
  }

  inline
  af::shared<std::complex<double> >
  object::valid_complex_anomalous(
    const char* column_label_ampl_plus,
    const char* column_label_phi_plus,
    const char* column_label_ampl_minus,
    const char* column_label_phi_minus)
  {
    column values_ampl_plus(lookup_column(column_label_ampl_plus));
    column values_phi_plus(lookup_column(column_label_phi_plus));
    column values_ampl_minus(lookup_column(column_label_ampl_minus));
    column values_phi_minus(lookup_column(column_label_phi_minus));
    af::shared<std::complex<double> > result((af::reserve(2*n_reflections())));
    for(int i=0;i<n_reflections();i++) {
      if (values_ampl_plus.is_ccp4_nan(i) != values_phi_plus.is_ccp4_nan(i)) {
        throw cctbx::error(std::string(
          "Unexpected NAN while extracting complex array from columns: ")
          + column_label_ampl_plus + ", " + column_label_phi_plus);
      }
      if (!values_ampl_plus.is_ccp4_nan(i)) {
        result.push_back(detail::polar_deg(
          values_ampl_plus.float_datum(i), values_phi_plus.float_datum(i)));
      }
      if (values_ampl_minus.is_ccp4_nan(i) != values_phi_minus.is_ccp4_nan(i)){
        throw cctbx::error(std::string(
          "Unexpected NAN while extracting complex array from columns: ")
          + column_label_ampl_minus + ", " + column_label_phi_minus);
      }
      if (!values_ampl_minus.is_ccp4_nan(i)) {
        result.push_back(detail::polar_deg(
          values_ampl_minus.float_datum(i), values_phi_minus.float_datum(i)));
      }
    }
    return result;
  }

  /*! http://www.ccp4.ac.uk/dist/html/mtzMADmod.html
        F(+) = F + 0.5*D
        F(-) = F - 0.5*D
        SIGF(+) = sqrt( SIGF**2 + 0.25*SIGD**2 )
        SIGF(-) = SIGF(+)
   */
  inline
  array_group
  object::valid_delta_anomalous(
    const char* column_label_f_data,
    const char* column_label_f_sigmas,
    const char* column_label_d_data,
    const char* column_label_d_sigmas)
  {
    column h(lookup_column("H"));
    column k(lookup_column("K"));
    column l(lookup_column("L"));
    column values_fd(lookup_column(column_label_f_data));
    column values_fs(lookup_column(column_label_f_sigmas));
    column values_dd(lookup_column(column_label_d_data));
    column values_ds(lookup_column(column_label_d_sigmas));
    array_group result(2*n_reflections());
    for(int i=0;i<n_reflections();i++) {
      if (   values_fs.is_ccp4_nan(i) != values_fd.is_ccp4_nan(i)
          || values_ds.is_ccp4_nan(i) != values_dd.is_ccp4_nan(i)
          || (values_fd.is_ccp4_nan(i) && !values_dd.is_ccp4_nan(i))) {
        throw cctbx::error(std::string(
          "Unexpected NAN while extracting complex array from columns: ")
          + column_label_f_data   + ", "
          + column_label_f_sigmas + ", "
          + column_label_d_data   + ", "
          + column_label_d_sigmas);
      }
      if (!values_fd.is_ccp4_nan(i)) {
        result.indices.push_back(cctbx::miller::index<>(
           h.int_datum(i),  k.int_datum(i),  l.int_datum(i)));
        double fd = values_fd.float_datum(i);
        double fs = values_fs.float_datum(i);
        if (values_dd.is_ccp4_nan(i)) {
          result.data.push_back(fd);
          result.sigmas.push_back(fs);
        }
        else {
          double ddh = values_dd.float_datum(i) * .5;
          double dsh = values_ds.float_datum(i) * .5;
          double s = std::sqrt(fs*fs + dsh*dsh);
          result.indices.push_back(cctbx::miller::index<>(
            -h.int_datum(i), -k.int_datum(i), -l.int_datum(i)));
          result.data.push_back(fd + ddh);
          result.data.push_back(fd - ddh);
          result.sigmas.push_back(std::sqrt(s));
          result.sigmas.push_back(std::sqrt(s));
        }
      }
    }
    return result;
  }

  inline
  af::shared<cctbx::hendrickson_lattman<> >
  object::valid_hl(
    const char* column_label_a,
    const char* column_label_b,
    const char* column_label_c,
    const char* column_label_d)
  {
    column values_a(lookup_column(column_label_a));
    column values_b(lookup_column(column_label_b));
    column values_c(lookup_column(column_label_c));
    column values_d(lookup_column(column_label_d));
    af::shared<cctbx::hendrickson_lattman<> >
      result((af::reserve(n_reflections())));
    for(int i=0;i<n_reflections();i++) {
      if (   values_b.is_ccp4_nan(i) != values_a.is_ccp4_nan(i)
          || values_c.is_ccp4_nan(i) != values_a.is_ccp4_nan(i)
          || values_d.is_ccp4_nan(i) != values_a.is_ccp4_nan(i)) {
        throw cctbx::error(std::string(
          "Unexpected NAN while extracting complex array from columns: ")
          + column_label_a + ", "
          + column_label_b + ", "
          + column_label_c + ", "
          + column_label_d);
      }
      if (!values_a.is_ccp4_nan(i)) {
        result.push_back(cctbx::hendrickson_lattman<>(
          values_a.float_datum(i),
          values_b.float_datum(i),
          values_c.float_datum(i),
          values_d.float_datum(i)));
      }
    }
    return result;
  }

}} // namespace iotbx::mtz

#endif // IOTBX_MTZ_COLUMN_H
