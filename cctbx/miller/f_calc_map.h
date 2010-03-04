#ifndef CCTBX_MILLER_MAP_H
#define CCTBX_MILLER_MAP_H

#include <cctbx/miller/index_span.h>
#include <cctbx/maptbx/utils.h>
#include <scitbx/array_family/versa_plain.h>

namespace cctbx { namespace miller {

  /// The accessor doing the magic behind the class f_calc_map
  class hermitian_accessor
  {
    public:
      typedef miller::index<> index_type;
      typedef miller::index<>::value_type index_value_type;

      struct index_1d_flag_conj
      {
        long index_1d;
        bool is_conj;
      };

      hermitian_accessor() {}

      hermitian_accessor(bool anomalous_flag,
                         miller::index<> const& n,
                         bool two_n_minus_one_flag)
      :
        anomalous_flag_(anomalous_flag),
        range_(two_n_minus_one_flag ? miller::index<>(2*n-1) : n)
      {
        if (!anomalous_flag_) range_[2] = n[2];
      }

      std::size_t
      size_1d() const { return range_.product(); }

      index_value_type
      operator[](std::size_t i) const { return range_[i]; }

      bool
      anomalous_flag() const { return anomalous_flag_; }

      index_1d_flag_conj
      operator()(miller::index<> h) const
      {
        miller::index<> i3d;
        index_1d_flag_conj result;
        result.index_1d = -1;
        result.is_conj = false;
        if (!anomalous_flag_) {
          if (h[2] < 0) {
            h = -h;
            result.is_conj = true;
          }
          for(std::size_t i=0;i<2;i++) {
            i3d[i] = maptbx::h_as_ih_exact(h[i], range_[i], false);
          }
          i3d[2] = maptbx::h_as_ih_exact(h[2], range_[2], true);
        }
        else {
          for(std::size_t i=0;i<3;i++) {
            i3d[i] = maptbx::h_as_ih_exact(h[i], range_[i], false);
          }
        }
        if (i3d.min() < 0) return result;
        // Manually optimized for best performance.
        result.index_1d = (i3d[0] * range_[1] + i3d[1]) * range_[2] + i3d[2];
        return result;
      }

    protected:
      bool anomalous_flag_;
      miller::index<> range_;
  };

  /// A mapping of miller indices to f_calc values with fast readonly access
  template <typename FloatType>
  class f_calc_map
  {
    public:
      typedef std::complex<FloatType> complex_type;
      typedef af::versa_plain<complex_type, hermitian_accessor>
                data_array_type;

      /// Construct an empty mapping
      f_calc_map() {}

      /// Construct a mapping of the given indices to the given f_calc
      f_calc_map(af::const_ref<miller::index<> > const &indices,
                 af::const_ref<complex_type> const &f_calc,
                 bool anomalous_flag)
      :
        data_(hermitian_accessor(anomalous_flag,
                                 miller::index_span(indices).abs_range(),
                                 true))
      {
        import(indices, f_calc);
      }

      /// The argument of same name passed to the constructor
      bool
      anomalous_flag() const { return data_.accessor().anomalous_flag(); }

      void import(af::const_ref<miller::index<> > const& miller_indices,
                  af::const_ref<complex_type> const& f_calc)
      {
        CCTBX_ASSERT(miller_indices.size() == f_calc.size());
        for(std::size_t i=0;i<f_calc.size();i++) {
          set(miller_indices[i], f_calc[i]);
        }
      }

      /// The value of f_calc for the given miller index
      /** or 0 if there is no mapping for h. Exception:
          if anomalous_flag is true, f[-h] is the complex conjugate of f[h].
       */
      complex_type
      operator[](miller::index<> const& h) const
      {
        hermitian_accessor::index_1d_flag_conj ic = data_.accessor()(h);
        if (ic.index_1d < 0) return complex_type(0);
        if (ic.is_conj) return std::conj(data_[ic.index_1d]);
        return data_[ic.index_1d];
      }

    protected:
      void
      set(miller::index<> const& h, complex_type const& val)
      {
        hermitian_accessor::index_1d_flag_conj ic = data_.accessor()(h);
        CCTBX_ASSERT(ic.index_1d >= 0);
        if (ic.is_conj) data_[ic.index_1d] = std::conj(val);
        else            data_[ic.index_1d] =           val;
        if (anomalous_flag() || h[2] != 0) return;
        ic = data_.accessor()(-h);
        CCTBX_ASSERT(ic.index_1d >= 0);
        if (ic.is_conj) data_[ic.index_1d] =           val;
        else            data_[ic.index_1d] = std::conj(val);
      }

      data_array_type data_;
  };


}}

#endif
