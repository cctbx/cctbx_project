#ifndef CCTBX_MAPTBX_MAP_ACCUMULATOR_H
#define CCTBX_MAPTBX_MAP_ACCUMULATOR_H

#include <scitbx/array_family/accessors/c_grid.h>
#include <cctbx/xray/sampling_base.h>

#if defined(_MSC_VER) && _MSC_VER < 1600
typedef unsigned char     uint8_t;
#endif


namespace cctbx { namespace maptbx {

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

// Compute median
double median(std::vector<double> data) {
    size_t n = data.size();
    std::nth_element(data.begin(), data.begin() + n / 2, data.end());
    double med = data[n / 2];
    if (n % 2 == 0) {
        std::nth_element(data.begin(), data.begin() + n / 2 - 1, data.end());
        med = 0.5 * (med + data[n / 2 - 1]);
    }
    return med;
}

// Compute IQR (Interquartile Range)
double iqr(std::vector<double> data) {
    size_t n = data.size();
    std::sort(data.begin(), data.end());
    double q1 = data[n / 4];
    double q3 = data[(3 * n) / 4];
    return q3 - q1;
}

// Auto bin count selection using Freedman-Diaconis
int auto_num_bins(const std::vector<double>& data) {
    double min_val = *std::min_element(data.begin(), data.end());
    double max_val = *std::max_element(data.begin(), data.end());
    double range = max_val - min_val;
    if (range == 0) return 1; // all values same

    double IQR = iqr(data);
    if (IQR == 0) {
        // Fallback to Sturges' rule
        return std::max(8, (int)std::ceil(std::log2(data.size()) + 1));
    }

    double bin_width = 2.0 * IQR / std::cbrt((double)data.size());
    int bins = (int)std::ceil(range / bin_width);

    // Clamp bins to reasonable range
    if (bins < 8) bins = 8;
    if (bins > 128) bins = 128;

    return bins;
}

// Fast robust mode or median with auto bin count and auto threshold
double fast_mode_or_median_auto_bins(const std::vector<double>& data) {
    if (data.empty()) return std::numeric_limits<double>::quiet_NaN();

    int num_bins = auto_num_bins(data);

    double min_val = *std::min_element(data.begin(), data.end());
    double max_val = *std::max_element(data.begin(), data.end());
    if (min_val == max_val) return min_val;

    double bin_width = (max_val - min_val) / num_bins;
    std::vector<int> counts(num_bins, 0);

    // Fill histogram
    for (size_t i = 0; i < data.size(); ++i) {
        int bin = (int)((data[i] - min_val) / bin_width);
        if (bin == num_bins) bin--; // edge case
        counts[bin]++;
    }

    // Find max bin
    int max_bin_idx = std::max_element(counts.begin(), counts.end()) - counts.begin();
    int max_count = counts[max_bin_idx];
    double avg_count = (double)data.size() / num_bins;

    // Auto threshold based on statistical noise
    double p = 1.0 / num_bins;
    double noise_level = std::sqrt((1.0 - p) / (data.size() * p));
    double prominence_threshold = noise_level * 2.0; // 2sigma rule

    double prominence = (max_count - avg_count) / (double)max_count;

    if (prominence < prominence_threshold) {
        return median(data); // flat -- median
    } else {
        return min_val + (max_bin_idx + 0.5) * bin_width; // mode midpoint
    }
}

template <typename FloatType, typename GridType>
class map_accumulator2 {
public:
  af::versa<af::shared<FloatType>, GridType> map_new;
  af::int3 n_real;
  cctbx::xray::detail::exponent_table<FloatType> exp_table_;


  map_accumulator2(
    af::int3 const& n_real_)
  :
  n_real(n_real_),
  exp_table_(-100)
  {
    map_new.resize(GridType(n_real));
    for(std::size_t i=0;i<map_new.size(); i++) {
      map_new[i] = af::shared<FloatType>();
    }
  }

  void add(af::const_ref<FloatType, GridType> const& map_data)
  {
    GridType a = map_data.accessor();
    for(int i = 0; i < 3; i++) CCTBX_ASSERT(a[i]==n_real[i]);
    for(std::size_t i=0;i<map_new.size(); i++)
      map_new[i].push_back( map_data[i] );
  }


  af::versa<FloatType, GridType>
  as_median_map()
  {
    af::versa<FloatType, GridType> result;
    result.resize(GridType(n_real), 0.0);
    for(int i = 0; i < n_real[0]; i++) {
      std::cout<<i<<std::endl;
      for(int j = 0; j < n_real[1]; j++) {
        for(int k = 0; k < n_real[2]; k++) {

          std::vector<FloatType> data;
          af::shared<FloatType> dataf = map_new(i,j,k);
          for(int N = 0; N < map_new(i,j,k).size(); N++) {
            FloatType d = dataf[N];
            if(std::abs(d)>0.001) data.push_back(d);
          }

          double representative = fast_mode_or_median_auto_bins(data);
          result(i,j,k) = representative;

//            std::size_t N = data.size();
//            double max = data[0];
//            // ---- Compute mean ----
//            FloatType sum = 0;
//            for (size_t i = 0; i < N; ++i) {
//                sum += data[i];
//                if(max < data[i]) max = data[i];
//            }
//            FloatType mean = sum / static_cast<FloatType>(N);
//
//            result(i,j,k) = mean;


//
//          if(std::abs(mean) < 0.2) result(i,j,k) = 0.0;
//          else {
//
//            // ---- Compute variance ----
//            FloatType var_sum = 0;
//            for (size_t i = 0; i < N; ++i) {
//                FloatType diff = data[i] - mean;
//                var_sum += diff * diff;
//            }
//
//            // Sample standard deviation (N-1 in denominator)
//            FloatType stdev = 0.;
//            if (N > 1) {
//                stdev = std::sqrt(var_sum / static_cast<FloatType>(N - 1));
//            }
//
//            double CV = stdev / std::abs(mean);
//            if(CV > 2.0) {
//              result(i,j,k) = 0;
//              cntr_kde += 1;
//            }
//            else {
//              ModeResult modeEst = estimateModeHistogram(data);
//              cntr_mode += 1;
//              if (modeEst.hasStrongMode) result(i,j,k) = modeEst.modeValue;
//              else                       result(i,j,k) = mean;
//            }
//
//          }

    }}}

    return result;
  }

};


template <typename FloatType, typename GridType>
class map_accumulator {
public:
  af::versa<af::shared<uint8_t>, GridType> map_new;
  af::shared<FloatType> v_values_;
  af::int3 n_real;
  cctbx::xray::detail::exponent_table<FloatType> exp_table_;
  FloatType smearing_b;
  FloatType max_peak_scale;
  int smearing_span;
  bool use_exp_table;
  bool use_max_map;

  map_accumulator(
    af::int3 const& n_real_,
    FloatType const& smearing_b_,
    FloatType const& max_peak_scale_,
    int const& smearing_span_,
    bool use_exp_table_,
    bool use_max_map_)
  :
  n_real(n_real_), exp_table_(-100), smearing_b(smearing_b_),
  max_peak_scale(max_peak_scale_), smearing_span(smearing_span_),
  use_exp_table(use_exp_table_), use_max_map(use_max_map_)
  {
    map_new.resize(GridType(n_real));
    for(std::size_t i=0;i<map_new.size(); i++) map_new[i]=af::shared<uint8_t>();
  }

  void add(af::const_ref<FloatType, GridType> const& map_data)
  {
    GridType a = map_data.accessor();
    for(int i = 0; i < 3; i++) CCTBX_ASSERT(a[i]==n_real[i]);
    FloatType map_min = 0;//af::min(map_data);
    for(std::size_t i=0;i<map_new.size(); i++)
      map_new[i].push_back((uint8_t)to_int(map_data[i], map_min));
  }

  uint8_t to_int(FloatType x, FloatType p0)
  {
    CCTBX_ASSERT(x>=0 && x<=1);
    if(x<=p0) return 0;
    return (uint8_t)std::min(int(256*(x-p0)/(1.-p0))+1, 255);
  }

  af::shared<int> at_index(af::int3 const& n)
  {
    af::shared<int> result;
    for(int i = 0; i < map_new(n).size(); i++) result.push_back(map_new(n)[i]);
    return result;
  }

  inline FloatType smear(FloatType x_mins_a, FloatType two_b_sq)
  {
    if(!this->use_exp_table) { return std::exp(-x_mins_a*x_mins_a/two_b_sq); }
    else { return exp_table_(-x_mins_a*x_mins_a/two_b_sq); }
  }

  af::shared<FloatType> int_to_float_at_index(af::int3 const& n)
  {
    af::shared<uint8_t> as = map_new(n);
    af::shared<FloatType> result;
    FloatType b = this->smearing_b;
    FloatType two_b_sq = 2 * b * b;
    int ss = this->smearing_span;
    result.resize(256, 0);
    for(int i = 0; i < as.size(); i++) {
      int a = (int)as[i];
      for(int j = -ss; j <=ss; j++) {
        int x = a + j;
        if(x>=0 && x<=255) {
          FloatType x_mins_a = (FloatType)x - (FloatType)a;
          result[x] += smear(x_mins_a, two_b_sq);
        }
    }}
    return result;
  }

  FloatType quadratic_approximation(FloatType x1,FloatType x2,FloatType x3,
    FloatType f1, FloatType f2, FloatType f3) {
    if(x1<x2 && x2<x3) {
      FloatType s21 = (f2-f1)/(x2-x1);
      FloatType s32 = (f3-f2)/(x3-x2);
      return (x1+x2)/2-s21*(x3-x1)/2./(s32-s21);
    }
    else return x2;
  }

  FloatType find_peaks(af::const_ref<FloatType> const& f)
  {
    CCTBX_ASSERT(f.size()==256);
    FloatType result = 0.;
    af::shared<FloatType> peaks;
    af::shared<int> peak_args;
    FloatType lv=0, rv=0, eps=1.e-3;
    // find peaks; does not handle flat peaks like 0 1 222 1 0
    FloatType p_min = 0, p_min_=1.e+9, p_max_=-1.e+9;
    FloatType p_max = 0;
    for(int i = 0; i < 256; i++) {
      FloatType v = f[i];
      if(std::abs(v-1.)>eps && v>1.) {
        if(i==0) {
          rv = f[i+1];
          if(v>rv) {
            peaks.push_back(v);
            peak_args.push_back(i);
            if(v<p_min_) p_min_ = v;
            if(v>p_max_) p_max_ = v;
          }
        }
        else if(i==255) {
          lv = f[i-1];
          if(v>lv) {
            peaks.push_back(v);
            peak_args.push_back(i);
            if(v<p_min_) p_min_ = v;
            if(v>p_max_) p_max_ = v;
          }
        }
        else {
          lv = f[i-1];
          rv = f[i+1];
          if(v>lv && v>rv) {
            peaks.push_back(v);
            peak_args.push_back(i);
            if(v<p_min_) p_min_ = v;
            if(v>p_max_) p_max_ = v;
          }
        }
      }
    }
    p_min = p_min_;
    p_max = p_max_;
    FloatType p_max_over_2 = p_max/this->max_peak_scale;
    // analyze peaks
    if(peaks.size()==0) return 0;
    int i_result;
    // only one peak
    if(peaks.size()==1) {
      CCTBX_ASSERT(peak_args.size()==1);
      i_result = peak_args[0];
    }
    // several similar peaks
    else {
      int i_result_ = 1.e+9; // << BUG ?
      for(int i = 0; i < peaks.size(); i++) {
        FloatType peak = peaks[i];
        if(use_max_map) {
          if(peak>=p_max) { // Max map
            FloatType pa = peak_args[i];
            if(pa<i_result_) i_result_ = pa;
          }
        }
        else {
          if(peak<=p_max && peak>=p_max_over_2) { // Min map
            FloatType pa = peak_args[i];
            if(pa<i_result_) i_result_ = pa;
          }
        }
      }
      i_result = i_result_;
    }
    FloatType i_result_f = (FloatType)i_result;
    if(i_result>0 && i_result<255) {
      i_result_f = quadratic_approximation(
        i_result-1,
        i_result,
        i_result+1,
        f[i_result-1],
        f[i_result],
        f[i_result+1]);
    }
    return i_result_f;
  }

  af::versa<FloatType, GridType>
  as_median_map()
  {
    af::versa<FloatType, GridType> result;
    result.resize(GridType(n_real), 0.0);
    for(int i = 0; i < n_real[0]; i++) {
      for(int j = 0; j < n_real[1]; j++) {
        for(int k = 0; k < n_real[2]; k++) {
          result(i,j,k) = find_peaks(
            int_to_float_at_index(af::int3(i,j,k)).ref());
    }}}
    return result;
  }

};

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_MAP_ACCUMULATOR_H
