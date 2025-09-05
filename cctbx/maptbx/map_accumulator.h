#ifndef CCTBX_MAPTBX_MAP_ACCUMULATOR_H
#define CCTBX_MAPTBX_MAP_ACCUMULATOR_H

#include <scitbx/array_family/accessors/c_grid.h>
#include <cctbx/xray/sampling_base.h>

#if defined(_MSC_VER) && _MSC_VER < 1600
typedef unsigned char     uint8_t;
#endif


namespace cctbx { namespace maptbx {

/*

double stdev(const std::vector<double>& data) {
    double m = 0.0;
    for (std::vector<double>::const_iterator it = data.begin(); it != data.end(); ++it) {
      m += *it;
    }
    m = m / data.size();
    double sumSq = 0.0;
    for (double v : data) sumSq += (v - m) * (v - m);
    return std::sqrt(sumSq / (data.size() - 1));
}


// Sorts the data so we can compute quartiles.
// Calculates IQR and bin width using Freedman–Diaconis.
// If Freedman–Diaconis produces too few bins, falls back to Sturges’ formula.
// Bins the data and finds the most populated bin.
// Returns the mean of all values in that bin as the mode estimate.
// Calculates mode prominence — how much higher the most populated bin is compared to the second-most populated bin, relative to the average bin count.
// Flags “no strong mode” if the prominence is below a threshold (default 15%).

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>

double computeQuantile(const std::vector<double>& data, double q) {
    // Assumes data is sorted
    double pos = (data.size() - 1) * q;
    size_t idx = (size_t)pos;
    double frac = pos - idx;
    if (idx + 1 < data.size())
        return data[idx] * (1 - frac) + data[idx + 1] * frac;
    else
        return data[idx];
}

struct ModeResult {
    double modeValue;
    double prominence;
    bool hasStrongMode;
};

ModeResult estimateModeHistogram(const std::vector<double>& data, double prominenceThreshold = 0.1) {
    ModeResult result{NAN, 0.0, false};

    if (data.empty()) return result;
    if (data.size() == 1) {
        result.modeValue = data[0];
        result.hasStrongMode = true;
        return result;
    }

    // Sort data for quantile computation
    std::vector<double> sorted = data;
    std::sort(sorted.begin(), sorted.end());

    size_t n = sorted.size();
    double minVal = sorted.front();
    double maxVal = sorted.back();

    // Compute IQR for Freedman–Diaconis
    double q1 = computeQuantile(sorted, 0.25);
    double q3 = computeQuantile(sorted, 0.75);
    double iqr = q3 - q1;

    // Freedman–Diaconis bin width
    double h = 2.0 * iqr / std::cbrt((double)n);
    int binCount;
    if (h > 0) {
        binCount = (int)std::ceil((maxVal - minVal) / h);
    } else {
        binCount = 1; // all values same
    }

    // If FD rule gives too few bins, use Sturges' formula
    if (binCount < 1) {
        binCount = (int)std::ceil(std::log2((double)n) + 1);
    }

    // Avoid too many bins
    if (binCount > (int)n) binCount = (int)n;

    // Bin edges
    double binWidth = (maxVal - minVal) / binCount;

    // Count values in each bin
    std::vector<int> counts(binCount, 0);
    for (double v : sorted) {
        int binIdx = (int)((v - minVal) / binWidth);
        if (binIdx == binCount) binIdx = binCount - 1; // edge case
        counts[binIdx]++;
    }

    // Find most populated bin
    int maxBinIdx = std::distance(counts.begin(),
                                  std::max_element(counts.begin(), counts.end()));
    int maxCount = counts[maxBinIdx];

    // Find second max count
    int secondMaxCount = 0;
    for (int c : counts) {
        if (c > secondMaxCount && c < maxCount) {
            secondMaxCount = c;
        }
    }

    double meanCount = 0.0;
    for (std::vector<double>::const_iterator it = counts.begin(); it != counts.end(); ++it) {
      meanCount += *it;
    }
    meanCount = meanCount / counts.size();

    double prominence = (double)(maxCount - secondMaxCount) / meanCount;
    result.prominence = prominence;

    // Check if mode is strong enough
    if (prominence < prominenceThreshold) {
        result.hasStrongMode = false;
        return result; // No strong mode
    }

    // Gather values in the most populated bin
    double binStart = minVal + maxBinIdx * binWidth;
    double binEnd = binStart + binWidth;
    std::vector<double> binValues;
    for (double v : sorted) {
        if ((v >= binStart && v < binEnd) ||
            (maxBinIdx == binCount - 1 && v == maxVal)) {
            binValues.push_back(v);
        }
    }

    // Compute mean of values in most populated bin
    double sum = 0.0;
    for (std::vector<double>::const_iterator it = binValues.begin(); it != binValues.end(); ++it) {
      sum += *it;
    }

    result.modeValue = sum / binValues.size();
    result.hasStrongMode = true;
    return result;
}

//int main() {
//    std::vector<double> data = {
//        3.1, 3.5, 3.6, 4.0, 4.1, 4.2, 5.5, 6.0, 6.1, 6.2
//    };
//
//    ModeResult modeEst = estimateModeHistogram(data);
//
//    if (modeEst.hasStrongMode) {
//        std::cout << "Estimated mode: " << modeEst.modeValue << "\n";
//        std::cout << "Prominence: " << modeEst.prominence * 100 << "% higher than average bin count\n";
//    } else {
//        std::cout << "No strong mode detected (prominence = "
//                  << modeEst.prominence * 100 << "%)\n";
//    }
//
//    return 0;
//}




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
    int cntr_mode = 0;
    int cntr_kde  = 0;

    af::versa<FloatType, GridType> result;
    result.resize(GridType(n_real), 0.0);
    for(int i = 0; i < n_real[0]; i++) {
      std::cout<<i<<std::endl;
      for(int j = 0; j < n_real[1]; j++) {
        for(int k = 0; k < n_real[2]; k++) {

          std::vector<FloatType> data;
          af::shared<FloatType> dataf = map_new(i,j,k);
          for(int N = 0; N < map_new(i,j,k).size(); N++) {
            data.push_back(dataf[N]);
          }

          ModeResult modeEst = estimateModeHistogram(data);
          //ModeResult modeEst = findModeWithKDEFallback(data, exp_table_);

          //if (modeEst.hasStrongMode) {
          //  result(i,j,k) =  modeEst.modeValue;
          //  cntr_mode += 1;
          //}
          //else {
          //  cntr_kde += 1;
          //  result(i,j,k) = 0;
          //}

          if (modeEst.hasStrongMode) {
            result(i,j,k) = modeEst.modeValue;
            cntr_mode += 1;
          }
          else {
            //FloatType sum = std::accumulate(data.begin(), data.end(), 0.0);
            //FloatType mean = sum / data.size();
            result(i,j,k) = 0.; //mean;
            //result(i,j,k) = estimateModeKDE(data, exp_table_); //0.; //mean;
            cntr_kde += 1;
          }
    }}}

    std::cout<<" cntr_mode "<< cntr_mode << std::endl;
    std::cout<<" cntr_kde  "<< cntr_kde << std::endl;

    return result;
  }

};



*/

























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
