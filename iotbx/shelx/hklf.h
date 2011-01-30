#ifndef IOTBX_SHELX_HKLF_H
#define IOTBX_SHELX_HKLF_H

#include <cctbx/miller.h>
#include <scitbx/array_family/shared.h>
#include <fable/fem/read.hpp>

namespace iotbx { namespace shelx {

namespace af = scitbx::af;

class hklf_reader
{
  public:
    typedef char char_t;
    typedef cctbx::miller::index<> miller_t;

  protected:
    af::shared<miller_t> indices_;
    af::shared<double> data_;
    af::shared<double> sigmas_;
    af::shared<int> batch_numbers_or_phases_;
    af::shared<double> wavelengths_;

    static
    void
    prepare_for_read(
      std::string& line,
      std::size_t target_size)
    {
      std::size_t initial_size = line.size();
      for(std::size_t i=0;i<initial_size;i++) {
        if (line[i] < ' ') line[i] = ' '; // emulates shelxl
      }
      if (initial_size < target_size) {
        line.append(target_size-initial_size, ' ');
      }
    }

    static
    bool
    substr_is_whitespace_only(
      std::string const& line,
      std::size_t start,
      std::size_t stop)
    {
      for(std::size_t i=start;i<stop;i++) {
        if (!fem::utils::is_whitespace(line[i])) return false;
      }
      return true;
    }

  public:
    hklf_reader(
      af::const_ref<std::string> const& lines)
    {
      std::size_t n = lines.size();
      indices_.reserve(n);
      data_.reserve(n);
      sigmas_.reserve(n);
      batch_numbers_or_phases_.reserve(n);
      wavelengths_.reserve(n);
      bool have_bp = false;
      bool have_wa = false;
      for(std::size_t i=0;i<n;i++) {
        std::string line = lines[i];
        miller_t h;
        double datum, sigma, wa;
        int bp;
        prepare_for_read(line, 40);
        fem::read_from_string(line, "(3i4,2f8.0,i4,f8.4)"),
          h[0], h[1], h[2], datum, sigma, bp, wa;
        if (h.is_zero()) break;
        indices_.push_back(h);
        data_.push_back(datum);
        sigmas_.push_back(sigma);
        batch_numbers_or_phases_.push_back(bp);
        wavelengths_.push_back(wa);
        if (!have_bp) have_bp = !substr_is_whitespace_only(line, 28, 32);
        if (!have_wa) have_wa = !substr_is_whitespace_only(line, 32, 40);
      }
      if (indices_.size() == 0) {
        throw std::runtime_error("No data in SHELX hklf file.");
      }
      if (!have_bp) batch_numbers_or_phases_ = af::shared<int>();// free memory
      if (!have_wa) wavelengths_ = af::shared<double>();// free memory
    }

    af::shared<miller_t> indices() { return indices_; };

    af::shared<double> data() { return data_; }

    af::shared<double> sigmas() { return sigmas_; }

    af::shared<int> alphas() { return batch_numbers_or_phases_; }

    af::shared<int> batch_numbers() { return batch_numbers_or_phases_; }

    af::shared<double> wavelengths() { return wavelengths_; }
};

}} // iotbx::shelx

#endif // GUARD
