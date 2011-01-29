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

    hklf_reader(
      af::const_ref<std::string> const& lines,
      bool strict=true)
    {
      for(std::size_t i_line=0;i_line<lines.size();i_line++) {
        std::string line = lines[i_line];
        std::size_t i_trailing = 3*4+2*8;
        miller_t h;
        double datum, sigma;
        int extra = 0;
        bool have_extra = false;
        if (substr_is_whitespace_only(
              line, i_trailing, std::min(i_trailing+4, line.size()))) {
          prepare_for_read(line, 28);
          fem::read_from_string(line, "(3i4,2f8.0)"),
            h[0], h[1], h[2], datum, sigma;
        }
        else {
          i_trailing += 4;
          have_extra = true;
          prepare_for_read(line, 32);
          fem::read_from_string(line, "(3i4,2f8.0,i4)"),
            h[0], h[1], h[2], datum, sigma, extra;
        }
        if (h.is_zero()) break;
        if (strict
              && !substr_is_whitespace_only(line, i_trailing, line.size())) {
          throw std::runtime_error("Not a SHELX hklf file.");
        }
        indices_.push_back(h);
        data_.push_back(datum);
        sigmas_.push_back(sigma);
        if (have_extra) {
          extra_.push_back(extra);
        }
      }
      if (indices_.size() == 0) {
        throw std::runtime_error("No data in SHELX hklf file.");
      }
    }

    af::shared<miller_t> indices() { return indices_; };

    af::shared<double> data() { return data_; }

    af::shared<double> sigmas() { return sigmas_; }

    af::shared<int> alphas() { return extra_; }

    af::shared<int> batch_numbers() { return extra_; }

  private:
    af::shared<miller_t> indices_;
    af::shared<double> data_, sigmas_;
    af::shared<int> extra_; // batch numbers or phases
};

}} // iotbx::shelx

#endif // GUARD
