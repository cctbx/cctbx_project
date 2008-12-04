#ifndef IOTBX_SHELX_HKLF_H
#define IOTBX_SHELX_HKLF_H

#include <cctbx/miller.h>
#include <scitbx/array_family/shared.h>
#include <boost/algorithm/string.hpp>
#include <string>
#include <fstream>
#include <cstdlib>

namespace iotbx { namespace shelx {

  class fast_hklf_reader
  {
    public:
      typedef cctbx::miller::index<> miller_t;

      fast_hklf_reader(std::string const &filename, bool strict=true)
      {
        std::ifstream hkl_file(filename.c_str());
        std::string line;
        const miller_t end_hkl(0, 0, 0);
        typedef std::string::size_type index_t;
        std::runtime_error not_hklf_error("Not a SHELX hklf file.");
        while (std::getline(hkl_file, line)) {
          miller_t hkl;
          boost::trim_right(line);
          if (line.size() == 0) break;
          if (line.size() < 28) throw not_hklf_error;
          if (strict && line.size() > 32) throw not_hklf_error;
          /*if (!boost::all(line, boost::is_any_of("+-.")
               || boost::is_digit()
               || boost::is_space()))
          {
            throw not_hklf_error;
          }*/
          for (index_t i = 0; i != 3; i++) {
            hkl[i] = std::atoi(line.substr(i*4, 4).c_str());
          }
          double datum = std::atof(line.substr(12, 8).c_str());
          double sigma = std::atof(line.substr(20, 8).c_str());
          if (hkl == end_hkl) break;
          indices_.push_back(hkl);
          data_.push_back(datum);
          sigmas_.push_back(sigma);
          if (indices_.size() == 1) {
            try {
              int extra = std::atoi(line.substr(28, 4).c_str());
              extras_.push_back(extra);
            }
            catch(std::out_of_range) {}
          }
          else if (extras_.size()) {
            int extra = std::atoi(line.substr(28, 4).c_str());
            extras_.push_back(extra);
          }
        }
        if (indices_.size() == 0) throw std::runtime_error(
          "No data in SHELX hklf file.");
      }

      scitbx::af::shared<miller_t> indices() { return indices_; }
      scitbx::af::shared<double> data() { return data_; }
      scitbx::af::shared<double> sigmas() { return sigmas_; }
      scitbx::af::shared<int> alphas() { return extras_; }
      scitbx::af::shared<int> batch_numbers() { return extras_; }

    private:
      scitbx::af::shared<miller_t> indices_;
      scitbx::af::shared<double> data_, sigmas_;
      scitbx::af::shared<int> extras_;
  };

}} //iotbx::shelx

#endif
