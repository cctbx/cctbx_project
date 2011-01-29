#ifndef IOTBX_SHELX_HKLF_SIMPLE_H
#define IOTBX_SHELX_HKLF_SIMPLE_H

#include <cctbx/miller.h>
#include <scitbx/array_family/shared.h>
#include <boost/algorithm/string.hpp>
#include <string>
#include <cstdlib>

namespace iotbx { namespace shelx {

  class simple_hklf_reader
  {
    public:
      typedef cctbx::miller::index<> miller_t;

      simple_hklf_reader(std::string const &content, bool strict=true)
      {
        miller_t const end_hkl(0, 0, 0);
        std::runtime_error not_hklf_error("Not a SHELX hklf file.");
        typedef boost::split_iterator<std::string::const_iterator> split_iter_t;
        for(split_iter_t i(content,
                           boost::token_finder(boost::is_any_of("\r\n"),
                                               boost::token_compress_on));
            i != split_iter_t();
            ++i)
        {
          std::string line = boost::copy_range<std::string>(*i);
          miller_t hkl;
          boost::trim_right(line);
          if (line.size() == 0) break;
          if (line.size() < 12) throw not_hklf_error;
          if (strict && line.size() > 32) throw not_hklf_error;
          /* That would make a strong check but it kills performance
          if (!boost::all(line, boost::is_any_of("+-.")
               || boost::is_digit()
               || boost::is_space()))
          {
            throw not_hklf_error;
          }*/
          for (std::string::size_type i=0; i != 3; i++) {
            hkl[i] = std::atoi(line.substr(i*4, 4).c_str());
          }
          if (hkl == end_hkl) break;
          if (line.size() < 28) throw not_hklf_error;
          double datum = std::atof(line.substr(12, 8).c_str());
          double sigma = std::atof(line.substr(20, 8).c_str());
          indices_.push_back(hkl);
          data_.push_back(datum);
          sigmas_.push_back(sigma);
          if (line.size() > 28) {
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
