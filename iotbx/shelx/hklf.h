#ifndef IOTBX_SHELX_HKLF_H
#define IOTBX_SHELX_HKLF_H

#include <scitbx/array_family/shared.h>
#include <cctbx/miller.h>
#include <scitbx/fortran_io/numeric_manipulators.h>
#include <scitbx/misc/file_utils.h>
#include <istream>
#include <iostream>

namespace iotbx { namespace shelx {

class hklf_reader
{
  public:
    typedef char char_t;
    typedef cctbx::miller::index<> miller_t;

    hklf_reader(std::istream &input, bool strict=true)
    {
      using namespace scitbx::fortran_io::manipulators;
      fortran_int i4(4, /*strict=*/true);
      fortran_real f8_2(8, 2, /*strict=*/true);

      std::runtime_error not_hklf("Not a SHELX hklf file.");
      while(!input.eof()) {
        miller_t h;
        double datum, sigma;
        int extra;
        std::string trailing;
        input >> sticky(i4) >> h[0] >> h[1] >> h[2];
        if (h.is_zero()) break;
        input >> sticky(f8_2) >> datum >> sigma;
        if (!input.good()) throw not_hklf;
        indices_.push_back(h);
        data_.push_back(datum);
        sigmas_.push_back(sigma);
        scitbx::misc::end_of_line eol(input);
        if (eol) continue;
        input >> i4 >> extra;
        std::getline(input, trailing);
        if (!input.good()) throw not_hklf;
        if (strict && trailing.size() && trailing != "\r") throw not_hklf;
        extra_.push_back(extra);
      }
      std::runtime_error empty_hklf("No data in SHELX hklf file.");
      if (indices_.size() == 0) throw empty_hklf;
    }

    scitbx::af::shared<miller_t> indices() { return indices_; };

    scitbx::af::shared<double> data() { return data_; }

    scitbx::af::shared<double> sigmas() { return sigmas_; }

    scitbx::af::shared<int> alphas() { return extra_; }

    scitbx::af::shared<int> batch_numbers() { return extra_; }

  private:
    scitbx::af::shared<miller_t> indices_;
    scitbx::af::shared<double> data_, sigmas_;
    scitbx::af::shared<int> extra_; // batch numbers or phases
};

}} // iotbx::shelx

#endif // GUARD
