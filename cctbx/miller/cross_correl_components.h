#ifndef CCTBX_MILLER_CROSS_CORREL_COMPONENTS_H
#define CCTBX_MILLER_CROSS_CORREL_COMPONENTS_H

#include <cctbx/miller.h>
#include <cmath>
#include <scitbx/array_family/shared.h>
#include <cctbx/miller/match_multi_indices.h>
#include <boost/python/object.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/extract.hpp>
#include <chrono>

namespace cctbx { namespace miller {

  struct cross_correlation_components {

    af::shared<int> count;
    af::shared<double> sum_xx, sum_yy, sum_xy, sum_x, sum_y;
    cross_correlation_components(af::shared<index<> > const &indices1,
                                 af::shared<index<> > const &indices2,
                                 af::shared<double> const &data1,
                                 af::shared<double> const &data2,
                                 boost::python::object hkl_bin_pydict,
                                 int n_bins) {
      auto t0 = std::chrono::high_resolution_clock::now();
      boost::python::extract<boost::python::dict> hkl_bin_dict_(hkl_bin_pydict);
      boost::python::dict hkl_bin_dict = hkl_bin_dict_;
      auto t_pydict = std::chrono::high_resolution_clock::now()-t0;
      count = af::shared<int>(n_bins, 0);
      sum_xx = af::shared<double>(n_bins, 0.0);
      sum_yy = af::shared<double>(n_bins, 0.0);
      sum_xy = af::shared<double>(n_bins, 0.0);
      sum_x = af::shared<double>(n_bins, 0.0);
      sum_y = af::shared<double>(n_bins, 0.0);
      match_multi_indices matching_indices = 
        match_multi_indices(indices1, indices2);
      for (int i=0; i<matching_indices.pairs().size(); i++) {
        pair_type pair = matching_indices.pairs()[i];
        int i1 = pair[0];
        int i2 = pair[1];
        index<> hkl = indices1[i1];
        CCTBX_ASSERT(indices2[i2] == hkl);
        //int i_bin = hkl_bin_dict[hkl];
        auto t1 = std::chrono::high_resolution_clock::now();
        if (! hkl_bin_dict.has_key(hkl)) {
          t_pydict += std::chrono::high_resolution_clock::now()-t1;
          continue;
        }
        int i_bin = boost::python::extract<int>(hkl_bin_dict[hkl]);
        t_pydict += std::chrono::high_resolution_clock::now()-t1;
        double I_x = data1[i1];
        double I_y = data2[i2];
        count[i_bin] += 1;
        sum_xx[i_bin] += std::pow(I_x, 2);
        sum_yy[i_bin] += std::pow(I_y, 2);
        sum_xy[i_bin] += I_x * I_y;
        sum_x[i_bin] += I_x;
        sum_y[i_bin] += I_y;
      }
      auto t_total = std::chrono::high_resolution_clock::now()-t0;
      std::cout<<"pydict: "<<std::chrono::duration_cast<std::chrono::nanoseconds>(t_pydict).count()<<", total: "<<std::chrono::duration_cast<std::chrono::nanoseconds>(t_total).count()<<std::endl;
    }

    af::shared<int> get_count() {
      return count;
    }
    af::shared<double> get_sum_xx() {
      return sum_xx;
    }
    af::shared<double> get_sum_xy() {
      return sum_xy;
    }
    af::shared<double> get_sum_yy() {
      return sum_yy;
    }
    af::shared<double> get_sum_x() {
      return sum_x;
    }
    af::shared<double> get_sum_y() {
      return sum_y;
    }

        
      

  };

}} // namespace cctbx::miller

#endif
