#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <cctbx/miller.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <cctbx/miller/match_indices.h>
#include <scitbx/math/linear_correlation.h>

namespace af = scitbx::af;

using namespace boost::python;

namespace cctbx { namespace merging {

  void update_wij_rij(int const& i, int const& j,
       af::shared<cctbx::miller::index<int> > millers_i, af::shared<cctbx::miller::index<int> > millers_j,
       af::shared<double> data_i, af::shared<double> data_j,
       int const& i_start, int const& j_start,
       af::versa< double, af::flex_grid<> > wij_, af::versa< double, af::flex_grid<> > rij_,
       double const& sign, bool const& use_weights
                      ){
    cctbx::miller::match_indices matches(millers_i, millers_j);
    af::shared<double> intensities_i, intensities_j;
    for (int ipair = 0; ipair < matches.pairs().size(); ++ipair){
      intensities_i.push_back(data_i[i_start + matches.pairs()[ipair][0]]);
      intensities_j.push_back(data_j[j_start + matches.pairs()[ipair][1]]);
    }
    scitbx::math::linear_correlation<> corr(intensities_i.const_ref(), intensities_j.const_ref());
    //printf("%5d <--> %5d, matches=%4d, corr=%7.4f,%4d\n",i,j,matches.pairs().size(),corr.coefficient(),
      //wij_.accessor().focus()[0]);
    if (corr.is_well_defined()){
      int size_1d = wij_.accessor().focus()[0];
      double* wij_begin = wij_.begin();
      double* rij_begin = rij_.begin();
      if (use_weights) {
        wij_begin[i*size_1d + j] += corr.n();
        wij_begin[j*size_1d + i] += corr.n();
      } else {
        wij_begin[i*size_1d + j] += 1.;
        wij_begin[j*size_1d + i] += 1.;
      }
      rij_begin[i*size_1d + j] += sign * corr.coefficient();
      rij_begin[j*size_1d + i] += sign * corr.coefficient();
    }
  };

}}

BOOST_PYTHON_MODULE(cctbx_merging_ext)
{
  def ("update_wij_rij",&cctbx::merging::update_wij_rij);
}
