#include <cctbx/boost_python/flex_fwd.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <cfloat>
#include <cctbx/uctbx/determine_unit_cell/NCDist.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace af = scitbx::af;

using namespace boost::python;

namespace cctbx { namespace uctbx {

  double NCDist_wrapper(af::tiny<double,6> mm1,af::tiny<double,6> mm2){
    return NCDist(&mm1[0],&mm2[0]);
  }

  scitbx::af::versa<double, scitbx::af::c_grid<2> > NCDist_matrix(scitbx::af::shared<double> MM){
    int NN = MM.size()/6;
    af::versa<double, af::c_grid<2> > result( af::c_grid<2> (NN,NN));
    double* MM_ptr = MM.begin();
    double* result_ptr = result.begin();

    # pragma omp parallel
    {
    # pragma omp for
    for (int i = 0; i < NN; ++i) {
      af::tiny<double,6> mm1(&(MM_ptr[i*6]), &(MM_ptr[(i+1)*6]));
      for (int j = i+1; j < NN; ++j) {
        af::tiny<double,6> mm2(&(MM_ptr[j*6]), &(MM_ptr[(j+1)*6]));
        double metric = NCDist(&mm1[0],&mm2[0]);
        result_ptr[i*NN + j] = metric;
        result_ptr[j*NN + i] = metric;
      }
    }
    }
    return result;
  }


}}

BOOST_PYTHON_MODULE(determine_unit_cell_ext)
{
  def ("NCDist",&cctbx::uctbx::NCDist_wrapper);
  def ("NCDist_matrix",&cctbx::uctbx::NCDist_matrix);
}
