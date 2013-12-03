#ifndef RSTBX_INDEXING_API_H
#define RSTBX_INDEXING_API_H
#include <boost/shared_ptr.hpp>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <cctbx/miller.h>
#include <cctbx/uctbx.h>
#include <cctbx/crystal_orientation.h>
#include <rstbx/dps_core/dps_core.h>
#include <annlib_adaptbx/ann_adaptor.h>
#include <dxtbx/model/detector.h>

namespace af = scitbx::af;
namespace rstbx {

  typedef af::const_ref<scitbx::vec3<double> > pointlist;

namespace indexing_api {

af::shared<int> cpp_absence_test(af::shared<cctbx::miller::index<> >,
                                 const int&, cctbx::miller::index<>);

struct dps_extended: public rstbx::dps_core {
  pointlistmm rawdata;// original spots, film coords x(mm),y(mm),phi(degrees)

  dps_extended();
  dps_extended(const dps_extended&){} // copy constructor (gives
    //a duplicate engine with parameters but no data)


  // XXX rename these methods (like "set/get raw_input_positions_mm")
  void setData(const pointlist &);
  inline pointlistmm getData(){return rawdata;}

  // XXX not used now, but use this for tweaking of the basis set
  Direction refine_direction(const Direction& candidate,
                             const double& current_grid,
                             const double& target_grid) const;

  // both of these unused?  XXX deprecate???
  double high()const;          //high-resolution limit defined by the provided rawdata
  double model_likelihood(double)const;//same as above except the random residual (mm)
                            //is replaced by the input of choice, usually the
                            //closest spot-to-spot separation.
};

af::shared< scitbx::vec3<double> > raw_spot_positions_mm_to_reciprocal_space_xyz(
  pointlist,
  dxtbx::model::Detector const&, double const&,
  scitbx::vec3<double> const& , scitbx::vec3<double> const&, af::shared<int>
);
af::shared< scitbx::vec3<double> > raw_spot_positions_mm_to_reciprocal_space_xyz(
  pointlist,
  dxtbx::model::Detector const&, double const&,
  scitbx::vec3<double> const& , af::shared<int>
);

} //namespace

} //namespace

#endif //RSTBX_INDEXING_API_H
