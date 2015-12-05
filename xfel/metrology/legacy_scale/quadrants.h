#ifndef XFEL_LEGACY_QUADRANTS_H
#define XFEL_LEGACY_QUADRANTS_H

#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/vec2.h>
#include <scitbx/mat2.h>
#include <vector>
#include <map>
#include <string>
#include <cstdio>
#include <boost/python/tuple.hpp>

namespace xfel_legacy {
typedef scitbx::af::shared<double> farray;
typedef scitbx::af::shared<int>    iarray;
typedef scitbx::af::shared<scitbx::vec3<double> > vec3array;
typedef scitbx::vec2<double>                      vec2;
typedef scitbx::vec2<int>                         vec2i;
typedef scitbx::mat2<double>                      mat2;

inline int iround(double const& x) {
  if (x < 0) {return static_cast<int>(x-0.5);}
  return static_cast<int>(x+.5);
}

//quadrant self correlation
boost::python::tuple
qsc(scitbx::af::flex_double asic,
    vec2 const& asic_origin,
    vec2 const& beam_center,
    mat2 const& rotmat,
    double const min_px_value){

  int F0 = asic.accessor().focus()[0];
  int F1 = asic.accessor().focus()[1];
  farray ref_data, rot_data;

  vec2 constant = rotmat*(asic_origin - beam_center) + beam_center - asic_origin;

  for (int xcoord = 0; xcoord < F0; ++xcoord){
    for (int ycoord = 0; ycoord < F1; ++ycoord){

      vec2 acoord(xcoord,ycoord);
      vec2 prime = rotmat*acoord + constant;
      vec2i primei(iround(prime[0]),iround(prime[1]));
      if (0<=primei[0] && primei[0]<F0 && 0<=primei[1] && primei[1]<F1) {
        if (asic[F1*xcoord+ycoord] >= min_px_value && asic[F1*primei[0] + primei[1]] >= min_px_value) {
          ref_data.push_back(asic[F1*xcoord+ycoord]);
          rot_data.push_back(asic[F1*primei[0] + primei[1]]);
        }
      }
    }
  }
  return boost::python::make_tuple(ref_data,rot_data);

}


} //namespace xfel
#endif// XFEL_LEGACY_QUADRANTS_H
