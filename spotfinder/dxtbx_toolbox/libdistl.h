#ifndef DXTBX_WRAP_LIBDISTL_H
#define DXTBX_WRAP_LIBDISTL_H

#include <spotfinder/core_toolbox/libdistl.h>
#include <dxtbx/model/panel.h>
#include <dxtbx/model/beam.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>

using namespace std;

namespace Distl {

class dxtbx_diffimage: public diffimage {

  ::dxtbx::model::Panel panel;
  ::dxtbx::model::Beam beam;

public:
  inline dxtbx_diffimage(){}
  inline ~dxtbx_diffimage(){}

  inline void set_panel(::dxtbx::model::Panel other){
    panel = other;}

  inline void set_beam(::dxtbx::model::Beam& other) {
    beam = ::dxtbx::model::Beam(other);
  }

  inline double xy2resol(const double x, const double y) const {
    double resolution = panel.get_resolution_at_pixel(beam.get_s0(), scitbx::vec2<double>(x,y));
    return resolution;
  }

  void search_icerings();
};

}


#endif
