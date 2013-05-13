#ifndef XFEL_PARAMETER_H
#define XFEL_PARAMETER_H

#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>
#include <vector>
#include <map>
#include <string>
#include <cstdio>
#include <cctbx/miller.h>
#include <xfel/mono_simulation/vector_collection.h>

namespace xfel {
typedef scitbx::af::shared<double> farray;
typedef scitbx::af::shared<int>    iarray;
typedef scitbx::af::shared<scitbx::vec3<double> > vec3array;
typedef scitbx::vec3<double>                      vec3;
typedef scitbx::af::shared<cctbx::miller::index<> > marray;

namespace parameter {

struct parameter_array {
  farray parameters;
  farray gradients;
  farray curvatures;
  int ndata;
  bool local_flag;

  inline parameter_array(){}
  parameter_array(int const&, int const&, farray);
  inline int size()const{return parameters.size();}
};

struct organizer_base {
  farray
  as_x_array() const;

  farray
  get_gradient_array() const;

  void set_gradient_array(std::string const&,farray);

  farray
  get_curvature_array() const;

  void
  from_x_array(farray const&);

  std::map<std::string,parameter_array> _P;

  parameter_array
  register_array
  (std::string const&, int const&, int const&, farray);

  parameter_array
  register_local_array
  (std::string const&, int const&, int const&, farray);

  void initialize_gradients_curvatures();
  void rezero_gradients_curvatures();

};
}
namespace algorithm {

struct mark5_iteration: public parameter::organizer_base{
  std::vector<double> sine,cosine;
  farray model_calcx, model_calcy;
  double calc_minus_To_x, calc_minus_To_y;
  double rotated_o_x, rotated_o_y;
  double partial_partial_theta_x, partial_partial_theta_y;
  double partial_sq_theta_x, partial_sq_theta_y;
  double functional;
  vec3array frame_origins;

  void set_refined_origins_to_c(vec3array);

  double
  compute_target(
    farray tox, farray toy, farray spotcx, farray spotcy,
    farray spotfx, farray spotfy,
    iarray master_tiles,iarray frames,
    vec3array partial_r_partial_distance);
  double
  compute_functional_only(
    farray tox, farray toy, farray spotcx, farray spotcy,
    farray spotfx, farray spotfy,
    iarray master_tiles,iarray frames,
    vec3array partial_r_partial_distance);
  inline mark5_iteration(){}

  xfel::parameter::vector_collection vecc;
  inline void set_vector_collection(xfel::parameter::vector_collection& v){
    vecc=v;
  }

  vec3array
  uncorrected_detector_to_laboratory_frame(
    farray tox, farray toy,
    farray spotfx, farray spotfy,
    iarray master_tiles) const;

};

};

} //namespace xfel
#endif// XFEL_PARAMETER_H
