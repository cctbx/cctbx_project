#ifndef RSTBX_DETECTOR_MODEL_SENSOR_H
#define RSTBX_DETECTOR_MODEL_SENSOR_H
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>

namespace rstbx { namespace detector_model {

class sensor {

 public:

  sensor(){}

  sensor(const scitbx::vec3<double>& origin,
         const scitbx::vec3<double>& dir1,
         const scitbx::vec3<double>& dir2,
         const scitbx::vec2<double>& lim1,
         const scitbx::vec2<double>& lim2);

  //getters
  double get_distance() const {return distance;}
  scitbx::vec3<double> get_origin() const {return origin;}
  scitbx::vec3<double> get_normal() const {return normal;}
  scitbx::vec3<double> get_dir1() const {return dir1;}
  scitbx::vec3<double> get_dir2() const {return dir2;}
  scitbx::vec2<double> get_lim1() const {return lim1;}
  scitbx::vec2<double> get_lim2() const {return lim2;}
  scitbx::mat3<double> get_d() const {return d;}
  scitbx::mat3<double> get_D() const {return D;}
  bool get_d_is_invertible() const {return d_is_invertible;}

  //setters

  /* Each must call update(). set all three
     properties of the sensor at once (i.e. dir1, dir2, origin) as we
     want to maintain the same "shape" coordinate frame. */

  void set_origin(const scitbx::vec3<double>& origin);

  void set_frame(const scitbx::vec3<double>& origin,
                 const scitbx::vec3<double>& dir1,
                 const scitbx::vec3<double>& dir2);

 private:

  //members
  scitbx::vec3<double> origin;
  scitbx::vec3<double> dir1;
  scitbx::vec3<double> dir2;
  scitbx::vec2<double> lim1;
  scitbx::vec2<double> lim2;

  scitbx::vec3<double> normal;
  scitbx::mat3<double> d;
  scitbx::mat3<double> D;
  bool d_is_invertible;
  double distance;

  //Update the d matrix, according to Lure I-II Contribution 8
  //(D. Thomas) amongst other things. Call this from any set method.
  void update();

};

}} //namespace rstbx::detector_model
#endif //RSTBX_DETECTOR_MODEL_SENSOR_H
