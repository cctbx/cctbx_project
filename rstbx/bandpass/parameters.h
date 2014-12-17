#ifndef RSTBX_BANDPASS_PARAMETERS_H
#define RSTBX_BANDPASS_PARAMETERS_H
#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <cctbx/miller.h>
#include <cctbx/crystal_orientation.h>
#include <annlib_adaptbx/ann_adaptor.h>

namespace rstbx { namespace bandpass {
  typedef scitbx::vec2<double> vec2;
  typedef scitbx::vec3<double> vec3;
  typedef vec3 const& vec3ref;

  struct parameters_bp3 {
      // full sphere Miller indices already computed to resolution limit
    scitbx::af::shared<cctbx::miller::index<> > indices;
    cctbx::crystal_orientation orientation;
      // direction of incident radiation;not necessarily normalized
    vec3 incident_beam;
    double half_mosaicity_rad;
    double wavelengthHE;
    double wavelengthLE;
    vec3 detector_normal;//unit normal to detector surface toward incident beam
    vec3 detector_fast;  // detector fast direction, length 1
    vec3 detector_slow;  // detector slow direction, length 1
    vec3 pixel_size;     // in (fast,slow,dummy) directions, mm
    vec3 pixel_offset;   // to render in (fast,slow,dummy) directions, pixels
    double distance;     // detector distance, mm
    vec3 detector_origin;// mm coordinates of the origin pixel (fast,slow,dummy)

    inline
    parameters_bp3 (scitbx::af::shared<cctbx::miller::index<> > indices,
                    cctbx::crystal_orientation const& orientation,
                    vec3ref incident_beam,
                    vec3ref packed_tophat,
                    vec3ref detector_normal,
                    vec3ref detector_fast,
                    vec3ref detector_slow,
                    vec3ref pixel_size,
                    vec3ref pixel_offset,
                    double const& distance,
                    vec3ref detector_origin):
         indices(indices), orientation(orientation),
         incident_beam(incident_beam), half_mosaicity_rad(packed_tophat[2]),
         wavelengthHE(packed_tophat[0]),
         wavelengthLE (packed_tophat[1]),
         detector_normal(detector_normal),
         detector_fast(detector_fast),detector_slow(detector_slow),
         pixel_size(pixel_size), pixel_offset(pixel_offset),
         distance(distance),detector_origin(detector_origin) {
          SCITBX_ASSERT (wavelengthHE <= wavelengthLE);
          SCITBX_ASSERT (wavelengthHE>0.);
         }
  };

  struct pad_sensor_model {
    double thickness_mm;
    double mu_rho;      //product of mass attenuation * amorphous silicon density in units of 1/mm
    inline
    pad_sensor_model():
    thickness_mm(0.5),mu_rho(8.36644){} //ersatz values for CS-PAD detector at 1.3 Angstrom

    inline
    vec3
    sensor_coords_in_pixels(double const& signal_penetration, parameters_bp3 const& camera,
                            vec3 const& q_unit, double const& q_dot_n){
      // signal penetration is the fraction of signal height lost at the returned position;
      // sensor top==0; sensor bottom approaches 1 at low Energy
      double path_thru_sensor = -thickness_mm/q_dot_n;
      double target_distance_thru_sensor = -std::log( 1. - signal_penetration )/mu_rho;
      double observed_thru_sensor_path = std::min(path_thru_sensor,target_distance_thru_sensor);
      vec3 diffraction_vec = (q_unit * camera.distance / q_dot_n);
      diffraction_vec = diffraction_vec * (1. + observed_thru_sensor_path/diffraction_vec.length());
      vec3 r_mod = diffraction_vec - camera.detector_origin;
      double x_mod = r_mod * camera.detector_fast;
      double y_mod = r_mod * camera.detector_slow;
      return vec3 ( (x_mod/camera.pixel_size[0])+camera.pixel_offset[0],
                    (y_mod/camera.pixel_size[1])+camera.pixel_offset[1],0. );
    }
  };

  struct active_area_filter{
    int NEAR;
    scitbx::af::shared<int> tiles;
    annlib_adaptbx::AnnAdaptor adapt;
    int tile_id; /*looks like an unsafe practice. Could be accessed before being
                   set.  Temporary workaround by setting it to zero first,
                   however this doesn't really fix the problem; tile_id could contain
                   a stale value but still be accessed. XXX A better fix might involve
                   the active_area_filter() function returning a struct containing
                   a bool is_active_area and a tile tile_id. */
    scitbx::af::shared<vec2> centers;
    inline
    active_area_filter(){}
    inline
    active_area_filter(scitbx::af::shared<int> IT): NEAR(2),tiles(IT),tile_id(0){

      af::shared<ANNcoord> reference;
      for (int i=0; i < tiles.size()/4; ++i){
        vec2 UL (static_cast<double>(tiles[4*i]),static_cast<double>(tiles[4*i+1]));
        vec2 LR (static_cast<double>(tiles[4*i+2]),static_cast<double>(tiles[4*i+3]));
        vec2 center = (UL+LR)/2.;
        reference.push_back(center[0]);
        reference.push_back(center[1]);
        centers.push_back(center);
      }
      adapt = annlib_adaptbx::AnnAdaptor(reference,2,NEAR);
    }
    inline
    bool operator()(vec3 prediction, int const& TOLERANCE=0){
      // if no tiles are set, assume all pixels are active
      if (tiles.size() == 0)
        return true;

      af::flex_int nearest_neighbours;
      if (tiles.size() == 4){
        // We have only one tile, AnnAdaptor chokes in this case but then there is
        // only one choice of nearest neighbour anyway!
        nearest_neighbours = af::flex_int (NEAR, 0);
      } else {
        af::shared<ANNcoord> query;
        query.push_back(prediction[0]);
        query.push_back(prediction[1]);
        adapt.query(query);
        SCITBX_ASSERT( adapt.nn.size()== NEAR);
        nearest_neighbours = adapt.nn;
      }

      bool is_in_active_area = false;
      for (int n = 0; n< NEAR;++n){
        int itile = nearest_neighbours[n];
        if (tiles[4*itile]-TOLERANCE <= prediction[0] && prediction[0] <= TOLERANCE + tiles[4*itile+2] &&
            tiles[4*itile+1]-TOLERANCE <= prediction[1] &&prediction[1] <= TOLERANCE + tiles[4*itile+3]){
          is_in_active_area = true;tile_id=itile;break;
          }
      }
      return is_in_active_area;
    }
  };


}} // namespace rstbx::bandpass

#endif// RSTBX_BANDPASS_PARAMETERS_H
