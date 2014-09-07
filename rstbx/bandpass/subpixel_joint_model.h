#ifndef RSTBX_BANDPASS_SUBPIXEL_JOINT_MODEL_H
#define RSTBX_BANDPASS_SUBPIXEL_JOINT_MODEL_H
#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <scitbx/mat2.h>

namespace rstbx { namespace bandpass {
  typedef scitbx::vec2<double> vec2;
  typedef scitbx::vec3<double> vec3;
  typedef scitbx::mat2<double> mat2;
  typedef vec3 const& vec3ref;

  struct subpixel_joint_model {
      // a description of CSPAD subpixel translations and rotations
      // along with math to convert back and forth between
      // fictitious space where pixels are grid-aligned and
      // laboratory space where each asic is individually placed

    scitbx::af::shared<double> subpixel;
    scitbx::af::shared<double> rotations_rad;
    scitbx::af::shared<mat2> L2F_R;
    void set_subpixel(scitbx::af::shared<double> s, scitbx::af::shared<double> rotations_deg){
      subpixel=s;
      rotations_rad=scitbx::af::shared<double>();
      for (int ixx=0; ixx< rotations_deg.size(); ++ixx){
        rotations_rad.push_back(scitbx::constants::pi_180*rotations_deg[ixx]);
      }
      SCITBX_ASSERT( s.size() == 2 * rotations_rad.size());
      // Now figure out the rotation matrices
      for (int ixx=0; ixx< rotations_deg.size(); ++ixx){
        L2F_R.push_back(mat2 ( std::cos(-rotations_rad[ixx]), -std::sin(-rotations_rad[ixx]),
                               std::sin(-rotations_rad[ixx]), std::cos(-rotations_rad[ixx])));
      }
    }

    subpixel_joint_model(scitbx::af::shared<double> s, scitbx::af::shared<double> rotations_deg){
      set_subpixel(s,rotations_deg);
    }
    subpixel_joint_model(){}

    vec3
    laboratory_to_fictitious(vec3ref rhs, int const& tile_id, vec2 const& tile_center){
      vec3 subpixel_trans(subpixel[2*tile_id],subpixel[1+2*tile_id],0.0);
      vec3 translated = rhs + subpixel_trans;

      //rotate the spot around tile center.  But these tile centers have xy coords swapped
      vec2 spot_in_tile_frame = vec2(translated[0]-tile_center[1],translated[1]-tile_center[0]);
      vec2 rotated = L2F_R[tile_id] * spot_in_tile_frame;
      return vec3( rotated[0] + tile_center[1], rotated[1] + tile_center[0], translated[2] );

    }

  };

}} // namespace rstbx::bandpass

#endif// RSTBX_BANDPASS_SUBPIXEL_JOINT_MODEL_H
