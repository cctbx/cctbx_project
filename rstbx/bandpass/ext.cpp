#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>

#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>
#include <cctbx/miller.h>
#include <cctbx/crystal_orientation.h>

namespace rstbx { namespace bandpass { namespace ext {

  static scitbx::af::shared<scitbx::vec3<double> >
  use_case_bp2_picture_fast_slow(
    scitbx::af::shared<cctbx::miller::index<> > indices, // full sphere Miller indices already computed to resolution limit
    cctbx::crystal_orientation const& orientation,
    scitbx::vec3<double> const& incident_beam,   // direction of incident radiation; not necessarily normalized
    double const& wavelength,                    // in Angstroms
    scitbx::vec3<double> const& detector_normal, // normal to the detector surface (toward incident beam), length 1
    scitbx::vec3<double> const& detector_fast,   // detector fast direction, length 1
    scitbx::vec3<double> const& detector_slow,   // detector slow direction, length 1
    scitbx::vec3<double> const& pixel_size,      // in (fast,slow,dummy) directions, mm
    scitbx::vec3<double> const& pixel_offset,    // for rendering in (fast,slow,dummy) directions, pixels
    double const& distance,                      // detector distance, mm
    scitbx::vec3<double> const& detector_origin, // coordinates of the origin pixel (fast,slow,dummy), mm
    double const& half_mosaicity_rad             // top-hat half mosaicity, radians
         ) {
    scitbx::af::shared<scitbx::vec3<double> > Z;
    scitbx::mat3<double> A = orientation.reciprocal_matrix();

    //s0:  parallel to the direction of incident radiation
    scitbx::vec3<double> s0 = (1./wavelength) * incident_beam;
    double s0_length = s0.length();
    scitbx::vec3<double> s0_unit = s0.normalize();

    //  Cn, the circular section through the Ewald sphere.
    for (int idx = 0; idx < indices.size(); ++idx){

        scitbx::vec3<double> H(indices[idx][0],indices[idx][1], indices[idx][2]); // the Miller index
        scitbx::vec3<double> s = A * H; //s, the reciprocal space coordinates, lab frame, of the oriented Miller index
        scitbx::vec3<double> rotax = s.normalize().cross(s0_unit); //The axis that most directly brings the Bragg spot onto Ewald sphere
        double s_rad_sq = s.length_sq();

        // take a page from ewald_sphere.cpp, determine intersection of two coplanar circles
        //  Co, the circle centered on reciprocal origin and containing the point s,
        //  Cn, the circle centered on -s0 (ewald sphere center) of radius (1/lambda) with normal rotax.
        // Consider the intersection of two circles:
        //  Co, the circle of rotation of H.
        // Cocenter = 0; so it falls out of the equations

        //   The chord of intersection between Co and Cn lies a
        //   distance x along the (Cocenter - Cncenter) vector
        scitbx::vec3<double> chord_direction =      (rotax.cross(s0)).normalize();

        double a = s.length_sq()/(2.*s0_length); // see diagram
        double b = std::sqrt(s.length_sq() - (a*a));   //  Calculate half-length of the chord of intersection
        //  Two intersection points
        scitbx::vec3<double> intersections_0p = -a * s0_unit+ b*chord_direction;
        scitbx::vec3<double> intersections_1p = -a * s0_unit- b*chord_direction;
        double iangle_0= std::acos ( (intersections_0p * s) / (s_rad_sq));
        double iangle_1= std::acos ( (intersections_1p * s) / (s_rad_sq));

        // assert approx_equal((intersections_0p+s0).length()-s0_length,0. )

        scitbx::vec3<double> intersection;
        if (iangle_0 < half_mosaicity_rad) {
          intersection = intersections_0p;
        } else if (iangle_1 < half_mosaicity_rad) {
          intersection = intersections_1p;
        } else {continue;}

        scitbx::vec3<double> q = (intersection + s0);
        scitbx::vec3<double> q_unit = q.normalize();

        // check if diffracted ray parallel to detector face

        double q_dot_n = q_unit * detector_normal;

        if (q_dot_n >= 0) { continue; }

        scitbx::vec3<double> r = (q_unit * distance / q_dot_n) - detector_origin;

        double x = r * detector_fast;
        double y = r * detector_slow;

        Z.push_back( scitbx::vec3<double> ( (x/pixel_size[0])+pixel_offset[0],(y/pixel_size[1])+pixel_offset[1],0. ));
    }
    return Z;
  }

  static boost::python::tuple
  use_case_bp3_picture_fast_slow(
    scitbx::af::shared<cctbx::miller::index<> > indices, // full sphere Miller indices already computed to resolution limit
    cctbx::crystal_orientation const& orientation,
    scitbx::vec3<double> const& incident_beam,   // direction of incident radiation; not necessarily normalized
    scitbx::vec3<double> const& tophat,          // top hat parameters. (high-energy wavelength in Angstroms,
                                                 // low-energy wavelength in Angstroms, top-hat half mosaicity, radians)
    scitbx::vec3<double> const& detector_normal, // normal to the detector surface (toward incident beam), length 1
    scitbx::vec3<double> const& detector_fast,   // detector fast direction, length 1
    scitbx::vec3<double> const& detector_slow,   // detector slow direction, length 1
    scitbx::vec3<double> const& pixel_size,      // in (fast,slow,dummy) directions, mm
    scitbx::vec3<double> const& pixel_offset,    // for rendering in (fast,slow,dummy) directions, pixels
    double const& distance,                      // detector distance, mm
    scitbx::vec3<double> const& detector_origin  // coordinates of the origin pixel (fast,slow,dummy), mm
         ) {

    double half_mosaicity_rad = tophat[2];
    double wavelengthHE = tophat[0];
    double wavelengthLE = tophat[1];
    SCITBX_ASSERT (wavelengthHE <= wavelengthLE);
    scitbx::af::shared<scitbx::vec3<double> > hi_E_limit(indices.size());
    scitbx::af::shared<scitbx::vec3<double> > lo_E_limit(indices.size());
    scitbx::af::shared<char > limit_types(indices.size());
    scitbx::af::shared<bool > observed_flag(indices.size());
    scitbx::mat3<double> A = orientation.reciprocal_matrix();

    //s0:  parallel to the direction of incident radiation
    scitbx::vec3<double> s0 = (1./wavelengthHE) * incident_beam;
    double s0_length = s0.length();
    scitbx::vec3<double> s0_unit = s0.normalize();
    scitbx::vec3<double> s1 = (1./wavelengthLE) * incident_beam;
    double s1_length = s1.length();
    scitbx::vec3<double> s1_unit = s1.normalize();

    //  Cn, the circular section through the Ewald sphere.
    for (int idx = 0; idx < indices.size(); ++idx){

        scitbx::vec3<double> H(indices[idx][0],indices[idx][1], indices[idx][2]); // the Miller index
        scitbx::vec3<double> s = A * H; //s, the reciprocal space coordinates, lab frame, of the oriented Miller index
        double s_rad_sq = s.length_sq();
        scitbx::vec3<double> rotax = s.normalize().cross(s0_unit); //The axis that most directly brings the Bragg spot onto Ewald sphere
        scitbx::vec3<double> chord_direction =      (rotax.cross(s0)).normalize();

       //  ###########  Look at the high-energy wavelength boundary

        double a = s.length_sq()/(2.*s0_length); // see diagram
        double b = std::sqrt(s.length_sq() - (a*a));   //  Calculate half-length of the chord of intersection

        scitbx::vec3<double> intersection = -a * s0_unit- b*chord_direction;
        double iangle_1= std::acos ( (intersection * s) / (s_rad_sq));

        // assert approx_equal((intersection+s0).length()-s0_length,0. )

        if (iangle_1 < half_mosaicity_rad) {

        scitbx::vec3<double> q = (intersection + s0);
        scitbx::vec3<double> q_unit = q.normalize();

        // check if diffracted ray parallel to detector face

        double q_dot_n = q_unit * detector_normal;

        scitbx::vec3<double> r = (q_unit * distance / q_dot_n) - detector_origin;

        double x = r * detector_fast;
        double y = r * detector_slow;

        limit_types[idx] += 1; // indicate that the high-energy boundary is found
        observed_flag[idx] = true;
        hi_E_limit[idx] = scitbx::vec3<double> ( (x/pixel_size[0])+pixel_offset[0],(y/pixel_size[1])+pixel_offset[1],0. );
        }
       //  ###########  Look at the low-energy wavelength boundary

        double alow = s.length_sq()/(2.*s1_length); // see diagram
        double blow = std::sqrt(s.length_sq() - (alow*alow));   //  Calculate half-length of the chord of intersection

        scitbx::vec3<double> intersectionlow = -alow * s0_unit- blow*chord_direction;
        double iangle_1low= std::acos ( (intersectionlow * s) / (s_rad_sq));

        // assert approx_equal((intersection+s0).length()-s0_length,0. )

        if (iangle_1low < half_mosaicity_rad) {

        scitbx::vec3<double> q = (intersectionlow + s1);
        scitbx::vec3<double> q_unit = q.normalize();

        // check if diffracted ray parallel to detector face

        double q_dot_n = q_unit * detector_normal;

        scitbx::vec3<double> r = (q_unit * distance / q_dot_n) - detector_origin;

        double x = r * detector_fast;
        double y = r * detector_slow;

        limit_types[idx] += 2; // indicate that the low-energy boundary is found
        observed_flag[idx] = true;
        lo_E_limit[idx] = scitbx::vec3<double> ( (x/pixel_size[0])+pixel_offset[0],(y/pixel_size[1])+pixel_offset[1],0. );
        }

       //  ###########  Look at rocking the crystal along rotax toward hiE reflection condition
        if (limit_types[idx]%2 == 0) { // ==3 or ==1 means that hiE test is unnecessary
          scitbx::vec3<double> s_rot_hi = s.rotate_around_origin(rotax,half_mosaicity_rad);
          double a_hi = -s_rot_hi * s0_unit;
          double r_n_hi = s_rad_sq/(2.*a_hi);
          double wavelength_hi = 1./r_n_hi;
          if (wavelengthHE < wavelength_hi && wavelength_hi < wavelengthLE) {
            scitbx::vec3<double> s0_hi = r_n_hi * incident_beam;
            scitbx::vec3<double> q = (s_rot_hi + s0_hi);
            scitbx::vec3<double> q_unit = q.normalize();

            double q_dot_n = q_unit * detector_normal;

            scitbx::vec3<double> r = (q_unit * distance / q_dot_n) - detector_origin;

            double x = r * detector_fast;
            double y = r * detector_slow;
            limit_types[idx] += 1; // indicate that the hi-energy boundary is found
            observed_flag[idx] = true;
            hi_E_limit[idx] = scitbx::vec3<double> ( (x/pixel_size[0])+pixel_offset[0],(y/pixel_size[1])+pixel_offset[1],0. );
          }
        }
       //  ###########  Look at rocking the crystal along rotax toward loE reflection condition
        if (limit_types[idx] < 2) { // >=2 means that loE test is unnecessary
          scitbx::vec3<double> s_rot_lo = s.rotate_around_origin(rotax,-half_mosaicity_rad);
          double a_lo = -s_rot_lo * s0_unit;
          double r_n_lo = s_rad_sq/(2.*a_lo);
          double wavelength_lo = 1./r_n_lo;
          if (wavelengthHE < wavelength_lo && wavelength_lo < wavelengthLE) {
            scitbx::vec3<double> s0_lo = r_n_lo * incident_beam;
            scitbx::vec3<double> q = (s_rot_lo + s0_lo);
            scitbx::vec3<double> q_unit = q.normalize();

            double q_dot_n = q_unit * detector_normal;

            scitbx::vec3<double> r = (q_unit * distance / q_dot_n) - detector_origin;

            double x = r * detector_fast;
            double y = r * detector_slow;
            limit_types[idx] += 2; // indicate that the hi-energy boundary is found
            observed_flag[idx] = true;
            lo_E_limit[idx] = scitbx::vec3<double> ( (x/pixel_size[0])+pixel_offset[0],(y/pixel_size[1])+pixel_offset[1],0. );
          }
        }
    }
    return boost::python::make_tuple(hi_E_limit,lo_E_limit,observed_flag);
  }

  struct bandpass_wrappers
  {

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef return_internal_reference<> rir;
      def("use_case_bp2_picture_fast_slow",&use_case_bp2_picture_fast_slow,
        (arg("indices"), arg("orientation"), arg("incident_beam"), arg("wavelength"),
         arg("detector_normal"), arg("detector_fast"),arg("detector_slow"),
         arg("pixel_size"), arg("pixel_offset"), arg("distance"),arg("detector_origin"),
         arg("half_mosaicity_rad")
         )
         );
      def("use_case_bp3_picture_fast_slow",&use_case_bp3_picture_fast_slow,
        (arg("indices"), arg("orientation"), arg("incident_beam"), arg("tophat"),
         arg("detector_normal"), arg("detector_fast"),arg("detector_slow"),
         arg("pixel_size"), arg("pixel_offset"), arg("distance"),arg("detector_origin")
         )
         );

    }
  };


  void init_module()
  {
    using namespace boost::python;
    bandpass_wrappers::wrap();
  }

}}} // namespace rstbx::bandpass::ext

BOOST_PYTHON_MODULE(rstbx_bandpass_ext)
{
  rstbx::bandpass::ext::init_module();
}
