#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>

#include <algorithm>
#include <rstbx/bandpass/parameters.h>
#include <rstbx/bandpass/subpixel_joint_model.h>

namespace rstbx { namespace bandpass {

  struct use_case_bp3 {
    parameters_bp3 P;
    pad_sensor_model sensor;
    double signal_penetration;
    scitbx::af::shared<vec3 > hi_E_limit;
    scitbx::af::shared<vec3 > lo_E_limit;
    scitbx::af::shared<bool > observed_flag;
    scitbx::af::shared<vec3 > spot_rectangle_vertices;
    scitbx::af::shared<vec3 > spot_rectregion_vertices;
    scitbx::af::shared<vec3 > part_distance;
    use_case_bp3 (parameters_bp3 const& P):P(P),subpixel_translations_set(false){set_ellipse_model();}
    active_area_filter aaf;
    void set_active_areas(scitbx::af::shared<int> IT){
      aaf = active_area_filter(IT);
    }
    void set_sensor_model(double const& thickness_mm, double const& mu_rho, double const& sp){
      sensor.thickness_mm = thickness_mm;
      sensor.mu_rho = mu_rho;
      signal_penetration = sp;
    }
    void
    prescreen_indices(double const& wavelength){
      // The effect of this function is to eliminate most of the HKLs in P.indices
      // based on Ewald proximity, to make repeat calculations much faster.
      scitbx::af::shared<cctbx::miller::index<> > indices_subset;
      scitbx::mat3<double> A = P.orientation.reciprocal_matrix();

      //s0:  parallel to the direction of incident radiation
      scitbx::vec3<double> s0 = (1./wavelength) * P.incident_beam;
      double s0_length = s0.length();

      for (int idx = 0; idx < P.indices.size(); ++idx){
          // the Miller index
          scitbx::vec3<double> H(P.indices[idx][0],P.indices[idx][1], P.indices[idx][2]);
          //s, the reciprocal space coordinates, lab frame, of the oriented Miller index
          scitbx::vec3<double> s = A * H;

          scitbx::vec3<double> q = (s + s0);
          double q_len = q.length();
          double ratio = q_len/s0_length;
          if (ratio > 0.96 && ratio < 1.04) indices_subset.push_back(P.indices[idx]);
      }
      //SCITBX_EXAMINE(P.indices.size());
      //SCITBX_EXAMINE(indices_subset.size());
      P.indices = indices_subset;
    }

    int ellip_gran;
    double ellip_vol;
    void set_ellipse_model(){
      // ad hoc initialization to calculate partialities
      ellip_gran = 20;
      ellip_vol = 0.;
      for (double x = -1. + 0.5 * (2./ellip_gran); x < 1. ; x += 2./ellip_gran){
        for (double y = -1. + 0.5 * (2./ellip_gran); y < 1. ; y += 2./ellip_gran){
          for (double z = -1. + 0.5 * (2./ellip_gran); z < 1. ; z += 2./ellip_gran){
            if (x*x + y*y + z*z <= 1.0){
              ellip_vol += 1.;
            }}}}
      /* Crude model for partialities, treating Bragg spots as ellipsoids of rotation.
         This initialization normalizes the volume of a sphere when pixelated within
         a box of dimensions ellip_gran x ellip_gran x ellip_gran */
    }

    scitbx::af::shared<double >
    selected_partialities()const{
      scitbx::af::shared<double > data;
      scitbx::mat3<double> A = P.orientation.reciprocal_matrix();

      //s0:  parallel to the direction of incident radiation
      scitbx::vec3<double> s0 = (1./P.wavelengthHE) * P.incident_beam.normalize();
      double s0_length = s0.length();
      scitbx::vec3<double> s0_unit = s0.normalize();
      scitbx::vec3<double> s1 = (1./P.wavelengthLE) * P.incident_beam.normalize();
      double s1_length = s1.length();
      SCITBX_ASSERT (s0_length > 0.);
      SCITBX_ASSERT (s1_length > 0.);

      for (int idx = 0; idx < lo_E_limit.size(); ++idx){
        if (!observed_flag[idx]) {continue;}
        scitbx::vec3<double> H(P.indices[idx][0],P.indices[idx][1], P.indices[idx][2]); // the Miller index
        double spot_resolution_ang = P.orientation.unit_cell().d(P.indices[idx]);

        // constructing the reciprocal space model of the spot
        // Use new notation, reciprocal vector is qvec
        scitbx::vec3<double> qvec = A * H; //qvec, the reciprocal space coordinates, lab frame,
                                           // of the oriented Miller index
        scitbx::vec3<double> qnorm= qvec.normalize();
        scitbx::vec3<double> rotax = qnorm.cross(s0_unit); //Tangent to Ewald sphere at the Bragg spot
        scitbx::vec3<double> chord_direction =      (rotax.cross(s0)).normalize();
        // These last three vectors are taken as principle unit axes of an ellipsoid of rotation
        // used to model the spot.  The semimajor axes in these directions are respecively,
        // half_domain, half_domain+half_mosaic, half_domain+half_mosaic
        double a_semi = 0.5/p_domain_size_ang;
        double b_semi = a_semi + P.half_mosaicity_rad/spot_resolution_ang;
        double test_volume = 0;
        for (double x = -1. + 0.5 * (2./ellip_gran); x < 1. ; x += 2./ellip_gran){
          for (double y = -1. + 0.5 * (2./ellip_gran); y < 1. ; y += 2./ellip_gran){
            for (double z = -1. + 0.5 * (2./ellip_gran); z < 1. ; z += 2./ellip_gran){
              if (x*x + y*y + z*z <= 1.0){
                scitbx::vec3<double> test_vector = qvec+ x * a_semi * qnorm +
                                                   y * b_semi * rotax + z * b_semi * chord_direction;
                // if vector is in between the two limiting spheres, increment test_volume
                scitbx::vec3<double> test_scatterHE = test_vector + s0;
                scitbx::vec3<double> test_scatterLE = test_vector + s1;

                if (test_scatterHE.length() < (1./P.wavelengthHE) &&
                    test_scatterLE.length() > (1./P.wavelengthLE)){
                  test_volume += 1.;
                }
              }}}}

        data.push_back( test_volume/ellip_vol );
      }
      return data;
    }

    void
    picture_fast_slow(){
      //The effect of this function is to set the following three state vectors:
      hi_E_limit = scitbx::af::shared<scitbx::vec3<double> >(P.indices.size());
      lo_E_limit = scitbx::af::shared<scitbx::vec3<double> >(P.indices.size());
      observed_flag = scitbx::af::shared<bool >(P.indices.size());

      scitbx::af::shared<char > limit_types(P.indices.size());
      scitbx::mat3<double> A = P.orientation.reciprocal_matrix();

      //s0:  parallel to the direction of incident radiation
      scitbx::vec3<double> s0 = (1./P.wavelengthHE) * P.incident_beam;
      double s0_length = s0.length();
      scitbx::vec3<double> s0_unit = s0.normalize();
      scitbx::vec3<double> s1 = (1./P.wavelengthLE) * P.incident_beam;
      double s1_length = s1.length();
      SCITBX_ASSERT (s0_length > 0.);
      SCITBX_ASSERT (s1_length > 0.);

      //  Cn, the circular section through the Ewald sphere.
      for (int idx = 0; idx < P.indices.size(); ++idx){
          double spot_resolution_ang = P.orientation.unit_cell().d(P.indices[idx]);
          double effective_half_mosaicity_rad = (p_domain_size_ang > 0.)? P.half_mosaicity_rad + spot_resolution_ang / (2. * p_domain_size_ang) : P.half_mosaicity_rad;

          scitbx::vec3<double> H(P.indices[idx][0],P.indices[idx][1], P.indices[idx][2]); // the Miller index
          scitbx::vec3<double> s = A * H; //s, the reciprocal space coordinates, lab frame, of the oriented Miller index
          double s_rad_sq = s.length_sq();
          SCITBX_ASSERT(s_rad_sq > 0.);
          scitbx::vec3<double> rotax = s.normalize().cross(s0_unit); //The axis that most directly brings the Bragg spot onto Ewald sphere
          scitbx::vec3<double> chord_direction =      (rotax.cross(s0)).normalize();

         //  ###########  Look at the high-energy wavelength boundary

          double a = s.length_sq()/(2.*s0_length); // see diagram
          double b = std::sqrt(s.length_sq() - (a*a));   //  Calculate half-length of the chord of intersection

          scitbx::vec3<double> intersection = -a * s0_unit- b*chord_direction;

          double acos_argument = (intersection * s) / (s_rad_sq);
          double iangle_1= std::acos ( std::min(1.0,acos_argument) );//avoid math domain error
          // assert approx_equal((intersection+s0).length()-s0_length,0. )

          if (iangle_1 < effective_half_mosaicity_rad) {

          scitbx::vec3<double> q = (intersection + s0);
          scitbx::vec3<double> q_unit = q.normalize();

          // check if diffracted ray parallel to detector face

          double q_dot_n = q_unit * P.detector_normal;
          SCITBX_ASSERT(q_dot_n != 0.);
          scitbx::vec3<double> r = (q_unit * P.distance / q_dot_n) - P.detector_origin;

          double x = r * P.detector_fast;
          double y = r * P.detector_slow;

          limit_types[idx] += 1; // indicate that the high-energy boundary is found
          observed_flag[idx] = true;
          hi_E_limit[idx] = scitbx::vec3<double> ( (x/P.pixel_size[0])+P.pixel_offset[0],(y/P.pixel_size[1])+P.pixel_offset[1],0. );
          }
         //  ###########  Look at the low-energy wavelength boundary

          double alow = s.length_sq()/(2.*s1_length); // see diagram
          double blow = std::sqrt(s.length_sq() - (alow*alow));   //  Calculate half-length of the chord of intersection

          scitbx::vec3<double> intersectionlow = -alow * s0_unit- blow*chord_direction;
          acos_argument = (intersectionlow * s) / (s_rad_sq);
          double iangle_1low= std::acos ( std::min(1.0,acos_argument) );//avoid math domain error

          // assert approx_equal((intersection+s0).length()-s0_length,0. )

          if (iangle_1low < effective_half_mosaicity_rad) {

          scitbx::vec3<double> q = (intersectionlow + s1);
          scitbx::vec3<double> q_unit = q.normalize();

          // check if diffracted ray parallel to detector face

          double q_dot_n = q_unit * P.detector_normal;
          SCITBX_ASSERT(q_dot_n != 0.);

          limit_types[idx] += 2; // indicate that the low-energy boundary is found
          observed_flag[idx] = true;
          lo_E_limit[idx] = sensor.sensor_coords_in_pixels(signal_penetration, P, q_unit, q_dot_n);
          }

         //  ###########  Look at rocking the crystal along rotax toward hiE reflection condition
          if (limit_types[idx]%2 == 0) { // ==3 or ==1 means that hiE test is unnecessary
            scitbx::vec3<double> s_rot_hi = s.rotate_around_origin(rotax,effective_half_mosaicity_rad);
            double a_hi = -s_rot_hi * s0_unit;
            SCITBX_ASSERT(a_hi != 0.);
            double r_n_hi = s_rad_sq/(2.*a_hi);
            double wavelength_hi = 1./r_n_hi;
            if (P.wavelengthHE < wavelength_hi && wavelength_hi < P.wavelengthLE) {
              scitbx::vec3<double> s0_hi = r_n_hi * P.incident_beam;
              scitbx::vec3<double> q = (s_rot_hi + s0_hi);
              scitbx::vec3<double> q_unit = q.normalize();

              double q_dot_n = q_unit * P.detector_normal;
              SCITBX_ASSERT(q_dot_n != 0.);
              scitbx::vec3<double> r = (q_unit * P.distance / q_dot_n) - P.detector_origin;

              double x = r * P.detector_fast;
              double y = r * P.detector_slow;
              limit_types[idx] += 1; // indicate that the hi-energy boundary is found
              observed_flag[idx] = true;
              hi_E_limit[idx] = scitbx::vec3<double> ( (x/P.pixel_size[0])+P.pixel_offset[0],(y/P.pixel_size[1])+P.pixel_offset[1],0. );
            }
          }
         //  ###########  Look at rocking the crystal along rotax toward loE reflection condition
          if (limit_types[idx] < 2) { // >=2 means that loE test is unnecessary
            scitbx::vec3<double> s_rot_lo = s.rotate_around_origin(rotax,-effective_half_mosaicity_rad);
            double a_lo = -s_rot_lo * s0_unit;
            SCITBX_ASSERT(a_lo != 0.);
            double r_n_lo = s_rad_sq/(2.*a_lo);
            double wavelength_lo = 1./r_n_lo;
            if (P.wavelengthHE < wavelength_lo && wavelength_lo < P.wavelengthLE) {
              scitbx::vec3<double> s0_lo = r_n_lo * P.incident_beam;
              scitbx::vec3<double> q = (s_rot_lo + s0_lo);
              scitbx::vec3<double> q_unit = q.normalize();

              double q_dot_n = q_unit * P.detector_normal;
              SCITBX_ASSERT(q_dot_n != 0.);

              limit_types[idx] += 2; // indicate that the hi-energy boundary is found
              observed_flag[idx] = true;
              lo_E_limit[idx] = sensor.sensor_coords_in_pixels(signal_penetration, P, q_unit, q_dot_n);
            }
          }
          if (observed_flag[idx]) {
            //Do some further tests to determine if the spot is within the flagged active area with peripheral margin
            vec3 central_position = ((lo_E_limit[idx]+hi_E_limit[idx])/2.);//already in pixel units
            if (!aaf(vec3(central_position[1],central_position[0],0.))){
              observed_flag[idx]=false;
            } else if (subpixel_translations_set) {
              lo_E_limit[idx] = sjm.laboratory_to_fictitious(lo_E_limit[idx], aaf.tile_id, aaf.centers[aaf.tile_id]);
              hi_E_limit[idx] = sjm.laboratory_to_fictitious(hi_E_limit[idx], aaf.tile_id, aaf.centers[aaf.tile_id]);
            }
          }
      }
    }
    void
    picture_fast_slow_force(){ //DEPRECATED; obsolete; sets translations but not rotations
      //The effect of this function is to set the following three state vectors:
      hi_E_limit = scitbx::af::shared<scitbx::vec3<double> >(P.indices.size(),scitbx::af::init_functor_null<scitbx::vec3<double> >());
      lo_E_limit = scitbx::af::shared<scitbx::vec3<double> >(P.indices.size(),scitbx::af::init_functor_null<scitbx::vec3<double> >());
      part_distance = scitbx::af::shared<scitbx::vec3<double> >(P.indices.size(),scitbx::af::init_functor_null<scitbx::vec3<double> >());
      observed_flag = scitbx::af::shared<bool >(P.indices.size());

      scitbx::af::shared<char > limit_types(P.indices.size());
      scitbx::mat3<double> A = P.orientation.reciprocal_matrix();

      //s0:  parallel to the direction of incident radiation
      scitbx::vec3<double> s0 = (1./P.wavelengthHE) * P.incident_beam;
      double s0_length = s0.length();
      scitbx::vec3<double> s0_unit = s0.normalize();
      scitbx::vec3<double> s1 = (1./P.wavelengthLE) * P.incident_beam;
      double s1_length = s1.length();
      SCITBX_ASSERT (s0_length > 0.);
      SCITBX_ASSERT (s1_length > 0.);
      scitbx::vec3<double> hi_E_part_r_part_distance = s0;//initialize only; values not used
      scitbx::vec3<double> lo_E_part_r_part_distance = s1;

      //  Cn, the circular section through the Ewald sphere.
      for (int idx = 0; idx < P.indices.size(); ++idx){
          double spot_resolution_ang = P.orientation.unit_cell().d(P.indices[idx]);
          double effective_half_mosaicity_rad = (p_domain_size_ang > 0.)? P.half_mosaicity_rad + spot_resolution_ang / (2. * p_domain_size_ang) : P.half_mosaicity_rad;

          scitbx::vec3<double> H(P.indices[idx][0],P.indices[idx][1], P.indices[idx][2]); // the Miller index
          scitbx::vec3<double> s = A * H; //s, the reciprocal space coordinates, lab frame, of the oriented Miller index
          double s_rad_sq = s.length_sq();
          SCITBX_ASSERT(s_rad_sq > 0.);
          scitbx::vec3<double> rotax = s.normalize().cross(s0_unit); //The axis that most directly brings the Bragg spot onto Ewald sphere
          scitbx::vec3<double> chord_direction =      (rotax.cross(s0)).normalize();

         //  ###########  Look at the high-energy wavelength boundary

          double a = s.length_sq()/(2.*s0_length); // see diagram
          double b = std::sqrt(s.length_sq() - (a*a));   //  Calculate half-length of the chord of intersection

          scitbx::vec3<double> intersection = -a * s0_unit- b*chord_direction;

          double acos_argument = (intersection * s) / (s_rad_sq);
          double iangle_1= std::acos ( std::min(1.0,acos_argument) );//avoid math domain error
          // assert approx_equal((intersection+s0).length()-s0_length,0. )

          if (iangle_1 < effective_half_mosaicity_rad) {

          scitbx::vec3<double> q = (intersection + s0);
          scitbx::vec3<double> q_unit = q.normalize();

          // check if diffracted ray parallel to detector face

          double q_dot_n = q_unit * P.detector_normal;
          SCITBX_ASSERT(q_dot_n != 0.);
          scitbx::vec3<double> r = (q_unit * P.distance / q_dot_n) - P.detector_origin;

          double x = r * P.detector_fast;
          double y = r * P.detector_slow;

          limit_types[idx] += 1; // indicate that the high-energy boundary is found
          observed_flag[idx] = true;
          hi_E_limit[idx] = scitbx::vec3<double> ( (x/P.pixel_size[0])+P.pixel_offset[0],(y/P.pixel_size[1])+P.pixel_offset[1],0. );
          scitbx::vec3<double> part_r_part_d( q_unit / q_dot_n );
          hi_E_part_r_part_distance = scitbx::vec3<double>( (part_r_part_d * P.detector_fast) / P.pixel_size[0],
                                                            (part_r_part_d * P.detector_slow) / P.pixel_size[1],0.);
          }
         //  ###########  Look at the low-energy wavelength boundary

          double alow = s.length_sq()/(2.*s1_length); // see diagram
          double blow = std::sqrt(s.length_sq() - (alow*alow));   //  Calculate half-length of the chord of intersection

          scitbx::vec3<double> intersectionlow = -alow * s0_unit- blow*chord_direction;
          acos_argument = (intersectionlow * s) / (s_rad_sq);
          double iangle_1low= std::acos ( std::min(1.0,acos_argument) );//avoid math domain error

          // assert approx_equal((intersection+s0).length()-s0_length,0. )

          if (iangle_1low < effective_half_mosaicity_rad) {

          scitbx::vec3<double> q = (intersectionlow + s1);
          scitbx::vec3<double> q_unit = q.normalize();

          // check if diffracted ray parallel to detector face

          double q_dot_n = q_unit * P.detector_normal;
          SCITBX_ASSERT(q_dot_n != 0.);

          limit_types[idx] += 2; // indicate that the low-energy boundary is found
          observed_flag[idx] = true;
          lo_E_limit[idx] = sensor.sensor_coords_in_pixels(signal_penetration, P, q_unit, q_dot_n);
          scitbx::vec3<double> part_r_part_d( q_unit / q_dot_n );
          lo_E_part_r_part_distance = scitbx::vec3<double>( (part_r_part_d * P.detector_fast) / P.pixel_size[0],
                                                            (part_r_part_d * P.detector_slow) / P.pixel_size[1],0.);
          }

         //  ###########  Look at rocking the crystal along rotax toward hiE reflection condition
          if (limit_types[idx]%2 == 0) { // ==3 or ==1 means that hiE test is unnecessary
            scitbx::vec3<double> s_rot_hi = s.rotate_around_origin(rotax,effective_half_mosaicity_rad);
            double a_hi = -s_rot_hi * s0_unit;
            SCITBX_ASSERT(a_hi != 0.);
            double r_n_hi = s_rad_sq/(2.*a_hi);
            double wavelength_hi = 1./r_n_hi;
            if (P.wavelengthHE < wavelength_hi && wavelength_hi < P.wavelengthLE) {
              scitbx::vec3<double> s0_hi = r_n_hi * P.incident_beam;
              scitbx::vec3<double> q = (s_rot_hi + s0_hi);
              scitbx::vec3<double> q_unit = q.normalize();

              double q_dot_n = q_unit * P.detector_normal;
              SCITBX_ASSERT(q_dot_n != 0.);
              scitbx::vec3<double> r = (q_unit * P.distance / q_dot_n) - P.detector_origin;

              double x = r * P.detector_fast;
              double y = r * P.detector_slow;
              limit_types[idx] += 1; // indicate that the hi-energy boundary is found
              observed_flag[idx] = true;
              hi_E_limit[idx] = scitbx::vec3<double> ( (x/P.pixel_size[0])+P.pixel_offset[0],(y/P.pixel_size[1])+P.pixel_offset[1],0. );
              scitbx::vec3<double> part_r_part_d( q_unit / q_dot_n );
              hi_E_part_r_part_distance = scitbx::vec3<double>( (part_r_part_d * P.detector_fast) / P.pixel_size[0],
                                                                (part_r_part_d * P.detector_slow) / P.pixel_size[1],0.);
            }
          }
         //  ###########  Look at rocking the crystal along rotax toward loE reflection condition
          if (limit_types[idx] < 2) { // >=2 means that loE test is unnecessary
            scitbx::vec3<double> s_rot_lo = s.rotate_around_origin(rotax,-effective_half_mosaicity_rad);
            double a_lo = -s_rot_lo * s0_unit;
            SCITBX_ASSERT(a_lo != 0.);
            double r_n_lo = s_rad_sq/(2.*a_lo);
            double wavelength_lo = 1./r_n_lo;
            if (P.wavelengthHE < wavelength_lo && wavelength_lo < P.wavelengthLE) {
              scitbx::vec3<double> s0_lo = r_n_lo * P.incident_beam;
              scitbx::vec3<double> q = (s_rot_lo + s0_lo);
              scitbx::vec3<double> q_unit = q.normalize();

              double q_dot_n = q_unit * P.detector_normal;
              SCITBX_ASSERT(q_dot_n != 0.);

              limit_types[idx] += 2; // indicate that the hi-energy boundary is found
              observed_flag[idx] = true;
              lo_E_limit[idx] = sensor.sensor_coords_in_pixels(signal_penetration, P, q_unit, q_dot_n);
              scitbx::vec3<double> part_r_part_d( q_unit / q_dot_n );
              lo_E_part_r_part_distance = scitbx::vec3<double>( (part_r_part_d * P.detector_fast) / P.pixel_size[0],
                                                                (part_r_part_d * P.detector_slow) / P.pixel_size[1],0.);
            }
          }
          if (!observed_flag[idx]) {
            //bogus position in the case where reflection is out of bounds

            //s0:  parallel to the direction of incident radiation
            double meanwave = (P.wavelengthHE+P.wavelengthLE)/2.;
            scitbx::vec3<double> s0_mean = (1./meanwave) * P.incident_beam;
            double s0_mean_length = s0_mean.length();
            scitbx::vec3<double> s0_mean_unit = s0_mean.normalize();


            double a = s.length_sq()/(2.*s0_mean_length); // see diagram
            double b = std::sqrt(s.length_sq() - (a*a));   //  Calculate half-length of the chord of intersection

            scitbx::vec3<double> intersection = -a * s0_mean_unit- b*chord_direction;

            scitbx::vec3<double> q = (intersection + s0_mean);
            scitbx::vec3<double> q_unit = q.normalize();

            // check if diffracted ray parallel to detector face

            double q_dot_n = q_unit * P.detector_normal;
            SCITBX_ASSERT(q_dot_n != 0.);
            scitbx::vec3<double> r = (q_unit * P.distance / q_dot_n) - P.detector_origin;

            double x = r * P.detector_fast;
            double y = r * P.detector_slow;

            observed_flag[idx] = true;
            hi_E_limit[idx] = scitbx::vec3<double> ( (x/P.pixel_size[0])+P.pixel_offset[0],(y/P.pixel_size[1])+P.pixel_offset[1],0. );
            lo_E_limit[idx] = hi_E_limit[idx];
            scitbx::vec3<double> part_r_part_d( q_unit / q_dot_n );
            lo_E_part_r_part_distance = scitbx::vec3<double>( (part_r_part_d * P.detector_fast) / P.pixel_size[0],
                                                              (part_r_part_d * P.detector_slow) / P.pixel_size[1],0.);
            hi_E_part_r_part_distance = lo_E_part_r_part_distance;
          }
          if (observed_flag[idx]) {
            if (subpixel_translations_set) {
              SCITBX_EXAMINE("NOT EXPECTED TO WORK BECAUSE aaf.tile_id NOT SET FOR THIS REFLECTION");
              vec3 subpixel_trans(subpixel[2*aaf.tile_id],subpixel[1+2*aaf.tile_id],0.0);
              lo_E_limit[idx] += subpixel_trans;
              hi_E_limit[idx] += subpixel_trans;
            }
            part_distance[idx] = (lo_E_part_r_part_distance + hi_E_part_r_part_distance)/2.;
          }
      }
      for ( int idx = 0; idx < P.indices.size(); ++idx){
        SCITBX_ASSERT (observed_flag[idx]);
      }
    }
    scitbx::af::shared<vec3 >
    spot_rectangles(vec3ref beam_coor){
      vec3 beam_pos(
        beam_coor[0]/P.pixel_size[0]+P.pixel_offset[0],
         beam_coor[1]/P.pixel_size[1]+P.pixel_offset[1],0.);
      scitbx::af::shared<vec3 > polydata;
      vec3 crystal_to_detector = -P.distance * P.detector_normal;
      SCITBX_ASSERT (P.pixel_size[0] == P.pixel_size[1]);
      // need to assert this because all calculations including detector distance
      // are done in units of pixels; if not equal then formulae will have to
      // be reimplemented in units of mm.
      for (int idx = 0; idx < lo_E_limit.size(); ++idx){
        if (!observed_flag[idx]) {continue;}
        vec3 hi_pos = hi_E_limit[idx];
        vec3 lo_pos = lo_E_limit[idx];
        vec3 radial_vector = (hi_pos-beam_pos);
        vec3 radial_unit_vec = radial_vector.normalize();
        double radius = radial_vector.length();
        vec3 tangential_unit_vec(-radial_unit_vec[1],radial_unit_vec[0],0.); // 90-degree rotation
        vec3 tangential_excursion = tangential_unit_vec * radius * P.half_mosaicity_rad;
        if (p_domain_size_ang > 0.) {
          double half_scherrer_broadening = ((crystal_to_detector.length())/P.pixel_size[0]) * P.wavelengthHE / (2.*p_domain_size_ang);
          tangential_excursion += tangential_unit_vec * half_scherrer_broadening;
          hi_pos -= half_scherrer_broadening*radial_unit_vec;
          lo_pos += half_scherrer_broadening*radial_unit_vec;
        }
        polydata.push_back( (hi_pos + tangential_excursion)-P.pixel_offset);
        polydata.push_back( (hi_pos - tangential_excursion)-P.pixel_offset);
        polydata.push_back( (lo_pos - tangential_excursion)-P.pixel_offset);
        polydata.push_back( (lo_pos + tangential_excursion)-P.pixel_offset);
        polydata.push_back( (hi_pos + tangential_excursion)-P.pixel_offset);
      }
      spot_rectangle_vertices = polydata;
      return polydata;
    }
    double margin;
    scitbx::af::shared<vec3 >
    spot_rectregions(vec3ref beam_coor,double const& MARGIN){
      margin = MARGIN;
      vec3 beam_pos(
        beam_coor[0]/P.pixel_size[0]+P.pixel_offset[0],
         beam_coor[1]/P.pixel_size[1]+P.pixel_offset[1],0.);
      scitbx::af::shared<vec3 > polydata;
      vec3 crystal_to_detector = -P.distance * P.detector_normal;

      for (int idx = 0; idx < lo_E_limit.size(); ++idx){
        if (!observed_flag[idx]) {continue;}
        vec3 hi_pos = hi_E_limit[idx];
        vec3 lo_pos = lo_E_limit[idx];
        vec3 radial_vector = (hi_pos-beam_pos);
        vec3 radial_unit_vec = radial_vector.normalize();
        double radius = radial_vector.length();
        hi_pos -= MARGIN*radial_unit_vec;
        lo_pos += MARGIN*radial_unit_vec;
        vec3 tangential_unit_vec(-radial_unit_vec[1],radial_unit_vec[0],0.); // 90-degree rotation
        vec3 tangential_excursion = tangential_unit_vec * (radius * P.half_mosaicity_rad + MARGIN);
        if (p_domain_size_ang > 0.) {
          double half_scherrer_broadening = (crystal_to_detector + radius).length() * P.wavelengthHE / (2.*p_domain_size_ang);
          tangential_excursion += tangential_unit_vec * half_scherrer_broadening;
          hi_pos -= half_scherrer_broadening*radial_unit_vec;
          lo_pos += half_scherrer_broadening*radial_unit_vec;
        }
        polydata.push_back( hi_pos + tangential_excursion);
        polydata.push_back( hi_pos - tangential_excursion);
        polydata.push_back( lo_pos - tangential_excursion);
        polydata.push_back( lo_pos + tangential_excursion);
        polydata.push_back( hi_pos + tangential_excursion);
      }
      spot_rectregion_vertices = polydata;
      return polydata;
    }

    scitbx::af::shared<vec3 >
    enclosed_pixels()const{
      // calculation is done in picture_fast_slow coordinates (units of pixels)
      scitbx::af::shared<vec3 > points;

      for (int idx = 0; idx < spot_rectangle_vertices.size(); idx+=5){
        //unroll everthing to set boundary box.  Must be a more concise expression of this?
        double slow_min = std::min<double>(spot_rectangle_vertices[idx][0],spot_rectangle_vertices[idx+1][0]);
               slow_min = std::min<double>(spot_rectangle_vertices[idx+2][0],slow_min);
               slow_min = std::min<double>(spot_rectangle_vertices[idx+3][0],slow_min);
        double slow_max = std::max<double>(spot_rectangle_vertices[idx][0],spot_rectangle_vertices[idx+1][0]);
               slow_max = std::max<double>(spot_rectangle_vertices[idx+2][0],slow_max);
               slow_max = std::max<double>(spot_rectangle_vertices[idx+3][0],slow_max);
        double fast_min = std::min<double>(spot_rectangle_vertices[idx][1],spot_rectangle_vertices[idx+1][1]);
               fast_min = std::min<double>(spot_rectangle_vertices[idx+2][1],fast_min);
               fast_min = std::min<double>(spot_rectangle_vertices[idx+3][1],fast_min);
        double fast_max = std::max<double>(spot_rectangle_vertices[idx][1],spot_rectangle_vertices[idx+1][1]);
               fast_max = std::max<double>(spot_rectangle_vertices[idx+2][1],fast_max);
               fast_max = std::max<double>(spot_rectangle_vertices[idx+3][1],fast_max);
        for (int islow = (int)std::floor(slow_min); islow < (int)std::ceil(slow_max); ++islow){
          for (int ifast = (int)std::floor(fast_min); ifast < (int)std::ceil(fast_max); ++ifast){
            vec3 candidate_pixel(islow+0.5,ifast+0.5,0.0);
            //algorithm to determine if the pixel is enclosed in the rectangle
            bool inside = false;
            double p1x = spot_rectangle_vertices[idx][0];
            double p1y = spot_rectangle_vertices[idx][1];
            double xinters= -10000000.;

            for (int iv = 0; iv < 5; ++iv){
              double p2x = spot_rectangle_vertices[idx+iv][0];
              double p2y = spot_rectangle_vertices[idx+iv][1];
              if (candidate_pixel[1] > std::min(p1y, p2y)){
                  if (candidate_pixel[1] <= std::max(p1y, p2y)){
                      if (candidate_pixel[0] <= std::max(p1x, p2x)){
                          if (p1y != p2y){
                              xinters = (candidate_pixel[1]-p1y)*(p2x-p1x)/(p2y-p1y) + p1x;
                          }
                          SCITBX_ASSERT(xinters!=-10000000.);
                          if (p1x == p2x || candidate_pixel[0] <= xinters){
                              inside = !inside;
                          }
                      }
                  }
              }
              p1x = p2x;
              p1y = p2y;
            }
            if (!inside){continue;}
            points.push_back(vec3(islow+0.5,ifast+0.5,0.0));
          }
        }
      }
      return points;
    }

    scitbx::af::shared<vec3 > margin_px;
    scitbx::af::shared<double > margin_distances;
    scitbx::af::shared<vec3 > enclosed_px;
    struct unlimited_pixel_accept_policy {
      static bool accept_pixel(vec3 const& vec,use_case_bp3* const thisuc){
        return true;
      }
    };
    struct resolution_limited_pixel_accept_policy {
      static bool accept_pixel(vec3 const& vec,use_case_bp3* const thisuc){
        vec3 pixel_center = vec*thisuc->P.pixel_size;
        vec3 offs = thisuc->P.detector_origin+pixel_center;
        double radius_mm = std::sqrt(offs[0]*offs[0]+offs[1]*offs[1]);
        double pixel_two_theta_rad = std::atan(radius_mm/thisuc->P.distance);
        double pixel_d_ang =   thisuc->P.wavelengthLE / (2.*std::sin (pixel_two_theta_rad/2.)) ;
        return (pixel_d_ang > thisuc->model_refinement_limiting_resolution);
      }
    };
    double model_refinement_limiting_resolution;
    void
    enclosed_pixels_and_margin_pixels(double const& model_refinement_limiting_resolution){
      //use this pointer to override the automatic variable
      this->model_refinement_limiting_resolution = model_refinement_limiting_resolution;
      enclosed_pixels_and_margin_pixels_detail<resolution_limited_pixel_accept_policy>();
    }
    void
    enclosed_pixels_and_margin_pixels(){
      enclosed_pixels_and_margin_pixels_detail<unlimited_pixel_accept_policy>();
    }
    template <typename pixel_accept_policy>
    void
    enclosed_pixels_and_margin_pixels_detail(){
      // The effect of this function is to set these three state vectors:
      margin_px = scitbx::af::shared<vec3 > ();
      margin_distances = scitbx::af::shared<double > ();
      enclosed_px = scitbx::af::shared<vec3 > ();
      // calculation is done in picture_fast_slow coordinates (units of pixels)

      for (int idx = 0; idx < spot_rectangle_vertices.size(); idx+=5){
        //unroll everthing to set boundary box.  Must be a more concise expression of this?
        double slow_min = std::min<double>(spot_rectregion_vertices[idx][0],spot_rectregion_vertices[idx+1][0]);
               slow_min = std::min<double>(spot_rectregion_vertices[idx+2][0],slow_min);
               slow_min = std::min<double>(spot_rectregion_vertices[idx+3][0],slow_min);
        double slow_max = std::max<double>(spot_rectregion_vertices[idx][0],spot_rectregion_vertices[idx+1][0]);
               slow_max = std::max<double>(spot_rectregion_vertices[idx+2][0],slow_max);
               slow_max = std::max<double>(spot_rectregion_vertices[idx+3][0],slow_max);
        double fast_min = std::min<double>(spot_rectregion_vertices[idx][1],spot_rectregion_vertices[idx+1][1]);
               fast_min = std::min<double>(spot_rectregion_vertices[idx+2][1],fast_min);
               fast_min = std::min<double>(spot_rectregion_vertices[idx+3][1],fast_min);
        double fast_max = std::max<double>(spot_rectregion_vertices[idx][1],spot_rectregion_vertices[idx+1][1]);
               fast_max = std::max<double>(spot_rectregion_vertices[idx+2][1],fast_max);
               fast_max = std::max<double>(spot_rectregion_vertices[idx+3][1],fast_max);
        for (int islow = (int)std::floor(slow_min); islow < (int)std::ceil(slow_max); ++islow){
          for (int ifast = (int)std::floor(fast_min); ifast < (int)std::ceil(fast_max); ++ifast){
            vec3 candidate_pixel(islow+0.5,ifast+0.5,0.0);
            //algorithm to determine if the pixel is enclosed in the rectangle
            bool inside_margin = false;
            double p1x = spot_rectregion_vertices[idx][0];
            double p1y = spot_rectregion_vertices[idx][1];
            double xinters= -10000000.;

            for (int iv = 0; iv < 5; ++iv){
              double p2x = spot_rectregion_vertices[idx+iv][0];
              double p2y = spot_rectregion_vertices[idx+iv][1];
              if (candidate_pixel[1] > std::min(p1y, p2y)){
                  if (candidate_pixel[1] <= std::max(p1y, p2y)){
                      if (candidate_pixel[0] <= std::max(p1x, p2x)){
                          if (p1y != p2y){
                              xinters = (candidate_pixel[1]-p1y)*(p2x-p1x)/(p2y-p1y) + p1x;
                          }
                          SCITBX_ASSERT(xinters!=-10000000.);
                          if (p1x == p2x || candidate_pixel[0] <= xinters){
                              inside_margin = !inside_margin;
                          }
                      }
                  }
              }
              p1x = p2x;
              p1y = p2y;
            }
            if (!inside_margin){continue;}
            // Don't bother with this pixel if it is not within active area
            if (!aaf(vec3(candidate_pixel[1],candidate_pixel[0],0.))){continue;}

            //algorithm to determine if the pixel is enclosed in the rectangle
            bool inside_rectangle = false;
            {
            double p1x = spot_rectangle_vertices[idx][0];
            double p1y = spot_rectangle_vertices[idx][1];
            double xinters= -10000000.;

            for (int iv = 0; iv < 5; ++iv){
              double p2x = spot_rectangle_vertices[idx+iv][0];
              double p2y = spot_rectangle_vertices[idx+iv][1];
              if (candidate_pixel[1] > std::min(p1y, p2y)){
                  if (candidate_pixel[1] <= std::max(p1y, p2y)){
                      if (candidate_pixel[0] <= std::max(p1x, p2x)){
                          if (p1y != p2y){
                              xinters = (candidate_pixel[1]-p1y)*(p2x-p1x)/(p2y-p1y) + p1x;
                          }
                          SCITBX_ASSERT(xinters!=-10000000.);
                          if (p1x == p2x || candidate_pixel[0] <= xinters){
                              inside_rectangle = !inside_rectangle;
                          }
                      }
                  }
              }
              p1x = p2x;
              p1y = p2y;
            }
            }
            vec3 candidate_point(islow+0.5,ifast+0.5,0.0);
            if (inside_rectangle) {
              if (pixel_accept_policy::accept_pixel(candidate_point,this)) {
                enclosed_px.push_back(candidate_point);
              }
            } else {
              //margin_px.push_back(candidate_point);
              //This is a margin pixel; but let's determine the distance-to-rectangle for the scoring function
              int is_positive = 0;
              double distance = 0;
              for (int iv = 1; iv < 5; ++iv){
                vec3 side = spot_rectangle_vertices[idx+iv]-spot_rectangle_vertices[idx+iv-1];
                vec3 to_point = candidate_point - spot_rectangle_vertices[idx+iv];
                double potential_distance = (to_point[0]*side[1]-side[0]*to_point[1])/side.length();
                if (potential_distance >= 0) {
                  is_positive += 1;
                  distance = potential_distance;
                }
              }
              if (is_positive == 1) {
                if (pixel_accept_policy::accept_pixel(candidate_point,this)) {
                  //we got our distance.
                  margin_px.push_back(candidate_point);
                  margin_distances.push_back(distance);
                }
              } else {
                //pixel is in the corner.  Need to find out which corner is closest
                for (int iv = 1; iv < 5; ++iv){
                  vec3 to_point = candidate_point - spot_rectangle_vertices[idx+iv];
                  double potential_distance = to_point.length();
                  if (potential_distance < margin){ // this test makes the margin a rounded rectangle rather than rectangle proper.
                    if (pixel_accept_policy::accept_pixel(candidate_point,this)) {
                      margin_px.push_back(candidate_point);
                      margin_distances.push_back(potential_distance);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    scitbx::af::shared<vec3 >
    selected_hi_predictions()const{
      scitbx::af::shared<vec3 > data;

      for (int idx = 0; idx < lo_E_limit.size(); ++idx){
        if (!observed_flag[idx]) {continue;}
        data.push_back(hi_E_limit[idx]);
      }
      return data;
    }

    scitbx::af::shared<vec3 >
    selected_lo_predictions()const{
      scitbx::af::shared<vec3 > data;

      for (int idx = 0; idx < lo_E_limit.size(); ++idx){
        if (!observed_flag[idx]) {continue;}
        data.push_back(lo_E_limit[idx]);
      }
      return data;
    }

    scitbx::af::shared<vec3 >
    selected_predictions()const{
      scitbx::af::shared<vec3 > data;

      for (int idx = 0; idx < lo_E_limit.size(); ++idx){
        if (!observed_flag[idx]) {continue;}
        vec3 hi_pos = hi_E_limit[idx];
        vec3 lo_pos = lo_E_limit[idx];
        data.push_back( (hi_pos + lo_pos)/2.);
      }
      return data;
    }

    scitbx::af::shared<vec3 >
    selected_predictions_labelit_format()const{
      scitbx::af::shared<vec3 > data;

      for (int idx = 0; idx < lo_E_limit.size(); ++idx){
        if (!observed_flag[idx]) {continue;}
        vec3 hi_pos = hi_E_limit[idx];
        vec3 lo_pos = lo_E_limit[idx];
        vec3 temp = (hi_pos + lo_pos)/2.;
        data.push_back(
        vec3((temp[1]-P.pixel_offset[1])*P.pixel_size[1],(temp[0]-P.pixel_offset[0])*P.pixel_size[0],temp[2]) );
      }
      return data;
    }

    scitbx::af::shared<cctbx::miller::index<> >
    selected_hkls()const{
      scitbx::af::shared<cctbx::miller::index<> > data;

      for (int idx = 0; idx < lo_E_limit.size(); ++idx){
        if (!observed_flag[idx]) {continue;}
        data.push_back( P.indices[idx] );
      }
      return data;
    }

    scitbx::af::shared<vec3 >
    restricted_to_active_areas(scitbx::af::shared<int> active_areas, scitbx::af::shared<vec3 > candidates)const{
      // calculation is done in picture_fast_slow coordinates (units of pixels)
      scitbx::af::shared<vec3 > active_area_traps;
      ;

      for (int idx = 0; idx < active_areas.size(); idx+=4){
      }
      return active_area_traps;
    }
    bool subpixel_translations_set;
    subpixel_joint_model sjm;
    scitbx::af::shared<double> subpixel; // DEPRECATED; delete with picture_fast_slow_force()
    scitbx::af::shared<double> rotations_rad; // DEPRECATED; same message
    void set_subpixel(scitbx::af::shared<double> s, scitbx::af::shared<double> rotations_deg){
      subpixel_translations_set=true;
      sjm = subpixel_joint_model(s,rotations_deg);
      subpixel=s; //DEPRECATED
      rotations_rad=scitbx::af::shared<double>(); //DEPRECATED
      for (int ixx=0; ixx< rotations_deg.size(); ++ixx){ //DEPRECATED
        rotations_rad.push_back(scitbx::constants::pi_180*rotations_deg[ixx]);
      }
      SCITBX_ASSERT( s.size() == 2 * rotations_rad.size()); //DEPRECATED
    }
    void set_mosaicity(double const& half_mosaicity_rad){
      P.half_mosaicity_rad=half_mosaicity_rad;}
    double get_mosaicity()const { return P.half_mosaicity_rad; }
    double p_domain_size_ang;
    void set_domain_size(double const& value){
      p_domain_size_ang=value;}
    double get_domain_size()const{ return p_domain_size_ang; }
    void set_bandpass(double const& wave_HI,double const& wave_LO){
      P.wavelengthHE = wave_HI;
      P.wavelengthLE = wave_LO;
      SCITBX_ASSERT (P.wavelengthHE <= P.wavelengthLE); SCITBX_ASSERT (P.wavelengthHE > 0.);   }
    double get_wave_HI()const{ return P.wavelengthHE; }
    double get_wave_LO()const{ return P.wavelengthLE; }
    void set_detector_origin(vec3 const& detector_origin){
      P.detector_origin = detector_origin;
    }
    void set_distance(double const& dist){
      P.distance = dist;
    }
    void set_orientation(cctbx::crystal_orientation const& orientation){
      P.orientation = orientation; }
    annlib_adaptbx::AnnAdaptorSelfInclude adapt;
    int N_bodypix;
    void set_adaptor(af::shared<ANNcoord> body_pixel_reference){
          N_bodypix = body_pixel_reference.size()/2;
          adapt = annlib_adaptbx::AnnAdaptorSelfInclude(body_pixel_reference,2,1);
    }
    double
    score_only_detail(double const& WEIGHT){
        // second set:  prediction box
        int N_enclosed = enclosed_px.size();
        int N_enclosed_body_pixels = 0;
        af::shared<ANNcoord> query;
        for (int idx = 0; idx < N_enclosed; ++idx){
          query.push_back(enclosed_px[idx][0]);
          query.push_back(enclosed_px[idx][1]);
        }
        adapt.query(query);

        for (int p = 0; p < N_enclosed; ++p){
          if (std::sqrt(adapt.distances[p]) < 0.1){
            N_enclosed_body_pixels += 1;
          }
        }

        // third set:  marginal
        int N_marginal = margin_px.size();
        int N_marginal_body_pixels = 0;
        double marginal_body = 0;
        double marginal_nonbody = 0;
        query = af::shared<ANNcoord>();
        for (int idx = 0; idx < N_marginal; ++idx){
          query.push_back(margin_px[idx][0]);
          query.push_back(margin_px[idx][1]);
        }
        adapt.query(query);
        for (int p = 0; p < N_marginal; ++p){
          if (std::sqrt(adapt.distances[p]) < 0.1){
            N_marginal_body_pixels += 1;
            marginal_body += 0.5 + 0.5 * std::cos (-scitbx::constants::pi * margin_distances[p]);//taking MARGIN==1
          } else {
            marginal_nonbody += 0.5 + 0.5 * std::cos (scitbx::constants::pi * margin_distances[p]);
          }
        }
        marginal_body *=WEIGHT;
        double Score = 0;
        Score += WEIGHT*(N_bodypix-N_enclosed_body_pixels-N_marginal_body_pixels);
        Score += marginal_body + marginal_nonbody;
        Score += N_enclosed-N_enclosed_body_pixels;
        return Score;
    }
  };


namespace ext {

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
    parameters_bp3 params( indices, orientation, incident_beam, tophat, detector_normal, detector_fast, detector_slow,
       pixel_size, pixel_offset, distance, detector_origin );
    use_case_bp3 use(params);
    use.picture_fast_slow();
    return boost::python::make_tuple(use.hi_E_limit,use.lo_E_limit,use.observed_flag);
  }

  struct bandpass_wrappers
  {

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<parameters_bp3>("parameters_bp3",
                                  init<scitbx::af::shared<cctbx::miller::index<> >,
                                  cctbx::crystal_orientation const&,
                                  vec3ref,vec3ref,vec3ref, vec3ref,vec3ref,
                                  vec3ref , vec3ref , double const& ,  vec3ref>
      ((arg("indices"), arg("orientation"), arg("incident_beam"),
        arg("packed_tophat"), arg("detector_normal"), arg("detector_fast"),
        arg("detector_slow"), arg("pixel_size"),
        arg("pixel_offset"), arg("distance"), arg("detector_origin"))))
        .add_property("distance",make_getter(&parameters_bp3::distance, rbv()))
      ;
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
      class_<use_case_bp3>("use_case_bp3",
        init<parameters_bp3 const&>(arg("parameters")))
        .def("set_active_areas", &use_case_bp3::set_active_areas)
        .def("set_sensor_model", &use_case_bp3::set_sensor_model,(
           arg("thickness_mm"), arg("mu_rho"), arg("signal_penetration")))
        .def("prescreen_indices", &use_case_bp3::prescreen_indices)
        .def("picture_fast_slow", &use_case_bp3::picture_fast_slow)
        .def("picture_fast_slow_force", &use_case_bp3::picture_fast_slow_force)//DEPRECATED
        .add_property("hi_E_limit",make_getter(&use_case_bp3::hi_E_limit, rbv()))
        .add_property("lo_E_limit",make_getter(&use_case_bp3::lo_E_limit, rbv()))
        .add_property("part_distance",make_getter(&use_case_bp3::part_distance, rbv()))
        .add_property("observed_flag",make_getter(&use_case_bp3::observed_flag, rbv()))
        .def("spot_rectangles", &use_case_bp3::spot_rectangles)
        .def("spot_rectregions", &use_case_bp3::spot_rectregions)
        .def("enclosed_pixels_and_margin_pixels",
              (void(use_case_bp3::*)()) &use_case_bp3::enclosed_pixels_and_margin_pixels)
        .def("enclosed_pixels_and_margin_pixels",
              (void(use_case_bp3::*)(double const&)) &use_case_bp3::enclosed_pixels_and_margin_pixels)
        .add_property("margin_px",make_getter(&use_case_bp3::margin_px, rbv()))
        .add_property("enclosed_px",make_getter(&use_case_bp3::enclosed_px, rbv()))
        .add_property("margin_distances",make_getter(&use_case_bp3::margin_distances, rbv()))
        .def("enclosed_pixels", &use_case_bp3::enclosed_pixels)
        .def("selected_partialities", &use_case_bp3::selected_partialities)
        .def("selected_hi_predictions",&use_case_bp3::selected_hi_predictions)
        .def("selected_lo_predictions",&use_case_bp3::selected_lo_predictions)
        .def("selected_predictions", &use_case_bp3::selected_predictions)
        .def("selected_predictions_labelit_format",
              &use_case_bp3::selected_predictions_labelit_format)
        .def("selected_hkls", &use_case_bp3::selected_hkls)
        .def("restricted_to_active_areas", &use_case_bp3::restricted_to_active_areas)
        .def("set_subpixel", &use_case_bp3::set_subpixel,(
           arg("translations"), arg("rotations_deg")))
        .def("set_mosaicity", &use_case_bp3::set_mosaicity)
        .def("get_mosaicity", &use_case_bp3::get_mosaicity)
        .def("set_domain_size", &use_case_bp3::set_domain_size)
        .def("get_domain_size", &use_case_bp3::get_domain_size)
        .def("set_bandpass", &use_case_bp3::set_bandpass,(
           arg("wave_HI"), arg("wave_LO")))
        .def("get_wave_HI", &use_case_bp3::get_wave_HI)
        .def("get_wave_LO", &use_case_bp3::get_wave_LO)
        .def("set_orientation", &use_case_bp3::set_orientation)
        .def("set_adaptor", &use_case_bp3::set_adaptor)
        .def("set_detector_origin", &use_case_bp3::set_detector_origin)
        .def("set_distance", &use_case_bp3::set_distance)
        .def("score_only_detail", &use_case_bp3::score_only_detail,(arg("weight")))
      ;

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
