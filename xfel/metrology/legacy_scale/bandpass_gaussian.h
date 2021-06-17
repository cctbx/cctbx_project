#ifndef LEGACY_SCALE_BANDPASS_GAUSSIAN_H
#define LEGACY_SCALE_BANDPASS_GAUSSIAN_H
#include <cmath>
#include <rstbx/bandpass/parameters.h>
#include <xfel/metrology/legacy_scale/vector_collection.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/lbfgs.h>

namespace xfel_legacy {

namespace parameter {

  typedef scitbx::vec2<double> vec2;
  typedef scitbx::vec3<double> vec3;
  typedef vec3 const& vec3ref;

  double best_fit_limit( farray const&, farray const&);

  struct streak_parameters{
    double wavelength;
    vec3 rotax;
    double rotax_excursion_rad;
    double domain_size_inv_ang;
    vec3 position;
    vec3 part_position_part_wavelength;
  };


  struct bandpass_gaussian {
    rstbx::bandpass::parameters_bp3 P;
    rstbx::bandpass::pad_sensor_model sensor;
    double signal_penetration;
    scitbx::af::shared<vec3 > hi_E_limit;
    scitbx::af::shared<vec3 > lo_E_limit;
    scitbx::af::shared<vec3 > mean_position;
    scitbx::af::shared<double > calc_radial_length;
    scitbx::af::shared<bool > observed_flag;
    scitbx::af::shared<vec3 > spot_rectangle_vertices;
    scitbx::af::shared<vec3 > spot_rectregion_vertices;
    scitbx::af::shared<vec3 > part_distance;
    bandpass_gaussian (rstbx::bandpass::parameters_bp3 const& P):P(P),subpixel_translations_set(false){}
    rstbx::bandpass::active_area_filter aaf;
    void set_active_areas(scitbx::af::shared<int> IT){
      aaf = rstbx::bandpass::active_area_filter(IT);
    }
    void set_sensor_model(double const& thickness_mm, double const& mu_rho, double const& sp){
      sensor.thickness_mm = thickness_mm;
      sensor.mu_rho = mu_rho;
      signal_penetration = sp;
    }

    void
    picture_fast_slow_force(){
      //The effect of this function is to set the following three state vectors:
      hi_E_limit = scitbx::af::shared<scitbx::vec3<double> >(P.indices.size(),scitbx::af::init_functor_null<scitbx::vec3<double> >());
      lo_E_limit = scitbx::af::shared<scitbx::vec3<double> >(P.indices.size(),scitbx::af::init_functor_null<scitbx::vec3<double> >());
      mean_position = scitbx::af::shared<scitbx::vec3<double> >(P.indices.size(),scitbx::af::init_functor_null<scitbx::vec3<double> >());
      calc_radial_length = scitbx::af::shared<double >(P.indices.size(),scitbx::af::init_functor_null<double >());
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
            calc_radial_length[idx] = (lo_E_limit[idx] - hi_E_limit[idx]).length();
            if (subpixel_translations_set) {
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

    void
    gaussian_fast_slow(){
      //The effect of this function is to set the following three state vectors:
      mean_position = scitbx::af::shared<scitbx::vec3<double> >(P.indices.size(),scitbx::af::init_functor_null<scitbx::vec3<double> >());
      part_distance = scitbx::af::shared<scitbx::vec3<double> >(P.indices.size(),scitbx::af::init_functor_null<scitbx::vec3<double> >());
      observed_flag = scitbx::af::shared<bool >(P.indices.size());
      calc_radial_length = scitbx::af::shared<double >(P.indices.size(),scitbx::af::init_functor_null<double >());

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
      scitbx::vec3<double> mean_position_part_r_part_distance = s1;
      scitbx::vec3<double> mean_position_part_r_part_half_mos = s1;

      vec3* part_half_mos_ptr = gradient_mapping["half_mosaicity_rad"];

      //  Cn, the circular section through the Ewald sphere.
      for (int idx = 0; idx < P.indices.size(); ++idx){
          double spot_resolution_ang = P.orientation.unit_cell().d(P.indices[idx]);
          double effective_half_mosaicity_rad = (p_domain_size_ang > 0.)? P.half_mosaicity_rad + spot_resolution_ang / (2. * p_domain_size_ang) : P.half_mosaicity_rad;

          scitbx::vec3<double> H(P.indices[idx][0],P.indices[idx][1], P.indices[idx][2]); // the Miller index
          scitbx::vec3<double> s = A * H; //s, the reciprocal space coordinates, lab frame, of the oriented Miller index
          double s_rad_sq = s.length_sq();
          scitbx::vec3<double> s_unit = s.normalize();
          SCITBX_ASSERT(s_rad_sq > 0.);
          scitbx::vec3<double> rotax = s.normalize().cross(s0_unit); //The axis that most directly brings the Bragg spot onto Ewald sphere
          scitbx::vec3<double> sxrx = s_unit.cross(rotax);
          scitbx::vec3<double> chord_direction =      (rotax.cross(s0)).normalize();

         //  ###########  Look at the high-energy wavelength boundary

          double a = s.length_sq()/(2.*s0_length); // see diagram
          double b = std::sqrt(s.length_sq() - (a*a));   //  Calculate half-length of the chord of intersection

          scitbx::vec3<double> intersection_hi = -a * s0_unit- b*chord_direction;

          double hi_energy_theta_rad = std::atan( (intersection_hi * sxrx)/(intersection_hi * s_unit));

         //  ###########  Look at the low-energy wavelength boundary

          double alow = s.length_sq()/(2.*s1_length); // see diagram
          double blow = std::sqrt(s.length_sq() - (alow*alow));   //  Calculate half-length of the chord of intersection

          scitbx::vec3<double> intersection_lo = -alow * s0_unit- blow*chord_direction;

          double lo_energy_theta_rad = std::atan( (intersection_lo * sxrx)/(intersection_lo * s_unit));

          //printf("Comparative rotations (deg) hi_e %7.4f lo_e %7.4f half mos %7.4f\n",
          //  hi_energy_theta_deg, lo_energy_theta_deg, effective_half_mosaicity_rad/scitbx::constants::pi_180);

          double sigma_energy_rad = (hi_energy_theta_rad - lo_energy_theta_rad)/2.;
          double variance_energy_rad =  sigma_energy_rad * sigma_energy_rad;
          double variance_mosaicity_rad = effective_half_mosaicity_rad*effective_half_mosaicity_rad;
          double mean_energy_rad = (hi_energy_theta_rad + lo_energy_theta_rad)/2.;
          double mean_product = (mean_energy_rad*variance_mosaicity_rad)/(variance_energy_rad+variance_mosaicity_rad);
          double sigma_product = std::sqrt((variance_energy_rad*variance_mosaicity_rad)/(variance_energy_rad+variance_mosaicity_rad));
          //-----------------
          scitbx::vec3<double> s_rot_mean_product = s.rotate_around_origin(rotax,-mean_product);
          double a_mean_product = -s_rot_mean_product * s0_unit;
          SCITBX_ASSERT(a_mean_product != 0.);
          double r_n_mean_product = s_rad_sq/(2.*a_mean_product);
          //double wavelength_mean_product = 1./r_n_mean_product;
          scitbx::vec3<double> s0_mean_product = r_n_mean_product * P.incident_beam;
          scitbx::vec3<double> q = (s_rot_mean_product + s0_mean_product);
          //-----------------
          scitbx::vec3<double> q_unit = q.normalize();

          double q_dot_n = q_unit * P.detector_normal;
          SCITBX_ASSERT(q_dot_n != 0.);
          scitbx::vec3<double> r = (q_unit * P.distance / q_dot_n) - P.detector_origin;

          double x = r * P.detector_fast;
          double y = r * P.detector_slow;
          observed_flag[idx] = true;
          mean_position[idx] = scitbx::vec3<double> ( (x/P.pixel_size[0])+P.pixel_offset[0],(y/P.pixel_size[1])+P.pixel_offset[1],0. );
          scitbx::vec3<double> part_r_part_d( q_unit / q_dot_n );
          mean_position_part_r_part_distance = scitbx::vec3<double>( (part_r_part_d * P.detector_fast) / P.pixel_size[0],
                                                                     (part_r_part_d * P.detector_slow) / P.pixel_size[1],0.);

          {//Special scope to calculate the half-radial elongation
            scitbx::vec3<double> s_rot_mean_product = s.rotate_around_origin(rotax,-mean_product+sigma_product);
            double a_mean_product = -s_rot_mean_product * s0_unit;
            double r_n_mean_product = s_rad_sq/(2.*a_mean_product);
            //double wavelength_mean_product = 1./r_n_mean_product;
            scitbx::vec3<double> s0_mean_product = r_n_mean_product * P.incident_beam;
            scitbx::vec3<double> q = (s_rot_mean_product + s0_mean_product);
            //-----------------
            scitbx::vec3<double> q_unit = q.normalize();

            double q_dot_n = q_unit * P.detector_normal;
            scitbx::vec3<double> r = (q_unit * P.distance / q_dot_n) - P.detector_origin;

            double x = r * P.detector_fast;
            double y = r * P.detector_slow;
            scitbx::vec3<double> sigma_position( (x/P.pixel_size[0])+P.pixel_offset[0],(y/P.pixel_size[1])+P.pixel_offset[1],0. );
            calc_radial_length[idx] = 2. * (mean_position[idx] - sigma_position).length();
          }

          if (refine_vectors["half_mosaicity_rad"]) {
            //derivatives of mosaicity.  qvec result from Mathematica
            double temp_denom = (variance_energy_rad + variance_mosaicity_rad)*(variance_energy_rad + variance_mosaicity_rad);
            double part_mu_prod_part_meff = (2.* mean_energy_rad * effective_half_mosaicity_rad)*
                                            (variance_energy_rad + variance_mosaicity_rad) -
                                            (mean_energy_rad * variance_mosaicity_rad) *
                                            (2 * effective_half_mosaicity_rad);
            part_mu_prod_part_meff /= temp_denom;

            // !**********************************************
            // differ from Mathematica worksheet by SM *= -1, and part_q_part_half_mos *= -1
            // in other words, handedness of rotation is opposite in cpp code vs Mathematica
            double CM = std::cos(mean_product); double SM = -std::sin(mean_product);
            double rds = s * rotax;
            scitbx::vec3<double> rs = rotax.cross(s);
            double pqpm_denom = 2. * std::pow((-rds * rotax + (rotax * rds - s) * CM - rs * SM) * s0_unit, 2);
            double pqpm_num_fac = (rs * CM + (rotax * rds - s) * SM) * s0_unit;

            scitbx::vec3<double> mathematica_part_qvec_part_mu_prod (
              rs[0] * CM + rotax[0] * rds * SM - s[0] * SM + (P.incident_beam[0]* pqpm_num_fac * s_rad_sq)/
                pqpm_denom,

              rs[1] * CM + rotax[1] * rds * SM - s[1] * SM + (P.incident_beam[1]* pqpm_num_fac * s_rad_sq)/
                pqpm_denom,

              rs[2] * CM + rotax[2] * rds * SM - s[2] * SM + (P.incident_beam[2]* pqpm_num_fac * s_rad_sq)/
                pqpm_denom
              );

            scitbx::vec3<double> part_q_part_half_mos = - mathematica_part_qvec_part_mu_prod * part_mu_prod_part_meff;

            // !**********************************************

            scitbx::mat3<double>jacobian_rvec_wrt_qvec( -P.distance/q[2], 0,              q[0]*P.distance/(q[2]*q[2]),
                                                                       0,-P.distance/q[2],q[1]*P.distance/(q[2]*q[2]),
                                                                       0,               0,                         0);
            scitbx::vec3<double>part_r_part_half_mos = jacobian_rvec_wrt_qvec * part_q_part_half_mos;

            mean_position_part_r_part_half_mos = scitbx::vec3<double>( (part_r_part_half_mos * P.detector_fast) / P.pixel_size[0],
                                                                       (part_r_part_half_mos * P.detector_slow) / P.pixel_size[1],0.);

            *part_half_mos_ptr++ = mean_position_part_r_part_half_mos;

            // end of mosaicity
          }

          if (observed_flag[idx]) {
            if (subpixel_translations_set) {
              vec3 subpixel_trans(subpixel[2*aaf.tile_id],subpixel[1+2*aaf.tile_id],0.0);
              mean_position[idx] += subpixel_trans;
            }
            part_distance[idx] = mean_position_part_r_part_distance;
          }
      }
      for ( int idx = 0; idx < P.indices.size(); ++idx){
        SCITBX_ASSERT (observed_flag[idx]);
      }
    }

    streak_parameters simple_forward_calculation_spot_position (
      double const& wavelength, int const& observation_no) const {
      streak_parameters result;

      scitbx::mat3<double> A = P.orientation.reciprocal_matrix();

      //s0:  parallel to the direction of incident radiation
      scitbx::vec3<double> s0 = (1./wavelength) * P.incident_beam;
      double s0_length = s0.length();
      scitbx::vec3<double> s0_hat = s0.normalize();

      scitbx::vec3<double> H(P.indices[observation_no][0],
                             P.indices[observation_no][1],
                             P.indices[observation_no][2]); // the Miller index
      scitbx::vec3<double> q = A * H; //q, recip space coord, lab frame, of oriented Milleridx
      double q_sq = q.length_sq();
      //The axis that most directly brings the Bragg spot onto Ewald sphere
      scitbx::vec3<double> rotax = q.normalize().cross(s0_hat);
      scitbx::vec3<double> chord_direction =      (rotax.cross(s0)).normalize();

      double a = q_sq/(2.*s0_length); // see diagram
      double b = std::sqrt(q_sq - (a*a));   //  half-length of the chord of intersection

      scitbx::vec3<double> intersection = -a * s0_hat - b*chord_direction;

      scitbx::vec3<double> S = (intersection + s0); // the reflection vector
      double S_sq = S * S;
      double S_length = std::sqrt(S_sq);
      scitbx::vec3<double> S_unit = S/S_length;

      double S_dot_n = S_unit * P.detector_normal;
      scitbx::vec3<double> r = (S_unit * P.distance / S_dot_n) - P.detector_origin;

      double x = r * P.detector_fast;
      double y = r * P.detector_slow;
      result.position = scitbx::vec3<double> ((x/P.pixel_size[0])+P.pixel_offset[0],
                                              (y/P.pixel_size[1])+P.pixel_offset[1],0. );

      //get some partial derivatives with respect to wavelength
      scitbx::vec3<double> s0_prime = -P.incident_beam/(wavelength*wavelength);
      scitbx::vec3<double> intersection_prime = (-q_sq/2.)*s0_hat + (a/b)*chord_direction;
      scitbx::vec3<double> S_prime = intersection_prime + s0_prime;
      scitbx::vec3<double> S_hat_prime = ((S_length*S_prime)-(S_unit*(S*S_prime)))/S_sq;
      double Sn_prime = S_hat_prime * P.detector_normal;
      scitbx::vec3<double> r_prime = (P.distance/(S_dot_n*S_dot_n)) *
       (S_hat_prime * S_dot_n - S_unit*Sn_prime);
      result.part_position_part_wavelength = scitbx::vec3<double>(
          r_prime * P.detector_fast/P.pixel_size[0],
          r_prime * P.detector_slow/P.pixel_size[1],
          0.0
        );

      //now some extra work to determine angular excursion about rotax to get from q to r.
      double angle_Q = std::atan2 ( q * s0_hat, q * chord_direction);
      double angle_R = std::atan2 ( intersection * s0_hat, intersection * chord_direction);
      result.rotax_excursion_rad = angle_Q - angle_R;
      /* WARNING! cases have been found with Q +179 and R -179; leading to unphysical Q-R difference.
         This was observed with simulated psI data, with zone nearly aligned on beam,
         and severe indexing misalignment.
         Caused an endless loop with memory increase
         Guard against this in the Python layer by checking for excursion values
           that are significantly different from zero.
       */
      result.rotax = rotax;
      result.wavelength = wavelength;
      return result;
    }

    double simple_part_excursion_part_rotxy (
      double const& wavelength,
      int const& observation_no,
      scitbx::mat3<double> dA_drotxy) const {

      scitbx::mat3<double> A = P.orientation.reciprocal_matrix();

      //s0:  parallel to the direction of incident radiation
      scitbx::vec3<double> s0 = (1./wavelength) * P.incident_beam;
      double s0_length = s0.length();
      scitbx::vec3<double> s0_hat = s0.normalize();

      scitbx::vec3<double> H(P.indices[observation_no][0],
                             P.indices[observation_no][1],
                             P.indices[observation_no][2]); // the Miller index
      scitbx::vec3<double> qvec = A * H; //q, recip space coord, lab frame, of oriented Milleridx
      scitbx::vec3<double> dqvec_drot = dA_drotxy * H;

      double q_sq = qvec.length_sq();
      double qlen = std::sqrt(q_sq);

      double dq_drot = (qvec * dqvec_drot)/qlen;

      scitbx::vec3<double> qhat = qvec/qlen;
      scitbx::vec3<double> dqhat_drot = (qlen*dqvec_drot - dq_drot*qvec)/q_sq;

      scitbx::vec3<double> e1 = qhat.cross(s0_hat);
      scitbx::vec3<double> de1_drot = dqhat_drot.cross(s0_hat);

      scitbx::vec3<double> chord_direction = e1.cross(s0_hat);
      scitbx::vec3<double> dc_drot = de1_drot.cross(s0_hat);

      double a = q_sq/(2.*s0_length);
      double da_drot = wavelength*(qvec * dqvec_drot);

      double b = std::sqrt(q_sq - (a*a));
      double db_drot = (s0_length - a) * da_drot / b;

      scitbx::vec3<double> intersection = -a * s0_hat - b*chord_direction;
      scitbx::vec3<double> di_drot = -(da_drot * s0_hat + b * dc_drot + db_drot * chord_direction);

      //double angle_Q = std::atan2 ( q * s0_hat, q * chord_direction);
      double x,y;
      y = qvec * s0_hat; x = qvec * chord_direction;
      double datan2_dx = -y/(x*x+y*y);
      double datan2_dy = x/(x*x+y*y);
      double dangle_Q_drot = datan2_dx * (dqvec_drot * chord_direction + qvec * dc_drot)
                            + datan2_dy * (dqvec_drot * s0_hat);

      //double angle_R = std::atan2 ( intersection * s0_hat, intersection * chord_direction);
      y = intersection * s0_hat; x = intersection * chord_direction;
      datan2_dx = -y/(x*x+y*y);
      datan2_dy = x/(x*x+y*y);
      double dangle_R_drot = datan2_dx * (di_drot * chord_direction + intersection * dc_drot)
                            + datan2_dy * (di_drot * s0_hat);
      return dangle_Q_drot - dangle_R_drot;
    }

    scitbx::vec3<double>
    detector_origin_as_pixels() const{
      scitbx::vec3<double> r = - P.detector_origin;
      double x = r * P.detector_fast;
      double y = r * P.detector_slow;
      return scitbx::vec3<double> ((x/P.pixel_size[0])+P.pixel_offset[0],
                                   (y/P.pixel_size[1])+P.pixel_offset[1],0. );
    }

    bool subpixel_translations_set;
    scitbx::af::shared<double> subpixel;
    void set_subpixel(scitbx::af::shared<double> s){
      subpixel_translations_set=true;
      subpixel=s;}
    void set_mosaicity(double const& half_mosaicity_rad){
      P.half_mosaicity_rad=half_mosaicity_rad;}
    double p_domain_size_ang;
    void set_domain_size(double const& value){
      p_domain_size_ang=value;}
    void set_bandpass(double const& wave_HI,double const& wave_LO){
      P.wavelengthHE = wave_HI;
      P.wavelengthLE = wave_LO;
      SCITBX_ASSERT (P.wavelengthHE <= P.wavelengthLE); SCITBX_ASSERT (P.wavelengthHE > 0.);   }
    void set_detector_origin(vec3 const& detector_origin){
      P.detector_origin = detector_origin;
    }
    void set_distance(double const& dist){
      P.distance = dist;
    }
    void set_orientation(cctbx::crystal_orientation const& orientation){
      P.orientation = orientation; }

    std::map<std::string, vec3*> gradient_mapping;
    std::map<std::string, vec3*> curvature_mapping;
    std::map<std::string, bool> refine_vectors;
    void set_vector_output_pointers(xfel_legacy::parameter::vector_collection& vectordata,
                                    int const& frame_id, bool const& refine){
      if (!refine) {return;}
      for (std::map<std::string,vector_array>::const_iterator tag=vectordata._V.begin(); tag!=vectordata._V.end(); ++tag){
        std::string field = tag->first;
        {vec3* first = &(vectordata._V.find(field)->second.gradients.ref()[vectordata.frame_first_index.find(frame_id)->second]);
        gradient_mapping[field]=first;
        }
        {vec3* first = &(vectordata._V.find(field)->second.curvatures.ref()[vectordata.frame_first_index.find(frame_id)->second]);
        curvature_mapping[field]=first;
        }
        refine_vectors[field]=refine;
      }
    }

    scitbx::af::flex_double wavelength_fit_ang;
    scitbx::af::flex_double mosaicity_fit_rad;
    scitbx::af::flex_double wavelength_fit_ang_sigma;
    scitbx::af::flex_double mosaicity_fit_rad_sigma;
    void
    measure_bandpass_and_mosaic_parameters(
            scitbx::af::shared<double > radial,
            scitbx::af::shared<double > azimut, double const& domain_sz_inv_ang,
            scitbx::af::shared<scitbx::vec3<double> > lab_frame_obs
            ){
      int datasz = radial.size();

      // row 0 == high energy fit, row 1 == low energy fit
      wavelength_fit_ang = af::flex_double(af::flex_grid<>(2,datasz));
      mosaicity_fit_rad = af::flex_double(af::flex_grid<>(2,datasz));
      wavelength_fit_ang_sigma = af::flex_double(af::flex_grid<>(2,datasz));
      mosaicity_fit_rad_sigma = af::flex_double(af::flex_grid<>(2,datasz));

      double* beginw = wavelength_fit_ang.begin();
      double* beginm = mosaicity_fit_rad.begin();
      double* beginws = wavelength_fit_ang_sigma.begin();
      double* beginms = mosaicity_fit_rad_sigma.begin();

      vec3 origin = detector_origin_as_pixels();
      for (int i = 0; i<datasz; ++i){
      // figure out positions of radial streak limits
        vec3 outward_radial_direction((lab_frame_obs[i] - origin).normalize());
        for (int row = 0; row <2; ++row){
          int SIGN=(row?-1:1); // -1, high energy; +1 low energy
          //add or subtract the radial half-width
          vec3 streak_limit = lab_frame_obs[i] + (SIGN * radial[i]/2.)* outward_radial_direction ;
          //printf("new point row %2d limitx %6.1f y %6.1f\n",row, streak_limit[0],streak_limit[1]);
          streak_parameters SLP = refine_streak_limit_one_obs(streak_limit,P.wavelengthHE,i);
          *(beginw+row*datasz) = SLP.wavelength;
          *(beginm+row*datasz) = SLP.rotax_excursion_rad;

          // add half a pixel to give an approximation of the error standard deviation
          streak_limit = lab_frame_obs[i] + (SIGN * ((radial[i]/2.) + 0.5))* outward_radial_direction ;
          SLP = refine_streak_limit_one_obs(streak_limit,P.wavelengthHE,i);
          *(beginws+row*datasz) = SLP.wavelength;
          *(beginms+row*datasz) = SLP.rotax_excursion_rad;
        }
        ++beginw;++beginws;
        ++beginm;++beginms;
      }
    }

    streak_parameters
    refine_streak_limit_one_obs (vec3 const& streak_limit,
                                 double const& W,
                                 int const& observation_no) const{
      bool verbose(false);
      bool verbose_pos(false);
      //First want to modify the streak limit position with rotz rotation
      //such that calc and obs positions are related by simple rotax
      streak_parameters provisional = simple_forward_calculation_spot_position (
            W, observation_no);
      vec3 origin = detector_origin_as_pixels();
      vec3 prov_calc = provisional.position - origin;
      vec3 prov_obs  = streak_limit - origin;
      double angle_calc = std::atan2 ( prov_calc * P.detector_fast, prov_calc * P.detector_slow);
      double angle_obs = std::atan2 ( prov_obs * P.detector_fast, prov_obs * P.detector_slow);
      double rotz_angle = angle_calc - angle_obs;

      vec3 rotz_obs = prov_obs.unit_rotate_around_origin( vec3(0.,0.,1.), rotz_angle);
      if (verbose) {printf("Rotz angle %7.4f after ward %12.9f\n",rotz_angle,
        ((rotz_obs*prov_calc)/(rotz_obs.length()*prov_calc.length())));}
      vec3 rotz_streak_limit = rotz_obs + origin;

      streak_parameters PAR;
      PAR.wavelength = W;
      double* xb = &(PAR.wavelength);
      double gradient;
      const int nparam = 1; // number of parameters
      double convergence_test_=0.001; //convergence threshhold in units of pixels

      scitbx::lbfgs::minimizer<double, int> minimizer(nparam);
      scitbx::lbfgs::traditional_convergence_test<double, int> is_converged(nparam);
      if (verbose) {
        std::cout << "n: " << minimizer.n() << std::endl;
        std::cout << "m: " << minimizer.m() << std::endl;
        std::cout << "xtol: " << minimizer.xtol() << std::endl;
      }

      for(;;) {
        double f = 0.;
        if (verbose) printf("input wavelength %8.6f\n",*xb);
        streak_parameters calc = simple_forward_calculation_spot_position (
          *xb, observation_no);
        PAR.rotax_excursion_rad = calc.rotax_excursion_rad;
        if (verbose_pos) {
          printf("calculated position x %12.7f y %12.7f\n", calc.position[0], calc.position[1]);
          printf("streak position obs x %12.7f y %12.7f\n", rotz_streak_limit[0], rotz_streak_limit[1]);
          printf("beam spot origin    x %12.7f y %12.7f\n", origin[0], origin[1]);
          double dotproduct = (rotz_streak_limit - origin)*(calc.position-origin)/ (
          (rotz_streak_limit - origin).length()*(calc.position-origin).length()
          );
          printf ("dot product is %15.10f\n",dotproduct);
        }
        vec3 delta = calc.position - rotz_streak_limit;
        f = delta.length_sq();
        if (verbose) printf("delta in pixels %7.4f\n",std::sqrt(f));
        gradient = 2.*delta*calc.part_position_part_wavelength;

        if (verbose) std::cout << "f: " << f << "x="<<*xb<< " gradient="<<gradient<<std::endl;
        if (std::sqrt(f) < convergence_test_) {
          if (verbose) printf("SEEMS TO HAVE CONVERGED\n");
          break;}
        if (minimizer.run(xb, f, &gradient)) continue;
        if (verbose_pos) {
          std::cout << "f: " << f << "x="<<*xb<< " gradient="<<gradient;
          std::cout << std::endl;
          std::cout << " " << minimizer.iter();
          std::cout << " " << minimizer.nfun();
          std::cout << " " << minimizer.stp();
          std::cout << std::endl;
        }
        if (is_converged(xb, &gradient)) break;
        if (minimizer.nfun() > 2000) break;
        minimizer.run(xb, f, &gradient);
      }
      return PAR;
    }

  };


}

} //namespace xfel
#endif// LEGACY_SCALE_BANDPASS_GAUSSIAN_H
