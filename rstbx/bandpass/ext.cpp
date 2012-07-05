#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>

#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <cctbx/miller.h>
#include <cctbx/crystal_orientation.h>
#include <algorithm>
#include <annlib_adaptbx/ann_adaptor.h>

namespace rstbx { namespace bandpass {
  typedef scitbx::vec2<double> vec2;
  typedef scitbx::vec3<double> vec3;
  typedef vec3 const& vec3ref;

  struct active_area_filter{
    int NEAR;
    scitbx::af::shared<int> tiles;
    annlib_adaptbx::AnnAdaptor adapt;
    int tile_id;
    active_area_filter(){}
    active_area_filter(scitbx::af::shared<int> IT): NEAR(2),tiles(IT){

      af::shared<ANNcoord> reference;
      for (int i=0; i < tiles.size()/4; ++i){
        vec2 UL (static_cast<double>(tiles[4*i]),static_cast<double>(tiles[4*i+1]));
        vec2 LR (static_cast<double>(tiles[4*i+2]),static_cast<double>(tiles[4*i+3]));
        vec2 center = (UL+LR)/2.;
        reference.push_back(center[0]);
        reference.push_back(center[1]);
      }
      adapt = annlib_adaptbx::AnnAdaptor(reference,2,NEAR);
    }
    bool operator()(vec3 prediction){
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
        if (tiles[4*itile] <= prediction[0] && prediction[0] <= tiles[4*itile+2] &&
            tiles[4*itile+1] <= prediction[1] &&prediction[1] <= tiles[4*itile+3]){
          is_in_active_area = true;tile_id=itile;break;
          }
      }
      return is_in_active_area;
    }
  };

  struct parameters_bp3 {
    scitbx::af::shared<cctbx::miller::index<> > indices; // full sphere Miller indices already computed to resolution limit
    cctbx::crystal_orientation orientation;
    vec3 incident_beam;                               // direction of incident radiation; not necessarily normalized
    double half_mosaicity_rad;
    double wavelengthHE;
    double wavelengthLE;
    vec3 detector_normal;                             // normal to the detector surface (toward incident beam), length 1
    vec3 detector_fast;                               // detector fast direction, length 1
    vec3 detector_slow;                               // detector slow direction, length 1
    vec3 pixel_size;                                  // in (fast,slow,dummy) directions, mm
    vec3 pixel_offset;                                // for rendering in (fast,slow,dummy) directions, pixels
    double const distance;                              // detector distance, mm
    vec3 detector_origin;                              // coordinates of the origin pixel (fast,slow,dummy), mm

    parameters_bp3 (scitbx::af::shared<cctbx::miller::index<> > indices, cctbx::crystal_orientation const& orientation,
      vec3ref incident_beam,vec3ref packed_tophat,vec3ref detector_normal, vec3ref detector_fast,vec3ref detector_slow,
      vec3ref pixel_size, vec3ref pixel_offset, double const& distance,  vec3ref detector_origin):
         indices(indices), orientation(orientation), incident_beam(incident_beam),
         half_mosaicity_rad(packed_tophat[2]), wavelengthHE(packed_tophat[0]), wavelengthLE (packed_tophat[1]),
         detector_normal(detector_normal), detector_fast(detector_fast),detector_slow(detector_slow),
         pixel_size(pixel_size), pixel_offset(pixel_offset), distance(distance),detector_origin(detector_origin) {
          SCITBX_ASSERT (wavelengthHE <= wavelengthLE);
          SCITBX_ASSERT (wavelengthHE>0.);
         }
  };

  struct use_case_bp3 {
    parameters_bp3 P;
    scitbx::af::shared<vec3 > hi_E_limit;
    scitbx::af::shared<vec3 > lo_E_limit;
    scitbx::af::shared<bool > observed_flag;
    scitbx::af::shared<vec3 > spot_rectangle_vertices;
    scitbx::af::shared<vec3 > spot_rectregion_vertices;
    use_case_bp3 (parameters_bp3 const& P):P(P),subpixel_translations_set(false){}
    active_area_filter aaf;
    void set_active_areas(scitbx::af::shared<int> IT){
      aaf = active_area_filter(IT);
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
      scitbx::vec3<double> s0_unit = s0.normalize();

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
      scitbx::vec3<double> s1_unit = s1.normalize();
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
          scitbx::vec3<double> r = (q_unit * P.distance / q_dot_n) - P.detector_origin;

          double x = r * P.detector_fast;
          double y = r * P.detector_slow;

          limit_types[idx] += 2; // indicate that the low-energy boundary is found
          observed_flag[idx] = true;
          lo_E_limit[idx] = scitbx::vec3<double> ( (x/P.pixel_size[0])+P.pixel_offset[0],(y/P.pixel_size[1])+P.pixel_offset[1],0. );
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
              scitbx::vec3<double> r = (q_unit * P.distance / q_dot_n) - P.detector_origin;

              double x = r * P.detector_fast;
              double y = r * P.detector_slow;
              limit_types[idx] += 2; // indicate that the hi-energy boundary is found
              observed_flag[idx] = true;
              lo_E_limit[idx] = scitbx::vec3<double> ( (x/P.pixel_size[0])+P.pixel_offset[0],(y/P.pixel_size[1])+P.pixel_offset[1],0. );
            }
          }
          if (observed_flag[idx]) {
            //Do some further tests to determine if the spot is within the flagged active area with peripheral margin
            vec3 central_position = ((lo_E_limit[idx]+hi_E_limit[idx])/2.);//already in pixel units
            if (!aaf(vec3(central_position[1],central_position[0],0.))){
              observed_flag[idx]=false;
            } else if (subpixel_translations_set) {
              vec3 subpixel_trans(subpixel[2*aaf.tile_id],subpixel[1+2*aaf.tile_id],0.0);
              lo_E_limit[idx] += subpixel_trans;
              hi_E_limit[idx] += subpixel_trans;
            }
          }
      }
    }
    scitbx::af::shared<vec3 >
    spot_rectangles(vec3ref beam_coor){
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
        vec3 tangential_unit_vec(-radial_unit_vec[1],radial_unit_vec[0],0.); // 90-degree rotation
        vec3 tangential_excursion = tangential_unit_vec * radius * P.half_mosaicity_rad;
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
      typedef return_internal_reference<> rir;
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
        .def("prescreen_indices", &use_case_bp3::prescreen_indices)
        .def("picture_fast_slow", &use_case_bp3::picture_fast_slow)
        .add_property("hi_E_limit",make_getter(&use_case_bp3::hi_E_limit, rbv()))
        .add_property("lo_E_limit",make_getter(&use_case_bp3::lo_E_limit, rbv()))
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
        .def("selected_predictions", &use_case_bp3::selected_predictions)
        .def("selected_predictions_labelit_format",
        &use_case_bp3::selected_predictions_labelit_format)
        .def("selected_hkls", &use_case_bp3::selected_hkls)
        .def("restricted_to_active_areas", &use_case_bp3::restricted_to_active_areas)
        .def("set_subpixel", &use_case_bp3::set_subpixel)
        .def("set_mosaicity", &use_case_bp3::set_mosaicity)
        .def("set_domain_size", &use_case_bp3::set_domain_size)
        .def("set_bandpass", &use_case_bp3::set_bandpass)
        .def("set_orientation", &use_case_bp3::set_orientation)
        .def("set_adaptor", &use_case_bp3::set_adaptor)
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
