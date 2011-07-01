#ifndef RSTBX_SIM_XFEL1_H
#define RSTBX_SIM_XFEL1_H

#include <scitbx/array_family/flex_types.h>
#include <rstbx/dps_core/direction.h>
#include <rstbx/diffraction/ewald_sphere.h>
#include <rstbx/diffraction/partial_spot_position_partial_H.h>
#include <cctbx/crystal_orientation.h>
#include <scitbx/random.h>
#include <scitbx/constants.h>

namespace rstbx {
//! XFEL Simulation #1
/*!
   Simulated X-ray Free Electron Laser Image -- Model 1

   Simulate a serial microcrystallography experiment.
     -- Assume that the crystalline sample has sufficient unit cells such that
        the shape transform is not observed.  Each Bragg spot is recorded
        as a single point.

   Parameters that are assumed to be constant:
     Unit cell
     Structure factors
     Size of the detector in pixels and pixel size
     Sample-to-detector distance

   Each instantiation of the xfel1 class produces a single image output,
     corresponding to a single shot from the XFEL.

   Work Flow:
     instantiate the xfel1 class

 */
struct xfel1 {
  //! Default constructor
  xfel1(){}

  inline void
  set_indices(af::shared<cctbx::miller::index<> > idx){
    indices_all = idx;
  }

  inline void
  set_intensities(af::shared<double > data){
    intensities_all = data;
  }

  af::shared<std::size_t>
  select_proximal_indices(int const& half_edge,
                          double const& detector_distance_m,
                          double const& pixel_size_m,
                          cctbx::crystal_orientation const& orientation,
                          double const& mosaicity_full_width,
                          double const& bandpass_full_width,
                          double const& wavelength,
                          double const& limiting_resolution){
    mosaicity_full_width_radians = mosaicity_full_width;
    bandpass_full_width_frac = bandpass_full_width;
    wavelength_m = wavelength;
    limiting_resolution_Ang = limiting_resolution;
    shot_orientation = orientation;

    af::shared<std::size_t> selected;
    spots.resize(0);

    /* Compute the maximum separation (proximity) between mean Ewald sphere
    and bragg spot of interest, given the full bandpass and full mosaicity.
    Use a fudge factor somewhere between 1 and 2 just to be sure of getting
    them all.*/

    /* Two effects:  bandpass & mosaicity; consider bandpass first.
    Consider two Ewald spheres that intersect at the reciprocal space origin.
    Diameters are r1 = 2/lambda and r2 = 2/(lambda[1 + B/2]), where B is the
    fractional bandpass_full_width.  So the max separation is 2B/[lambda(2+B)];
    this maximum is only realized for backscattered reflections.  So for the
    actual limiting resolution L, the effect is reduced by a factor of roughly
    lambda/2L.  Again, this is a rough estimate.  So the maximum proximity
    in reciprocal space is about 2B/[(2+B)L].

    Now consider mosaicity, assume isotropic mosaicity with angular full
    width M.  This is just an angular rotation; proximity range is
    approx. (sin M/2L), or M/2L since the angles are small.

    Just add up these two effects for worst case separation & multiply
    by fudge factor 1.5 */
    double proximity_cutoff = 1.5 *(
      (mosaicity_full_width/(2.* limiting_resolution_Ang))+
      (2.*bandpass_full_width/
      ((2.+bandpass_full_width)*limiting_resolution_Ang))
    ); //inverse Angstroms

    scitbx::vec3<double> XTD(0.,0.,detector_distance_m);//crystal to detector vector
    scitbx::vec3<double> DX = XTD + scitbx::vec3<double>(1.,0.,0.); //detector x
    scitbx::vec3<double> DY = XTD + scitbx::vec3<double>(0.,1.,0.); //detector y
    double numerator = scitbx::mat3<double>(DX[0],DY[0],0.,
                                            DX[1],DY[1],0.,
                                            DX[2],DY[2],0.).determinant() -
                       scitbx::mat3<double>(XTD[0],DY[0],0.,
                                            XTD[1],DY[1],0.,
                                            XTD[2],DY[2],0.).determinant() +
                       scitbx::mat3<double>(XTD[0],DX[0],0.,
                                            XTD[1],DX[1],0.,
                                            XTD[2],DX[2],0.).determinant() -
                       scitbx::mat3<double>(XTD[0],DX[0],DY[0],
                                            XTD[1],DX[1],DY[1],
                                            XTD[2],DX[2],DY[2]).determinant();
    scitbx::mat3<double> orientation_matrix = orientation.reciprocal_matrix();
    scitbx::vec3<double> beam_vector_B(0.,0.,1./(wavelength_m*1.E10));
    scitbx::vec2<double> full_pass( beam_vector_B[2]+proximity_cutoff,
                                    beam_vector_B[2]-proximity_cutoff);

    for (int x = 0; x < indices_all.size(); ++x){
      cctbx::miller::index<> hkl = indices_all[x];
      scitbx::vec3<double> hkld (hkl[0],hkl[1],hkl[2]);
      scitbx::vec3<double> H = orientation_matrix*hkld;
      if (H.length()==0.0) { continue; }
      if (1./H.length() < limiting_resolution_Ang) { continue; }//resol. cutoff

      double t1 = 0.5 * (H*H) / (-beam_vector_B*H);
      if (t1 <= 0) { continue; }

      //actual vector to center Ewald Sphere
      scitbx::vec3<double>C = t1 * -beam_vector_B;

      double Clen = C.length();
      if (Clen < full_pass[0] && Clen > full_pass[1]) {
        selected.push_back( x );
        scitbx::vec3<double> H1 =  H - C;

        double denominator = scitbx::mat3<double>(DX[0],DY[0],H1[0],
                                                  DX[1],DY[1],H1[1],
                                                  DX[2],DY[2],H1[2]).determinant() -
                             scitbx::mat3<double>(XTD[0],DY[0],H1[0],
                                                  XTD[1],DY[1],H1[1],
                                                  XTD[2],DY[2],H1[2]).determinant() +
                             scitbx::mat3<double>(XTD[0],DX[0],H1[0],
                                                  XTD[1],DX[1],H1[1],
                                                  XTD[2],DX[2],H1[2]).determinant();

        double t = numerator/denominator;
        scitbx::vec3<double> meter_coordinates = -t*H1 ;
        //conversion of meters to pixels.
        double x = half_edge+(meter_coordinates[0]/pixel_size_m);
        double y = half_edge+(meter_coordinates[1]/pixel_size_m);
        spots.push_back(scitbx::vec3<double>(x,y,0.));

      }
    }

    return selected;
  }

  inline
  af::versa<int, af::flex_grid<> >
  raw_diffraction(af::shared<std::size_t> selected,
                  af::versa<double, af::flex_grid<> > pix,
                  size_t const& mosaic_impacts,
                  double const& detector_distance_m,
                  double const& pixel_size_m,
                  double const& darwin_factor){
    pixels.resize(af::flex_grid<>(pix.accessor().focus()), 0);
    selection_raw_counts.resize(selected.size(), 0.);
    selection_partiality.resize(selected.size(), 0.);

    //Detector must have an even number of pixels on edge;
    // assumption that direct beam is at a pixel corner.
    int half_edge = pix.accessor().focus()[0]/2;
    SCITBX_ASSERT( pix.accessor().focus()[0]%2==0 );
    SCITBX_ASSERT( pix.accessor().focus()[0]==pix.accessor().focus()[1] );
    int full_edge = half_edge*2;

    // a lot of work to generate the mosaic ensemble of orientations:
    mosaic_perturbations.resize(0);
    scitbx::boost_random::mt19937 generator = scitbx::boost_random::mt19937();
    af::shared<double> random_u = scitbx::random::mersenne_twister(
      generator).random_double(4 * mosaic_impacts);
    for (int mcount=0; mcount<mosaic_impacts; ++mcount){
          double u1=random_u[4 * mcount];
          double u2=random_u[4 * mcount + 1];
          double u3=random_u[4 * mcount + 2];
          double u4=random_u[4 * mcount + 3];
          static const double two_pi=scitbx::constants::two_pi;

          //random unit quaternion, Steven M. LaValle, Planning Algorithms
          //Cambridge University Press, 2006, p. 198.
          //http://planning.cs.uiuc.edu
          double q[4];
          q[0] = std::sqrt(1.-u1)*std::sin(two_pi*u2);
          q[1] = std::sqrt(1.-u1)*std::cos(two_pi*u2);
          q[2] = std::sqrt(u1)*std::sin(two_pi*u3);
          q[3] = std::sqrt(u1)*std::cos(two_pi*u3);

          // discard the angle from LaValle approach; set angle to a random
          // value within the range determined by the mosaicity.
          double LaValle_angle =2.* std::acos(q[0]);

          double angle = (u4-0.5)  * mosaicity_full_width_radians;

          //find vector of rotation
          double denom = std::sin(LaValle_angle/2.);
          scitbx::vec3<double> vector(q[1]/denom, q[2]/denom, q[3]/denom);

          mosaic_perturbations.push_back(
            shot_orientation.rotate_thru(vector,angle));

    }
    // finished creating the mosaic orientation ensemble

    // Now loop over all selected reflections and all mosaic orientations
    scitbx::vec3<double> XTD(0.,0.,detector_distance_m);//crystal to detector vector
    scitbx::vec3<double> DX = XTD + scitbx::vec3<double>(1.,0.,0.); //detector x
    scitbx::vec3<double> DY = XTD + scitbx::vec3<double>(0.,1.,0.); //detector y
    double numerator = scitbx::mat3<double>(DX[0],DY[0],0.,
                                            DX[1],DY[1],0.,
                                            DX[2],DY[2],0.).determinant() -
                       scitbx::mat3<double>(XTD[0],DY[0],0.,
                                            XTD[1],DY[1],0.,
                                            XTD[2],DY[2],0.).determinant() +
                       scitbx::mat3<double>(XTD[0],DX[0],0.,
                                            XTD[1],DX[1],0.,
                                            XTD[2],DX[2],0.).determinant() -
                       scitbx::mat3<double>(XTD[0],DX[0],DY[0],
                                            XTD[1],DX[1],DY[1],
                                            XTD[2],DX[2],DY[2]).determinant();
    scitbx::vec3<double> beam_vector_B(0.,0.,1./(wavelength_m*1.E10));
    scitbx::vec2<double> full_pass( beam_vector_B[2]*(1+(bandpass_full_width_frac/2.)),
                                    beam_vector_B[2]*(1-(bandpass_full_width_frac/2.)));

    for (int isel = 0; isel < selected.size(); ++isel){
      for (int imos = 0; imos < mosaic_perturbations.size(); ++imos){
        cctbx::miller::index<> hkl = indices_all[selected[isel]];
        scitbx::vec3<double> hkld (hkl[0],hkl[1],hkl[2]);
        scitbx::vec3<double> H = mosaic_perturbations[imos].reciprocal_matrix()*hkld;
        //if (H.length()==0.0) { continue; }
        //if (1./H.length() < resolution) { continue; }//resolution cutoff

        double t1 = 0.5 * (H*H) / (-beam_vector_B*H);
        if (t1 <= 0) { continue; }

        //actual vector to center Ewald Sphere
        scitbx::vec3<double>C = t1 * -beam_vector_B;

        double Clen = C.length();
        if (Clen < full_pass[0] && Clen > full_pass[1]) {

          scitbx::vec3<double> H1 =  H - C;

          double denominator = scitbx::mat3<double>(DX[0],DY[0],H1[0],
                                                    DX[1],DY[1],H1[1],
                                                    DX[2],DY[2],H1[2]).determinant() -
                               scitbx::mat3<double>(XTD[0],DY[0],H1[0],
                                                    XTD[1],DY[1],H1[1],
                                                    XTD[2],DY[2],H1[2]).determinant() +
                               scitbx::mat3<double>(XTD[0],DX[0],H1[0],
                                                    XTD[1],DX[1],H1[1],
                                                    XTD[2],DX[2],H1[2]).determinant();

          double t = numerator/denominator;
          scitbx::vec3<double> meter_coordinates = -t*H1 ;
          //conversion of meters to pixels.
          int x = half_edge+int(std::floor(meter_coordinates[0]/pixel_size_m));
          int y = half_edge+int(std::floor(meter_coordinates[1]/pixel_size_m));

          if (x>=0 && x < full_edge && y>=0 && y < full_edge) {

            //two_theta = uc.two_theta(h_k_l,self.sim.xray_wavelength)
            // **** Not taking into account the polarization of incident radiation ****
            // polarization_factor = (1. + math.cos(two_theta)**2) / 2.

            double raw_counts = intensities_all[selected[isel]] *
                                darwin_factor/mosaic_impacts;
            pixels(x,y)+=int( raw_counts );
            selection_partiality[isel] += 1./mosaic_impacts;
            selection_raw_counts[isel] += raw_counts;
          }
        }
      }
    }

    return pixels;
  }

  af::shared<cctbx::miller::index<> > indices_all;
  af::shared<double > intensities_all;
  af::shared<cctbx::crystal_orientation > mosaic_perturbations;
  af::versa<int, af::flex_grid<> > pixels;
  af::shared<scitbx::vec3<double> > spots;
  af::shared<double> selection_raw_counts;
  af::shared<double> selection_partiality;
  double mosaicity_full_width_radians;
  double bandpass_full_width_frac;
  double wavelength_m;
  double limiting_resolution_Ang;
  cctbx::crystal_orientation shot_orientation;

};
}
#endif //RSTBX_SIM_XFEL1_H
