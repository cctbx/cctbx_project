#include <simtbx/nanoBragg/nanoBragg.h>

//Contributed by Nicholas Sauter,LBNL.

namespace simtbx {
namespace nanoBragg {

struct source_const{
  source_const(int const& source, const double* diffracted, const nanoBragg* n){
    /* retrieve stuff from cache */
    incident[1] = -n->source_X[source];
    incident[2] = -n->source_Y[source];
    incident[3] = -n->source_Z[source];
    double lambda = n->source_lambda[source];

    /* construct the incident beam unit vector while recovering source distance */
    source_path = unitize(incident,incident);

    /* construct the scattering vector for this pixel */
    scattering[1] = (diffracted[1]-incident[1])/lambda;
    scattering[2] = (diffracted[2]-incident[2])/lambda;
    scattering[3] = (diffracted[3]-incident[3])/lambda;

    /* sin(theta)/lambda is half the scattering vector length */
    stol = 0.5*magnitude(scattering);
  }
  double incident[4],scattering[4];
  double stol,source_path;
};

struct phitic_const{
  phitic_const(int const& mostic, double const& phi, const nanoBragg* n){
    for (int icopy=0; icopy<4; ++icopy){
        a0[icopy] = n->a0[icopy];
        b0[icopy] = n->b0[icopy];
        c0[icopy] = n->c0[icopy];
        spindle_vector[icopy] = n->spindle_vector[icopy];
    }
    if( phi != 0.0 ) {
      /* rotate about spindle if neccesary */
      rotate_axis(a0,ap,spindle_vector,phi);
      rotate_axis(b0,bp,spindle_vector,phi);
      rotate_axis(c0,cp,spindle_vector,phi);
    } else {
      for (int icopy=0; icopy<4; ++icopy){
        ap[icopy] = n->a0[icopy];
        bp[icopy] = n->b0[icopy];
        cp[icopy] = n->c0[icopy];
      }
    }
  }
  double ap[4],bp[4],cp[4],a0[4],b0[4],c0[4],spindle_vector[4];
};

struct mostic_const{
  mostic_const( int const& mostic, double* n_diffracted,
      double const& capture_fraction, int const& source,
      double* n_ap, double* n_bp, double* n_cp, double* n_scattering,
      double* n_incident, const nanoBragg* n):
      F_cell(n->default_F),polar(1.0),I_increment(0){
      for (int icopy=0; icopy<4; ++icopy){
        incident[icopy] = n_incident[icopy];
        diffracted[icopy] = n_diffracted[icopy];
        axis[icopy] = n->polar_vector[icopy];
      }

    vec3 a,b,c;
    vec3 ap(n_ap[1],n_ap[2],n_ap[3]);
    vec3 bp(n_bp[1],n_bp[2],n_bp[3]);
    vec3 cp(n_cp[1],n_cp[2],n_cp[3]);
    mat3 u_mat = mat3(
      n->mosaic_umats[mostic*9],n->mosaic_umats[mostic*9+1],n->mosaic_umats[mostic*9+2],
      n->mosaic_umats[mostic*9+3],n->mosaic_umats[mostic*9+4],n->mosaic_umats[mostic*9+5],
      n->mosaic_umats[mostic*9+6],n->mosaic_umats[mostic*9+7],n->mosaic_umats[mostic*9+8]);

    /* apply mosaic rotation after phi rotation */
    if( n->mosaic_spread > 0.0 ){ a = u_mat * ap; b = u_mat * bp; c = u_mat * cp;
    } else { a = ap; b = bp; c = cp;}
    vec3 scattering(n_scattering[1],n_scattering[2],n_scattering[3]);

    /* construct fractional Miller indicies */
    h = a * scattering;k = b * scattering;l = c * scattering;
    /* round off to nearest whole index */
    h0 = static_cast<int>(ceil(h-0.5));
    k0 = static_cast<int>(ceil(k-0.5));
    l0 = static_cast<int>(ceil(l-0.5));

    // structure factor of the lattice (paralelpiped crystal)
    F_latt = 1.0;
    hrad_sqr = 0.;
    if(n->xtal_shape == SQUARE){ /* xtal is a paralelpiped */
      if(n->Na>1){ F_latt *= sincg(M_PI*h,n->Na); }
      if(n->Nb>1){ F_latt *= sincg(M_PI*k,n->Nb); }
      if(n->Nc>1){ F_latt *= sincg(M_PI*l,n->Nc); }
    } else { /* handy radius in reciprocal space, squared */
      hrad_sqr = (h-h0)*(h-h0)*n->Na*n->Na + (k-k0)*(k-k0)*n->Nb*n->Nb + (l-l0)*(l-l0)*n->Nc*n->Nc ;
    }
    if(n->xtal_shape == ROUND){ /* use sinc3 for elliptical xtal shape,
                                correcting for sqrt of volume ratio between cube and sphere */
      F_latt = n->Na*n->Nb*n->Nc*0.723601254558268*sinc3(M_PI*sqrt( hrad_sqr * n->fudge ) );
    }
    if(n->xtal_shape == GAUSS){
               /* fudge the radius so that volume and FWHM are similar to square_xtal spots */
      F_latt = n->Na*n->Nb*n->Nc*exp(-( hrad_sqr / 0.63 * n->fudge ));
    }
    if(n->xtal_shape == TOPHAT) {
               /* make a flat-top spot of same height and volume as square_xtal spots */
      F_latt = n->Na*n->Nb*n->Nc*(hrad_sqr*n->fudge < 0.3969 );
    }
    if(F_latt == 0.0) {return;}

    /* structure factor of the unit cell */
    if ( (h0<=n->h_max) && (h0>=n->h_min) &&
         (k0<=n->k_max) && (k0>=n->k_min) &&
         (l0<=n->l_max) && (l0>=n->l_min)  ) {
      /* just take nearest-neighbor */
      F_cell = n->Fhkl[h0-n->h_min][k0-n->k_min][l0-n->l_min];
    }else{
      F_cell = n->default_F; // usually zero
    }
    /* now we have the structure factor for this pixel */
    /* polarization factor */
    if(! n->nopolar){
      /* need to compute polarization factor */
      polar = polarization_factor( n->polarization, incident,
                                   diffracted,axis);
    } else {
      polar = 1.0;
    }
    I_increment = F_cell*F_cell*F_latt*F_latt*n->source_I[source]*capture_fraction;
  };
  double h,k,l;
  int h0,k0,l0;
  double F_latt, F_cell, hrad_sqr, polar, I_increment;
  double diffracted[4], incident[4], axis[4];
  vec3 as_vec3(const double* vec4){ return vec3(vec4[1],vec4[2],vec4[3]); }
};

/* add spots from nanocrystal simulation */
void
nanoBragg::add_nanoBragg_spots_nks(boost_adaptbx::python::streambuf & output)
{
  boost_adaptbx::python::streambuf::ostream os(output);
  floatimage = raw_pixels.begin();

  if(verbose) {printf("TESTING sincg(1,1)= %f\n",sincg(1,1));}

  /* make sure we are normalizing with the right number of sub-steps */
  steps = phisteps*mosaic_domains*oversample*oversample;
  subpixel_size = pixel_size/oversample;

  double sum = 0.0; //reduction variable
  double sumsqr = 0.0; //reduction variable

  int maxidx = spixels * fpixels;
  # pragma omp parallel for reduction (+:sum, sumsqr)
  for(int imgidx = 0; imgidx < maxidx; ++imgidx){ //single loop over pixels
            int fpixel = imgidx % fpixels;
            int spixel = imgidx / fpixels;
            double polar=0, omega_pixel=1; // per-pixel locality needed for parallelism
            /* allow for just one part of detector to be rendered */
            if(fpixel < roi_xmin || fpixel > roi_xmax || spixel < roi_ymin || spixel > roi_ymax)
            {
                /*++imgidx*/; continue;
            }
            /* allow for the use of a mask */
            if(maskimage != NULL)
            {
                /* skip any flagged pixels in the mask */
                if(maskimage[imgidx] == 0)
                {
                    /*++imgidx*/; continue;
                }
            }

            /* reset photon count for this pixel */
            double I = 0;

            /* loop over sub-pixels */
            for(int subS=0;subS<oversample;++subS)
            {
                for(int subF=0;subF<oversample;++subF)
                {

                    /* absolute mm position on detector (relative to its origin) */
                    double Fdet = subpixel_size*(fpixel*oversample + subF ) + subpixel_size/2.0;
                    double Sdet = subpixel_size*(spixel*oversample + subS ) + subpixel_size/2.0;
//                  Fdet = pixel_size*fpixel;
//                  Sdet = pixel_size*spixel;

                    for(int thick_tic=0;thick_tic<detector_thicksteps;++thick_tic)
                    {
                        /* assume "distance" is to the front of the detector sensor layer */
                        double Odet = thick_tic*detector_thickstep;

                        /* construct detector subpixel position in 3D space */
//                      pixel_X = distance;
//                      pixel_Y = Sdet-Ybeam;
//                      pixel_Z = Fdet-Xbeam;

                        vec3 pixel_pos(Fdet*fdet_vector[1]+Sdet*sdet_vector[1]+Odet*odet_vector[1]+pix0_vector[1],
                                       Fdet*fdet_vector[2]+Sdet*sdet_vector[2]+Odet*odet_vector[2]+pix0_vector[2],
                                       Fdet*fdet_vector[3]+Sdet*sdet_vector[3]+Odet*odet_vector[3]+pix0_vector[3]
                                      );
                        SCITBX_ASSERT(!curved_detector);
                        /* construct the diffracted-beam unit vector to this sub-pixel */
                        double airpath = pixel_pos.length();
                        vec3 diffracted_v;
                        if (airpath != 0.0) {diffracted_v = pixel_pos.normalize();}
double diffracted[4];
diffracted[0]=diffracted_v.length();
diffracted[1]=diffracted_v[0];
diffracted[2]=diffracted_v[1];
diffracted[3]=diffracted_v[2];
                        if (subS == oversample-1 && subF == oversample-1 && thick_tic==detector_thicksteps-1){
                          /* solid angle subtended by a pixel: (pix/airpath)^2*cos(2theta) */
                          omega_pixel = pixel_size*pixel_size/airpath/airpath*close_distance/airpath;
                          /* option to turn off obliquity effect, inverse-square-law only */
                          if(point_pixel) omega_pixel = 1.0/airpath/airpath;
                        }

                        /* now calculate detector thickness effects */
                        double capture_fraction = 1.;
                        if(detector_thick > 0.0 && detector_attnlen > 0.0) {
                            /* inverse of effective thickness increase */
                            parallax = dot_product(diffracted,odet_vector);
                            capture_fraction = exp(-thick_tic*detector_thickstep/detector_attnlen/parallax)
                                              -exp(-(thick_tic+1)*detector_thickstep/detector_attnlen/parallax);
                        }

                    /* loop over sources now */
                    for(int source=0;source<sources;++source){

                        source_const SC(source, diffracted, this);
                        //source_path = SC.source_path; //breaks const correctness

                        /* rough cut to speed things up when we aren't using whole detector */
                        if(dmin > 0.0 && SC.stol > 0.0)
                        {
                            if(dmin > 0.5/SC.stol)
                            {
                                continue;
                            }
                        }

                        /* sweep over phi angles */
                        for(int phi_tic = 0; phi_tic < phisteps; ++phi_tic)
                        {
                            double phi = phi0 + phistep*phi_tic;
                            phitic_const PC(phi_tic, phi, this);

                            /* enumerate mosaic domains */
                            double I_reduction=0;
                            //# pragma omp parallel for reduction(+:I_reduction)
                            for(int mos_tic=0;mos_tic<mosaic_domains;++mos_tic)
                            {

                                mostic_const MC(mos_tic,diffracted,capture_fraction, source,
                                  PC.ap, PC.bp, PC.cp, SC.scattering, SC.incident, this);
                                /* polarization factor */
                                if (subS == oversample-1 && subF == oversample-1 && thick_tic==detector_thicksteps-1){
                                if (source == sources-1 && phi_tic == phisteps-1 && mos_tic == mosaic_domains-1) {
                                  polar = MC.polar;  //breaks const correctness
if (printout){
  h0 = MC.h0;
  k0 = MC.k0;
  l0 = MC.l0;
  F_cell = MC.F_cell;
}
                                }
                                }
                                /* convert amplitudes into intensity (photons per steradian) */
                                I_reduction += MC.I_increment; //breaks const correctness
                            }
                            /* end of mosaic loop */
                            I+=I_reduction;
                        }
                        /* end of phi loop */
                    }
                    /* end of source loop */
                }
                    /* end of detector thickness loop */
                }
                /* end of sub-pixel y loop */
            }
            /* end of sub-pixel x loop */
            floatimage[imgidx] += r_e_sqr*fluence*polar*I/steps*omega_pixel;
            sum += floatimage[imgidx];
            sumsqr += floatimage[imgidx]*floatimage[imgidx];
            SCITBX_ASSERT(!progress_meter);
            if( printout ){
                if((fpixel==printout_fpixel && spixel==printout_spixel) || printout_fpixel < 0){
                    os << "NanoBragg Structure factor of " << h0 <<" "<< k0<<" "<<l0<<" ";
                    os.precision(17);
                    os << "intensity is " << std::fixed << (F_cell*F_cell) <<"\n";
                    /*
                    twotheta = atan2(sqrt(pixel_pos[1]*pixel_pos[1]+pixel_pos[2]*pixel_pos[2]),pixel_pos[0]);
                    test = sin(twotheta/2.0)/(lambda0*1e10);
                    printf("%4d %4d : stol = %g or %g\n", fpixel,spixel,stol,test);
                    printf("at %g %g %g\n", pixel_pos[0],pixel_pos[1],pixel_pos[2]);
                    printf("hkl= %f %f %f  hkl0= %d %d %d\n", h,k,l,h0,k0,l0);
                    printf(" F_cell=%g  F_latt=%g   I = %g\n", F_cell,F_latt,I);
                    printf("I/steps %15.10g\n", I/steps);
                    printf("polar   %15.10g\n", polar);
                    printf("omega   %15.10g\n", omega_pixel);
                    printf("pixel   %15.10g\n", floatimage[i]);
                    printf("real-space cell vectors (Angstrom):\n");
                    printf("     %-10s  %-10s  %-10s\n","a","b","c");
                    printf("X: %11.8f %11.8f %11.8f\n",a[1]*1e10,b[1]*1e10,c[1]*1e10);
                    printf("Y: %11.8f %11.8f %11.8f\n",a[2]*1e10,b[2]*1e10,c[2]*1e10);
                    printf("Z: %11.8f %11.8f %11.8f\n",a[3]*1e10,b[3]*1e10,c[3]*1e10);
                    */
                }
            }
            /*++imgidx*/;

  } // for loop over pixels
}
// end of add_nanoBragg_spots()

// private copy of add_noise() to enforce const correctness
af::flex_double
nanoBragg::add_noise(af::flex_double panel_pixels) const
{
    SCITBX_ASSERT (panel_pixels.accessor().focus()[0] == raw_pixels.accessor().focus()[0]);
    SCITBX_ASSERT (panel_pixels.accessor().focus()[1] == raw_pixels.accessor().focus()[1]);
    SCITBX_ASSERT (psf_type == UNKNOWN ); // cannot call apply_psf as it is not const correct
    double *floatimage;
    double sum = 0.0;
    int sumn = 0;
    double max_I = 0.0;
    double max_I_x = 0.0; double max_I_y = 0.0; // = 0.0; location of max pixel value
    int fpixel,spixel;
    encapsulated_twodev image_deviates;
    encapsulated_twodev pixel_deviates;
    int i = 0;
    long cseed;

    double expected_photons,observed_photons,adu;
    /* refer to panel_pixels pixel data */
    floatimage = panel_pixels.begin();

    /* don't bother with this loop if calibration is perfect
     NOTE: applying calibration before Poisson noise simulates loss of photons before the detector
     NOTE: applying calibration after Poisson noise simulates systematics in read-out electronics
     here we do the latter */

    /* re-start the RNG */
    long localseed = -labs(seed);

    if(verbose) printf("applying calibration at %g%%, flicker noise at %g%%\n",calibration_noise*100.,flicker_noise*100.);
    i = 0;
    for(spixel=0;spixel<spixels;++spixel)
    {
        for(fpixel=0;fpixel<fpixels;++fpixel)
        {
            /* allow for just one part of detector to be rendered */
            if(fpixel < roi_xmin || fpixel > roi_xmax || spixel < roi_ymin || spixel > roi_ymax)
            {
                ++i; continue;
            }
            /* allow for the use of a mask */
            if(maskimage != NULL)
            {
                /* skip any flagged pixels in the mask */
                if(maskimage[i] == 0)
                {
                    ++i; continue;
                }
            }

            /* take input image to be ideal photons/pixel */
            expected_photons = floatimage[i];

            /* negative photons should be taken as invalid? */
            if(expected_photons < 0.0)
            {
                ++i; continue;
            }

            /* simulate 1/f noise in source */
            if(flicker_noise > 0.0){
                expected_photons *= ( 1.0 + flicker_noise * image_deviates.gaussdev( &localseed ) );
            }
            /* simulate photon-counting error */
            observed_photons = image_deviates.poidev( expected_photons, &localseed );

            /* now we overwrite the flex array, it is now observed, rather than expected photons */
            floatimage[i] = observed_photons;

            /* accumulate number of photons, and keep track of max */
            if(floatimage[i] > max_I) {
                max_I = floatimage[i];
                max_I_x = fpixel;
                max_I_y = spixel;
            }
            sum += observed_photons;
            ++sumn;

            ++i;
        }
    }
    if(verbose) printf("%.0f photons generated on noise image, max= %f at ( %.0f, %.0f )\n",sum,max_I,max_I_x,max_I_y);

    if(calibration_noise > 0.0)
    {
        /* calibration is same from shot to shot, so use well-known seed */
        cseed = -labs(calib_seed);
        sum = max_I = 0.0;
        i = sumn = 0;
        for(spixel=0;spixel<spixels;++spixel)
        {
            for(fpixel=0;fpixel<fpixels;++fpixel)
            {
                /* allow for just one part of detector to be rendered */
                if(fpixel < roi_xmin || fpixel > roi_xmax || spixel < roi_ymin || spixel > roi_ymax)
                {
                    ++i; continue;
                }
                /* allow for the use of a mask */
                if(maskimage != NULL)
                {
                    /* skip any flagged pixels in the mask */
                    if(maskimage[i] == 0)
                    {
                        ++i; continue;
                    }
                }

                /* calibration is same from shot to shot, but varies from pixel to pixel */
                floatimage[i] *= ( 1.0 + calibration_noise * pixel_deviates.gaussdev( &cseed ) );

                /* accumulate number of photons, and keep track of max */
                if(floatimage[i] > max_I) {
                    max_I = floatimage[i];
                    max_I_x = fpixel;
                    max_I_y = spixel;
                }
                sum += floatimage[i];
                ++sumn;

                ++i;
            }
        }
    }
    if(verbose) printf("%.0f photons after calibration error, max= %f at ( %.0f, %.0f )\n",sum,max_I,max_I_x,max_I_y);

    /* now would be a good time to implement PSF?  before we add read-out noise */

    /* now that we have photon count at each point, implement any PSF */
    //if(psf_type != UNKNOWN && psf_fwhm > 0.0)
    //{
    //    /* report on sum before the PSF is applied */
    //    if(verbose) printf("%.0f photons on noise image before PSF\n",sum);
    //    /* start with a clean slate */
    //    if(verbose) printf("  applying PSF width = %g um\n",psf_fwhm*1e6);

    //    apply_psf(psf_type, psf_fwhm/pixel_size, 0);

    //    /* the flex array is now the blurred version of itself, ready for read-out noise */
    //}


    if(verbose) printf("adu = quantum_gain= %g * observed_photons + offset= %g + readout_noise= %g\n",quantum_gain,adc_offset,readout_noise);
    sum = max_I = 0.0;
    i = sumn = 0;
    for(spixel=0;spixel<spixels;++spixel)
    {
        for(fpixel=0;fpixel<fpixels;++fpixel)
        {
            /* allow for just one part of detector to be rendered */
            if(fpixel < roi_xmin || fpixel > roi_xmax || spixel < roi_ymin || spixel > roi_ymax)
            {
                ++i; continue;
            }
            /* allow for the use of a mask */
            if(maskimage != NULL)
            {
                /* skip any flagged pixels in the mask */
                if(maskimage[i] == 0)
                {
                    ++i; continue;
                }
            }

                /* convert photon signal to pixel units */
                adu = floatimage[i]*quantum_gain + adc_offset;

                /* readout noise is in pixel units (adu) */
                if(readout_noise > 0.0){
                    adu += readout_noise * image_deviates.gaussdev( &localseed );
            }

            /* once again, overwriting flex array, this time in ADU units */
            floatimage[i] = adu;

            if(adu > max_I) {
                max_I = adu;
                max_I_x = fpixel;
                max_I_y = spixel;
            }
            sum += adu;
            ++sumn;
            ++i;
        }
    }
    if(verbose) printf("%.0f net adu generated on final image, max= %f at ( %.0f, %.0f )\n",sum-adc_offset*sumn,max_I,max_I_x,max_I_y);
    return panel_pixels;
}
// end of add_noise()


}}// namespace simtbx::nanoBragg
