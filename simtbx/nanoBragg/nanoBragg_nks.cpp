#include <simtbx/nanoBragg/nanoBragg.h>

//Contributed by Nicholas Sauter,LBNL.

namespace simtbx {
namespace nanoBragg {

/* add spots from nanocrystal simulation */
void
nanoBragg::add_nanoBragg_spots_nks()
{
  int imgidx = 0;
  floatimage = raw_pixels.begin();

  if(verbose) {printf("TESTING sincg(1,1)= %f\n",sincg(1,1));}

  /* make sure we are normalizing with the right number of sub-steps */
  steps = phisteps*mosaic_domains*oversample*oversample;
  subpixel_size = pixel_size/oversample;

  double sum = 0.0; //reduction variable
  double sumsqr = 0.0; //reduction variable

  for(int spixel=0;spixel<spixels;++spixel){
    for(int fpixel=0;fpixel<fpixels;++fpixel){
            /* allow for just one part of detector to be rendered */
            if(fpixel < roi_xmin || fpixel > roi_xmax || spixel < roi_ymin || spixel > roi_ymax)
            {
                ++imgidx; continue;
            }
            /* allow for the use of a mask */
            if(maskimage != NULL)
            {
                /* skip any flagged pixels in the mask */
                if(maskimage[imgidx] == 0)
                {
                    ++imgidx; continue;
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
                        /* solid angle subtended by a pixel: (pix/airpath)^2*cos(2theta) */
                        omega_pixel = pixel_size*pixel_size/airpath/airpath*close_distance/airpath;
                        /* option to turn off obliquity effect, inverse-square-law only */
                        if(point_pixel) omega_pixel = 1.0/airpath/airpath;

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

                        /* retrieve stuff from cache */
                        incident[1] = -source_X[source];
                        incident[2] = -source_Y[source];
                        incident[3] = -source_Z[source];
                        lambda = source_lambda[source];

                        /* construct the incident beam unit vector while recovering source distance */
                        source_path = unitize(incident,incident);

                        /* construct the scattering vector for this pixel */
                        scattering[1] = (diffracted[1]-incident[1])/lambda;
                        scattering[2] = (diffracted[2]-incident[2])/lambda;
                        scattering[3] = (diffracted[3]-incident[3])/lambda;

                        /* sin(theta)/lambda is half the scattering vector length */
                        stol = 0.5*magnitude(scattering);

                        /* rough cut to speed things up when we aren't using whole detector */
                        if(dmin > 0.0 && stol > 0.0)
                        {
                            if(dmin > 0.5/stol)
                            {
                                continue;
                            }
                        }

                        /* sweep over phi angles */
                        for(int phi_tic = 0; phi_tic < phisteps; ++phi_tic)
                        {
                            phi = phi0 + phistep*phi_tic;

                            if( phi != 0.0 )
                            {
                                /* rotate about spindle if neccesary */
                                rotate_axis(a0,ap,spindle_vector,phi);
                                rotate_axis(b0,bp,spindle_vector,phi);
                                rotate_axis(c0,cp,spindle_vector,phi);
                            }

                            /* enumerate mosaic domains */
                            for(mos_tic=0;mos_tic<mosaic_domains;++mos_tic)
                            {
                                /* apply mosaic rotation after phi rotation */
                                if( mosaic_spread > 0.0 )
                                {
                                    rotate_umat(ap,a,&mosaic_umats[mos_tic*9]);
                                    rotate_umat(bp,b,&mosaic_umats[mos_tic*9]);
                                    rotate_umat(cp,c,&mosaic_umats[mos_tic*9]);
                                }
                                else
                                {
                                    a[1]=ap[1];a[2]=ap[2];a[3]=ap[3];
                                    b[1]=bp[1];b[2]=bp[2];b[3]=bp[3];
                                    c[1]=cp[1];c[2]=cp[2];c[3]=cp[3];
                                }
//                              printf("%d %f %f %f\n",mos_tic,mosaic_umats[mos_tic*9+0],mosaic_umats[mos_tic*9+1],mosaic_umats[mos_tic*9+2]);
//                              printf("%d %f %f %f\n",mos_tic,mosaic_umats[mos_tic*9+3],mosaic_umats[mos_tic*9+4],mosaic_umats[mos_tic*9+5]);
//                              printf("%d %f %f %f\n",mos_tic,mosaic_umats[mos_tic*9+6],mosaic_umats[mos_tic*9+7],mosaic_umats[mos_tic*9+8]);

                                /* construct fractional Miller indicies */
                                double h = dot_product(a,scattering);
                                double k = dot_product(b,scattering);
                                double l = dot_product(c,scattering);

                                /* round off to nearest whole index */
                                int h0 = static_cast<int>(ceil(h-0.5));
                                int k0 = static_cast<int>(ceil(k-0.5));
                                int l0 = static_cast<int>(ceil(l-0.5));


                                // structure factor of the lattice (paralelpiped crystal)
                                double F_latt = 1.0;
                                double hrad_sqr = 0.;
                                if(xtal_shape == SQUARE){ /* xtal is a paralelpiped */
                                  if(Na>1){ F_latt *= sincg(M_PI*h,Na); }
                                  if(Nb>1){ F_latt *= sincg(M_PI*k,Nb); }
                                  if(Nc>1){ F_latt *= sincg(M_PI*l,Nc); }
                                } else { /* handy radius in reciprocal space, squared */
                                  hrad_sqr = (h-h0)*(h-h0)*Na*Na + (k-k0)*(k-k0)*Nb*Nb + (l-l0)*(l-l0)*Nc*Nc ;
                                }
                                if(xtal_shape == ROUND){ /* use sinc3 for elliptical xtal shape,
                                           correcting for sqrt of volume ratio between cube and sphere */
                                  F_latt = Na*Nb*Nc*0.723601254558268*sinc3(M_PI*sqrt( hrad_sqr * fudge ) );
                                }
                                if(xtal_shape == GAUSS){
                                        /* fudge the radius so that volume and FWHM are similar to square_xtal spots */
                                  F_latt = Na*Nb*Nc*exp(-( hrad_sqr / 0.63 * fudge ));
                                }
                                if(xtal_shape == TOPHAT) {
                                        /* make a flat-top spot of same height and volume as square_xtal spots */
                                  F_latt = Na*Nb*Nc*(hrad_sqr*fudge < 0.3969 );
                                }
                                /* no need to go further if result will be zero */
                                if(F_latt == 0.0) continue;

                                /* find nearest point on Ewald sphere surface? */
                                SCITBX_ASSERT(!integral_form);
                                SCITBX_ASSERT(!interpolate);
                                /* structure factor of the unit cell */
                                double F_cell;
                                if ( (h0<=h_max) && (h0>=h_min) && (k0<=k_max) && (k0>=k_min) && (l0<=l_max) && (l0>=l_min)  ) {
                                        /* just take nearest-neighbor */
                                        F_cell = Fhkl[h0-h_min][k0-k_min][l0-l_min];
                                }else{
                                        F_cell = default_F; // usually zero
                                }

                                /* now we have the structure factor for this pixel */

                                /* polarization factor */
                                if(! nopolar){
                                    /* need to compute polarization factor */
                                    polar = polarization_factor(polarization,incident,diffracted,polar_vector);
                                }
                                else
                                {
                                    polar = 1.0;
                                }

                                /* convert amplitudes into intensity (photons per steradian) */
                                I += F_cell*F_cell*F_latt*F_latt*source_I[source]*capture_fraction;
                            }
                            /* end of mosaic loop */
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
                }
            }
            ++imgidx;
    }  // for loop over fast pixels
  } // for loop over slow pixels
}
// end of add_nanoBragg_spots()



}}// namespace simtbx::nanoBragg
