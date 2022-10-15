#include <simtbx/kokkos/kokkos_types.h>
#include <simtbx/nanoBragg/nanotypes.h>
#include <simtbx/kokkos/kernel_math.h>
using simtbx::nanoBragg::shapetype;
using simtbx::nanoBragg::hklParams;
using simtbx::nanoBragg::SQUARE;
using simtbx::nanoBragg::ROUND;
using simtbx::nanoBragg::GAUSS;
using simtbx::nanoBragg::GAUSS_ARGCHK;
using simtbx::nanoBragg::TOPHAT;

void kokkosSpotsKernel(int spixels, int fpixels, int roi_xmin, int roi_xmax,
    int roi_ymin, int roi_ymax, int oversample, int point_pixel,
    CUDAREAL pixel_size, CUDAREAL subpixel_size, int steps, CUDAREAL detector_thickstep,
    int detector_thicksteps, CUDAREAL detector_thick, CUDAREAL detector_mu,
    const vector_cudareal_t sdet_vector, const vector_cudareal_t fdet_vector,
    const vector_cudareal_t odet_vector, const vector_cudareal_t pix0_vector,
    int curved_detector, CUDAREAL distance, CUDAREAL close_distance,
     const vector_cudareal_t beam_vector,
    CUDAREAL Xbeam, CUDAREAL Ybeam, CUDAREAL dmin, CUDAREAL phi0, CUDAREAL phistep,
    int phisteps, const vector_cudareal_t spindle_vector, int sources,
    const vector_cudareal_t source_X, const vector_cudareal_t source_Y,
    const vector_cudareal_t source_Z,
    const vector_cudareal_t source_I, const vector_cudareal_t source_lambda,
    const vector_cudareal_t a0, const vector_cudareal_t b0,
    const vector_cudareal_t c0, shapetype xtal_shape, CUDAREAL mosaic_spread,
    int mosaic_domains, const vector_cudareal_t mosaic_umats,
    CUDAREAL Na, CUDAREAL Nb,
    CUDAREAL Nc, CUDAREAL V_cell,
    CUDAREAL water_size, CUDAREAL water_F, CUDAREAL water_MW, CUDAREAL r_e_sqr,
    CUDAREAL fluence, CUDAREAL Avogadro, CUDAREAL spot_scale, int integral_form,
    CUDAREAL default_F,
    int interpolate, const vector_cudareal_t Fhkl,
    const hklParams FhklParams, int nopolar,
    const vector_cudareal_t polar_vector, CUDAREAL polarization, CUDAREAL fudge,
    const vector_ushort_t * maskimage, vector_float_t floatimage /*out*/,
    vector_float_t omega_reduction /*out*/, vector_float_t max_I_x_reduction/*out*/,
    vector_float_t max_I_y_reduction /*out*/, vector_bool_t rangemap) {

        const int s_h_min = FhklParams.h_min;
        const int s_k_min = FhklParams.k_min;
        const int s_l_min = FhklParams.l_min;
        const int s_h_range = FhklParams.h_range;
        const int s_k_range = FhklParams.k_range;
        const int s_l_range = FhklParams.l_range;
        const int s_h_max = s_h_min + s_h_range - 1;
        const int s_k_max = s_k_min + s_k_range - 1;
        const int s_l_max = s_l_min + s_l_range - 1;



        const int total_pixels = spixels * fpixels;

        // add background from something amorphous
        CUDAREAL F_bg = water_F;
        CUDAREAL I_bg = F_bg * F_bg * r_e_sqr * fluence * water_size * water_size * water_size * 1e6 * Avogadro / water_MW;

       Kokkos::parallel_for("kokkosSpotsKernel", total_pixels, KOKKOS_LAMBDA(const int& pixIdx) {

                const int fpixel = pixIdx % fpixels;
                const int spixel = pixIdx / fpixels;
                // allow for just a region of interest (roi) on detector to be rendered
                if (fpixel < roi_xmin || fpixel > roi_xmax || spixel < roi_ymin || spixel > roi_ymax) {
                        return;
                }

                // allow for the use of a mask
                if (maskimage != NULL) {
                        // skip any flagged pixels in the mask
                        if ((*maskimage)(pixIdx) == 0) {
                                return;
                        }
                }

                // reset photon count for this pixel
                CUDAREAL I = I_bg;
                CUDAREAL omega_sub_reduction = 0.0;
                CUDAREAL max_I_x_sub_reduction = 0.0;
                CUDAREAL max_I_y_sub_reduction = 0.0;
                CUDAREAL polar = 0.0;
                if (nopolar) {
                        polar = 1.0;
                }

                // add this now to avoid problems with skipping later
                // move this to the bottom to avoid accessing global device memory. floatimage[j] = I_bg;
                // loop over sub-pixels
                int subS, subF;
                for (subS = 0; subS < oversample; ++subS) { // Y voxel
                        for (subF = 0; subF < oversample; ++subF) { // X voxel
                                // absolute mm position on detector (relative to its origin)
                                CUDAREAL Fdet = subpixel_size * (fpixel * oversample + subF) + subpixel_size / 2.0; // X voxel
                                CUDAREAL Sdet = subpixel_size * (spixel * oversample + subS) + subpixel_size / 2.0; // Y voxel
                                // Fdet = pixel_size*fpixel;
                                // Sdet = pixel_size*spixel;

                                max_I_x_sub_reduction = Fdet;
                                max_I_y_sub_reduction = Sdet;

                                int thick_tic;
                                for (thick_tic = 0; thick_tic < detector_thicksteps; ++thick_tic) {
                                        // assume "distance" is to the front of the detector sensor layer
                                        CUDAREAL Odet = thick_tic * detector_thickstep; // Z Orthagonal voxel.

                                        // construct detector subpixel position in 3D space
                                        //                      pixel_X = distance;
                                        //                      pixel_Y = Sdet-Ybeam;
                                        //                      pixel_Z = Fdet-Xbeam;
                                        //CUDAREAL * pixel_pos = tmpVector1;
                                        CUDAREAL pixel_pos[4];
                                        pixel_pos[1] = Fdet * fdet_vector(1)
                                                     + Sdet * sdet_vector(1)
                                                     + Odet * odet_vector(1)
                                                            + pix0_vector(1); // X
                                        pixel_pos[2] = Fdet * fdet_vector(2)
                                                     + Sdet * sdet_vector(2)
                                                     + Odet * odet_vector(2)
                                                            + pix0_vector(2); // Y
                                        pixel_pos[3] = Fdet * fdet_vector(3)
                                                     + Sdet * sdet_vector(3)
                                                     + Odet * odet_vector(3)
                                                            + pix0_vector(3); // Z

                                        if (curved_detector) {
                                                // construct detector pixel that is always "distance" from the sample
                                                CUDAREAL dbvector[] = { 0.0, 0.0, 0.0, 0.0 };
                                                dbvector[1] = distance * beam_vector(1);
                                                dbvector[2] = distance * beam_vector(2);
                                                dbvector[3] = distance * beam_vector(3);
                                                // treat detector pixel coordinates as radians
                                                CUDAREAL newvector[] = { 0.0, 0.0, 0.0, 0.0 };
                                                rotate_axis(dbvector, newvector, sdet_vector, pixel_pos[2] / distance);
                                                rotate_axis(newvector, pixel_pos, fdet_vector, pixel_pos[3] / distance);
                                                // rotate(vector,pixel_pos,0,pixel_pos[3]/distance,pixel_pos[2]/distance);
                                        }

                                        // construct the diffracted-beam unit vector to this sub-pixel
                                        //CUDAREAL * diffracted = tmpVector2;
                                        CUDAREAL diffracted[4];
                                        CUDAREAL airpath = unitize(pixel_pos, diffracted);

                                        // solid angle subtended by a pixel: (pix/airpath)^2*cos(2theta)
                                        CUDAREAL omega_pixel = pixel_size * pixel_size / airpath / airpath * close_distance / airpath;
                                        // option to turn off obliquity effect, inverse-square-law only
                                        if (point_pixel) {
                                                omega_pixel = 1.0 / airpath / airpath;
                                        }

                                        // now calculate detector thickness effects
                                        CUDAREAL capture_fraction = 1.0;
                                        if (detector_thick > 0.0 && detector_mu> 0.0) {
                                                // inverse of effective thickness increase
                                                CUDAREAL odet[4];
                                                odet[1] = odet_vector(1);
                                                odet[2] = odet_vector(2);
                                                odet[3] = odet_vector(3);
                                                CUDAREAL parallax = dot_product(odet, diffracted);
                                                capture_fraction = exp(-thick_tic * detector_thickstep / detector_mu / parallax)
                                                                - exp(-(thick_tic + 1) * detector_thickstep / detector_mu / parallax);
                                        }

                                        // loop over sources now
                                        int source;
                                        for (source = 0; source < sources; ++source) {

                                                // retrieve stuff from cache
                                                CUDAREAL incident[4];
                                                incident[1] = -source_X(source);
                                                incident[2] = -source_Y(source);
                                                incident[3] = -source_Z(source);
                                                CUDAREAL lambda = source_lambda(source);
                                                CUDAREAL source_fraction = source_I(source);

                                                // construct the incident beam unit vector while recovering source distance
                                                // TODO[Giles]: Optimization! We can unitize the source vectors before passing them in.
                                                unitize(incident, incident);

                                                // construct the scattering vector for this pixel
                                                CUDAREAL scattering[4];
                                                scattering[1] = (diffracted[1] - incident[1]) / lambda;
                                                scattering[2] = (diffracted[2] - incident[2]) / lambda;
                                                scattering[3] = (diffracted[3] - incident[3]) / lambda;

                                                #ifdef __CUDA_ARCH__
                                                CUDAREAL stol = 0.5 * norm3d(scattering[1], scattering[2], scattering[3]);
                                                #else
                                                CUDAREAL stol = 0.5 * sqrt(scattering[1]*scattering[1] + scattering[2]*scattering[2] + scattering[3]*scattering[3]);
                                                #endif

                                                // rough cut to speed things up when we aren't using whole detector
                                                if (dmin > 0.0 && stol > 0.0) {
                                                        if (dmin > 0.5 / stol) {
                                                                continue;
                                                        }
                                                }

                                                // polarization factor
                                                if (!nopolar) {
                                                        // need to compute polarization factor
                                                        polar = polarization_factor(polarization, incident, diffracted, polar_vector);
                                                } else {
                                                        polar = 1.0;
                                                }

                                                // sweep over phi angles
                                                for (int phi_tic = 0; phi_tic < phisteps; ++phi_tic) {
                                                        CUDAREAL phi = phistep * phi_tic + phi0;

                                                        CUDAREAL ap[4];
                                                        CUDAREAL bp[4];
                                                        CUDAREAL cp[4];

                                                        // rotate about spindle if necessary
                                                        rotate_axis(a0, ap, spindle_vector, phi);
                                                        rotate_axis(b0, bp, spindle_vector, phi);
                                                        rotate_axis(c0, cp, spindle_vector, phi);

                                                        // enumerate mosaic domains
                                                        for (int mos_tic = 0; mos_tic < mosaic_domains; ++mos_tic) {
                                                                // apply mosaic rotation after phi rotation
                                                                CUDAREAL a[4];
                                                                CUDAREAL b[4];
                                                                CUDAREAL c[4];

                                                                if (mosaic_spread > 0.0) {
                                                                        CUDAREAL umat[] = {mosaic_umats(mos_tic * 9 + 0),
                                                                                           mosaic_umats(mos_tic * 9 + 1),
                                                                                           mosaic_umats(mos_tic * 9 + 2),
                                                                                           mosaic_umats(mos_tic * 9 + 3),
                                                                                           mosaic_umats(mos_tic * 9 + 4),
                                                                                           mosaic_umats(mos_tic * 9 + 5),
                                                                                           mosaic_umats(mos_tic * 9 + 6),
                                                                                           mosaic_umats(mos_tic * 9 + 7),
                                                                                           mosaic_umats(mos_tic * 9 + 8)};

                                                                        rotate_umat(ap, a, umat);
                                                                        rotate_umat(bp, b, umat);
                                                                        rotate_umat(cp, c, umat);
                                                                } else {
                                                                        a[1] = ap[1];
                                                                        a[2] = ap[2];
                                                                        a[3] = ap[3];
                                                                        b[1] = bp[1];
                                                                        b[2] = bp[2];
                                                                        b[3] = bp[3];
                                                                        c[1] = cp[1];
                                                                        c[2] = cp[2];
                                                                        c[3] = cp[3];
                                                                }


                                                                // construct fractional Miller indicies

                                                                CUDAREAL h = dot_product(a, scattering);
                                                                CUDAREAL k = dot_product(b, scattering);
                                                                CUDAREAL l = dot_product(c, scattering);

                                                                // round off to nearest whole index
                                                                int h0 = ceil(h - 0.5);
                                                                int k0 = ceil(k - 0.5);
                                                                int l0 = ceil(l - 0.5);

                                                                // structure factor of the lattice (paralelpiped crystal)
                                                                // F_latt = sin(M_PI*s_Na*h)*sin(M_PI*s_Nb*k)*sin(M_PI*s_Nc*l)/sin(M_PI*h)/sin(M_PI*k)/sin(M_PI*l);
                                                                CUDAREAL F_latt = 1.0; // Shape transform for the crystal.
                                                                CUDAREAL hrad_sqr = 0.0;

                                                                if (xtal_shape == SQUARE) {
                                                                        // xtal is a paralelpiped
                                                                        if (Na > 1) {
                                                                                // F_latt *= sincgrad(h, s_Na);
                                                                                F_latt *= sincg(M_PI * h, Na);
                                                                        }
                                                                        if (Nb > 1) {
                                                                                // F_latt *= sincgrad(k, s_Nb);
                                                                                F_latt *= sincg(M_PI * k, Nb);
                                                                        }
                                                                        if (Nc > 1) {
                                                                                // F_latt *= sincgrad(l, s_Nc);
                                                                                F_latt *= sincg(M_PI * l, Nc);
                                                                        }
                                                                } else {
                                                                        // handy radius in reciprocal space, squared
                                                                        hrad_sqr = (h - h0) * (h - h0) * Na * Na + (k - k0) * (k - k0) * Nb * Nb + (l - l0) * (l - l0) * Nc * Nc;
                                                                }
                                                                if (xtal_shape == ROUND) {
                                                                        // use sinc3 for elliptical xtal shape,
                                                                        // correcting for sqrt of volume ratio between cube and sphere
                                                                        F_latt = Na * Nb * Nc * 0.723601254558268 * sinc3(M_PI * sqrt(hrad_sqr * fudge));
                                                                }
                                                                if (xtal_shape == GAUSS) {
                                                                        // fudge the radius so that volume and FWHM are similar to square_xtal spots
                                                                        F_latt = Na * Nb * Nc * exp(-(hrad_sqr / 0.63 * fudge));
                                                                }
                                                                if (xtal_shape == GAUSS_ARGCHK) {
                                                                        // fudge the radius so that volume and FWHM are similar to square_xtal spots
                                                                        double my_arg = hrad_sqr / 0.63 * fudge;
                                                                        if (my_arg<35.){ F_latt = Na * Nb * Nc * exp(-(my_arg));
                                                                        } else { F_latt = 0.; } // warps coalesce when blocks of 32 pixels have no Bragg signal
                                                                }
                                                                if (xtal_shape == TOPHAT) {
                                                                        // make a flat-top spot of same height and volume as square_xtal spots
                                                                        F_latt = Na * Nb * Nc * (hrad_sqr * fudge < 0.3969);
                                                                }
                                                                // no need to go further if result will be zero?
                                                                if (F_latt == 0.0 && water_size == 0.0)
                                                                        continue;

                                                                // structure factor of the unit cell
                                                                CUDAREAL F_cell = default_F;
                                                                //F_cell = quickFcell_ldg(s_hkls, s_h_max, s_h_min, s_k_max, s_k_min, s_l_max, s_l_min, h0, k0, l0, s_h_range, s_k_range, s_l_range, default_F, Fhkl);
                                                                if (
                                                                        h0 < s_h_min ||
                                                                        k0 < s_k_min ||
                                                                        l0 < s_l_min ||
                                                                        h0 > s_h_max ||
                                                                        k0 > s_k_max ||
                                                                        l0 > s_l_max
                                                                        ) {
                                                                        F_cell = 0.;
                                                                } else {
                                                                        const int hkl_index = (h0-s_h_min)*s_k_range*s_l_range + (k0-s_k_min)*s_l_range + (l0-s_l_min);
                                                                        F_cell = Fhkl[hkl_index];
                                                                }

                                                                // now we have the structure factor for this pixel

                                                                // convert amplitudes into intensity (photons per steradian)
                                                                I += F_cell * F_cell * F_latt * F_latt * source_fraction * capture_fraction * omega_pixel;
                                                                omega_sub_reduction += omega_pixel;
                                                        }
                                                        // end of mosaic loop
                                                }
                                                // end of phi loop
                                        }
                                        // end of source loop
                                }
                                // end of detector thickness loop
                        }
                        // end of sub-pixel y loop
                }
                // end of sub-pixel x loop
                const double photons = I_bg + (r_e_sqr * spot_scale * fluence * polar * I) / steps;
                floatimage( pixIdx ) = photons;
                omega_reduction( pixIdx ) = omega_sub_reduction; // shared contention
                max_I_x_reduction( pixIdx ) = max_I_x_sub_reduction;
                max_I_y_reduction( pixIdx ) = max_I_y_sub_reduction;
                rangemap( pixIdx ) = true;
        });
    }

void debranch_maskall_Kernel(int npanels, int spixels, int fpixels, int total_pixels,
    int oversample, int point_pixel,
    CUDAREAL pixel_size, CUDAREAL subpixel_size, int steps,
    CUDAREAL detector_thickstep, int detector_thicksteps, CUDAREAL detector_thick, CUDAREAL detector_mu,
    const int vec_len,
    const vector_cudareal_t sdet_vector, const vector_cudareal_t fdet_vector,
    const vector_cudareal_t odet_vector, const vector_cudareal_t pix0_vector,
    const vector_cudareal_t distance, const vector_cudareal_t close_distance,
    const vector_cudareal_t beam_vector,
    const vector_cudareal_t Xbeam, const vector_cudareal_t Ybeam, // not even used, after all the work
    CUDAREAL dmin, CUDAREAL phi0, CUDAREAL phistep, int phisteps,
    const vector_cudareal_t spindle_vector, int sources,
    const vector_cudareal_t source_X, const vector_cudareal_t source_Y,
    const vector_cudareal_t source_Z,
    const vector_cudareal_t source_I, const vector_cudareal_t source_lambda,
    const vector_cudareal_t a0, const vector_cudareal_t b0,
    const vector_cudareal_t c0, shapetype xtal_shape,
    int mosaic_domains, const vector_cudareal_t mosaic_umats,
    CUDAREAL Na, CUDAREAL Nb, CUDAREAL Nc, CUDAREAL V_cell, CUDAREAL water_size, CUDAREAL water_F, CUDAREAL water_MW,
    CUDAREAL r_e_sqr, CUDAREAL fluence,
    CUDAREAL Avogadro, CUDAREAL spot_scale, int integral_form, CUDAREAL default_F,
    const vector_cudareal_t Fhkl, const hklParams FhklParams,
    int nopolar, const vector_cudareal_t polar_vector,
    CUDAREAL polarization, CUDAREAL fudge,
    const vector_size_t pixel_lookup,
    vector_float_t floatimage /*out*/, vector_float_t omega_reduction/*out*/,
    vector_float_t max_I_x_reduction/*out*/, vector_float_t max_I_y_reduction /*out*/, vector_bool_t rangemap) {


                const int s_h_min = FhklParams.h_min;
                const int s_k_min = FhklParams.k_min;
                const int s_l_min = FhklParams.l_min;
                const int s_h_range = FhklParams.h_range;
                const int s_k_range = FhklParams.k_range;
                const int s_l_range = FhklParams.l_range;
                const int s_h_max = s_h_min + s_h_range - 1;
                const int s_k_max = s_k_min + s_k_range - 1;
                const int s_l_max = s_l_min + s_l_range - 1;

// Implementation notes.  This kernel is aggressively debranched, therefore the assumptions are:
// 1) mosaicity non-zero positive
// 2) xtal shape is "Gauss" i.e. 3D spheroid.
// 3) No bounds check for access to the structure factor array.
// 4) No check for Flatt=0.


        // add background from something amorphous
        CUDAREAL F_bg = water_F;
        CUDAREAL I_bg = F_bg * F_bg * r_e_sqr * fluence * water_size * water_size * water_size * 1e6 * Avogadro / water_MW;

        Kokkos::parallel_for("debranch_maskall", total_pixels, KOKKOS_LAMBDA(const int& pixIdx) {
                // position in pixel array
                const int j = pixel_lookup(pixIdx);//pixIdx: index into pixel subset; j: index into the data.
                const int i_panel = j / (fpixels*spixels); // the panel number
                const int j_panel = j % (fpixels*spixels); // the pixel number within the panel
                const int fpixel = j_panel % fpixels;
                const int spixel = j_panel / fpixels;

                // reset photon count for this pixel
                CUDAREAL I = I_bg;
                CUDAREAL omega_sub_reduction = 0.0;
                CUDAREAL max_I_x_sub_reduction = 0.0;
                CUDAREAL max_I_y_sub_reduction = 0.0;
                CUDAREAL polar = 0.0;
                if (nopolar) {
                        polar = 1.0;
                }

                // add this now to avoid problems with skipping later
                // move this to the bottom to avoid accessing global device memory. floatimage[j] = I_bg;
                // loop over sub-pixels
                int subS, subF;
                for (subS = 0; subS < oversample; ++subS) { // Y voxel
                        for (subF = 0; subF < oversample; ++subF) { // X voxel
                                // absolute mm position on detector (relative to its origin)
                                CUDAREAL Fdet = subpixel_size * (fpixel * oversample + subF) + subpixel_size / 2.0; // X voxel
                                CUDAREAL Sdet = subpixel_size * (spixel * oversample + subS) + subpixel_size / 2.0; // Y voxel
                                // Fdet = pixel_size*fpixel;
                                // Sdet = pixel_size*spixel;

                                max_I_x_sub_reduction = Fdet;
                                max_I_y_sub_reduction = Sdet;

                                int thick_tic;
                                for (thick_tic = 0; thick_tic < detector_thicksteps; ++thick_tic) {
                                        // assume "distance" is to the front of the detector sensor layer
                                        CUDAREAL Odet = thick_tic * detector_thickstep; // Z Orthagonal voxel.

                                        // construct detector subpixel position in 3D space
                                        //                      pixel_X = distance;
                                        //                      pixel_Y = Sdet-Ybeam;
                                        //                      pixel_Z = Fdet-Xbeam;
                                        //CUDAREAL * pixel_pos = tmpVector1;
                                        CUDAREAL pixel_pos[4];
                                        int iVL = vec_len * i_panel;
                                        pixel_pos[1] = Fdet * fdet_vector(iVL+1)
                                                     + Sdet * sdet_vector(iVL+1)
                                                     + Odet * odet_vector(iVL+1)
                                                                    + pix0_vector(iVL+1); // X
                                        pixel_pos[2] = Fdet * fdet_vector(iVL+2)
                                                     + Sdet * sdet_vector(iVL+2)
                                                     + Odet * odet_vector(iVL+2)
                                                            + pix0_vector(iVL+2); // Y
                                        pixel_pos[3] = Fdet * fdet_vector(iVL+3)
                                                     + Sdet * sdet_vector(iVL+3)
                                                     + Odet * odet_vector(iVL+3)
                                                            + pix0_vector(iVL+3); // Z

                                        // construct the diffracted-beam unit vector to this sub-pixel
                                        //CUDAREAL * diffracted = tmpVector2;
                                        CUDAREAL diffracted[4];
                                        CUDAREAL airpath = unitize(pixel_pos, diffracted);

                                        // solid angle subtended by a pixel: (pix/airpath)^2*cos(2theta)
                                        CUDAREAL omega_pixel = pixel_size * pixel_size / airpath / airpath * close_distance(i_panel) / airpath;
                                        // option to turn off obliquity effect, inverse-square-law only
                                        if (point_pixel) {
                                                omega_pixel = 1.0 / airpath / airpath;
                                        }

                                        // now calculate detector thickness effects
                                        CUDAREAL capture_fraction = 1.0;
                                        if (detector_thick > 0.0 && detector_mu> 0.0) {
                                                // inverse of effective thickness increase
                                                CUDAREAL odet[4];
                                                odet[1] = odet_vector(iVL+1);
                                                odet[2] = odet_vector(iVL+2);
                                                odet[3] = odet_vector(iVL+3);
                                                CUDAREAL parallax = dot_product(odet, diffracted);
                                                capture_fraction = exp(-thick_tic * detector_thickstep / detector_mu / parallax)
                                                                - exp(-(thick_tic + 1) * detector_thickstep / detector_mu / parallax);
                                        }

                                        // loop over sources now
                                        int source;
                                        for (source = 0; source < sources; ++source) {

                                                // retrieve stuff from cache
                                                CUDAREAL incident[4];
                                                incident[1] = -source_X(source);
                                                incident[2] = -source_Y(source);
                                                incident[3] = -source_Z(source);
                                                CUDAREAL lambda = source_lambda(source);
                                                CUDAREAL source_fraction = source_I(source);

                                                // construct the incident beam unit vector while recovering source distance
                                                // TODO[Giles]: Optimization! We can unitize the source vectors before passing them in.
                                                unitize(incident, incident);

                                                // construct the scattering vector for this pixel
                                                CUDAREAL scattering[4];
                                                scattering[1] = (diffracted[1] - incident[1]) / lambda;
                                                scattering[2] = (diffracted[2] - incident[2]) / lambda;
                                                scattering[3] = (diffracted[3] - incident[3]) / lambda;

                                                #ifdef __CUDA_ARCH__
                                                CUDAREAL stol = 0.5 * norm3d(scattering[1], scattering[2], scattering[3]);
                                                #else
                                                CUDAREAL stol = 0.5 * sqrt(scattering[1]*scattering[1] + scattering[2]*scattering[2] + scattering[3]*scattering[3]);
                                                #endif

                                                // rough cut to speed things up when we aren't using whole detector
                                                if (dmin > 0.0 && stol > 0.0) {
                                                        if (dmin > 0.5 / stol) {
                                                                continue;
                                                        }
                                                }

                                                // polarization factor
                                                if (!nopolar) {
                                                        // need to compute polarization factor
                                                        polar = polarization_factor(polarization, incident, diffracted, polar_vector);
                                                } else {
                                                        polar = 1.0;
                                                }

                                                // sweep over phi angles
                                                for (int phi_tic = 0; phi_tic < phisteps; ++phi_tic) {
                                                        CUDAREAL phi = phistep * phi_tic + phi0;

                                                        CUDAREAL ap[4];
                                                        CUDAREAL bp[4];
                                                        CUDAREAL cp[4];

                                                        // rotate about spindle if necessary
                                                        rotate_axis(a0, ap, spindle_vector, phi);
                                                        rotate_axis(b0, bp, spindle_vector, phi);
                                                        rotate_axis(c0, cp, spindle_vector, phi);

                                                        // enumerate mosaic domains
                                                        for (int mos_tic = 0; mos_tic < mosaic_domains; ++mos_tic) {
                                                                // apply mosaic rotation after phi rotation
                                                                CUDAREAL a[4];
                                                                CUDAREAL b[4];
                                                                CUDAREAL c[4];

                                                                CUDAREAL umat[] = {mosaic_umats(mos_tic * 9 + 0),
                                                                                   mosaic_umats(mos_tic * 9 + 1),
                                                                                   mosaic_umats(mos_tic * 9 + 2),
                                                                                   mosaic_umats(mos_tic * 9 + 3),
                                                                                   mosaic_umats(mos_tic * 9 + 4),
                                                                                   mosaic_umats(mos_tic * 9 + 5),
                                                                                   mosaic_umats(mos_tic * 9 + 6),
                                                                                   mosaic_umats(mos_tic * 9 + 7),
                                                                                   mosaic_umats(mos_tic * 9 + 8)};

                                                                rotate_umat(ap, a, umat);
                                                                rotate_umat(bp, b, umat);
                                                                rotate_umat(cp, c, umat);

                                                                // construct fractional Miller indicies

                                                                CUDAREAL h = dot_product(a, scattering);
                                                                CUDAREAL k = dot_product(b, scattering);
                                                                CUDAREAL l = dot_product(c, scattering);

                                                                // round off to nearest whole index
                                                                int h0 = ceil(h - 0.5);
                                                                int k0 = ceil(k - 0.5);
                                                                int l0 = ceil(l - 0.5);

                                                                // structure factor of the lattice (paralelpiped crystal)
                                                                // F_latt = sin(M_PI*s_Na*h)*sin(M_PI*s_Nb*k)*sin(M_PI*s_Nc*l)/sin(M_PI*h)/sin(M_PI*k)/sin(M_PI*l);

                                                                CUDAREAL F_latt = 1.0; // Shape transform for the crystal.
                                                                CUDAREAL hrad_sqr = 0.0;
                                                                // handy radius in reciprocal space, squared
                                                                hrad_sqr = (h - h0) * (h - h0) * Na * Na + (k - k0) * (k - k0) * Nb * Nb + (l - l0) * (l - l0) * Nc * Nc;
                                                                // fudge the radius so that volume and FWHM are similar to square_xtal spots
                                                                double my_arg = hrad_sqr / 0.63 * fudge;
                                                                F_latt = Na * Nb * Nc * exp(-(my_arg));

                                                                // structure factor of the unit cell
                                                                CUDAREAL F_cell = default_F;
                                                                //F_cell = quickFcell_ldg(s_hkls, s_h_max, s_h_min, s_k_max, s_k_min, s_l_max, s_l_min, h0, k0, l0, s_h_range, s_k_range, s_l_range, default_F, Fhkl);
                                                                if (
                                                                        h0 < s_h_min ||
                                                                        k0 < s_k_min ||
                                                                        l0 < s_l_min ||
                                                                        h0 > s_h_max ||
                                                                        k0 > s_k_max ||
                                                                        l0 > s_l_max
                                                                        ) {
                                                                        F_cell = 0.;
                                                                } else {
                                                                        const int hkl_index = (h0-s_h_min)*s_k_range*s_l_range + (k0-s_k_min)*s_l_range + (l0-s_l_min);
                                                                        F_cell = Fhkl[hkl_index];
                                                                }

                                                                // now we have the structure factor for this pixel

                                                                // convert amplitudes into intensity (photons per steradian)
                                                                I += F_cell * F_cell * F_latt * F_latt * source_fraction * capture_fraction * omega_pixel;
                                                                omega_sub_reduction += omega_pixel;
                                                        }
                                                        // end of mosaic loop
                                                }
                                                // end of phi loop
                                        }
                                        // end of source loop
                                }
                                // end of detector thickness loop
                        }
                        // end of sub-pixel y loop
                }
                // end of sub-pixel x loop
                const double photons = I_bg + (r_e_sqr * spot_scale * fluence * polar * I) / steps;
                floatimage( j ) = photons;
                omega_reduction( j ) = omega_sub_reduction; // shared contention
                max_I_x_reduction( j ) = max_I_x_sub_reduction;
                max_I_y_reduction( j ) = max_I_y_sub_reduction;
                rangemap( j ) = true;
        });
}

// __global__ void nanoBraggSpotsInitCUDAKernel(int spixels, int fpixesl, float * floatimage, float * omega_reduction,
//                 float * max_I_x_reduction,
//                 float * max_I_y_reduction, bool * rangemap);


void add_background_kokkos_kernel(int sources, int nanoBragg_oversample,
    CUDAREAL pixel_size, int spixels, int fpixels, int detector_thicksteps,
    CUDAREAL detector_thickstep, CUDAREAL detector_attnlen,
    const vector_cudareal_t  sdet_vector, const vector_cudareal_t  fdet_vector,
    const vector_cudareal_t  odet_vector, const vector_cudareal_t  pix0_vector,
    CUDAREAL close_distance, int point_pixel, CUDAREAL detector_thick,
    const vector_cudareal_t  source_X, const vector_cudareal_t  source_Y,
    const vector_cudareal_t  source_Z,
    const vector_cudareal_t  source_lambda, const vector_cudareal_t  source_I,
    int stols, const vector_cudareal_t stol_of, const vector_cudareal_t Fbg_of,
    int nopolar, CUDAREAL polarization, const vector_cudareal_t  polar_vector,
    CUDAREAL r_e_sqr, CUDAREAL fluence, CUDAREAL amorphous_molecules,
    vector_float_t floatimage)
{

    int oversample=-1, override_source=-1; //override features that usually slow things down,
                                           //like oversampling pixels & multiple sources
    int source_start = 0;
    // allow user to override automated oversampling decision at call time with arguments
    if(oversample<=0) oversample = nanoBragg_oversample;
    if(oversample<=0) oversample = 1;
    if(override_source>=0) {
        // user-specified source in the argument
        source_start = override_source;
        sources = source_start +1;
    }
    // make sure we are normalizing with the right number of sub-steps
    int steps = oversample*oversample;
    CUDAREAL subpixel_size = pixel_size/oversample;

    // sweep over detector
    const int total_pixels = spixels * fpixels;
    // const int fstride = gridDim.x * blockDim.x;
    // const int sstride = gridDim.y * blockDim.y;
    // const int stride = fstride * sstride;
    Kokkos::parallel_for("add_background", total_pixels, KOKKOS_LAMBDA(const int& pixIdx) {

        const int fpixel = pixIdx % fpixels;
        const int spixel = pixIdx / fpixels;
        // reset background photon count for this pixel
        CUDAREAL Ibg = 0;
        int nearest = 0; // sort-stable alogorithm, instead of holding value over from previous pixel
        // loop over sub-pixels
        for(int subS=0; subS<oversample; ++subS) {
            for(int subF=0; subF<oversample; ++subF) {
                // absolute mm position on detector (relative to its origin)
                CUDAREAL Fdet = subpixel_size*(fpixel*oversample + subF ) + subpixel_size/2.0;
                CUDAREAL Sdet = subpixel_size*(spixel*oversample + subS ) + subpixel_size/2.0;

                for(int thick_tic=0; thick_tic<detector_thicksteps; ++thick_tic) {
                    // assume "distance" is to the front of the detector sensor layer
                    CUDAREAL Odet = thick_tic*detector_thickstep;
                    CUDAREAL pixel_pos[4];

                    pixel_pos[0] = 0.0;
                    pixel_pos[1] = Fdet * fdet_vector(1) + Sdet * sdet_vector(1) + Odet * odet_vector(1) + pix0_vector(1); // X
                    pixel_pos[2] = Fdet * fdet_vector(2) + Sdet * sdet_vector(2) + Odet * odet_vector(2) + pix0_vector(2); // Y
                    pixel_pos[3] = Fdet * fdet_vector(3) + Sdet * sdet_vector(3) + Odet * odet_vector(3) + pix0_vector(3); // Z

                    // no curved detector option (future implementation)
                    // construct the diffracted-beam unit vector to this pixel
                    CUDAREAL diffracted[4];
                    CUDAREAL airpath = unitize(pixel_pos, diffracted);

                    // solid angle subtended by a pixel: (pix/airpath)^2*cos(2theta)
                    CUDAREAL omega_pixel = pixel_size*pixel_size/airpath/airpath*close_distance/airpath;
                    // option to turn off obliquity effect, inverse-square-law only
                    if(point_pixel) omega_pixel = 1.0/airpath/airpath;

                    // now calculate detector thickness effects
                    CUDAREAL capture_fraction = 1.0;
                    if(detector_thick > 0.0){
                        // inverse of effective thickness increase
                        CUDAREAL parallax = diffracted[1] * odet_vector(1) + diffracted[2] * odet_vector(2) + diffracted[3] * odet_vector(3);
                        capture_fraction = exp(-thick_tic*detector_thickstep/detector_attnlen/parallax)
                                            -exp(-(thick_tic+1)*detector_thickstep/detector_attnlen/parallax);
                    }

                    // loop over sources now
                    for(int source=source_start; source<sources; ++source) {

                        // retrieve stuff from cache
                        CUDAREAL incident[4];
                        incident[1] = -source_X(source);
                        incident[2] = -source_Y(source);
                        incident[3] = -source_Z(source);
                        CUDAREAL lambda = source_lambda(source);
                        CUDAREAL source_fraction = source_I(source);
                        // construct the incident beam unit vector while recovering source distance
                        unitize(incident, incident);

                        // construct the scattering vector for this pixel
                        CUDAREAL scattering[4];
                        scattering[1] = (diffracted[1]-incident[1])/lambda;
                        scattering[2] = (diffracted[2]-incident[2])/lambda;
                        scattering[3] = (diffracted[3]-incident[3])/lambda;
                        magnitude(scattering);
                        // sin(theta)/lambda is half the scattering vector length
                        CUDAREAL stol = 0.5*scattering[0];

                        // now we need to find the nearest four "stol file" points
                        while(stol > stol_of(nearest) && nearest <= stols){ ++nearest; };
                        while(stol < stol_of(nearest) && nearest >= 2){ --nearest; };

                        // cubic spline interpolation
                        CUDAREAL Fbg;
                        CUDAREAL stol_points[4], Fbg_points[4];
                        stol_points[0] = stol_of(nearest-1);
                        stol_points[1] = stol_of(nearest-1+1);
                        stol_points[2] = stol_of(nearest-1+2);
                        stol_points[3] = stol_of(nearest-1+3);
                        Fbg_points[0] = Fbg_of(nearest-1);
                        Fbg_points[1] = Fbg_of(nearest-1+1);
                        Fbg_points[2] = Fbg_of(nearest-1+2);
                        Fbg_points[3] = Fbg_of(nearest-1+3);

                        polint(stol_points, Fbg_points, stol, &Fbg);

                        // allow negative F values to yield negative intensities
                        CUDAREAL sign=1.0;
                        if(Fbg<0.0) sign=-1.0;

                        // now we have the structure factor for this pixel

                        // polarization factor
                        CUDAREAL polar = 1.0;
                        if(! nopolar){
                            // need to compute polarization factor
                            CUDAREAL axis[] = {polar_vector(0), polar_vector(1), polar_vector(2), polar_vector(3)};
                            polar = polarization_factor(polarization, incident, diffracted, polar_vector);
                        }

                        // accumulate unscaled pixel intensity from this
                        Ibg += sign*Fbg*Fbg*polar*omega_pixel*source_fraction*capture_fraction;
                    } // end of source loop
                } // end of detector thickness loop
            } // end of sub-pixel y loop
        } // end of sub-pixel x loop
        // save photons/pixel (if fluence specified), or F^2/omega if no fluence given
        floatimage(pixIdx) += Ibg*r_e_sqr*fluence*amorphous_molecules/steps;
    }); // end of pixIdx loop
}
