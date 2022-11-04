#ifndef SIMTBX_GPU_SIMULATION_CUH
#define SIMTBX_GPU_SIMULATION_CUH

#include <simtbx/nanoBragg/nanotypes.h>
#include <simtbx/nanoBragg/nanoBraggCUDA.cuh>
using simtbx::nanoBragg::shapetype;
using simtbx::nanoBragg::hklParams;

__global__ void nanoBraggSpotsCUDAKernel(int spixels, int fpixels, int roi_xmin, int roi_xmax,
    int roi_ymin, int roi_ymax, int oversample, int point_pixel,
    CUDAREAL pixel_size, CUDAREAL subpixel_size, int steps, CUDAREAL detector_thickstep,
    int detector_thicksteps, CUDAREAL detector_thick, CUDAREAL detector_mu,
    const CUDAREAL * __restrict__ sdet_vector, const CUDAREAL * __restrict__ fdet_vector,
    const CUDAREAL * __restrict__ odet_vector,
    const CUDAREAL * __restrict__ pix0_vector, int curved_detector, CUDAREAL distance,
    CUDAREAL close_distance, const CUDAREAL * __restrict__ beam_vector,
    CUDAREAL Xbeam, CUDAREAL Ybeam, CUDAREAL dmin, CUDAREAL phi0, CUDAREAL phistep,
    int phisteps, const CUDAREAL * __restrict__ spindle_vector, int sources,
    const CUDAREAL * __restrict__ source_X, const CUDAREAL * __restrict__ source_Y,
    const CUDAREAL * __restrict__ source_Z,
    const CUDAREAL * __restrict__ source_I, const CUDAREAL * __restrict__ source_lambda,
    const CUDAREAL * __restrict__ a0, const CUDAREAL * __restrict__ b0,
    const CUDAREAL * __restrict c0, shapetype xtal_shape, CUDAREAL mosaic_spread,
    int mosaic_domains, const CUDAREAL * __restrict__ mosaic_umats,
    CUDAREAL Na, CUDAREAL Nb,
    CUDAREAL Nc, CUDAREAL V_cell,
    CUDAREAL water_size, CUDAREAL water_F, CUDAREAL water_MW, CUDAREAL r_e_sqr,
    CUDAREAL fluence, CUDAREAL Avogadro, CUDAREAL spot_scale, int integral_form,
    CUDAREAL default_F,
    int interpolate, const CUDAREAL * __restrict__ Fhkl,
    const hklParams * __restrict__ Fhklparams, int nopolar,
    const CUDAREAL * __restrict__ polar_vector, CUDAREAL polarization, CUDAREAL fudge,
    const int unsigned short * __restrict__ maskimage, float * floatimage /*out*/,
    float * omega_reduction/*out*/, float * max_I_x_reduction/*out*/,
    float * max_I_y_reduction /*out*/, bool * rangemap);

__global__ void debranch_maskall_CUDAKernel(int npanels, int spixels, int fpixels, int total_pixels,
    int oversample, int point_pixel,
    CUDAREAL pixel_size, CUDAREAL subpixel_size, int steps,
    CUDAREAL detector_thickstep, int detector_thicksteps, CUDAREAL detector_thick, CUDAREAL detector_mu,
    const int vec_len,
    const CUDAREAL * __restrict__ sdet_vector, const CUDAREAL * __restrict__ fdet_vector,
    const CUDAREAL * __restrict__ odet_vector,
    const CUDAREAL * __restrict__ pix0_vector,
    const CUDAREAL * __restrict__ distance, const CUDAREAL * __restrict__ close_distance,
    const CUDAREAL * __restrict__ beam_vector,
    const CUDAREAL * __restrict__ Xbeam, const CUDAREAL * __restrict__ Ybeam, // not even used, after all the work
    CUDAREAL dmin, CUDAREAL phi0, CUDAREAL phistep, int phisteps,
    const CUDAREAL * __restrict__ spindle_vector, int sources,
    const CUDAREAL * __restrict__ source_X, const CUDAREAL * __restrict__ source_Y,
    const CUDAREAL * __restrict__ source_Z,
    const CUDAREAL * __restrict__ source_I, const CUDAREAL * __restrict__ source_lambda,
    const CUDAREAL * __restrict__ a0, const CUDAREAL * __restrict__ b0,
    const CUDAREAL * __restrict__ c0, shapetype xtal_shape,
    int mosaic_domains, const CUDAREAL * __restrict__ mosaic_umats,
    CUDAREAL Na, CUDAREAL Nb, CUDAREAL Nc, CUDAREAL V_cell, CUDAREAL water_size, CUDAREAL water_F, CUDAREAL water_MW,
    CUDAREAL r_e_sqr, CUDAREAL fluence,
    CUDAREAL Avogadro, CUDAREAL spot_scale, int integral_form, CUDAREAL default_F,
    const CUDAREAL * __restrict__ Fhkl, const hklParams * __restrict__ FhklParams,
    int nopolar, const CUDAREAL * __restrict__ polar_vector,
    CUDAREAL polarization, CUDAREAL fudge,
    const std::size_t * __restrict__ pixel_lookup,
    float * floatimage /*out*/, float * omega_reduction/*out*/,
    float * max_I_x_reduction/*out*/, float * max_I_y_reduction /*out*/, bool * rangemap) {
        __shared__ int s_vec_len;
	__shared__ CUDAREAL s_dmin;

	__shared__ bool s_nopolar;

	__shared__ int s_phisteps;
	__shared__ CUDAREAL s_phi0, s_phistep;
	__shared__ int s_mosaic_domains;

	__shared__ CUDAREAL s_Na, s_Nb, s_Nc;
	__shared__ int s_h_min, s_k_min, s_l_min, s_h_range, s_k_range, s_l_range,
                       s_h_max, s_k_max, s_l_max;

	if (threadIdx.x == 0 && threadIdx.y == 0) {
                s_vec_len = vec_len;
		s_dmin = dmin;

		s_nopolar = nopolar;

		s_phisteps = phisteps;
		s_phi0 = phi0;
		s_phistep = phistep;

		s_mosaic_domains = mosaic_domains;

		s_Na = Na;
		s_Nb = Nb;
		s_Nc = Nc;

		s_h_min = FhklParams->h_min;
		s_k_min = FhklParams->k_min;
		s_l_min = FhklParams->l_min;
		s_h_range = FhklParams->h_range;
		s_k_range = FhklParams->k_range;
		s_l_range = FhklParams->l_range;
                s_h_max = s_h_min + s_h_range - 1;
                s_k_max = s_k_min + s_k_range - 1;
                s_l_max = s_l_min + s_l_range - 1;

	}
	__syncthreads();
/* Implementation notes.  This kernel is aggressively debranched, therefore the assumptions are:
1) mosaicity non-zero positive
2) xtal shape is "Gauss" i.e. 3D spheroid.
3) No bounds check for access to the structure factor array.
4) No check for Flatt=0.
*/
	//NKS new design, one-function call covers all panels with mask
	const int fstride = gridDim.x * blockDim.x;
	const int sstride = gridDim.y * blockDim.y;
	const int stride = fstride * sstride;

	/* add background from something amorphous */
	CUDAREAL F_bg = water_F;
	CUDAREAL I_bg = F_bg * F_bg * r_e_sqr * fluence * water_size * water_size * water_size * 1e6 * Avogadro / water_MW;

	for (int pixIdx = (blockDim.y * blockIdx.y + threadIdx.y) * fstride + blockDim.x * blockIdx.x + threadIdx.x;
             pixIdx < total_pixels;
             pixIdx += stride) {
		/* position in pixel array */
		const int j = pixel_lookup[pixIdx];//pixIdx: index into pixel subset; j: index into the data.
                const int i_panel = j / (fpixels*spixels); // the panel number
                const int j_panel = j % (fpixels*spixels); // the pixel number within the panel
		const int fpixel = j_panel % fpixels;
		const int spixel = j_panel / fpixels;

		/* reset photon count for this pixel */
		CUDAREAL I = I_bg;
		CUDAREAL omega_sub_reduction = 0.0;
		CUDAREAL max_I_x_sub_reduction = 0.0;
		CUDAREAL max_I_y_sub_reduction = 0.0;
		CUDAREAL polar = 0.0;
		if (s_nopolar) {
			polar = 1.0;
		}

		/* add this now to avoid problems with skipping later */
		// move this to the bottom to avoid accessing global device memory. floatimage[j] = I_bg;
		/* loop over sub-pixels */
		int subS, subF;
		for (subS = 0; subS < oversample; ++subS) { // Y voxel
			for (subF = 0; subF < oversample; ++subF) { // X voxel
				/* absolute mm position on detector (relative to its origin) */
				CUDAREAL Fdet = subpixel_size * (fpixel * oversample + subF) + subpixel_size / 2.0; // X voxel
				CUDAREAL Sdet = subpixel_size * (spixel * oversample + subS) + subpixel_size / 2.0; // Y voxel
				//                  Fdet = pixel_size*fpixel;
				//                  Sdet = pixel_size*spixel;

				max_I_x_sub_reduction = Fdet;
				max_I_y_sub_reduction = Sdet;

				int thick_tic;
				for (thick_tic = 0; thick_tic < detector_thicksteps; ++thick_tic) {
					/* assume "distance" is to the front of the detector sensor layer */
					CUDAREAL Odet = thick_tic * detector_thickstep; // Z Orthagonal voxel.

					/* construct detector subpixel position in 3D space */
					//                      pixel_X = distance;
					//                      pixel_Y = Sdet-Ybeam;
					//                      pixel_Z = Fdet-Xbeam;
					//CUDAREAL * pixel_pos = tmpVector1;
					CUDAREAL pixel_pos[4];
                                        int iVL = s_vec_len * i_panel;
                                        pixel_pos[1] = Fdet * __ldg(&fdet_vector[iVL+1]) +
                                                       Sdet * __ldg(&sdet_vector[iVL+1]) +
                                                       Odet * __ldg(&odet_vector[iVL+1]) +
                                                              __ldg(&pix0_vector[iVL+1]); // X
                                        pixel_pos[2] = Fdet * __ldg(&fdet_vector[iVL+2]) +
                                                       Sdet * __ldg(&sdet_vector[iVL+2]) +
                                                       Odet * __ldg(&odet_vector[iVL+2]) +
                                                              __ldg(&pix0_vector[iVL+2]); // X
                                        pixel_pos[3] = Fdet * __ldg(&fdet_vector[iVL+3]) +
                                                       Sdet * __ldg(&sdet_vector[iVL+3]) +
                                                       Odet * __ldg(&odet_vector[iVL+3]) +
                                                              __ldg(&pix0_vector[iVL+3]); // X

					/* construct the diffracted-beam unit vector to this sub-pixel */
					//CUDAREAL * diffracted = tmpVector2;
					CUDAREAL diffracted[4];
					CUDAREAL airpath = unitize(pixel_pos, diffracted);

					/* solid angle subtended by a pixel: (pix/airpath)^2*cos(2theta) */
					CUDAREAL omega_pixel = pixel_size * pixel_size / airpath / airpath * close_distance[i_panel] / airpath;
					/* option to turn off obliquity effect, inverse-square-law only */
					if (point_pixel) {
						omega_pixel = 1.0 / airpath / airpath;
					}

					/* now calculate detector thickness effects */
					CUDAREAL capture_fraction = 1.0;
					if (detector_thick > 0.0 && detector_mu> 0.0) {
						/* inverse of effective thickness increase */
						CUDAREAL parallax = dot_product_ldg(&(odet_vector[iVL]), diffracted);
						capture_fraction = exp(-thick_tic * detector_thickstep / detector_mu / parallax)
								- exp(-(thick_tic + 1) * detector_thickstep / detector_mu / parallax);
					}

					/* loop over sources now */
					int source;
					for (source = 0; source < sources; ++source) {

						/* retrieve stuff from cache */
						CUDAREAL incident[4];
						incident[1] = -__ldg(&source_X[source]);
						incident[2] = -__ldg(&source_Y[source]);
						incident[3] = -__ldg(&source_Z[source]);
						CUDAREAL lambda = __ldg(&source_lambda[source]);
						CUDAREAL source_fraction = __ldg(&source_I[source]);

						/* construct the incident beam unit vector while recovering source distance */
						// TODO[Giles]: Optimization! We can unitize the source vectors before passing them in.
						unitize(incident, incident);

						/* construct the scattering vector for this pixel */
						CUDAREAL scattering[4];
						scattering[1] = (diffracted[1] - incident[1]) / lambda;
						scattering[2] = (diffracted[2] - incident[2]) / lambda;
						scattering[3] = (diffracted[3] - incident[3]) / lambda;

						CUDAREAL stol = 0.5 * norm3d(scattering[1], scattering[2], scattering[3]);

						/* rough cut to speed things up when we aren't using whole detector */
						if (s_dmin > 0.0 && stol > 0.0) {
							if (s_dmin > 0.5 / stol) {
								continue;
							}
						}

						/* polarization factor */
						if (!s_nopolar) {
							/* need to compute polarization factor */
							polar = polarization_factor(polarization, incident, diffracted, polar_vector);
						} else {
							polar = 1.0;
						}

						/* sweep over phi angles */
						for (int phi_tic = 0; phi_tic < s_phisteps; ++phi_tic) {
							CUDAREAL phi = s_phistep * phi_tic + s_phi0;

							CUDAREAL ap[4];
							CUDAREAL bp[4];
							CUDAREAL cp[4];

							/* rotate about spindle if necessary */
							rotate_axis_ldg(a0, ap, spindle_vector, phi);
							rotate_axis_ldg(b0, bp, spindle_vector, phi);
							rotate_axis_ldg(c0, cp, spindle_vector, phi);

							/* enumerate mosaic domains */
							for (int mos_tic = 0; mos_tic < s_mosaic_domains; ++mos_tic) {
								/* apply mosaic rotation after phi rotation */
								CUDAREAL a[4];
								CUDAREAL b[4];
								CUDAREAL c[4];

								rotate_umat_ldg(ap, a, &mosaic_umats[mos_tic * 9]);
								rotate_umat_ldg(bp, b, &mosaic_umats[mos_tic * 9]);
								rotate_umat_ldg(cp, c, &mosaic_umats[mos_tic * 9]);

								/* construct fractional Miller indicies */

								CUDAREAL h = dot_product(a, scattering);
								CUDAREAL k = dot_product(b, scattering);
								CUDAREAL l = dot_product(c, scattering);

								/* round off to nearest whole index */
								int h0 = ceil(h - 0.5);
								int k0 = ceil(k - 0.5);
								int l0 = ceil(l - 0.5);

								/* structure factor of the lattice (paralelpiped crystal)
								 F_latt = sin(M_PI*s_Na*h)*sin(M_PI*s_Nb*k)*sin(M_PI*s_Nc*l)/sin(M_PI*h)/sin(M_PI*k)/sin(M_PI*l);
								 */
								CUDAREAL F_latt = 1.0; // Shape transform for the crystal.
								CUDAREAL hrad_sqr = 0.0;
								/* handy radius in reciprocal space, squared */
								hrad_sqr = (h - h0) * (h - h0) * s_Na * s_Na + (k - k0) * (k - k0) * s_Nb * s_Nb + (l - l0) * (l - l0) * s_Nc * s_Nc;
                                                                /* fudge the radius so that volume and FWHM are similar to square_xtal spots */
                                                                double my_arg = hrad_sqr / 0.63 * fudge;
                                                                F_latt = s_Na * s_Nb * s_Nc * exp(-(my_arg));

								/* structure factor of the unit cell */
								CUDAREAL F_cell = default_F;
								//F_cell = quickFcell_ldg(s_hkls, s_h_max, s_h_min, s_k_max, s_k_min, s_l_max, s_l_min, h0, k0, l0, s_h_range, s_k_range, s_l_range, default_F, Fhkl);
                                                                if (
                                                                    h0 < s_h_min ||
                                                                    k0 < s_k_min ||
                                                                    l0 < s_l_min ||
                                                                    h0 > s_h_max ||
                                                                    k0 > s_k_max ||
                                                                    l0 > s_l_max
                                                                   )
                                                                  F_cell = 0.;
                                                                else
                                                                  F_cell = __ldg(&Fhkl[(h0-s_h_min)*s_k_range*s_l_range + (k0-s_k_min)*s_l_range + (l0-s_l_min)]);

								/* now we have the structure factor for this pixel */

								/* convert amplitudes into intensity (photons per steradian) */
								I += F_cell * F_cell * F_latt * F_latt * source_fraction * capture_fraction * omega_pixel;
								omega_sub_reduction += omega_pixel;
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
		const double photons = I_bg + (r_e_sqr * spot_scale * fluence * polar * I) / steps;
		floatimage[j] = photons;
		omega_reduction[j] = omega_sub_reduction; // shared contention
		max_I_x_reduction[j] = max_I_x_sub_reduction;
		max_I_y_reduction[j] = max_I_y_sub_reduction;
		rangemap[j] = true;
	}
}

__global__ void nanoBraggSpotsInitCUDAKernel(int spixels, int fpixesl, float * floatimage, float * omega_reduction,
                float * max_I_x_reduction,
                float * max_I_y_reduction, bool * rangemap);

__global__ void add_array_CUDAKernel(double * lhs, float * rhs, int array_size){
  const int total_pixels = array_size;
  const int fstride = gridDim.x * blockDim.x;
  const int sstride = gridDim.y * blockDim.y;
  const int stride = fstride * sstride;
  for (int pixIdx =
      (blockDim.y * blockIdx.y + threadIdx.y) * fstride + blockDim.x * blockIdx.x + threadIdx.x;
      pixIdx < total_pixels; pixIdx += stride) {
    const int j = pixIdx; /* position in pixel array */
    lhs[j] = lhs[j] + (double)rhs[j]; // specifically add low precision to high precision array
    rhs[j] = 0; // reset deepest array for next calculation
    }
  }

__global__ void add_background_CUDAKernel(int sources, int nanoBragg_oversample,
    CUDAREAL pixel_size, int spixels, int fpixels, int detector_thicksteps,
    CUDAREAL detector_thickstep, CUDAREAL detector_attnlen,
    const CUDAREAL * __restrict__ sdet_vector, const CUDAREAL * __restrict__ fdet_vector,
    const CUDAREAL * __restrict__ odet_vector, const CUDAREAL * __restrict__ pix0_vector,
    CUDAREAL close_distance, int point_pixel, CUDAREAL detector_thick,
    const CUDAREAL * __restrict__ source_X, const CUDAREAL * __restrict__ source_Y,
    const CUDAREAL * __restrict__ source_Z,
    const CUDAREAL * __restrict__ source_lambda, const CUDAREAL * __restrict__ source_I,
    int stols, const CUDAREAL * stol_of, const CUDAREAL * Fbg_of,
    int nopolar, CUDAREAL polarization, const CUDAREAL * __restrict__ polar_vector,
    CUDAREAL r_e_sqr, CUDAREAL fluence, CUDAREAL amorphous_molecules,
    float * floatimage);

__global__ void add_background_CUDAKernel(int sources, int nanoBragg_oversample, int override_source,
    CUDAREAL pixel_size, int spixels, int fpixels, int detector_thicksteps,
    CUDAREAL detector_thickstep, CUDAREAL detector_attnlen,
    const CUDAREAL * __restrict__ sdet_vector, const CUDAREAL * __restrict__ fdet_vector,
    const CUDAREAL * __restrict__ odet_vector, const CUDAREAL * __restrict__ pix0_vector,
    CUDAREAL close_distance, int point_pixel, CUDAREAL detector_thick,
    const CUDAREAL * __restrict__ source_X, const CUDAREAL * __restrict__ source_Y,
    const CUDAREAL * __restrict__ source_Z,
    const CUDAREAL * __restrict__ source_lambda, const CUDAREAL * __restrict__ source_I,
    int stols, const CUDAREAL * stol_of, const CUDAREAL * Fbg_of,
    int nopolar, CUDAREAL polarization, const CUDAREAL * __restrict__ polar_vector,
    CUDAREAL r_e_sqr, CUDAREAL fluence, CUDAREAL amorphous_molecules,
    float * floatimage)
{
    int oversample=-1;                     //override features that usually slow things down,
                                           //like oversampling pixels & multiple sources
    int source_start = 0;
    int orig_sources = sources;
    int end_sources = sources;
    /* allow user to override automated oversampling decision at call time with arguments */
    if(oversample<=0) oversample = nanoBragg_oversample;
    if(oversample<=0) oversample = 1;
    bool have_single_source = false;
    if(override_source>=0) {
        /* user-specified source in the argument */
        source_start = override_source;
        end_sources = source_start +1;
        have_single_source = true;
    }
    /* make sure we are normalizing with the right number of sub-steps */
    int steps = oversample*oversample;
    CUDAREAL subpixel_size = pixel_size/oversample;

    /* sweep over detector */
    const int total_pixels = spixels * fpixels;
    const int fstride = gridDim.x * blockDim.x;
    const int sstride = gridDim.y * blockDim.y;
    const int stride = fstride * sstride;
    for (int pixIdx = (blockDim.y * blockIdx.y + threadIdx.y) * fstride + blockDim.x * blockIdx.x + threadIdx.x;
         pixIdx < total_pixels; pixIdx += stride) {
      const int fpixel = pixIdx % fpixels;
      const int spixel = pixIdx / fpixels;
      /* position in pixel array */
      const int j = pixIdx;
      /* reset background photon count for this pixel */
      CUDAREAL Ibg = 0;
      int nearest = 0; // sort-stable alogorithm, instead of holding value over from previous pixel
            /* loop over sub-pixels */
            for(int subS=0;subS<oversample;++subS){
                for(int subF=0;subF<oversample;++subF){
                    /* absolute mm position on detector (relative to its origin) */
                    CUDAREAL Fdet = subpixel_size*(fpixel*oversample + subF ) + subpixel_size/2.0;
                    CUDAREAL Sdet = subpixel_size*(spixel*oversample + subS ) + subpixel_size/2.0;

                    for(int thick_tic=0;thick_tic<detector_thicksteps;++thick_tic){
                        /* assume "distance" is to the front of the detector sensor layer */
                        CUDAREAL Odet = thick_tic*detector_thickstep;
                        CUDAREAL pixel_pos[4];

                        pixel_pos[1] = Fdet * __ldg(&fdet_vector[1]) + Sdet * __ldg(&sdet_vector[1]) + Odet * __ldg(&odet_vector[1]) + __ldg(&pix0_vector[1]); // X
                        pixel_pos[2] = Fdet * __ldg(&fdet_vector[2]) + Sdet * __ldg(&sdet_vector[2]) + Odet * __ldg(&odet_vector[2]) + __ldg(&pix0_vector[2]); // X
                        pixel_pos[3] = Fdet * __ldg(&fdet_vector[3]) + Sdet * __ldg(&sdet_vector[3]) + Odet * __ldg(&odet_vector[3]) + __ldg(&pix0_vector[3]); // X
                        pixel_pos[0] = 0.0;
                        /* no curved detector option (future implementation) */
                        /* construct the diffracted-beam unit vector to this pixel */
                        CUDAREAL diffracted[4];
                        CUDAREAL airpath = unitize(pixel_pos,diffracted);

                        /* solid angle subtended by a pixel: (pix/airpath)^2*cos(2theta) */
                        CUDAREAL omega_pixel = pixel_size*pixel_size/airpath/airpath*close_distance/airpath;
                        /* option to turn off obliquity effect, inverse-square-law only */
                        if(point_pixel) omega_pixel = 1.0/airpath/airpath;

                        /* now calculate detector thickness effects */
                        CUDAREAL capture_fraction = 1.0;
                        if(detector_thick > 0.0){
                            /* inverse of effective thickness increase */
                            CUDAREAL parallax = dot_product(diffracted,odet_vector);
                            capture_fraction = exp(-thick_tic*detector_thickstep/detector_attnlen/parallax)
                                              -exp(-(thick_tic+1)*detector_thickstep/detector_attnlen/parallax);
                        }

                        /* loop over sources now */
                        for(int source=source_start; source < end_sources; ++source){
                            CUDAREAL n_source_scale = (have_single_source) ? orig_sources : __ldg(&source_I[source]);

                            /* retrieve stuff from cache */
                            CUDAREAL incident[4];
                            incident[1] = -__ldg(&source_X[source]);
                            incident[2] = -__ldg(&source_Y[source]);
                            incident[3] = -__ldg(&source_Z[source]);
                            CUDAREAL lambda = __ldg(&source_lambda[source]);
                            /* construct the incident beam unit vector while recovering source distance */
                            unitize(incident,incident);

                            /* construct the scattering vector for this pixel */
                            CUDAREAL scattering[4];
                            scattering[1] = (diffracted[1]-incident[1])/lambda;
                            scattering[2] = (diffracted[2]-incident[2])/lambda;
                            scattering[3] = (diffracted[3]-incident[3])/lambda;
                            magnitude(scattering);
                            /* sin(theta)/lambda is half the scattering vector length */
                            CUDAREAL stol = 0.5*scattering[0];

                            /* now we need to find the nearest four "stol file" points */
                            while(stol > stol_of[nearest] && nearest <= stols){++nearest; };
                            while(stol < stol_of[nearest] && nearest >= 2){--nearest; };

                            /* cubic spline interpolation */
                            CUDAREAL Fbg;
                            polint(stol_of+nearest-1, Fbg_of+nearest-1, stol, &Fbg);

                            /* allow negative F values to yield negative intensities */
                            CUDAREAL sign=1.0;
                            if(Fbg<0.0) sign=-1.0;

                            /* now we have the structure factor for this pixel */

                            /* polarization factor */
                            CUDAREAL polar = 1.0;
                            if(! nopolar){
                                /* need to compute polarization factor */
                                polar = polarization_factor(polarization,incident,diffracted,polar_vector);
                            }

                            /* accumulate unscaled pixel intensity from this */
                            Ibg += sign*Fbg*Fbg*polar*omega_pixel*capture_fraction*n_source_scale;
                        } /* end of source loop */
                    } /* end of detector thickness loop */
                } /* end of sub-pixel y loop */
            } /* end of sub-pixel x loop */
            /* save photons/pixel (if fluence specified), or F^2/omega if no fluence given */
            floatimage[j] += Ibg*r_e_sqr*fluence*amorphous_molecules/steps;    } // end of pixIdx loop
}


#endif // SIMTBX_GPU_SIMULATION_CUH
