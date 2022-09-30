#include <cstdio>
#include "diffBraggKOKKOS.h"

// using ::Kokkos::Experimental::exp;
// using ::Kokkos::Experimental::sin;
// using ::Kokkos::Experimental::cos;
// using ::Kokkos::Experimental::sqrt;

void kokkos_sum_over_steps(
    int Npix_to_model,
    vector_uint_t panels_fasts_slows,
    vector_cudareal_t floatimage,
    vector_cudareal_t wavelenimage,
    vector_cudareal_t d_Umat_images,
    vector_cudareal_t d2_Umat_images,
    vector_cudareal_t d_Bmat_images,
    vector_cudareal_t d2_Bmat_images,
    vector_cudareal_t d_Ncells_images,
    vector_cudareal_t d2_Ncells_images,
    vector_cudareal_t d_fcell_images,
    vector_cudareal_t d2_fcell_images,
    vector_cudareal_t d_eta_images,
    vector_cudareal_t d2_eta_images,
    vector_cudareal_t d_lambda_images,
    vector_cudareal_t d2_lambda_images,
    vector_cudareal_t d_panel_rot_images,
    vector_cudareal_t d2_panel_rot_images,
    vector_cudareal_t d_panel_orig_images,
    vector_cudareal_t d2_panel_orig_images,
    vector_cudareal_t d_fp_fdp_images,
    const int Nsteps,
    int printout_fpixel,
    int printout_spixel,
    bool printout,
    CUDAREAL default_F,
    int oversample,
    bool oversample_omega,
    CUDAREAL subpixel_size,
    CUDAREAL pixel_size,
    CUDAREAL detector_thickstep,
    CUDAREAL detector_thick,
    const vector_cudareal_t close_distances,
    CUDAREAL detector_attnlen,
    int detector_thicksteps,
    int sources,
    int phisteps,
    int mosaic_domains,
    bool use_lambda_coefficients,
    CUDAREAL lambda0,
    CUDAREAL lambda1,
    MAT3 eig_U,
    MAT3 eig_O,
    MAT3 eig_B,
    MAT3 RXYZ,
    vector_vec3_t dF_vecs,
    vector_vec3_t dS_vecs,
    const vector_mat3_t UMATS_RXYZ,
    vector_mat3_t UMATS_RXYZ_prime,
    vector_mat3_t UMATS_RXYZ_dbl_prime,
    vector_mat3_t RotMats,
    vector_mat3_t dRotMats,
    vector_mat3_t d2RotMats,
    vector_mat3_t UMATS,
    vector_mat3_t dB_mats,
    vector_mat3_t dB2_mats,
    vector_mat3_t Amatrices,
    const vector_cudareal_t source_X,
    const vector_cudareal_t source_Y,
    const vector_cudareal_t source_Z,
    const vector_cudareal_t source_lambda,
    const vector_cudareal_t source_I,
    CUDAREAL kahn_factor,
    CUDAREAL Na,
    CUDAREAL Nb,
    CUDAREAL Nc,
    CUDAREAL Nd,
    CUDAREAL Ne,
    CUDAREAL Nf,
    CUDAREAL phi0,
    CUDAREAL phistep,
    VEC3 spindle_vec,
    VEC3 polarization_axis,
    int h_range,
    int k_range,
    int l_range,
    int h_max,
    int h_min,
    int k_max,
    int k_min,
    int l_max,
    int l_min,
    CUDAREAL dmin,
    CUDAREAL fudge,
    bool complex_miller,
    int verbose,
    bool only_save_omega_kahn,
    bool isotropic_ncells,
    bool compute_curvatures,
    const vector_cudareal_t FhklLinear,
    const vector_cudareal_t Fhkl2Linear,
    vector_bool_t refine_Bmat,
    vector_bool_t refine_Ncells,
    bool refine_Ncells_def,
    vector_bool_t refine_panel_origin,
    vector_bool_t refine_panel_rot,
    bool refine_fcell,
    vector_bool_t refine_lambda,
    bool refine_eta,
    vector_bool_t refine_Umat,
    const vector_cudareal_t fdet_vectors,
    const vector_cudareal_t sdet_vectors,
    const vector_cudareal_t odet_vectors,
    const vector_cudareal_t pix0_vectors,
    bool nopolar,
    bool point_pixel,
    CUDAREAL fluence,
    CUDAREAL r_e_sqr,
    CUDAREAL spot_scale,
    int Npanels,
    bool aniso_eta,
    bool no_Nabc_scale,
    const vector_cudareal_t fpfdp,
    const vector_cudareal_t fpfdp_derivs,
    const vector_cudareal_t atom_data,
    int num_atoms,
    bool refine_fp_fdp,
    const vector_int_t nominal_hkl,
    bool use_nominal_hkl,
    MAT3 anisoU,
    MAT3 anisoG,
    bool use_diffuse,
    vector_cudareal_t d_diffuse_gamma_images,
    vector_cudareal_t d_diffuse_sigma_images,
    bool refine_diffuse,
    bool gamma_miller_units,
    bool refine_Icell,
    bool save_wavelenimage) {  // BEGIN GPU kernel

    // int tid = blockIdx.x * blockDim.x + threadIdx.x;
    // int thread_stride = blockDim.x * gridDim.x;
    // __shared__ bool s_refine_Icell;
    // __shared__ bool s_use_diffuse;
    // __shared__ bool s_use_nominal_hkl;
    // __shared__ bool s_refine_fp_fdp;
    // __shared__ bool s_complex_miller;
    // __shared__ int s_num_atoms;
    // __shared__ bool s_aniso_eta;
    // __shared__ bool s_no_Nabc_scale;
    // __shared__ bool s_compute_curvatures;
    // __shared__ MAT3 s_Ot;
    // __shared__ bool s_refine_diffuse;
    // __shared__ bool s_gamma_miller_units;
    // __shared__ MAT3 _NABC;
    // __shared__ MAT3 s_dN;
    // __shared__ CUDAREAL C;
    // __shared__ CUDAREAL two_C;
    // __shared__ MAT3 Bmat_realspace;
    // __shared__ MAT3 Amat_init;
    //__shared__ CUDAREAL s_Na;
    //__shared__ CUDAREAL s_Nb;
    //__shared__ CUDAREAL s_Nc;
    // __shared__ CUDAREAL s_NaNbNc_squared; // unused?
    // __shared__ int s_h_max, s_k_max, s_l_max, s_h_min, s_k_min, s_l_min, s_k_range, s_l_range;
    // __shared__ int s_sources, s_mosaic_domains;
    // __shared__ CUDAREAL s_detector_attnlen, s_lambda0, s_lambda1;
    // __shared__ bool s_printout;
    // __shared__ VEC3 s_polarization_axis;

    // __shared__ bool s_refine_Umat[3];
    // __shared__ bool s_refine_panel_origin[3];
    // __shared__ bool s_refine_panel_rot[3];
    // __shared__ bool s_refine_Ncells[3];
    // __shared__ bool s_refine_eta;
    // __shared__ bool s_refine_Ncells_def;
    // __shared__ bool s_refine_fcell;
    // __shared__ bool s_refine_Bmat[6];
    // __shared__ bool s_refine_lambda[2];
    // __shared__ double s_NABC_det, s_NABC_det_sq;
    // extern __shared__ CUDAREAL det_vecs[];
    //__shared__ int det_stride;

    // TODO can we get speed gains by dividing up the following definitions over more threads ?
    // if (threadIdx.x == 0) {
    // for (int i = 0; i < 3; i++) {
    // s_refine_Ncells[i] = refine_Ncells[i];
    // s_refine_Umat[i] = refine_Umat[i];
    // s_refine_panel_origin[i] = refine_panel_origin[i];
    // s_refine_panel_rot[i] = refine_panel_rot[i];
    // }
    // s_refine_Icell = refine_Icell;
    // s_use_nominal_hkl = use_nominal_hkl;
    // s_aniso_eta = aniso_eta;
    // s_no_Nabc_scale = no_Nabc_scale;
    // s_complex_miller = complex_miller;
    // s_refine_lambda[0] = refine_lambda[0];
    // s_refine_lambda[1] = refine_lambda[1];
    // for (int i = 0; i < 6; i++) {
    //     s_refine_Bmat[i] = refine_Bmat[i];
    // }
    // s_use_diffuse = use_diffuse;
    // s_num_atoms = num_atoms;
    // s_refine_fcell = refine_fcell;
    // s_refine_eta = refine_eta;
    // s_refine_Ncells_def = refine_Ncells_def;
    // s_compute_curvatures = compute_curvatures;
    // s_refine_fp_fdp = refine_fp_fdp;
    // s_refine_diffuse = refine_diffuse;
    // s_gamma_miller_units = gamma_miller_units;

    const MAT3 Bmat_realspace = eig_B * 1e10;
    const MAT3 eig_Otranspose = eig_O.transpose();
    const MAT3 Amat_init = eig_U * Bmat_realspace * eig_Otranspose;
    const MAT3 _NABC {Na, Nd, Nf, Nd, Nb, Ne, Nf, Ne, Nc};
    const double NABC_det = _NABC.determinant();  // TODO is this slow ?
    const double NABC_det_sq = NABC_det * NABC_det;
    const CUDAREAL C = 2 / 0.63 * fudge;
    const CUDAREAL two_C = 2 * C;
    // s_Na = Na;
    // s_Nb = Nb;
    // s_Nc = Nc;
    // s_NaNbNc_squared = (Na * Nb * Nc);
    // s_NaNbNc_squared *= s_NaNbNc_squared;
    // s_h_max = h_max;
    // s_k_max = k_max;
    // s_l_max = l_max;
    // s_h_min = h_min;
    // s_k_min = k_min;
    // s_l_min = l_min;
    // s_h_range = h_range;
    // s_k_range = k_range;
    // s_l_range = l_range;

    // s_oversample = oversample;
    // s_detector_thicksteps = detector_thicksteps;
    // s_sources = sources;
    // s_mosaic_domains = mosaic_domains;
    // s_detector_thickstep = detector_thickstep;
    // s_detector_attnlen = detector_attnlen;
    // s_subpixel_size = subpixel_size;
    // s_pixel_size = pixel_size;
    // s_detector_thick = _detector_thick;
    // s_lambda0 = lambda0;
    // s_lambda1 = lambda1;
    // s_oversample_omega = _oversample_omega;
    // s_printout = printout;
    // s_printout_fpixel = printout_fpixel;
    // s_printout_spixel = printout_spixel;
    // s_default_F = default_F;
    // s_verbose = verbose;
    // s_polarization_axis = polarization_axis;
    // s_kahn_factor = kahn_factor;
    // s_nopolar = nopolar;
    // sX0 = source_X[0];
    // sY0 = source_Y[0];
    // sZ0 = source_Z[0];
    // s_Nsteps = Nsteps;

    // }

    const CUDAREAL overall_scale = r_e_sqr * spot_scale * fluence / Nsteps;

    Kokkos::parallel_for(
        "sum_over_steps", Npix_to_model, KOKKOS_LAMBDA(const int& pixIdx) {
            // __syncthreads();

            // for (int pixIdx=tid; pixIdx < Npix_to_model; pixIdx+= thread_stride){
            int _pid = panels_fasts_slows(pixIdx * 3);
            int _fpixel = panels_fasts_slows(pixIdx * 3 + 1);
            int _spixel = panels_fasts_slows(pixIdx * 3 + 2);

            // int fcell_idx=1;
            int nom_h, nom_k, nom_l;
            if (use_nominal_hkl) {
                nom_h = nominal_hkl(pixIdx * 3);
                nom_k = nominal_hkl(pixIdx * 3 + 1);
                nom_l = nominal_hkl(pixIdx * 3 + 2);
            }
            CUDAREAL close_distance = close_distances(_pid);

            // reset photon count for this pixel
            double _I = 0;
            double Ilambda = 0;

            // reset derivative photon counts for the various parameters
            double rot_manager_dI[3] = {0, 0, 0};
            double rot_manager_dI2[3] = {0, 0, 0};
            double ucell_manager_dI[6] = {0, 0, 0, 0, 0, 0};
            double ucell_manager_dI2[6] = {0, 0, 0, 0, 0, 0};
            double Ncells_manager_dI[6] = {0, 0, 0, 0, 0, 0};
            double Ncells_manager_dI2[6] = {0, 0, 0, 0, 0, 0};
            double pan_orig_manager_dI[3] = {0, 0, 0};
            double pan_orig_manager_dI2[3] = {0, 0, 0};
            double pan_rot_manager_dI[3] = {0, 0, 0};
            double pan_rot_manager_dI2[3] = {0, 0, 0};
            double fcell_manager_dI = 0;
            double fcell_manager_dI2 = 0;
            double eta_manager_dI[3] = {0, 0, 0};
            double eta_manager_dI2[3] = {0, 0, 0};
            double lambda_manager_dI[2] = {0, 0};
            double lambda_manager_dI2[2] = {0, 0};
            double fp_fdp_manager_dI[2] = {0, 0};
            double dI_diffuse[6] = {0, 0, 0, 0, 0, 0};

            for (int _subS = 0; _subS < oversample; ++_subS) {
                for (int _subF = 0; _subF < oversample; ++_subF) {
                    // absolute mm position on detector (relative to its origin)
                    CUDAREAL _Fdet =
                        subpixel_size * (_fpixel * oversample + _subF) + subpixel_size / 2.0;
                    CUDAREAL _Sdet =
                        subpixel_size * (_spixel * oversample + _subS) + subpixel_size / 2.0;

                    // assume "distance" is to the front of the detector sensor layer
                    int pid_x = _pid * 3;
                    int pid_y = _pid * 3 + 1;
                    int pid_z = _pid * 3 + 2;

                    CUDAREAL fx = fdet_vectors(pid_x);
                    CUDAREAL fy = fdet_vectors(pid_y);
                    CUDAREAL fz = fdet_vectors(pid_z);
                    CUDAREAL sx = sdet_vectors(pid_x);
                    CUDAREAL sy = sdet_vectors(pid_y);
                    CUDAREAL sz = sdet_vectors(pid_z);
                    CUDAREAL ox = odet_vectors(pid_x);
                    CUDAREAL oy = odet_vectors(pid_y);
                    CUDAREAL oz = odet_vectors(pid_z);
                    CUDAREAL px = pix0_vectors(pid_x);
                    CUDAREAL py = pix0_vectors(pid_y);
                    CUDAREAL pz = pix0_vectors(pid_z);

                    VEC3 _o_vec(ox, oy, oz);

                    for (int _thick_tic = 0; _thick_tic < detector_thicksteps; ++_thick_tic) {
                        CUDAREAL _Odet = _thick_tic * detector_thickstep;

                        CUDAREAL pixposX = _Fdet * fx + _Sdet * sx + _Odet * ox + px;
                        CUDAREAL pixposY = _Fdet * fy + _Sdet * sy + _Odet * oy + py;
                        CUDAREAL pixposZ = _Fdet * fz + _Sdet * sz + _Odet * oz + pz;
                        VEC3 _pixel_pos(pixposX, pixposY, pixposZ);

                        CUDAREAL _airpath = _pixel_pos.length();
                        VEC3 _diffracted = _pixel_pos.get_unit_vector();

                        // solid angle subtended by a pixel: (pix/airpath)^2*cos(2theta)
                        CUDAREAL _omega_pixel = pixel_size * pixel_size / _airpath / _airpath *
                                                close_distance / _airpath;

                        // option to turn off obliquity effect, inverse-square-law only
                        if (point_pixel)
                            _omega_pixel = 1.0 / _airpath / _airpath;

                        // now calculate detector thickness effects
                        CUDAREAL _capture_fraction = 1;

                        if (detector_thick > 0.0 && detector_attnlen > 0.0) {
                            // inverse of effective thickness increase
                            CUDAREAL _parallax = _diffracted.dot(_o_vec);
                            _capture_fraction = ::Kokkos::Experimental::exp(
                                                    -_thick_tic * detector_thickstep /
                                                    detector_attnlen / _parallax) -
                                                ::Kokkos::Experimental::exp(
                                                    -(_thick_tic + 1) * detector_thickstep /
                                                    detector_attnlen / _parallax);
                        }
                        CUDAREAL cap_frac_times_omega = _capture_fraction * _omega_pixel;

                        for (int _source = 0; _source < sources; ++_source) {
                            VEC3 _incident(
                                -source_X(_source), -source_Y(_source), -source_Z(_source));
                            CUDAREAL _lambda = source_lambda(_source);
                            CUDAREAL sI = source_I(_source);
                            CUDAREAL lambda_ang = _lambda * 1e10;
                            if (use_lambda_coefficients) {
                                lambda_ang = lambda0 + lambda1 * lambda_ang;
                                _lambda = lambda_ang * 1e-10;
                            }

                            VEC3 _scattering = (_diffracted - _incident) / _lambda;

                            VEC3 q_vec(_scattering[0], _scattering[1], _scattering[2]);
                            q_vec *= 1e-10;

                            // TODO rename
                            CUDAREAL texture_scale = 1;
                            texture_scale *= cap_frac_times_omega;
                            texture_scale *= sI;

                            for (int _mos_tic = 0; _mos_tic < mosaic_domains; ++_mos_tic) {
                                int amat_idx = _mos_tic;
                                MAT3 UBO = Amatrices(amat_idx);

                                VEC3 H_vec = UBO * q_vec;
                                CUDAREAL _h = H_vec[0];
                                CUDAREAL _k = H_vec[1];
                                CUDAREAL _l = H_vec[2];

                                int _h0 = ceil(_h - 0.5);
                                int _k0 = ceil(_k - 0.5);
                                int _l0 = ceil(_l - 0.5);

                                VEC3 H0(_h0, _k0, _l0);

                                VEC3 delta_H = H_vec - H0;
                                VEC3 V = _NABC * delta_H;
                                CUDAREAL _hrad_sqr = V.dot(V);
                                CUDAREAL exparg = _hrad_sqr * C / 2;
                                CUDAREAL I0 = 0;

                                if (exparg < 35)
                                    if (no_Nabc_scale)
                                        I0 = ::Kokkos::Experimental::exp(-2 * exparg);
                                    else
                                        I0 = (NABC_det_sq) *
                                             ::Kokkos::Experimental::exp(-2 * exparg);

                                // are we doing diffuse scattering
                                // CUDAREAL I_latt_diffuse = 0;
                                double step_diffuse_param[6] = {0, 0, 0, 0, 0, 0};
                                if (use_diffuse) {
                                    MAT3 Amat = UBO;
                                    MAT3 Ainv = UBO.inverse();
                                    MAT3 local_anisoG = anisoG;
                                    if (gamma_miller_units) {
                                        local_anisoG = local_anisoG.dot(Amat);
                                    }
                                    MAT3 Ginv = local_anisoG.inverse();
                                    CUDAREAL anisoG_determ = local_anisoG.determinant();
                                    for (int hh = 0; hh < 1; hh++) {
                                        for (int kk = 0; kk < 1; kk++) {
                                            for (int ll = 0; ll < 1; ll++) {
                                                VEC3 H0_offset(_h0 + hh, _k0 + kk, _l0 + ll);
                                                VEC3 Q0 = Ainv * H0_offset;
                                                CUDAREAL exparg =
                                                    4 * M_PI * M_PI * Q0.dot(anisoU.dot(Q0));
                                                VEC3 delta_H_offset = H_vec - H0_offset;
                                                VEC3 delta_Q = Ainv * delta_H_offset;
                                                VEC3 anisoG_q = local_anisoG * delta_Q;

                                                CUDAREAL V_dot_V = anisoG_q.dot(anisoG_q);
                                                CUDAREAL gamma_portion =
                                                    8. * M_PI * anisoG_determ /
                                                    pow((1. + V_dot_V * 4 * M_PI * M_PI), 2.);

                                                // if (exparg >= 0.5)
                                                //     exparg = 1;

                                                CUDAREAL this_I_latt_diffuse =
                                                    exparg * gamma_portion;

                                                I0 += this_I_latt_diffuse;
                                                if (refine_diffuse) {
                                                    for (int i_gam = 0; i_gam < 3; i_gam++) {
                                                        MAT3 dG_dgam;
                                                        dG_dgam(i_gam, i_gam) = 1;
                                                        if (gamma_miller_units) {
                                                            dG_dgam = dG_dgam.dot(Amat);
                                                        }
                                                        VEC3 dV = dG_dgam * delta_Q;
                                                        CUDAREAL V_dot_dV = anisoG_q.dot(dV);
                                                        CUDAREAL deriv =
                                                            (Ginv.dot(dG_dgam)).trace() -
                                                            16 * M_PI * M_PI * V_dot_dV /
                                                                (1 + 4 * M_PI * M_PI * V_dot_V);
                                                        step_diffuse_param[i_gam] +=
                                                            gamma_portion * deriv * exparg;
                                                    }
                                                    MAT3 dU_dsigma;
                                                    for (int i_sig = 0; i_sig < 3; i_sig++) {
                                                        dU_dsigma(i_sig, i_sig) =
                                                            2. * ::Kokkos::Experimental::sqrt(
                                                                     anisoU(i_sig, i_sig));
                                                        CUDAREAL dexparg =
                                                            4 * M_PI * M_PI *
                                                            Q0.dot(dU_dsigma.dot(Q0));
                                                        dU_dsigma(i_sig, i_sig) = 0.;
                                                        //    if (exparg  >= .5) // only valid up to
                                                        //    a point
                                                        //      dexparg = 0;

                                                        step_diffuse_param[i_sig + 3] +=
                                                            gamma_portion * dexparg;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                // CUDAREAL I_latt_diffuse = 0;
                                // if (use_diffuse){
                                //     MAT3 Ainv = UBO.inverse();
                                //     VEC3 Q0 = Ainv*H0;
                                //     CUDAREAL exparg_diffuse = 4*M_PI*M_PI*Q0.dot(anisoU*Q0);
                                //     //CUDAREAL dwf = exp(-exparg_diffuse);

                                //    VEC3 delta_Q = Ainv*delta_H;
                                //    VEC3 anisoG_q = anisoG*delta_Q;
                                //    I_latt_diffuse = 4.*M_PI*anisoG.determinant() /
                                //            (1.+ anisoG_q.dot(anisoG_q)* 4*M_PI*M_PI);
                                //    if (exparg_diffuse  < .5) // only valid up to a point
                                //        I_latt_diffuse *= (exparg_diffuse);
                                //    I0 += I_latt_diffuse;
                                //}

                                CUDAREAL _F_cell = default_F;
                                CUDAREAL _F_cell2 = 0;

                                if ((_h0 <= h_max) && (_h0 >= h_min) && (_k0 <= k_max) &&
                                    (_k0 >= k_min) && (_l0 <= l_max) && (_l0 >= l_min)) {
                                    int Fhkl_linear_index = (_h0 - h_min) * k_range * l_range +
                                                            (_k0 - k_min) * l_range + (_l0 - l_min);
                                    //_F_cell = __ldg(&FhklLinear[Fhkl_linear_index]);
                                    _F_cell = FhklLinear(Fhkl_linear_index);
                                    // if (complex_miller) _F_cell2 =
                                    // __ldg(&Fhkl2Linear[Fhkl_linear_index]);
                                    if (complex_miller)
                                        _F_cell2 = Fhkl2Linear(Fhkl_linear_index);
                                }

                                CUDAREAL c_deriv_Fcell = 0;
                                CUDAREAL d_deriv_Fcell = 0;
                                if (complex_miller) {
                                    CUDAREAL c_deriv_Fcell_real = 0;
                                    CUDAREAL c_deriv_Fcell_imag = 0;
                                    CUDAREAL d_deriv_Fcell_real = 0;
                                    CUDAREAL d_deriv_Fcell_imag = 0;
                                    if (num_atoms > 0) {
                                        CUDAREAL S_2 = 1.e-20 * (_scattering[0] * _scattering[0] +
                                                                 _scattering[1] * _scattering[1] +
                                                                 _scattering[2] * _scattering[2]);

                                        // fp is always followed by the fdp value
                                        CUDAREAL val_fp = fpfdp(2 * _source);
                                        CUDAREAL val_fdp = fpfdp(2 * _source + 1);

                                        CUDAREAL c_deriv_prime = 0;
                                        CUDAREAL c_deriv_dblprime = 0;
                                        CUDAREAL d_deriv_prime = 0;
                                        CUDAREAL d_deriv_dblprime = 0;
                                        if (refine_fp_fdp) {
                                            //   currently only supports two parameter model
                                            int d_idx = 2 * _source;
                                            c_deriv_prime = fpfdp_derivs(d_idx);
                                            c_deriv_dblprime = fpfdp_derivs(d_idx + 1);
                                            d_deriv_prime = fpfdp_derivs(d_idx + 2 * sources);
                                            d_deriv_dblprime =
                                                fpfdp_derivs(d_idx + 1 + 2 * sources);
                                        }

                                        for (int i_atom = 0; i_atom < num_atoms; i_atom++) {
                                            // fractional atomic coordinates
                                            CUDAREAL atom_x = atom_data(i_atom * 5);
                                            CUDAREAL atom_y = atom_data(i_atom * 5 + 1);
                                            CUDAREAL atom_z = atom_data(i_atom * 5 + 2);
                                            CUDAREAL B = atom_data(i_atom * 5 + 3);  // B factor
                                            B = ::Kokkos::Experimental::exp(
                                                -B * S_2 / 4.0);  // TODO: speed me up?
                                            CUDAREAL occ = atom_data(i_atom * 5 + 4);  // occupancy
                                            CUDAREAL r_dot_h =
                                                _h0 * atom_x + _k0 * atom_y + _l0 * atom_z;
                                            CUDAREAL phase = 2 * M_PI * r_dot_h;
                                            CUDAREAL s_rdoth = ::Kokkos::Experimental::sin(phase);
                                            CUDAREAL c_rdoth = ::Kokkos::Experimental::cos(phase);
                                            CUDAREAL Bocc = B * occ;
                                            CUDAREAL BC = B * c_rdoth;
                                            CUDAREAL BS = B * s_rdoth;
                                            CUDAREAL real_part = BC * val_fp - BS * val_fdp;
                                            CUDAREAL imag_part = BS * val_fp + BC * val_fdp;
                                            _F_cell += real_part;
                                            _F_cell2 += imag_part;
                                            if (refine_fp_fdp) {
                                                c_deriv_Fcell_real +=
                                                    BC * c_deriv_prime - BS * c_deriv_dblprime;
                                                c_deriv_Fcell_imag +=
                                                    BS * c_deriv_prime + BC * c_deriv_dblprime;

                                                d_deriv_Fcell_real +=
                                                    BC * d_deriv_prime - BS * d_deriv_dblprime;
                                                d_deriv_Fcell_imag +=
                                                    BS * d_deriv_prime + BC * d_deriv_dblprime;
                                            }
                                        }
                                    }
                                    CUDAREAL Freal = _F_cell;
                                    CUDAREAL Fimag = _F_cell2;
                                    _F_cell =
                                        ::Kokkos::Experimental::sqrt(Freal * Freal + Fimag * Fimag);
                                    if (refine_fp_fdp) {
                                        c_deriv_Fcell =
                                            Freal * c_deriv_Fcell_real + Fimag * c_deriv_Fcell_imag;
                                        d_deriv_Fcell =
                                            Freal * d_deriv_Fcell_real + Fimag * d_deriv_Fcell_imag;
                                    }
                                }
                                if (!oversample_omega)
                                    _omega_pixel = 1;

                                CUDAREAL _I_cell = _F_cell;
                                if (!refine_Icell)
                                    _I_cell *= _F_cell;
                                CUDAREAL _I_total = _I_cell * I0;
                                CUDAREAL Iincrement = _I_total * texture_scale;
                                _I += Iincrement;
                                if (save_wavelenimage)
                                    Ilambda += Iincrement * lambda_ang;

                                if (refine_diffuse) {
                                    CUDAREAL step_scale = texture_scale * _F_cell * _F_cell;
                                    for (int i_diff = 0; i_diff < 6; i_diff++) {
                                        dI_diffuse[i_diff] +=
                                            step_scale * step_diffuse_param[i_diff];
                                    }
                                }

                                if (refine_fp_fdp) {
                                    CUDAREAL I_noFcell = texture_scale * I0;
                                    fp_fdp_manager_dI[0] += 2 * I_noFcell * (c_deriv_Fcell);
                                    fp_fdp_manager_dI[1] += 2 * I_noFcell * (d_deriv_Fcell);
                                }

                                if (verbose > 3)
                                    printf(
                                        "hkl= %f %f %f  hkl1= %d %d %d  Fcell=%f\n", _h, _k, _l,
                                        _h0, _k0, _l0, _F_cell);

                                MAT3 UBOt;
                                if (refine_Umat(0) || refine_Umat(1) || refine_Umat(2) ||
                                    refine_eta) {
                                    UBOt = Amat_init;
                                }
                                if (refine_Umat(0)) {
                                    MAT3 RyRzUBOt = RotMats(1) * RotMats(2) * UBOt;
                                    VEC3 delta_H_prime =
                                        (UMATS(_mos_tic) * dRotMats(0) * RyRzUBOt)
                                            .transpose()
                                            .dot(q_vec);
                                    CUDAREAL V_dot_dV = V.dot(_NABC.dot(delta_H_prime));
                                    CUDAREAL value = -two_C * V_dot_dV * Iincrement;
                                    CUDAREAL value2 = 0;
                                    if (compute_curvatures) {
                                        VEC3 delta_H_dbl_prime =
                                            (UMATS(_mos_tic).dot(d2RotMats(0).dot(RyRzUBOt)))
                                                .transpose()
                                                .dot(q_vec);
                                        CUDAREAL dV_dot_dV = (_NABC.dot(delta_H_prime))
                                                                 .dot(_NABC.dot(delta_H_prime));
                                        CUDAREAL dV2_dot_V =
                                            (_NABC.dot(delta_H)).dot(_NABC.dot(delta_H_dbl_prime));
                                        value2 =
                                            two_C *
                                            (two_C * V_dot_dV * V_dot_dV - dV2_dot_V - dV_dot_dV) *
                                            Iincrement;
                                    }
                                    rot_manager_dI[0] += value;
                                    rot_manager_dI2[0] += value2;
                                }
                                if (refine_Umat(1)) {
                                    MAT3 UmosRx = UMATS(_mos_tic).dot(RotMats(0));
                                    MAT3 RzUBOt = RotMats(2).dot(UBOt);
                                    VEC3 delta_H_prime = (UmosRx.dot(dRotMats(1).dot(RzUBOt)))
                                                             .transpose()
                                                             .dot(q_vec);
                                    CUDAREAL V_dot_dV = V.dot(_NABC.dot(delta_H_prime));
                                    CUDAREAL value = -two_C * V_dot_dV * Iincrement;

                                    CUDAREAL value2 = 0;
                                    if (compute_curvatures) {
                                        VEC3 delta_H_dbl_prime =
                                            (UmosRx.dot(d2RotMats(1).dot(RzUBOt)))
                                                .transpose()
                                                .dot(q_vec);
                                        CUDAREAL dV_dot_dV = (_NABC.dot(delta_H_prime))
                                                                 .dot(_NABC.dot(delta_H_prime));
                                        CUDAREAL dV2_dot_V =
                                            (_NABC.dot(delta_H)).dot(_NABC.dot(delta_H_dbl_prime));
                                        value2 =
                                            two_C *
                                            (two_C * V_dot_dV * V_dot_dV - dV2_dot_V - dV_dot_dV) *
                                            Iincrement;
                                    }
                                    rot_manager_dI[1] += value;
                                    rot_manager_dI2[1] += value2;
                                }
                                if (refine_Umat(2)) {
                                    MAT3 UmosRxRy = UMATS(_mos_tic).dot(RotMats(0).dot(RotMats(1)));
                                    VEC3 delta_H_prime = (UmosRxRy.dot(dRotMats(2).dot(UBOt)))
                                                             .transpose()
                                                             .dot(q_vec);
                                    CUDAREAL V_dot_dV = V.dot(_NABC.dot(delta_H_prime));
                                    CUDAREAL value = -two_C * V_dot_dV * Iincrement;

                                    CUDAREAL value2 = 0;
                                    if (compute_curvatures) {
                                        VEC3 delta_H_dbl_prime =
                                            (UmosRxRy.dot(d2RotMats(2).dot(UBOt)))
                                                .transpose()
                                                .dot(q_vec);
                                        CUDAREAL dV_dot_dV = (_NABC.dot(delta_H_prime))
                                                                 .dot(_NABC.dot(delta_H_prime));
                                        CUDAREAL dV2_dot_V =
                                            (_NABC.dot(delta_H)).dot(_NABC.dot(delta_H_dbl_prime));
                                        value2 =
                                            two_C *
                                            (two_C * V_dot_dV * V_dot_dV - dV2_dot_V - dV_dot_dV) *
                                            Iincrement;
                                    }
                                    rot_manager_dI[2] += value;
                                    rot_manager_dI2[2] += value2;
                                }
                                // Checkpoint for unit cell derivatives
                                // MAT3 Ot = eig_O.transpose();
                                MAT3 UmosRxRyRzU;
                                VEC3 delta_H_prime;
                                for (int i_uc = 0; i_uc < 6; i_uc++) {
                                    if (refine_Bmat(i_uc)) {
                                        UmosRxRyRzU = UMATS_RXYZ(_mos_tic).dot(eig_U);
                                        delta_H_prime =
                                            (UmosRxRyRzU.dot(dB_mats(i_uc).dot(eig_Otranspose)))
                                                .transpose()
                                                .dot(q_vec);
                                        CUDAREAL V_dot_dV = V.dot(_NABC.dot(delta_H_prime));
                                        CUDAREAL value = -two_C * V_dot_dV * Iincrement;
                                        CUDAREAL value2 = 0;
                                        if (compute_curvatures) {
                                            VEC3 delta_H_dbl_prime =
                                                (UmosRxRyRzU.dot(
                                                     dB2_mats(i_uc).dot(eig_Otranspose)))
                                                    .transpose()
                                                    .dot(q_vec);
                                            CUDAREAL dV_dot_dV = (_NABC.dot(delta_H_prime))
                                                                     .dot(_NABC.dot(delta_H_prime));
                                            CUDAREAL dV2_dot_V =
                                                (_NABC.dot(delta_H))
                                                    .dot(_NABC.dot(delta_H_dbl_prime));
                                            value2 = two_C *
                                                     (two_C * V_dot_dV * V_dot_dV - dV2_dot_V -
                                                      dV_dot_dV) *
                                                     Iincrement;
                                        }
                                        ucell_manager_dI[i_uc] += value;
                                        ucell_manager_dI2[i_uc] += value2;
                                    }
                                }  // end ucell deriv

                                // Checkpoint for Ncells manager
                                if (refine_Ncells(0)) {
                                    int num_ncell_deriv = 1;
                                    if (!isotropic_ncells)
                                        num_ncell_deriv = 3;
                                    for (int i_nc = 0; i_nc < num_ncell_deriv; i_nc++) {
                                        MAT3 dN;
                                        dN(i_nc, i_nc) = 1;
                                        if (num_ncell_deriv == 1) {
                                            dN(0, 0) = 1;
                                            dN(1, 1) = 1;
                                            dN(2, 2) = 1;
                                        }
                                        CUDAREAL N_i = _NABC(i_nc, i_nc);
                                        VEC3 dV_dN = dN.dot(delta_H);
                                        // TODO speedops: precompute these, store shared var
                                        // _NABC.inverse
                                        CUDAREAL determ_deriv = (_NABC.inverse().dot(dN)).trace();
                                        CUDAREAL deriv_coef = determ_deriv - C * (dV_dN.dot(V));
                                        CUDAREAL value = 2 * Iincrement * deriv_coef;
                                        CUDAREAL value2 = 0;
                                        if (compute_curvatures) {
                                            dN(i_nc, i_nc) = 0;  // TODO check maths
                                            value2 = (-1 / N_i / N_i - C * (dV_dN.dot(dV_dN))) * 2 *
                                                     Iincrement;
                                            value2 += deriv_coef * 2 * value;
                                        }
                                        Ncells_manager_dI[i_nc] += value;
                                        Ncells_manager_dI2[i_nc] += value2;
                                    }
                                }  // end Ncells manager deriv

                                if (refine_Ncells_def) {
                                    for (int i_nc = 3; i_nc < 6; i_nc++) {
                                        MAT3 dN;
                                        if (i_nc == 3)
                                            dN = MAT3{0, 1, 0, 1, 0, 0, 0, 0, 0};
                                        else if (i_nc == 4)
                                            dN = MAT3{0, 0, 0, 0, 0, 1, 0, 1, 0};
                                        else
                                            dN = MAT3{0, 0, 1, 0, 0, 0, 1, 0, 0};
                                        VEC3 dV_dN = dN.dot(delta_H);
                                        // TODO speedops: precompute these
                                        CUDAREAL determ_deriv = (_NABC.inverse().dot(dN)).trace();
                                        CUDAREAL deriv_coef = determ_deriv - C * (dV_dN.dot(V));
                                        CUDAREAL value = 2 * Iincrement * deriv_coef;
                                        Ncells_manager_dI[i_nc] += value;
                                        CUDAREAL value2 = 0;
                                        if (compute_curvatures) {
                                            value2 = deriv_coef * value;
                                            value2 += -2 * C * Iincrement * (dV_dN.dot(dV_dN));
                                            Ncells_manager_dI2[i_nc] += value2;
                                        }
                                    }
                                }

                                // Checkpoint for Origin manager
                                for (int i_pan_orig = 0; i_pan_orig < 3; i_pan_orig++) {
                                    if (refine_panel_origin(i_pan_orig)) {
                                        CUDAREAL per_k = 1 / _airpath;
                                        CUDAREAL per_k3 = pow(per_k, 3.);
                                        CUDAREAL per_k5 = pow(per_k, 5.);
                                        CUDAREAL lambda_ang = _lambda * 1e10;

                                        MAT3 M = -two_C * (_NABC.dot(UBO)) / lambda_ang;
                                        VEC3 dk;
                                        if (i_pan_orig == 0)
                                            dk = VEC3{0, 0, 1};
                                        else if (i_pan_orig == 1)
                                            dk = VEC3{1, 0, 0};
                                        else
                                            dk = VEC3{0, 1, 0};

                                        CUDAREAL G = dk.dot(_pixel_pos);
                                        CUDAREAL pix2 = subpixel_size * subpixel_size;
                                        VEC3 dk_hat = -per_k3 * G * _pixel_pos + per_k * dk;
                                        CUDAREAL coef = (M.dot(dk_hat)).dot(V);
                                        CUDAREAL coef2 =
                                            -3 * pix2 * per_k5 * G * (_o_vec.dot(_pixel_pos));
                                        coef2 += pix2 * per_k3 * (_o_vec.dot(dk));
                                        CUDAREAL value =
                                            coef * Iincrement + coef2 * Iincrement / _omega_pixel;

                                        pan_orig_manager_dI[i_pan_orig] += value;
                                        pan_orig_manager_dI2[i_pan_orig] += 0;

                                    }  // end origin manager deriv
                                }

                                for (int i_pan_rot = 0; i_pan_rot < 3; i_pan_rot++) {
                                    if (refine_panel_rot(i_pan_rot)) {
                                        CUDAREAL per_k = 1 / _airpath;
                                        CUDAREAL per_k3 = pow(per_k, 3.);
                                        CUDAREAL per_k5 = pow(per_k, 5.);
                                        CUDAREAL lambda_ang = _lambda * 1e10;
                                        MAT3 M = -two_C * (_NABC.dot(UBO)) / lambda_ang;
                                        VEC3 dk = _Fdet * (dF_vecs(_pid * 3 + i_pan_rot)) +
                                                  _Sdet * (dS_vecs(_pid * 3 + i_pan_rot));
                                        CUDAREAL G = dk.dot(_pixel_pos);
                                        CUDAREAL pix2 = subpixel_size * subpixel_size;
                                        VEC3 dk_hat = -per_k3 * G * _pixel_pos + per_k * dk;
                                        CUDAREAL coef = (M.dot(dk_hat)).dot(V);
                                        CUDAREAL coef2 =
                                            -3 * pix2 * per_k5 * G * (_o_vec.dot(_pixel_pos));
                                        coef2 += pix2 * per_k3 * (_o_vec.dot(dk));
                                        CUDAREAL value =
                                            coef * Iincrement + coef2 * Iincrement / _omega_pixel;

                                        pan_rot_manager_dI[i_pan_rot] += value;
                                        pan_rot_manager_dI2[i_pan_rot] += 0;
                                    }
                                }

                                // checkpoint for Fcell manager
                                if (refine_fcell) {
                                    CUDAREAL value;
                                    if (refine_Icell)
                                        value = I0 * texture_scale;
                                    else
                                        value = 2 * I0 * _F_cell *
                                                texture_scale;  // Iincrement/_F_cell ;
                                    CUDAREAL value2 = 0;
                                    if (compute_curvatures) {
                                        //    NOTE if _Fcell >0
                                        value2 = 2 * I0 * texture_scale;
                                    }
                                    // if (fcell_idx >=0 && fcell_idx <=2){
                                    if (use_nominal_hkl) {
                                        if (_h0 == nom_h && _k0 == nom_k && _l0 == nom_l) {
                                            fcell_manager_dI += value;
                                            fcell_manager_dI2 += value2;
                                        }
                                    } else {
                                        fcell_manager_dI += value;
                                        fcell_manager_dI2 += value2;
                                    }
                                }  // end of fcell man deriv

                                // checkpoint for eta manager
                                if (refine_eta) {
                                    for (int i_eta = 0; i_eta < 3; i_eta++) {
                                        if (i_eta > 0 && !aniso_eta)
                                            continue;
                                        int mtic2 = _mos_tic + i_eta * mosaic_domains;
                                        VEC3 DeltaH_deriv = (UMATS_RXYZ_prime(mtic2).dot(UBOt))
                                                                .transpose()
                                                                .dot(q_vec);
                                        // vector V is _Nabc*Delta_H
                                        VEC3 dV = _NABC.dot(DeltaH_deriv);
                                        CUDAREAL V_dot_dV = V.dot(dV);
                                        CUDAREAL Iprime = -two_C * (V_dot_dV)*Iincrement;
                                        eta_manager_dI[i_eta] += Iprime;
                                        CUDAREAL Idbl_prime = 0;
                                        if (compute_curvatures) {
                                            VEC3 DeltaH_second_deriv =
                                                (UMATS_RXYZ_dbl_prime(mtic2).dot(UBOt))
                                                    .transpose()
                                                    .dot(q_vec);
                                            VEC3 dV2 = _NABC.dot(DeltaH_second_deriv);
                                            Idbl_prime =
                                                -two_C * (dV.dot(dV) + V.dot(dV2)) * Iincrement;
                                            Idbl_prime += -two_C * (V_dot_dV)*Iprime;
                                        }
                                        eta_manager_dI2[i_eta] += Idbl_prime;
                                    }
                                }  // end of eta man deriv

                                // checkpoint for lambda manager
                                for (int i_lam = 0; i_lam < 2; i_lam++) {
                                    if (refine_lambda(i_lam)) {
                                        CUDAREAL lambda_ang = _lambda * 1e10;
                                        CUDAREAL NH_dot_V = (_NABC.dot(H_vec)).dot(V);
                                        CUDAREAL dg_dlambda;
                                        if (i_lam == 0)
                                            dg_dlambda = 1;
                                        else  // i_lam==1
                                            dg_dlambda = lambda_ang;
                                        CUDAREAL coef =
                                            NH_dot_V * two_C * (dg_dlambda) / lambda_ang;
                                        CUDAREAL value = coef * Iincrement;
                                        CUDAREAL value2 = 0;
                                        lambda_manager_dI[i_lam] += value;
                                        lambda_manager_dI2[i_lam] += value2;
                                    }
                                }
                                // end of lambda deriv
                                if (printout) {
                                    if (_subS == 0 && _subF == 0 && _thick_tic == 0 &&
                                        _source == 0 && _mos_tic == 0) {
                                        if ((_fpixel == printout_fpixel &&
                                             _spixel == printout_spixel) ||
                                            printout_fpixel < 0) {
                                            printf(
                                                "%4d %4d :  lambda = %g\n", _fpixel, _spixel,
                                                _lambda);
                                            printf(
                                                "at %g %g %g\n", _pixel_pos[0], _pixel_pos[1],
                                                _pixel_pos[2]);
                                            printf(
                                                "Fdet= %g; Sdet= %g ; Odet= %g\n", _Fdet, _Sdet,
                                                _Odet);
                                            printf(
                                                "PIX0: %f %f %f\n", pix0_vectors(pid_x),
                                                pix0_vectors(pid_y), pix0_vectors(pid_z));
                                            printf(
                                                "F: %f %f %f\n", fdet_vectors(pid_x),
                                                fdet_vectors(pid_y), fdet_vectors(pid_z));
                                            printf(
                                                "S: %f %f %f\n", sdet_vectors(pid_x),
                                                sdet_vectors(pid_y), sdet_vectors(pid_z));
                                            printf(
                                                "O: %f %f %f\n", odet_vectors(pid_x),
                                                odet_vectors(pid_y), odet_vectors(pid_z));
                                            printf(
                                                "pid_x=%d, pid_y=%d; pid_z=%d\n", pid_x, pid_y,
                                                pid_z);

                                            printf(
                                                "QVECTOR: %f %f %f\n", q_vec[0], q_vec[1],
                                                q_vec[2]);
                                            MAT3 UU = UMATS_RXYZ(_mos_tic);
                                            printf(
                                                "UMAT_RXYZ :\n%f  %f  %f\n%f  %f  %f\n%f  %f  %f\n",
                                                UU(0, 0), UU(0, 1), UU(0, 2), UU(1, 0), UU(1, 1),
                                                UU(1, 2), UU(2, 0), UU(2, 1), UU(2, 2));
                                            UU = Bmat_realspace;
                                            printf(
                                                "Bmat_realspace :\n%f  %f  %f\n%f  %f  %f\n%f  %f  %f\n",
                                                UU(0, 0), UU(0, 1), UU(0, 2), UU(1, 0), UU(1, 1),
                                                UU(1, 2), UU(2, 0), UU(2, 1), UU(2, 2));
                                            UU = UBO;
                                            printf(
                                                "UBO :\n%f  %f  %f\n%f  %f  %f\n%f  %f  %f\n",
                                                UU(0, 0), UU(0, 1), UU(0, 2), UU(1, 0), UU(1, 1),
                                                UU(1, 2), UU(2, 0), UU(2, 1), UU(2, 2));

                                            UU = UBOt;
                                            printf(
                                                "UBOt :\n%f  %f  %f\n%f  %f  %f\n%f  %f  %f\n",
                                                UU(0, 0), UU(0, 1), UU(0, 2), UU(1, 0), UU(1, 1),
                                                UU(1, 2), UU(2, 0), UU(2, 1), UU(2, 2));

                                            UU = UmosRxRyRzU;
                                            printf(
                                                "UmosRxRyRzU :\n%f  %f  %f\n%f  %f  %f\n%f  %f  %f\n",
                                                UU(0, 0), UU(0, 1), UU(0, 2), UU(1, 0), UU(1, 1),
                                                UU(1, 2), UU(2, 0), UU(2, 1), UU(2, 2));
                                            VEC3 AA = delta_H_prime;
                                            printf(
                                                "delta_H_prime :\n%f  %f  %f\n", AA[0], AA[1],
                                                AA[2]);
                                            printf("Iincrement: %f\n", Iincrement);
                                            printf(
                                                "hkl= %f %f %f  hkl0= %d %d %d\n", _h, _k, _l, _h0,
                                                _k0, _l0);
                                            printf(
                                                " F_cell=%g  F_cell2=%g I_latt=%g   I = %g\n",
                                                _F_cell, _F_cell2, I0, _I);
                                            printf("I/steps %15.10g\n", _I / Nsteps);
                                            // printf("Ilatt diffuse %15.10g\n", I_latt_diffuse);
                                            printf("omega   %15.10g\n", _omega_pixel);
                                            printf("default_F= %f\n", default_F);
                                            printf(
                                                "Incident[0]=%g, Incident[1]=%g, Incident[2]=%g\n",
                                                _incident[0], _incident[1], _incident[2]);
                                            if (complex_miller)
                                                printf("COMPLEX MILLER!\n");
                                            if (no_Nabc_scale)
                                                printf("No Nabc scale!\n");
                                        }
                                    }
                                }

                            }  // end of mos_tic loop
                        }      // end of source loop
                    }          // end of thick step loop
                }              // end of fpos loop
            }                  // end of spos loop

            CUDAREAL _Fdet_ave = pixel_size * _fpixel + pixel_size / 2.0;
            CUDAREAL _Sdet_ave = pixel_size * _spixel + pixel_size / 2.0;
            CUDAREAL _Odet_ave = 0;  // Odet;
            // TODO maybe make this more general for thick detectors?

            VEC3 _pixel_pos_ave(0, 0, 0);
            int pid_x = _pid * 3;
            int pid_y = _pid * 3 + 1;
            int pid_z = _pid * 3 + 2;

            CUDAREAL fx = fdet_vectors(pid_x);
            CUDAREAL fy = fdet_vectors(pid_y);
            CUDAREAL fz = fdet_vectors(pid_z);

            CUDAREAL sx = sdet_vectors(pid_x);
            CUDAREAL sy = sdet_vectors(pid_y);
            CUDAREAL sz = sdet_vectors(pid_z);

            CUDAREAL ox = odet_vectors(pid_x);
            CUDAREAL oy = odet_vectors(pid_y);
            CUDAREAL oz = odet_vectors(pid_z);

            CUDAREAL px = pix0_vectors(pid_x);
            CUDAREAL py = pix0_vectors(pid_y);
            CUDAREAL pz = pix0_vectors(pid_z);

            _pixel_pos_ave[0] = _Fdet_ave * fx + _Sdet_ave * sx + _Odet_ave * ox + px;
            _pixel_pos_ave[1] = _Fdet_ave * fy + _Sdet_ave * sy + _Odet_ave * oy + py;
            _pixel_pos_ave[2] = _Fdet_ave * fz + _Sdet_ave * sz + _Odet_ave * oz + pz;

            CUDAREAL _airpath_ave = _pixel_pos_ave.length();
            VEC3 _diffracted_ave = _pixel_pos_ave.get_unit_vector();
            CUDAREAL _omega_pixel_ave = pixel_size * pixel_size / _airpath_ave / _airpath_ave *
                                        close_distance / _airpath_ave;

            CUDAREAL _polar = 1;
            if (!nopolar) {
                VEC3 _incident(-source_X(0), -source_Y(0), -source_Z(0));
                _incident.normalize();
                // component of diffracted unit vector along _incident beam unit vector
                CUDAREAL cos2theta = _incident.dot(_diffracted_ave);
                CUDAREAL cos2theta_sqr = cos2theta * cos2theta;
                CUDAREAL sin2theta_sqr = 1 - cos2theta_sqr;

                CUDAREAL _psi = 0;
                if (kahn_factor != 0.0) {
                    // cross product to get "vertical" axis that is orthogonal to the cannonical
                    // "polarization"
                    VEC3 B_in = polarization_axis.cross(_incident);
                    // cross product with _incident beam to get E-vector direction
                    VEC3 E_in = _incident.cross(B_in);
                    // get components of diffracted ray projected onto the E-B plane
                    CUDAREAL _kEi = _diffracted_ave.dot(E_in);
                    CUDAREAL _kBi = _diffracted_ave.dot(B_in);
                    // compute the angle of the diffracted ray projected onto the incident E-B plane
                    _psi = -atan2(_kBi, _kEi);
                }
                // correction for polarized _incident beam
                _polar =
                    0.5 * (1.0 + cos2theta_sqr -
                           kahn_factor * ::Kokkos::Experimental::cos(2 * _psi) * sin2theta_sqr);
            }

            CUDAREAL _om = 1;
            if (!oversample_omega)
                _om = _omega_pixel_ave;
            // final scale term to being everything to photon number units
            CUDAREAL _scale_term = _polar * _om * overall_scale;
            floatimage(pixIdx) = _scale_term * _I;
            if (save_wavelenimage)
                wavelenimage(pixIdx) = Ilambda / _I;

            // udpate the rotation derivative images*
            for (int i_rot = 0; i_rot < 3; i_rot++) {
                if (refine_Umat(i_rot)) {
                    CUDAREAL value = _scale_term * rot_manager_dI[i_rot];
                    CUDAREAL value2 = _scale_term * rot_manager_dI2[i_rot];
                    int idx = i_rot * Npix_to_model + pixIdx;
                    d_Umat_images(idx) = value;
                    d2_Umat_images(idx) = value2;
                }
            }  // end rot deriv image increment

            // update the ucell derivative images
            for (int i_uc = 0; i_uc < 6; i_uc++) {
                if (refine_Bmat(i_uc)) {
                    CUDAREAL value = _scale_term * ucell_manager_dI[i_uc];
                    CUDAREAL value2 = _scale_term * ucell_manager_dI2[i_uc];
                    int idx = i_uc * Npix_to_model + pixIdx;
                    d_Bmat_images(idx) = value;
                    d2_Bmat_images(idx) = value2;
                }
            }  // end ucell deriv image increment

            // update the Ncells derivative image
            if (refine_Ncells(0)) {
                CUDAREAL value = _scale_term * Ncells_manager_dI[0];
                CUDAREAL value2 = _scale_term * Ncells_manager_dI2[0];
                int idx = pixIdx;
                d_Ncells_images(idx) = value;
                d2_Ncells_images(idx) = value2;

                if (!isotropic_ncells) {
                    value = _scale_term * Ncells_manager_dI[1];
                    value2 = _scale_term * Ncells_manager_dI2[1];
                    idx = Npix_to_model + pixIdx;
                    d_Ncells_images(idx) = value;
                    d2_Ncells_images(idx) = value2;

                    value = _scale_term * Ncells_manager_dI[2];
                    value2 = _scale_term * Ncells_manager_dI2[2];
                    idx = Npix_to_model * 2 + pixIdx;
                    d_Ncells_images(idx) = value;
                    d2_Ncells_images(idx) = value2;
                }
            }  // end Ncells deriv image increment
            if (refine_Ncells_def) {
                for (int i_nc = 3; i_nc < 6; i_nc++) {
                    CUDAREAL value = _scale_term * Ncells_manager_dI[i_nc];
                    CUDAREAL value2 = _scale_term * Ncells_manager_dI2[i_nc];
                    int idx = i_nc * Npix_to_model + pixIdx;
                    d_Ncells_images(idx) = value;
                    d2_Ncells_images(idx) = value2;
                }
            }

            // update Fcell derivative image
            if (refine_fcell) {
                CUDAREAL value = _scale_term * fcell_manager_dI;
                CUDAREAL value2 = _scale_term * fcell_manager_dI2;
                d_fcell_images(pixIdx) = value;
                d2_fcell_images(pixIdx) = value2;
            }  // end Fcell deriv image increment

            if (refine_fp_fdp) {
                // c derivative
                CUDAREAL value = _scale_term * fp_fdp_manager_dI[0];
                d_fp_fdp_images(pixIdx) = value;
                // d derivative
                value = _scale_term * fp_fdp_manager_dI[1];
                d_fp_fdp_images(Npix_to_model + pixIdx) = value;
            }
            if (refine_diffuse) {
                for (int i_gam = 0; i_gam < 3; i_gam++) {
                    CUDAREAL val = dI_diffuse[i_gam] * _scale_term;
                    int img_idx = Npix_to_model * i_gam + pixIdx;
                    d_diffuse_gamma_images(img_idx) = val;
                }
                for (int i_sig = 0; i_sig < 3; i_sig++) {
                    CUDAREAL val = dI_diffuse[i_sig + 3] * _scale_term;
                    int img_idx = Npix_to_model * i_sig + pixIdx;
                    d_diffuse_sigma_images(img_idx) = val;
                }
            }

            // update eta derivative image
            if (refine_eta) {
                for (int i_eta = 0; i_eta < 3; i_eta++) {
                    if (i_eta > 0 && !aniso_eta)
                        continue;
                    int idx = pixIdx + Npix_to_model * i_eta;
                    CUDAREAL value = _scale_term * eta_manager_dI[i_eta];
                    CUDAREAL value2 = _scale_term * eta_manager_dI2[i_eta];
                    d_eta_images(idx) = value;
                    d2_eta_images(idx) = value2;
                }
            }  // end eta deriv image increment

            // update the lambda derivative images
            for (int i_lam = 0; i_lam < 2; i_lam++) {
                if (refine_lambda(i_lam)) {
                    CUDAREAL value = _scale_term * lambda_manager_dI[i_lam];
                    CUDAREAL value2 = _scale_term * lambda_manager_dI2[i_lam];
                    int idx = i_lam * Npix_to_model + pixIdx;
                    d_lambda_images(idx) = value;
                    // d2_lambda_images(idx) = value2;
                }
            }  // end lambda deriv image increment

            for (int i_pan_rot = 0; i_pan_rot < 3; i_pan_rot++) {
                if (refine_panel_rot(i_pan_rot)) {
                    CUDAREAL value = _scale_term * pan_rot_manager_dI[i_pan_rot];
                    CUDAREAL value2 = _scale_term * pan_rot_manager_dI2[i_pan_rot];
                    int idx = i_pan_rot * Npix_to_model + pixIdx;
                    d_panel_rot_images(idx) = value;
                    // d2_panel_rot_images(idx) = value2;
                }
            }  // end panel rot deriv image increment

            for (int i_pan_orig = 0; i_pan_orig < 3; i_pan_orig++) {
                if (refine_panel_origin(i_pan_orig)) {
                    CUDAREAL value = _scale_term * pan_orig_manager_dI[i_pan_orig];
                    CUDAREAL value2 = _scale_term * pan_orig_manager_dI2[i_pan_orig];
                    int idx = i_pan_orig * Npix_to_model + pixIdx;
                    d_panel_orig_images(idx) = value;
                    // d2_panel_orig_images(idx) = value2;
                }
            }  // end panel orig deriv image increment
        });    // end pixIdx loop

}  // END of GPU kernel
