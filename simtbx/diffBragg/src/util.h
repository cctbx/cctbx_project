
#ifndef SIMTBX_DIFFBRAGG_UTIL
#define SIMTBX_DIFFBRAGG_UTIL

#include <Eigen/Dense>
#include<Eigen/StdVector>
#include<vector>

typedef std::vector<double> image_type;
typedef Eigen::Matrix<double,3,1> VEC3;
typedef Eigen::Matrix<double,3,3> MAT3;
typedef std::vector<MAT3,Eigen::aligned_allocator<MAT3> > eigMat3_vec;
typedef std::vector<VEC3,Eigen::aligned_allocator<VEC3> > eigVec3_vec;

struct timer_variables{
    double add_spots_pre=0; // times the initializations for add spots kernel
    double add_spots_post=0; // times the copies that occur after add spots kernel
    double add_spots_kernel_wrapper=0; // times the add spots kernel overall, either CPU or GPU
    double cuda_alloc=0; // times the allocation of the device
    double cuda_copy_to_dev=0; // times the copying from host to device
    double cuda_copy_from_dev=0; // times the copying back from device to host
    double cuda_kernel=0; // times the GPU kernel
    int timings=0; // how many times these variables were incremented
    bool recording=true;
  };


// CONTAINERS
struct images{
    image_type wavelength; // image for storing mean wavelength of each pixel
    image_type Umat; // umatrix gradients
    image_type Bmat;  // Bmatrix gradients
    image_type Ncells; // mosaic domain size gradients
    image_type fcell; // structure factor gradients
    image_type eta; // mosaic spread gradients
    image_type lambda; // spectrum affine transform gradients
    image_type panel_rot; // panel rotation gradients
    image_type panel_orig; // panel translation gradients
    image_type fp_fdp;  // fprime and fdblprime gradients
    image_type diffuse_gamma; // diffuse gamma gradients
    image_type diffuse_sigma; // diffuse sigma gradients
};


struct step_arrays{
   int* subS_pos; // stepping through the slow-scan detector axis
   int* subF_pos; // ''   fast-scan ''
   int* thick_pos; // stepping through the detector thickness
   int* source_pos; // stepping through the beam wavelengths
   int* phi_pos;  // stepping through the gonio scan
   int* mos_pos;  // stepping through mosaic blocks (for mosaic spread)
   int Nsteps; // total number of steps
};


struct cuda_flags{
    int device_Id=0;  // gpu device id
    int Npix_to_allocate; // how much space to allocate for simulating forward model and gradients
    // these following flags indicate whether to update quantities on the GPU device prior to running the kernel
    // ( of course they are all set prior to running the kernel for the first time)
    bool update_step_positions;  // step arrays
    bool update_panels_fasts_slows; // pixels to simulatoe (panel id, fast scan, slow scan)
    bool update_sources;  // beam sources
    bool update_umats; // umatrices for mosaic blocks
    bool update_dB_mats; // derivative of the orthogonalization matrix (for unit cell derivatives)
    bool update_rotmats; // rotation matrices (for Umat derivatives)
    bool update_Fhkl; // structure factors
    bool update_detector; // detector vectors (origin, slow-axis, fast-axis, orth-axis)
    bool update_refine_flags;  // refinement flags (in case one is iteratively freezing parameters)
    bool update_panel_deriv_vecs; // if one is refining the detector vectors)
};

struct flags{
    bool wavelength_img=false;
    bool track_Fhkl; // for CPU kernel only, track the HKLS evaluated in the inner most loop
    bool printout; // whether to printout debug info for a pixel
    bool nopolar; // disable polarization effects
    bool point_pixel; // approximate solid angle effects
    bool only_save_omega_kahn; // only save the polarization and solid angle corrections (deprecated)
    bool compute_curvatures; // whether to compute the curvatures in addition to gradients
    bool isotropic_ncells; // one mosaic domain parameter
    bool complex_miller;  // is the miller array complex (such thet Fhkl_linear and Fhkl2_linear are both defined)
    bool no_Nabc_scale; // no Nabc prefactor
    bool refine_diffuse; // flag for computing diffuse gradients
    std::vector<bool> refine_Bmat;  //  Bmatrix
    std::vector<bool> refine_Ncells; // mosaic domain size
    bool refine_Ncells_def; // mosaic domain size off diag
    std::vector<bool> refine_panel_origin; // panel shift
    std::vector<bool> refine_panel_rot; // detector panel rotation
    bool refine_fcell; // structure factor
    std::vector<bool> refine_lambda; // spectrum affine correction
    bool refine_eta; // mosaic spread
    std::vector<bool> refine_Umat; // missetting angle umatrix
    bool refine_fp_fdp; // fprime and fbl prime
    bool use_lambda_coefficients; // affine correction lam0 , lam1
    bool oversample_omega; // omega is computed separately for each sub-pixel
    int printout_fpixel, printout_spixel; // debug printout pixel (fast scan, slow scan) // TODO add panel id
    int verbose; // nanoBragg verbosity flag
    bool use_diffuse = false; // model  diffuse
    bool only_diffuse = false; // model  diffuse scattering (experimental)
    bool refine_Icell = false; // option to refine the structure factor intensity directly (F_cell^2)
                                // The miller array used by nanoBragg/diffBragg is double precision, and hence
                                // allows for negative values. If refine_Icell=True, then the value stored in the data
                                // component of the miller array is assumed to be an intensity, as opposed to an amplitude.
                                // In such cases, the gradient of the scattering w.r.t. this parameter is modified accordingly
                                // such that one could use those gradients as part of a refinement protocol to optimize I_cell

    bool gamma_miller_units = false; // use Miller index units for diffuse gamma matrix
};

struct crystal{
    int laue_group_num=1;
    int stencil_size=0;
    Eigen::Matrix3d anisoG;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > dG_dgamma;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > dU_dsigma;
    Eigen::Matrix3d anisoU;
    int mosaic_domains; // number of mosaic domains to model
    double Na, Nb, Nc, Nd, Ne, Nf; // mosaic domain terms
    double phi0; // gonio
    double phistep;
    double phisteps;
    double fudge; // factor for Bragg peak exponential falloff adjustment
    double spot_scale; // factor applied to intensity
    int h_range, k_range, l_range;
    int h_max, h_min, k_max, k_min, l_max, l_min;
    double dmin; //res
    std::vector<double> FhklLinear, Fhkl2Linear; // structure factor amps magnitude (or real, image of complex)
    std::vector<double> fpfdp; // fprim fdblprime
    std::vector<double> fpfdp_derivs; // fprime fdblprime deriv
    std::vector<double> atom_data; // heavy atom data
    std::vector<int> nominal_hkl; // h,k,l of the pixel (expected)
    double default_F; // place holder amplitude
    double r_e_sqr; // electron rad

    Eigen::Matrix3d eig_U; // Umatrix
    Eigen::Matrix3d eig_O; // O-matrix
    Eigen::Matrix3d eig_B; // B matrix
    Eigen::Matrix3d RXYZ; // Rx*Ry*Rz misset perturtbation matrix (this is whats refined)
    Eigen::Vector3d spindle_vec; // gonio

    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > UMATS_RXYZ;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > UMATS_RXYZ_prime;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > UMATS_RXYZ_dbl_prime;

    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > RotMats;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > dRotMats;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > d2RotMats;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > UMATS;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > UMATS_prime;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > UMATS_dbl_prime;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > dB_Mats;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > dB2_Mats;

};

struct beam{
    Eigen::Vector3d polarization_axis;
    double fluence; // total fluence
    double kahn_factor; // polarization factor
    double *source_X, *source_Y, *source_Z, *source_lambda, *source_I;   // beam vectors, wavelenths, intensities
    double lambda0,lambda1; // affine correction to spectra
    int number_of_sources; // number of beams
};

struct detector{
    std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > dF_vecs; // derivative of the panel fast direction
    std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > dS_vecs; // derivative of the panel slow direction
    double detector_thickstep, detector_thicksteps, detector_thick, detector_attnlen;
    std::vector<double> close_distances; // offsets to the detector origins (Z direction)
    int oversample; // determines the pixel subsampling rate
    double subpixel_size, pixel_size;
    std::vector<double> fdet_vectors, sdet_vectors, odet_vectors, pix0_vectors; // these define the detector (fast, slow, orth, origin)
};

#endif
