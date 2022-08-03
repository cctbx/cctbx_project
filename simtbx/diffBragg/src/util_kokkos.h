#ifndef SIMTBX_DIFFBRAGG_UTIL_KOKKOS
#define SIMTBX_DIFFBRAGG_UTIL_KOKKOS

#include <vector>
#include "kokkostbx/kokkos_matrix3.h"
// #include "kokkostbx/kokkos_types.h"
#include "kokkostbx/kokkos_vector3.h"
#include "simtbx/diffBragg/src/util.h"

#ifndef CUDAREAL
    #define CUDAREAL double
#endif

using image_type = std::vector<CUDAREAL>;
using VEC3 = kokkostbx::vector3<CUDAREAL>;
using MAT3 = kokkostbx::matrix3<CUDAREAL>;
// using vector_vec3_t = view_1d_t<VEC3>;
// using vector_mat3_t = view_1d_t<MAT3>;

VEC3 to_vec3(const Eigen::Vector3d& v) {
    return VEC3(v[0], v[1], v[2]);
}

MAT3 to_mat3(const Eigen::Matrix3d& m) {
    // Eigen matrix is column-major!
    return MAT3(m(0, 0), m(0, 1), m(0, 2), m(1, 0), m(1, 1), m(1, 2), m(2, 0), m(2, 1), m(2, 2));
}

// CONTAINERS
struct kokkos_crystal {
    MAT3 anisoG;
    std::vector<MAT3> dG_dgamma;
    std::vector<MAT3> dU_dsigma;
    MAT3 anisoU;
    int mosaic_domains;               // number of mosaic domains to model
    CUDAREAL Na, Nb, Nc, Nd, Ne, Nf;  // mosaic domain terms
    CUDAREAL phi0;                    // gonio
    CUDAREAL phistep;
    int phisteps;
    CUDAREAL fudge;       // factor for Bragg peak exponential falloff adjustment
    CUDAREAL spot_scale;  // factor applied to intensity
    int h_range, k_range, l_range;
    int h_max, h_min, k_max, k_min, l_max, l_min;
    CUDAREAL dmin;  // res
    std::vector<CUDAREAL> FhklLinear,
        Fhkl2Linear;              // structure factor amps magnitude (or real, image of complex)
    std::vector<CUDAREAL> fpfdp;  // fprim fdblprime
    std::vector<CUDAREAL> fpfdp_derivs;  // fprime fdblprime deriv
    std::vector<CUDAREAL> atom_data;     // heavy atom data
    std::vector<int> nominal_hkl;        // h,k,l of the pixel (expected)
    CUDAREAL default_F;                  // place holder amplitude
    CUDAREAL r_e_sqr;                    // electron rad

    MAT3 eig_U;        // Umatrix
    MAT3 eig_O;        // O-matrix
    MAT3 eig_B;        // B matrix
    MAT3 RXYZ;         // Rx*Ry*Rz misset perturtbation matrix (this is whats refined)
    VEC3 spindle_vec;  // gonio

    std::vector<MAT3> UMATS_RXYZ;
    std::vector<MAT3> UMATS_RXYZ_prime;
    std::vector<MAT3> UMATS_RXYZ_dbl_prime;

    std::vector<MAT3> RotMats;
    std::vector<MAT3> dRotMats;
    std::vector<MAT3> d2RotMats;
    std::vector<MAT3> UMATS;
    std::vector<MAT3> UMATS_prime;
    std::vector<MAT3> UMATS_dbl_prime;
    std::vector<MAT3> dB_Mats;
    std::vector<MAT3> dB2_Mats;

    kokkos_crystal(crystal T) :
    anisoG(to_mat3(T.anisoG)),
    anisoU(to_mat3(T.anisoU)),
    mosaic_domains(T.mosaic_domains),
    Na(T.Na),
    Nb(T.Nb),
    Nc(T.Nc),
    Nd(T.Nd),
    Ne(T.Ne),
    Nf(T.Nf),
    phi0(T.phi0),
    phistep(T.phistep),
    phisteps(T.phisteps),
    fudge(T.fudge),
    spot_scale(T.spot_scale),
    h_range(T.h_range),
    k_range(T.k_range),
    l_range(T.l_range),
    h_max(T.h_max),
    h_min(T.h_min),
    k_max(T.k_max),
    k_min(T.k_min),
    l_max(T.l_max),
    l_min(T.l_min),
    dmin(T.dmin),
    FhklLinear(T.FhklLinear),
    Fhkl2Linear(T.Fhkl2Linear),
    fpfdp(T.fpfdp),
    fpfdp_derivs(T.fpfdp_derivs),
    atom_data(T.atom_data),
    nominal_hkl(T.nominal_hkl),
    default_F(T.default_F),
    r_e_sqr(T.r_e_sqr),
    eig_U(to_mat3(T.eig_U)),
    eig_O(to_mat3(T.eig_O)),
    eig_B(to_mat3(T.eig_B)),
    RXYZ(to_mat3(T.RXYZ)),
    spindle_vec(to_vec3(T.spindle_vec))    
    {
        for (auto& mat : T.dG_dgamma) {
            dG_dgamma.push_back(to_mat3(mat));
        }
        for (auto& mat : T.dU_dsigma) {
            dU_dsigma.push_back(to_mat3(mat));
        }
        for (auto& mat : T.UMATS_RXYZ) {
            UMATS_RXYZ.push_back(to_mat3(mat));
        }
        for (auto& mat : T.UMATS_RXYZ_prime) {
            UMATS_RXYZ_prime.push_back(to_mat3(mat));
        }
        for (auto& mat : T.UMATS_RXYZ_dbl_prime) {
            UMATS_RXYZ_dbl_prime.push_back(to_mat3(mat));
        }
        for (auto& mat : T.RotMats) {
            RotMats.push_back(to_mat3(mat));
        }    
        for (auto& mat : T.dRotMats) {
            dRotMats.push_back(to_mat3(mat));
        } 
        for (auto& mat : T.d2RotMats) {
            d2RotMats.push_back(to_mat3(mat));
        }     
        for (auto& mat : T.UMATS) {
            UMATS.push_back(to_mat3(mat));
        }  
        for (auto& mat : T.UMATS_prime) {
            UMATS_prime.push_back(to_mat3(mat));
        }  
        for (auto& mat : T.UMATS_dbl_prime) {
            UMATS_dbl_prime.push_back(to_mat3(mat));
        }     
        for (auto& mat : T.dB_Mats) {
            dB_Mats.push_back(to_mat3(mat));
        } 
        for (auto& mat : T.dB2_Mats) {
            dB2_Mats.push_back(to_mat3(mat));
        }                                                              
    }
};

struct kokkos_beam {
    VEC3 polarization_axis;
    CUDAREAL fluence;                          // total fluence
    CUDAREAL kahn_factor;                      // polarization factor
    CUDAREAL *source_X, *source_Y, *source_Z;  // beam vectors
    CUDAREAL* source_lambda;                   // wavelengths
    CUDAREAL* source_I;                        // intensities
    CUDAREAL lambda0, lambda1;                 // affine correction to spectra
    int number_of_sources;                     // number of beams

    kokkos_beam(beam T)
        : fluence(T.fluence),
          kahn_factor(T.kahn_factor),
          source_X(T.source_X),
          source_Y(T.source_Y),
          source_Z(T.source_Z),
          source_lambda(T.source_lambda),
          source_I(T.source_I),
          lambda0(T.lambda0),
          lambda1(T.lambda1),
          number_of_sources(T.number_of_sources),
          polarization_axis(to_vec3(T.polarization_axis)){};
};

struct kokkos_detector {
    std::vector<VEC3> dF_vecs;  // derivative of the panel fast direction
    std::vector<VEC3> dS_vecs;  // derivative of the panel slow direction
    CUDAREAL detector_thickstep, detector_thicksteps, detector_thick, detector_attnlen;
    std::vector<CUDAREAL> close_distances;  // offsets to the detector origins (Z direction)
    int oversample;                         // determines the pixel subsampling rate
    CUDAREAL subpixel_size, pixel_size;
    std::vector<CUDAREAL> fdet_vectors, sdet_vectors, odet_vectors,
        pix0_vectors;  // these define the detector (fast, slow, orth, origin)

    kokkos_detector(detector T)
        : detector_thickstep(T.detector_thickstep),
          detector_thicksteps(T.detector_thicksteps),
          detector_thick(T.detector_thick),
          detector_attnlen(T.detector_attnlen),
          close_distances(T.close_distances),
          oversample(T.oversample),
          subpixel_size(T.subpixel_size),
          pixel_size(T.pixel_size),
          fdet_vectors(T.fdet_vectors),
          sdet_vectors(T.sdet_vectors),
          odet_vectors(T.odet_vectors),
          pix0_vectors(T.pix0_vectors) {
        for (auto& vec : T.dF_vecs) {
            dF_vecs.push_back(to_vec3(vec));
        }
        for (auto& vec : T.dS_vecs) {
            dS_vecs.push_back(to_vec3(vec));
        }
    }
};

#endif
