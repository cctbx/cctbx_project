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
using KOKKOS_VEC3 = kokkostbx::vector3<CUDAREAL>;
using KOKKOS_MAT3 = kokkostbx::matrix3<CUDAREAL>;
// using vector_vec3_t = view_1d_t<KOKKOS_VEC3>;
// using vector_mat3_t = view_1d_t<KOKKOS_MAT3>;

inline KOKKOS_VEC3 to_vec3(const Eigen::Vector3d& v) {
    return KOKKOS_VEC3(v[0], v[1], v[2]);
}

inline KOKKOS_MAT3 to_mat3(const Eigen::Matrix3d& m) {
    // Eigen matrix is column-major!
    return KOKKOS_MAT3(m(0, 0), m(0, 1), m(0, 2),
                       m(1, 0), m(1, 1), m(1, 2),
                       m(2, 0), m(2, 1), m(2, 2));
}

template <class T, class U>
inline void transfer(T& dst, U& src, int size=-1) {
    const int length = size>=0 ? size : src.size();
    for (int i=0; i<length; ++i) {
        dst(i) = src[i];
    }
}

template <class T, class U>
inline void transfer_KOKKOS_VEC3(T& dst, U& src) {
    auto host_view = Kokkos::create_mirror_view(dst);
    auto src_ptr = src.begin();
    for (int i = 0; i < src.size(); ++i) {
        host_view(i) = to_vec3(src_ptr[i]);
    }
    Kokkos::deep_copy(dst, host_view);
}

template <class T, class U>
inline void transfer_KOKKOS_MAT3(T& dst, U& src) {
    auto host_view = Kokkos::create_mirror_view(dst);
    auto src_ptr = src.begin();
    for (int i = 0; i < src.size(); ++i) {
        host_view(i) = to_mat3(src_ptr[i]);
    }
    Kokkos::deep_copy(dst, host_view);
}

// CONTAINERS
struct kokkos_crystal {
    double Friedel_beta = 1e10; // restraint factor for Friedel pairs
    double Finit_beta = 1e10; // restraint factor for Friedel pairs
    std::vector<int> pos_inds; // indices of the positive Friedel mate
    std::vector<int> neg_inds; // indices of the negative Friedel mate
    double Fhkl_beta=1e10;
    bool use_geometric_mean=false;
    std::unordered_set<int> Fhkl_grad_idx_tracker;
    int num_Fhkl_channels=1;
    int laue_group_num=1;
    int stencil_size=0;
    KOKKOS_MAT3 anisoG;
    std::vector<KOKKOS_MAT3> dG_dgamma;
    std::vector<KOKKOS_MAT3> dU_dsigma;
    KOKKOS_MAT3 anisoU;
    KOKKOS_MAT3 rotate_principal_axes;


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
    std::vector<double> dspace_bins;
    std::vector<CUDAREAL> FhklLinear,
        Fhkl2Linear;              // structure factor amps magnitude (or real, image of complex)
    std::vector<double> ASU_dspace, ASU_Fcell;
    std::vector<int> FhklLinear_ASUid;
    std::unordered_map<std::string, int> ASUid_map;
    int Num_ASU;
    std::string hall_symbol =" P 4nw 2abw";
    std::vector<CUDAREAL> fpfdp;  // fprim fdblprime
    std::vector<CUDAREAL> fpfdp_derivs;  // fprime fdblprime deriv
    std::vector<CUDAREAL> atom_data;     // heavy atom data
    std::vector<int> nominal_hkl;        // h,k,l of the pixel (expected)
    CUDAREAL default_F;                  // place holder amplitude
    CUDAREAL r_e_sqr;                    // electron rad

    KOKKOS_MAT3 eig_U;        // Umatrix
    KOKKOS_MAT3 eig_O;        // O-matrix
    KOKKOS_MAT3 eig_B;        // B matrix
    KOKKOS_MAT3 RXYZ;         // Rx*Ry*Rz misset perturtbation matrix (this is whats refined)
    KOKKOS_VEC3 spindle_vec;  // gonio

    std::vector<KOKKOS_MAT3> UMATS_RXYZ;
    std::vector<KOKKOS_MAT3> UMATS_RXYZ_prime;
    std::vector<KOKKOS_MAT3> UMATS_RXYZ_dbl_prime;

    std::vector<KOKKOS_MAT3> RotMats;
    std::vector<KOKKOS_MAT3> dRotMats;
    std::vector<KOKKOS_MAT3> d2RotMats;
    std::vector<KOKKOS_MAT3> UMATS;
    std::vector<KOKKOS_MAT3> UMATS_prime;
    std::vector<KOKKOS_MAT3> UMATS_dbl_prime;
    std::vector<KOKKOS_MAT3> dB_Mats;
    std::vector<KOKKOS_MAT3> dB2_Mats;

    kokkos_crystal(crystal T)
        : Friedel_beta(T.Friedel_beta),
          Finit_beta(T.Finit_beta),
          pos_inds(T.pos_inds),
          neg_inds(T.neg_inds),
          Fhkl_beta(T.Fhkl_beta),
          use_geometric_mean(T.use_geometric_mean),
          Fhkl_grad_idx_tracker(T.Fhkl_grad_idx_tracker),
          num_Fhkl_channels(T.num_Fhkl_channels),
          laue_group_num(T.laue_group_num),
          stencil_size(T.stencil_size),
          anisoG(to_mat3(T.anisoG)),
          anisoU(to_mat3(T.anisoU)),
          rotate_principal_axes(to_mat3(T.rotate_principal_axes)),
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
          dspace_bins(T.dspace_bins),
          FhklLinear(T.FhklLinear),
          Fhkl2Linear(T.Fhkl2Linear),
          ASU_dspace(T.ASU_dspace),
          ASU_Fcell(T.ASU_Fcell),
          FhklLinear_ASUid(T.FhklLinear_ASUid),
          ASUid_map(T.ASUid_map),
          Num_ASU(T.Num_ASU),
          hall_symbol(T.hall_symbol),
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
          spindle_vec(to_vec3(T.spindle_vec)) {
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
    KOKKOS_VEC3 polarization_axis;
    std::vector<int> Fhkl_channels; // if refining scale factors for wavelength dependent structure factor intensities
    CUDAREAL fluence;                          // total fluence
    CUDAREAL kahn_factor;                      // polarization factor
    CUDAREAL *source_X, *source_Y, *source_Z;  // beam vectors
    CUDAREAL* source_lambda;                   // wavelengths
    CUDAREAL* source_I;                        // intensities
    CUDAREAL lambda0, lambda1;                 // affine correction to spectra
    int number_of_sources;                     // number of beams

    kokkos_beam(beam T)
        : polarization_axis(to_vec3(T.polarization_axis)),
          Fhkl_channels(T.Fhkl_channels),
          fluence(T.fluence),
          kahn_factor(T.kahn_factor),
          source_X(T.source_X),
          source_Y(T.source_Y),
          source_Z(T.source_Z),
          source_lambda(T.source_lambda),
          source_I(T.source_I),
          lambda0(T.lambda0),
          lambda1(T.lambda1),
          number_of_sources(T.number_of_sources){ };
};

struct kokkos_detector {
    std::vector<KOKKOS_VEC3> dF_vecs;  // derivative of the panel fast direction
    std::vector<KOKKOS_VEC3> dS_vecs;  // derivative of the panel slow direction
    CUDAREAL detector_thickstep;
    int detector_thicksteps;
    CUDAREAL detector_thick;
    CUDAREAL detector_attnlen;
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

struct kokkos_manager {
    KOKKOS_VEC3 rot;
    double ucell[6] = {0, 0, 0, 0, 0, 0};
    double Ncells[6] = {0, 0, 0, 0, 0, 0};
    KOKKOS_VEC3 pan_orig;
    KOKKOS_VEC3 pan_rot;
    double fcell = 0;
    KOKKOS_VEC3 eta;
    double lambda[2] = {0, 0};
    double fp_fdp[2] = {0, 0};
    double diffuse[6] = {0, 0, 0, 0, 0, 0};

    KOKKOS_INLINE_FUNCTION void reset() {
        for (int i=0; i<6; ++i) {
            ucell[i] = 0;
            Ncells[i] = 0;
            diffuse[i] = 0;
        }
        for (int i=0; i<2; ++i) {
            lambda[i] = 0;
            fp_fdp[i] = 0;
        }
        rot.zero();
        pan_orig.zero();
        pan_rot.zero();
        eta.zero();
        fcell = 0;
    }
};

#endif
