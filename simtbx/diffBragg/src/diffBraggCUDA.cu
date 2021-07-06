#include <sys/time.h>
#include "diffBraggCUDA.h"
#include "diffBragg_gpu_kernel.h"
#include <stdio.h>
//lkalskdlaksdlkalsd

//#define BLOCKSIZE 128
//#define NUMBLOCKS 128
//https://stackoverflow.com/a/14038590/2077270
#define gpuErr(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

void error_msg(cudaError_t err, char* msg){
    if (err != cudaSuccess){
        printf("%s: CUDA error message: %s\n", msg, cudaGetErrorString(err));
        exit(err);
    }
}

void diffBragg_loopy(
        int Npix_to_model, std::vector<unsigned int>& panels_fasts_slows,
        image_type& floatimage,
        image_type& d_Umat_images, image_type& d2_Umat_images,
        image_type& d_Bmat_images, image_type& d2_Bmat_images,
        image_type& d_Ncells_images, image_type& d2_Ncells_images,
        image_type& d_fcell_images, image_type& d2_fcell_images,
        image_type& d_eta_images, image_type& d2_eta_images,
        image_type& d_lambda_images, image_type& d2_lambda_images,
        image_type& d_panel_rot_images, image_type& d2_panel_rot_images,
        image_type& d_panel_orig_images, image_type& d2_panel_orig_images,
        image_type& d_sausage_XYZ_scale_images,
        image_type& d_fp_fdp_images,
        const int Nsteps, int _printout_fpixel, int _printout_spixel, bool _printout, CUDAREAL _default_F,
        int oversample, bool _oversample_omega, CUDAREAL subpixel_size, CUDAREAL pixel_size,
        CUDAREAL detector_thickstep, CUDAREAL _detector_thick, std::vector<CUDAREAL>& close_distances, CUDAREAL detector_attnlen,
        bool use_lambda_coefficients, CUDAREAL lambda0, CUDAREAL lambda1,
        MAT3& eig_U, MAT3& eig_O, MAT3& eig_B, MAT3& RXYZ,
        std::vector<VEC3,Eigen::aligned_allocator<VEC3> >& dF_vecs,
        std::vector<VEC3,Eigen::aligned_allocator<VEC3> >& dS_vecs,
        std::vector<MAT3,Eigen::aligned_allocator<MAT3> >& UMATS_RXYZ,
        std::vector<MAT3,Eigen::aligned_allocator<MAT3> >& UMATS_RXYZ_prime,
        std::vector<MAT3,Eigen::aligned_allocator<MAT3> >& UMATS_RXYZ_dbl_prime,
        std::vector<MAT3,Eigen::aligned_allocator<MAT3> >& RotMats,
        std::vector<MAT3,Eigen::aligned_allocator<MAT3> >& dRotMats,
        std::vector<MAT3,Eigen::aligned_allocator<MAT3> >& d2RotMats,
        std::vector<MAT3,Eigen::aligned_allocator<MAT3> >& UMATS,
        std::vector<MAT3,Eigen::aligned_allocator<MAT3> >& dB_Mats,
        std::vector<MAT3,Eigen::aligned_allocator<MAT3> >& dB2_Mats,
        std::vector<MAT3,Eigen::aligned_allocator<MAT3> >&sausages_RXYZ,
        std::vector<MAT3,Eigen::aligned_allocator<MAT3> >& d_sausages_RXYZ,
        std::vector<MAT3,Eigen::aligned_allocator<MAT3> >& sausages_U,
        image_type& sausages_scale, // TODO adjust sausages_scale type
        CUDAREAL* source_X, CUDAREAL* source_Y, CUDAREAL* source_Z, CUDAREAL* source_lambda, CUDAREAL* source_I,
        CUDAREAL kahn_factor,
        CUDAREAL Na, CUDAREAL Nb, CUDAREAL Nc,
        CUDAREAL Nd, CUDAREAL Ne, CUDAREAL Nf,
        CUDAREAL phi0, CUDAREAL phistep,
        VEC3& spindle_vec, VEC3 _polarization_axis,
        int h_range, int k_range, int l_range,
        int h_max, int h_min, int k_max, int k_min, int l_max, int l_min, CUDAREAL dmin,
        CUDAREAL fudge, bool complex_miller, int verbose, bool only_save_omega_kahn,
        bool isotropic_ncells, bool compute_curvatures,
        std::vector<CUDAREAL>& _FhklLinear, std::vector<CUDAREAL>& _Fhkl2Linear,
        std::vector<bool>& refine_Bmat, std::vector<bool>& refine_Ncells, bool refine_Ncells_def, std::vector<bool>& refine_panel_origin,
        std::vector<bool>& refine_panel_rot,
        bool refine_fcell, std::vector<bool>& refine_lambda, bool refine_eta, std::vector<bool>& refine_Umat,
        bool refine_sausages, int num_sausages,
        bool refine_fp_fdp,
        std::vector<CUDAREAL>& fdet_vectors, std::vector<CUDAREAL>& sdet_vectors,
        std::vector<CUDAREAL>& odet_vectors, std::vector<CUDAREAL>& pix0_vectors,
        bool _nopolar, bool _point_pixel, CUDAREAL _fluence, CUDAREAL _r_e_sqr, CUDAREAL _spot_scale,
        int number_of_sources, int device_Id,
        diffBragg_cudaPointers& cp,
        bool update_step_positions, bool update_panels_fasts_slows, bool update_sources, bool update_umats,
        bool update_dB_mats, bool update_rotmats, bool update_Fhkl, bool update_detector, bool update_refine_flags ,
        bool update_panel_deriv_vecs, bool update_sausages_on_device, int detector_thicksteps, int phisteps,
        int Npix_to_allocate, bool no_Nabc_scale,
        std::vector<CUDAREAL>& fpfdp,
        std::vector<CUDAREAL>& fpfdp_derivs,
        std::vector<CUDAREAL>&atom_data,
        std::vector<int>&nominal_hkl){ // diffBragg cuda loopy


    if (phi0 != 0 || phisteps > 1){
        printf("PHI (goniometer position) not supported in GPU code: phi0=%f phisteps=%d phistep=%f\n", phi0, phisteps, phistep);
        exit(-1);
    }

    int numblocks;
    int blocksize;
    char* diffBragg_blocks = getenv("DIFFBRAGG_NUM_BLOCKS");
    char* diffBragg_threads = getenv("DIFFBRAGG_THREADS_PER_BLOCK");
    if (diffBragg_threads==NULL)
        blocksize=128;
    else
        blocksize=atoi(diffBragg_threads);

    if (diffBragg_blocks==NULL)
        numblocks = (Npix_to_model+blocksize-1)/blocksize;
    else
        numblocks = atoi(diffBragg_blocks);

    int cuda_devices;
    cudaGetDeviceCount(&cuda_devices);

    if (num_sausages > 6){
        printf("Too many sausages! Should be less than 6 to run on GPU\n");
        exit(-1);
    }

    error_msg(cudaGetLastError(), "after device count");
    if (verbose > 1)
        printf("Found %d CUDA-capable devices\n", cuda_devices);

    //if (device_Id <= cuda_devices)
    gpuErr(cudaSetDevice(device_Id));

    double time;
    struct timeval t1, t2, t3 ,t4;
    gettimeofday(&t1, 0);

//  determine if we need to allocate pixels, and how many.
//  For best usage, one should use the diffBragg property (visible from Python) Npix_to_allocate
//  in order to just allocate to the GPU - this is useful for ensemble refinement, where each shot
//  can have a variable number of pixels being modeled, and ony only needs to allocate the device once
//  (with the largest expected number of pixels for a given shot)
   // TODO clean up this logic a bit
    if (cp.device_is_allocated && (cp.npix_allocated < Npix_to_model)){
        printf("Need to re-allocate pixels, currently have %d allocated, but trying to model %d\n",
            cp.npix_allocated, Npix_to_model);
        exit(-1);
    }
    else if (Npix_to_allocate==-1){
        Npix_to_allocate = Npix_to_model;
    }
    else if (Npix_to_model > Npix_to_allocate){
        printf("Npix to model=%d is greater than the number of pixel requested for allocation (%d)!\n",
            Npix_to_model, Npix_to_allocate);
        exit(-1);
    }

    if(cp.device_is_allocated){
        if (verbose){
           printf("Will model %d pixels (GPU has %d pre-allocated pix)\n", Npix_to_model, cp.npix_allocated);
        }
    }
    else{
        if (verbose){
           printf("Will model %d pixels and allocate %d pix\n", Npix_to_model, Npix_to_allocate);
        }
        gpuErr(cudaMallocManaged(&cp.cu_source_X, number_of_sources*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_source_Y, number_of_sources*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_source_Z, number_of_sources*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_source_I, number_of_sources*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_source_lambda, number_of_sources*sizeof(CUDAREAL)));

        gpuErr(cudaMallocManaged((void **)&cp.cu_UMATS, UMATS.size()*sizeof(MAT3)));
        gpuErr(cudaMallocManaged((void **)&cp.cu_UMATS_RXYZ, UMATS_RXYZ.size()*sizeof(MAT3)));
        gpuErr(cudaMallocManaged((void **)&cp.cu_AMATS, UMATS_RXYZ.size()*sausages_U.size()*sizeof(MAT3)));
        if (UMATS_RXYZ_prime.size()>0)
            gpuErr(cudaMallocManaged((void **)&cp.cu_UMATS_RXYZ_prime, UMATS_RXYZ_prime.size()*sizeof(MAT3)));
        if (UMATS_RXYZ_dbl_prime.size()>0)
            gpuErr(cudaMallocManaged((void **)&cp.cu_UMATS_RXYZ_dbl_prime, UMATS_RXYZ_dbl_prime.size()*sizeof(MAT3)));

        gpuErr(cudaMallocManaged((void **)&cp.cu_dB_Mats, dB_Mats.size()*sizeof(MAT3)));
        gpuErr(cudaMallocManaged((void **)&cp.cu_dB2_Mats, dB2_Mats.size()*sizeof(MAT3)));

        gpuErr(cudaMallocManaged((void **)&cp.cu_RotMats, RotMats.size()*sizeof(MAT3)));
        gpuErr(cudaMallocManaged((void **)&cp.cu_dRotMats, dRotMats.size()*sizeof(MAT3)));
        gpuErr(cudaMallocManaged((void **)&cp.cu_d2RotMats, d2RotMats.size()*sizeof(MAT3)));

        gpuErr(cudaMallocManaged(&cp.cu_fdet_vectors, fdet_vectors.size()*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_sdet_vectors, fdet_vectors.size()*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_odet_vectors, fdet_vectors.size()*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_pix0_vectors, fdet_vectors.size()*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_close_distances, close_distances.size()*sizeof(CUDAREAL)));

        if (fpfdp.size() > 0){
            gpuErr(cudaMallocManaged(&cp.cu_fpfdp, fpfdp.size()*sizeof(CUDAREAL)));
            gpuErr(cudaMallocManaged(&cp.cu_atom_data, atom_data.size()*sizeof(CUDAREAL)));
        }
        if(fpfdp_derivs.size() > 0)
            gpuErr(cudaMallocManaged(&cp.cu_fpfdp_derivs, fpfdp_derivs.size()*sizeof(CUDAREAL)));

        gpuErr(cudaMallocManaged(&cp.cu_refine_Bmat, 6*sizeof(bool)));
        gpuErr(cudaMallocManaged(&cp.cu_refine_Umat, 3*sizeof(bool)));
        gpuErr(cudaMallocManaged(&cp.cu_refine_Ncells, 3*sizeof(bool)));
        gpuErr(cudaMallocManaged(&cp.cu_refine_panel_origin, 3*sizeof(bool)));
        gpuErr(cudaMallocManaged(&cp.cu_refine_panel_rot, 3*sizeof(bool)));
        gpuErr(cudaMallocManaged(&cp.cu_refine_lambda, 2*sizeof(bool)));

        gpuErr(cudaMallocManaged(&cp.cu_Fhkl, _FhklLinear.size()*sizeof(CUDAREAL)));
        if (complex_miller)
            gpuErr(cudaMallocManaged(&cp.cu_Fhkl2, _FhklLinear.size()*sizeof(CUDAREAL)));

        gpuErr(cudaMallocManaged((void **)&cp.cu_dF_vecs, dF_vecs.size()*sizeof(VEC3)));
        gpuErr(cudaMallocManaged((void **)&cp.cu_dS_vecs, dF_vecs.size()*sizeof(VEC3)));

        gpuErr(cudaMallocManaged( (void**)&cp.cu_sausages_RXYZ, sausages_RXYZ.size()*sizeof(MAT3) ));
        gpuErr(cudaMallocManaged( (void**)&cp.cu_d_sausages_RXYZ, d_sausages_RXYZ.size()*sizeof(MAT3) ));
        gpuErr(cudaMallocManaged( (void**)&cp.cu_sausages_U, sausages_U.size()*sizeof(MAT3) ));
        gpuErr(cudaMallocManaged( &cp.cu_sausages_scale, sausages_scale.size()*sizeof(CUDAREAL) ));

        //gettimeofday(&t3, 0));
        gpuErr(cudaMallocManaged(&cp.cu_floatimage, Npix_to_allocate*sizeof(CUDAREAL) ));
        gpuErr(cudaMallocManaged(&cp.cu_d_fcell_images, Npix_to_allocate*3*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_d_eta_images, Npix_to_allocate*3*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_d2_eta_images, Npix_to_allocate*3*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_d_Umat_images, Npix_to_allocate*3*sizeof(CUDAREAL) ));
        gpuErr(cudaMallocManaged(&cp.cu_d_Ncells_images, Npix_to_allocate*6*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_d_panel_rot_images, Npix_to_allocate*3*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_d_panel_orig_images, Npix_to_allocate*3*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_d_lambda_images, Npix_to_allocate*2*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_d_Bmat_images, Npix_to_allocate*6*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_d_sausage_XYZ_scale_images, num_sausages*Npix_to_allocate*4*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_d_fp_fdp_images, Npix_to_allocate*2*sizeof(CUDAREAL)));

        // allocate curvatures
        //gpuErr(cudaMallocManaged(&cp.cu_d_eta_images, Npix_to_allocate*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_d2_Umat_images, Npix_to_allocate*3*sizeof(CUDAREAL) ));
        gpuErr(cudaMallocManaged(&cp.cu_d2_Ncells_images, Npix_to_allocate*6*sizeof(CUDAREAL)));
        //gpuErr(cudaMallocManaged(&cp.cu_d_panel_rot_images, Npix_to_allocate*3*sizeof(CUDAREAL)));
        //gpuErr(cudaMallocManaged(&cp.cu_d_panel_orig_images, Npix_to_allocate*3*sizeof(CUDAREAL)));
        //gpuErr(cudaMallocManaged(&cp.cu_d_lambda_images, Npix_to_allocate*2*sizeof(CUDAREAL)));
        gpuErr(cudaMallocManaged(&cp.cu_d2_Bmat_images, Npix_to_allocate*6*sizeof(CUDAREAL)));
        if(nominal_hkl.size() >0)
            gpuErr(cudaMallocManaged(&cp.cu_nominal_hkl, Npix_to_allocate*3*sizeof(int)));
        //gpuErr(cudaMallocManaged(&cp.cu_d_sausage_XYZ_scale_images, num_sausages*Npix_to_allocate*4*sizeof(CUDAREAL)));

        //gettimeofday(&t4, 0);
        //time = (1000000.0*(t4.tv_sec-t3.tv_sec) + t4.tv_usec-t3.tv_usec)/1000.0;
        //printf("TIME SPENT ALLOCATING (IMAGES ONLY):  %3.10f ms \n", time);
        gpuErr(cudaMallocManaged(&cp.cu_panels_fasts_slows, Npix_to_allocate*3*sizeof(panels_fasts_slows[0])));
        cp.npix_allocated = Npix_to_allocate;
    } // END of allocation

    bool ALLOC = !cp.device_is_allocated; // shortcut variable

    gettimeofday(&t2, 0);
    time = (1000000.0*(t2.tv_sec-t1.tv_sec) + t2.tv_usec-t1.tv_usec)/1000.0;
    if(verbose>1)
        printf("TIME SPENT ALLOCATING (TOTAL):  %3.10f ms \n", time);

    //ALLOC = false;
//  BEGIN COPYING DATA
    gettimeofday(&t1, 0);
    bool FORCE_COPY=true;

//  END step position
    int kaladin_stormblessed = 777;


//  BEGIN sources
    if (update_sources || ALLOC || FORCE_COPY){
        for (int i=0; i< number_of_sources; i++){
            VEC3 incident(source_X[i], source_Y[i], source_Z[i]);
            incident /= incident.norm();
            cp.cu_source_X[i] = incident[0];
            cp.cu_source_Y[i] = incident[1];
            cp.cu_source_Z[i] = incident[2];
            cp.cu_source_I[i] = source_I[i];
            cp.cu_source_lambda[i] = source_lambda[i];
        }
        if(verbose>1 )
          printf("H2D sources\n");
    }
//  END sources


//  UMATS
    if (update_umats || ALLOC||FORCE_COPY){
        for (int i=0; i< UMATS.size(); i++)
            cp.cu_UMATS[i] = UMATS[i];
        //int idx=0;
        //for (int i=0; i < UMATS_RXYZ.size(); i++){
        //    for (auto elem: UMATS_RXYZ[i].reshaped()){
        //        cp.cu_UMATS_RXYZ[idx] = elem;
        //        idx ++;
        //    }
        //}
        for (int i=0; i < UMATS_RXYZ.size(); i++)
            cp.cu_UMATS_RXYZ[i] = UMATS_RXYZ[i];
        for (int i=0; i < UMATS_RXYZ_prime.size(); i++)
            cp.cu_UMATS_RXYZ_prime[i] = UMATS_RXYZ_prime[i];
        for (int i=0; i < UMATS_RXYZ_dbl_prime.size(); i++)
            cp.cu_UMATS_RXYZ_dbl_prime[i] = UMATS_RXYZ_dbl_prime[i];
        if(verbose>1)
            printf("H2D Done copying Umats\n") ;
    }
//  END UMATS


    if (update_umats || update_sausages_on_device|| ALLOC||FORCE_COPY){
        MAT3 Amat_init = eig_U*eig_B*1e10*(eig_O.transpose());
        for (int i_sausage=0; i_sausage< sausages_U.size(); i_sausage++){
            for(int i_mos =0; i_mos< UMATS_RXYZ.size(); i_mos++){
                int idx = UMATS_RXYZ.size()*i_sausage + i_mos;
                cp.cu_AMATS[idx] = (UMATS_RXYZ[i_mos]*sausages_U[i_sausage]*Amat_init).transpose();
            }
        }
        if(verbose>1)
            printf("H2D Done copying Amats\n") ;
    }


//  BMATS
    if(update_dB_mats || ALLOC || FORCE_COPY){
        for (int i=0; i< dB_Mats.size(); i++)
            cp.cu_dB_Mats[i] = dB_Mats[i];
        for (int i=0; i< dB2_Mats.size(); i++)
            cp.cu_dB2_Mats[i] = dB2_Mats[i];
        if(verbose>1)
            printf("H2D Done copying dB_Mats\n") ;
    }
//  END BMATS


//  ROT MATS
    if(update_rotmats || ALLOC || FORCE_COPY){
        for (int i=0; i<RotMats.size(); i++)
            cp.cu_RotMats[i] = RotMats[i];
        for (int i=0; i<dRotMats.size(); i++)
            cp.cu_dRotMats[i] = dRotMats[i];
        for (int i=0; i<d2RotMats.size(); i++)
            cp.cu_d2RotMats[i] = d2RotMats[i];
        if (verbose>1)
          printf("H2D Done copying rotmats\n");
    }
//  END ROT MATS

//  sausages
    if(update_sausages_on_device || ALLOC || FORCE_COPY){
        for (int i=0; i<sausages_RXYZ.size(); i++)
            cp.cu_sausages_RXYZ[i] = sausages_RXYZ[i];
        for (int i=0; i<sausages_U.size(); i++)
            cp.cu_sausages_U[i] = sausages_U[i];
        for (int i=0; i<d_sausages_RXYZ.size(); i++)
            cp.cu_d_sausages_RXYZ[i] = d_sausages_RXYZ[i];
        for (int i=0; i< sausages_scale.size(); i++)
            cp.cu_sausages_scale[i] = sausages_scale[i];
        if (verbose>1)
          printf("H2D Done copying sausages\n");
    }
//  END ROT MATS


//  DETECTOR VECTORS
    if (update_detector || ALLOC || FORCE_COPY){
        for (int i=0; i<fdet_vectors.size(); i++){
            cp.cu_fdet_vectors[i] = fdet_vectors[i];
            cp.cu_sdet_vectors[i] = sdet_vectors[i];
            cp.cu_odet_vectors[i] = odet_vectors[i];
            cp.cu_pix0_vectors[i] = pix0_vectors[i];
        }
        for(int i=0; i < close_distances.size();i++){
            cp.cu_close_distances[i] = close_distances[i];
        }
        if (verbose>1)
          printf("H2D Done copying detector vectors\n");
    }
//  END  DETECTOR VECTORS

    if ( ALLOC || FORCE_COPY){
      for(int i=0; i< nominal_hkl.size(); i++){
        cp.cu_nominal_hkl[i] = nominal_hkl[i];
      }
      for (int i=0; i< atom_data.size(); i++){
        cp.cu_atom_data[i] = atom_data[i];
      }
      if (verbose>1)
        printf("H2D Done copying atom data\n");
      for(int i=0; i< fpfdp.size(); i++){
        cp.cu_fpfdp[i] = fpfdp[i];
      }
      for(int i=0; i< fpfdp_derivs.size(); i++){
        cp.cu_fpfdp_derivs[i] = fpfdp_derivs[i];
      }
      if (verbose>1)
        printf("H2D Done copying fprime and fdblprime\n");
    }


//  BEGIN REFINEMENT FLAGS
    if (update_refine_flags || ALLOC || FORCE_COPY){
        for (int i=0; i<3; i++){
            cp.cu_refine_Umat[i] = refine_Umat[i];
            cp.cu_refine_Ncells[i] = refine_Ncells[i];
            cp.cu_refine_panel_origin[i] = refine_panel_origin[i];
            cp.cu_refine_panel_rot[i] = refine_panel_rot[i];
        }
        for(int i=0; i<2; i++)
            cp.cu_refine_lambda[i] = refine_lambda[i];
        for(int i=0; i<6; i++)
            cp.cu_refine_Bmat[i] = refine_Bmat[i];
        if (verbose>1)
          printf("H2D Done copying refinement flags\n");
    }
//  END REFINEMENT FLAGS


//  BEGIN Fhkl
    if (update_Fhkl || ALLOC || FORCE_COPY){
        for(int i=0; i < _FhklLinear.size(); i++){
          cp.cu_Fhkl[i] = _FhklLinear[i];
          if (complex_miller)
              cp.cu_Fhkl2[i] = _Fhkl2Linear[i];
        }
        if (verbose>1)
            printf("H2D Done copying step Fhkl\n");
    }
//  END Fhkl

//  BEGIN panel derivative vecs
    if(update_panel_deriv_vecs || ALLOC || FORCE_COPY){
        for (int i=0; i<dF_vecs.size(); i++){
            cp.cu_dF_vecs[i] = dF_vecs[i];
            cp.cu_dS_vecs[i] = dS_vecs[i];
        }
        if (verbose>1)
            printf("H2D Done copying step panel derivative vectors\n");
    }
//  END panel derivative vecs

//  BEGIN panels fasts slows
    if (update_panels_fasts_slows || ALLOC || FORCE_COPY){
        for (int i=0; i< panels_fasts_slows.size(); i++)
            cp.cu_panels_fasts_slows[i] = panels_fasts_slows[i];
        if (verbose>1)
            printf("H2D Done copying panels_fasts_slows\n");
    }
//  END panels fasts slows


    gettimeofday(&t2, 0);
    time = (1000000.0*(t2.tv_sec-t1.tv_sec) + t2.tv_usec-t1.tv_usec)/1000.0;
    if(verbose>1)
        printf("TIME SPENT COPYING DATA HOST->DEV:  %3.10f ms \n", time);

    cp.device_is_allocated = true;
    error_msg(cudaGetLastError(), "after copy to device");

    gettimeofday(&t1, 0);

    int Npanels = fdet_vectors.size()/3;
    int num_atoms = atom_data.size()/5;
    if (fpfdp.size() == 0){ // note cannot use atom data if fpfdp is 0, make this cleaner
        num_atoms=0;
    }
    //int sm_size = number_of_sources*5*sizeof(CUDAREAL);
    //gpu_sum_over_steps<<<numblocks, blocksize, sm_size >>>(
    bool aniso_eta = UMATS_RXYZ.size() != UMATS_RXYZ_prime.size();
    bool use_nominal_hkl = !nominal_hkl.empty();
    gpu_sum_over_steps<<<numblocks, blocksize>>>(
        Npix_to_model, cp.cu_panels_fasts_slows,
        cp.cu_floatimage,
        cp.cu_d_Umat_images, cp.cu_d2_Umat_images,
        cp.cu_d_Bmat_images, cp.cu_d2_Bmat_images,
        cp.cu_d_Ncells_images, cp.cu_d2_Ncells_images,
        cp.cu_d_fcell_images, cp.cu_d2_fcell_images,
        cp.cu_d_eta_images, cp.cu_d2_eta_images,
        cp.cu_d_lambda_images, cp.cu_d2_lambda_images,
        cp.cu_d_panel_rot_images, cp.cu_d2_panel_rot_images,
        cp.cu_d_panel_orig_images, cp.cu_d2_panel_orig_images,
        cp.cu_d_sausage_XYZ_scale_images,
        cp.cu_d_fp_fdp_images,
        Nsteps, _printout_fpixel, _printout_spixel, _printout, _default_F,
        oversample,  _oversample_omega, subpixel_size, pixel_size,
        detector_thickstep, _detector_thick, cp.cu_close_distances, detector_attnlen,
        detector_thicksteps, number_of_sources, phisteps, UMATS.size(),
        use_lambda_coefficients, lambda0, lambda1,
        eig_U, eig_O, eig_B, RXYZ,
        cp.cu_dF_vecs,
        cp.cu_dS_vecs,
        cp.cu_UMATS_RXYZ,
        cp.cu_UMATS_RXYZ_prime,
        cp.cu_UMATS_RXYZ_dbl_prime,
        cp.cu_RotMats,
        cp.cu_dRotMats,
        cp.cu_d2RotMats,
        cp.cu_UMATS,
        cp.cu_dB_Mats,
        cp.cu_dB2_Mats,
        cp.cu_AMATS,
        cp.cu_sausages_RXYZ, cp.cu_d_sausages_RXYZ, cp.cu_sausages_U, cp.cu_sausages_scale,
        cp.cu_source_X, cp.cu_source_Y, cp.cu_source_Z, cp.cu_source_lambda, cp.cu_source_I,
        kahn_factor,
        Na, Nb, Nc,
        Nd, Ne, Nf,
        phi0, phistep,
        spindle_vec, _polarization_axis,
        h_range, k_range, l_range,
        h_max, h_min, k_max, k_min, l_max, l_min, dmin,
        fudge, complex_miller, verbose, only_save_omega_kahn,
        isotropic_ncells, compute_curvatures,
        cp.cu_Fhkl, cp.cu_Fhkl2,
        cp.cu_refine_Bmat, cp.cu_refine_Ncells, refine_Ncells_def, cp.cu_refine_panel_origin, cp.cu_refine_panel_rot,
        refine_fcell, cp.cu_refine_lambda, refine_eta, cp.cu_refine_Umat,
        refine_sausages, num_sausages,
        cp.cu_fdet_vectors, cp.cu_sdet_vectors,
        cp.cu_odet_vectors, cp.cu_pix0_vectors,
        _nopolar, _point_pixel, _fluence, _r_e_sqr, _spot_scale, Npanels, aniso_eta, no_Nabc_scale,
        cp.cu_fpfdp,  cp.cu_fpfdp_derivs, cp.cu_atom_data, num_atoms,
        refine_fp_fdp, cp.cu_nominal_hkl, use_nominal_hkl);

    error_msg(cudaGetLastError(), "after kernel call");

    cudaDeviceSynchronize();
    error_msg(cudaGetLastError(), "after kernel completion");

    if(verbose>1)
        printf("KERNEL_COMPLETE gpu_sum_over_steps\n");
    gettimeofday(&t2, 0);
    time = (1000000.0*(t2.tv_sec-t1.tv_sec) + t2.tv_usec-t1.tv_usec)/1000.0;
    if(verbose>1)
        printf("TIME SPENT(KERNEL):  %3.10f ms \n", time);

    gettimeofday(&t1, 0);
//  COPY BACK FROM DEVICE
    for (int i=0; i< Npix_to_model; i++){
        floatimage[i] = cp.cu_floatimage[i];
    }
    for (int i=0; i<3*Npix_to_model; i++){
        d_fcell_images[i] = cp.cu_d_fcell_images[i];
        d_Umat_images[i] = cp.cu_d_Umat_images[i];
        d2_Umat_images[i] = cp.cu_d2_Umat_images[i];
        d_panel_rot_images[i] = cp.cu_d_panel_rot_images[i];
        d_panel_orig_images[i] = cp.cu_d_panel_orig_images[i];
        d_eta_images[i] = cp.cu_d_eta_images[i];
        d2_eta_images[i] = cp.cu_d2_eta_images[i];
    }

    for(int i=0; i<6*Npix_to_model; i++){
        d_Ncells_images[i] = cp.cu_d_Ncells_images[i];
        d2_Ncells_images[i] = cp.cu_d2_Ncells_images[i];
        d_Bmat_images[i] = cp.cu_d_Bmat_images[i];
        d2_Bmat_images[i] = cp.cu_d2_Bmat_images[i];
    }
    for(int i=0; i<2*Npix_to_model; i++)
        d_lambda_images[i] = cp.cu_d_lambda_images[i];

    for (int i=0; i< num_sausages*4*Npix_to_model; i++)
        d_sausage_XYZ_scale_images[i] = cp.cu_d_sausage_XYZ_scale_images[i];

    for (int i=0; i< 2*Npix_to_model; i++)
        d_fp_fdp_images[i] = cp.cu_d_fp_fdp_images[i];

    gettimeofday(&t2, 0);
    time = (1000000.0*(t2.tv_sec-t1.tv_sec) + t2.tv_usec-t1.tv_usec)/1000.0;
    if(verbose>1)
        printf("TIME SPENT COPYING BACK :  %3.10f ms \n", time);
    error_msg(cudaGetLastError(), "After copy to host");
}


void freedom(diffBragg_cudaPointers& cp){

    if (cp.device_is_allocated){
        gpuErr(cudaFree( cp.cu_floatimage));
        gpuErr(cudaFree( cp.cu_d_Umat_images));
        gpuErr(cudaFree( cp.cu_d_Bmat_images));
        gpuErr(cudaFree( cp.cu_d_Ncells_images));
        gpuErr(cudaFree( cp.cu_d2_Umat_images));
        gpuErr(cudaFree( cp.cu_d2_Bmat_images));
        gpuErr(cudaFree( cp.cu_d2_Ncells_images));
        gpuErr(cudaFree( cp.cu_d_eta_images));
        gpuErr(cudaFree( cp.cu_d2_eta_images));
        gpuErr(cudaFree( cp.cu_d_fcell_images));
        gpuErr(cudaFree( cp.cu_d_lambda_images));
        gpuErr(cudaFree( cp.cu_d_panel_rot_images));
        gpuErr(cudaFree( cp.cu_d_panel_orig_images));
        gpuErr(cudaFree( cp.cu_d_sausage_XYZ_scale_images));
        gpuErr(cudaFree( cp.cu_d_fp_fdp_images));

        gpuErr(cudaFree(cp.cu_Fhkl));
        if (cp.cu_Fhkl2 != NULL)
            gpuErr(cudaFree(cp.cu_Fhkl2));

        gpuErr(cudaFree(cp.cu_fdet_vectors));
        gpuErr(cudaFree(cp.cu_sdet_vectors));
        gpuErr(cudaFree(cp.cu_odet_vectors));
        gpuErr(cudaFree(cp.cu_pix0_vectors));
        gpuErr(cudaFree(cp.cu_close_distances));
        if (cp.cu_nominal_hkl != NULL)
          gpuErr(cudaFree(cp.cu_nominal_hkl));
        if (cp.cu_atom_data != NULL)
            gpuErr(cudaFree(cp.cu_atom_data));
        if(cp.cu_fpfdp != NULL)
            gpuErr(cudaFree(cp.cu_fpfdp));
        if(cp.cu_fpfdp_derivs != NULL)
            gpuErr(cudaFree(cp.cu_fpfdp_derivs));

        gpuErr(cudaFree(cp.cu_source_X));
        gpuErr(cudaFree(cp.cu_source_Y));
        gpuErr(cudaFree(cp.cu_source_Z));
        gpuErr(cudaFree(cp.cu_source_I));
        gpuErr(cudaFree(cp.cu_source_lambda));

        gpuErr(cudaFree(cp.cu_UMATS));
        gpuErr(cudaFree(cp.cu_UMATS_RXYZ));
        gpuErr(cudaFree(cp.cu_AMATS));
        if(cp.cu_UMATS_RXYZ_prime != NULL)
            gpuErr(cudaFree(cp.cu_UMATS_RXYZ_prime));
        if(cp.cu_UMATS_RXYZ_dbl_prime != NULL)
            gpuErr(cudaFree(cp.cu_UMATS_RXYZ_dbl_prime));
        gpuErr(cudaFree(cp.cu_RotMats));
        gpuErr(cudaFree(cp.cu_dRotMats));
        gpuErr(cudaFree(cp.cu_d2RotMats));
        gpuErr(cudaFree(cp.cu_dB_Mats));
        gpuErr(cudaFree(cp.cu_dB2_Mats));
        gpuErr(cudaFree(cp.cu_sausages_RXYZ));
        gpuErr(cudaFree(cp.cu_d_sausages_RXYZ));
        gpuErr(cudaFree(cp.cu_sausages_U));
        gpuErr(cudaFree(cp.cu_sausages_scale));

        gpuErr(cudaFree(cp.cu_dF_vecs));
        gpuErr(cudaFree(cp.cu_dS_vecs));

        gpuErr(cudaFree(cp.cu_refine_Bmat));
        gpuErr(cudaFree(cp.cu_refine_Umat));
        gpuErr(cudaFree(cp.cu_refine_Ncells));
        gpuErr(cudaFree(cp.cu_refine_lambda));
        gpuErr(cudaFree(cp.cu_refine_panel_origin));
        gpuErr(cudaFree(cp.cu_refine_panel_rot));

        gpuErr(cudaFree(cp.cu_panels_fasts_slows));

        cp.device_is_allocated = false;
        cp.npix_allocated = 0;
    }
}



// Kernel function to add the elements of two arrays
__global__
void phat_add(int n, float *x, float *y)
{
  for (int i = 0; i < n; i++)
    y[i] = x[i] + y[i];
}

int phat_main(void)
{
  int N = 1<<20;
  float *x, *y;

  // Allocate Unified Memory  accessible from CPU or GPU
  cudaMallocManaged(&x, N*sizeof(float));
  cudaMallocManaged(&y, N*sizeof(float));

  // initialize x and y arrays on the host
  for (int i = 0; i < N; i++) {
    x[i] = 1.0f;
    y[i] = 2.0f;
  }

  // Run kernel on 1M elements on the GPU
  phat_add<<<1, 1>>>(N, x, y);

  // Wait for GPU to finish before accessing on host
  cudaDeviceSynchronize();

  // Check for errors (all values should be 3.0f)
  float maxError = 0.0f;
  for (int i = 0; i < N; i++)
    maxError = fmax(maxError, fabs(y[i]-3.0f));
  std::cout << "Max error: " << maxError << std::endl;

  // Free memory
  cudaFree(x);
  cudaFree(y);

  return 0;
}
