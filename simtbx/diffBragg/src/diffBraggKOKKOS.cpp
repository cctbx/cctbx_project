#include <stdio.h>
#include <sys/time.h>

#include "diffBraggKOKKOS.h"
#include "diffBragg_kokkos_kernel.h"

void diffBraggKOKKOS::diffBragg_sum_over_steps_kokkos(
    int Npix_to_model,
    std::vector<unsigned int>& panels_fasts_slows,
    image_type& floatimage,
    images& d_image,
    images& d2_image,
    step_arrays& db_steps,
    detector& db_det,
    beam& db_beam,
    crystal& db_cryst,
    flags& db_flags,
    cuda_flags& db_cu_flags,
    // diffBragg_kokkosPointers& kp,
    timer_variables& TIMERS) {
    if (db_cryst.phi0 != 0 || db_cryst.phisteps > 1) {
        printf(
            "PHI (goniometer position) not supported in GPU code: phi0=%f phisteps=%d, phistep=%d\n",
            db_cryst.phi0, db_cryst.phisteps, db_cryst.phistep);
        exit(-1);
    }

    // static bool m_device_is_allocated = false;
    // static int m_npix_allocated=0;
    // static int m_previous_nsource=0;

    // static vector_uint_t m_panels_fasts_slows =
    // vector_uint_t("m_panels_fasts_slows", 0);

    // CUDAREAL_VIEW(m_floatimage);
    // CUDAREAL_VIEW(m_wavelenimage);
    // CUDAREAL_VIEW(m_d_diffuse_sigma_images);
    // CUDAREAL_VIEW(m_d_diffuse_gamma_images);
    // CUDAREAL_VIEW(m_d_Umat_images);
    // CUDAREAL_VIEW(m_d_Bmat_images);
    // CUDAREAL_VIEW(m_d_Ncells_images);
    // CUDAREAL_VIEW(m_d_fcell_images);
    // CUDAREAL_VIEW(m_d_eta_images);
    // CUDAREAL_VIEW(m_d2_eta_images);
    // CUDAREAL_VIEW(m_d_lambda_images);
    // CUDAREAL_VIEW(m_d_panel_rot_images);
    // CUDAREAL_VIEW(m_d_panel_orig_images);

    // CUDAREAL_VIEW(m_d2_Umat_images);
    // CUDAREAL_VIEW(m_d2_Bmat_images);
    // CUDAREAL_VIEW(m_d2_Ncells_images);
    // CUDAREAL_VIEW(m_d2_fcell_images);
    // CUDAREAL_VIEW(m_d2_lambda_images);
    // CUDAREAL_VIEW(m_d2_panel_rot_images);
    // CUDAREAL_VIEW(m_d2_panel_orig_images);

    // CUDAREAL_VIEW(m_d_sausage_XYZ_scale_images);
    // CUDAREAL_VIEW(m_d_fp_fdp_images);

    // INTEGER_VIEW(m_subS_pos);
    // INTEGER_VIEW(m_subF_pos);
    // INTEGER_VIEW(m_thick_pos);
    // INTEGER_VIEW(m_source_pos);
    // INTEGER_VIEW(m_mos_pos);
    // INTEGER_VIEW(m_phi_pos);
    // INTEGER_VIEW(m_sausage_pos);

    // CUDAREAL_VIEW(m_Fhkl);
    // CUDAREAL_VIEW(m_Fhkl2);

    // CUDAREAL_VIEW(m_fdet_vectors);
    // CUDAREAL_VIEW(m_sdet_vectors);
    // CUDAREAL_VIEW(m_odet_vectors);
    // CUDAREAL_VIEW(m_pix0_vectors);
    // CUDAREAL_VIEW(m_close_distances);

    // INTEGER_VIEW(m_nominal_hkl);
    // CUDAREAL_VIEW(m_fpfdp);
    // CUDAREAL_VIEW(m_fpfdp_derivs);
    // CUDAREAL_VIEW(m_atom_data);

    // CUDAREAL_VIEW(m_source_X);
    // CUDAREAL_VIEW(m_source_Y);
    // CUDAREAL_VIEW(m_source_Z);
    // CUDAREAL_VIEW(m_source_I);
    // CUDAREAL_VIEW(m_source_lambda);
    // static int m_sources;
    // static bool m_sources_are_allocated = false;
    // static bool m_sources_recopy = false;

    // MATRIX3_VIEW(m_UMATS);
    // MATRIX3_VIEW(m_dB_Mats);
    // MATRIX3_VIEW(m_dB2_Mats);
    // MATRIX3_VIEW(m_UMATS_RXYZ);
    // MATRIX3_VIEW(m_UMATS_RXYZ_prime);
    // MATRIX3_VIEW(m_UMATS_RXYZ_dbl_prime);
    // MATRIX3_VIEW(m_RotMats);
    // MATRIX3_VIEW(m_dRotMats);
    // MATRIX3_VIEW(m_d2RotMats);

    // MATRIX3_VIEW(m_AMATS);

    // static vector_vec3_t m_dF_vecs = vector_vec3_t("m_dF_vecs", 0);
    // static vector_vec3_t m_dS_vecs = vector_vec3_t("m_dS_vecs", 0);

    // MATRIX3_VIEW(m_sausages_RXYZ);
    // MATRIX3_VIEW(m_d_sausages_RXYZ);
    // MATRIX3_VIEW(m_sausages_U);
    // CUDAREAL_VIEW(m_sausages_scale);

    // static vector_bool_t m_refine_Bmat = vector_bool_t("m_refine_Bmat", 6);
    // static vector_bool_t m_refine_Umat = vector_bool_t("m_refine_Umat", 3);
    // static vector_bool_t m_refine_Ncells = vector_bool_t("m_refine_Ncells", 3);
    // static vector_bool_t m_refine_panel_origin = vector_bool_t("m_refine_panel_origin", 3);
    // static vector_bool_t m_refine_panel_rot = vector_bool_t("m_refine_panel_rot", 3);
    // static vector_bool_t m_refine_lambda = vector_bool_t("m_refine_lambda", 2);

    // int numblocks;
    // int blocksize;
    // char* diffBragg_blocks = getenv("DIFFBRAGG_NUM_BLOCKS");
    // char* diffBragg_threads = getenv("DIFFBRAGG_THREADS_PER_BLOCK");
    // if (diffBragg_threads==NULL)
    //     blocksize=128;
    // else
    //     blocksize=atoi(diffBragg_threads);

    // if (diffBragg_blocks==NULL)
    //     numblocks = (Npix_to_model+blocksize-1)/blocksize;
    // else
    //     numblocks = atoi(diffBragg_blocks);

    // int cuda_devices;
    // cudaGetDeviceCount(&cuda_devices);

    // error_msg(cudaGetLastError(), "after device count");
    // if (db_flags.verbose > 1)
    //     printf("Found %d CUDA-capable devices\n", cuda_devices);

    // //if (device_Id <= cuda_devices)
    // gpuErr(cudaSetDevice(db_cu_flags.device_Id));

    double time;
    struct timeval t1, t2;  //, t3 ,t4;
    gettimeofday(&t1, 0);

    //  determine if we need to allocate pixels, and how many.
    //  For best usage, one should use the diffBragg property (visible from Python) Npix_to_allocate
    //  in order to just allocate to the GPU - this is useful for ensemble refinement, where each
    //  shot can have a variable number of pixels being modeled, and ony only needs to allocate the
    //  device once (with the largest expected number of pixels for a given shot)
    // TODO clean up this logic a bit
    if (m_device_is_allocated && (Npix_to_model > m_npix_allocated)) {
        printf(
            "Need to re-allocate pixels, currently have %d allocated, but trying to model %d\n",
            m_npix_allocated, Npix_to_model);
        exit(-1);
    } else if (db_cu_flags.Npix_to_allocate == -1) {
        db_cu_flags.Npix_to_allocate = Npix_to_model;
    } else if (Npix_to_model > db_cu_flags.Npix_to_allocate) {
        printf(
            "Npix to model=%d is greater than the number of pixel requested for allocation (%d)!\n",
            Npix_to_model, db_cu_flags.Npix_to_allocate);
        exit(-1);
    }

    //  support dynamic allocation for different numbers of sources
    if (m_previous_nsource != 0 && m_previous_nsource != db_beam.number_of_sources) {
        printf("Resizing for %d sources!:\n", db_beam.number_of_sources);
        resize(m_source_X, db_beam.number_of_sources);
        resize(m_source_Y, db_beam.number_of_sources);
        resize(m_source_Z, db_beam.number_of_sources);
        resize(m_source_I, db_beam.number_of_sources);
        resize(m_source_lambda, db_beam.number_of_sources);
        m_previous_nsource = db_beam.number_of_sources;
    }

    if (m_device_is_allocated) {
        if (db_flags.verbose) {
            printf(
                "Will model %d pixels (GPU has %d pre-allocated pix)\n", Npix_to_model,
                m_npix_allocated);
        }
    } else {
        if (db_flags.verbose) {
            printf(
                "Will model %d pixels and allocate %d pix\n", Npix_to_model,
                db_cu_flags.Npix_to_allocate);
        }
        resize(m_source_X, db_beam.number_of_sources);
        resize(m_source_Y, db_beam.number_of_sources);
        resize(m_source_Z, db_beam.number_of_sources);
        resize(m_source_I, db_beam.number_of_sources);
        resize(m_source_lambda, db_beam.number_of_sources);
        m_previous_nsource = db_beam.number_of_sources;

        resize(m_UMATS, db_cryst.UMATS.size());
        resize(m_UMATS_RXYZ, db_cryst.UMATS_RXYZ.size());
        resize(m_AMATS, db_cryst.UMATS_RXYZ.size());
        // gpuErr(cudaMallocManaged((void **)&m_cu_UMATS, db_cryst.UMATS.size()*sizeof(MAT3)));
        // gpuErr(cudaMallocManaged((void **)&m_cu_UMATS_RXYZ, db_cryst.UMATS_RXYZ.size()*sizeof(MAT3)));
        // gpuErr(cudaMallocManaged((void **)&m_cu_AMATS, db_cryst.UMATS_RXYZ.size()*sizeof(MAT3)));
        if (db_cryst.UMATS_RXYZ_prime.size() > 0)
            resize(m_UMATS_RXYZ_prime, db_cryst.UMATS_RXYZ_prime.size());
        // gpuErr(cudaMallocManaged((void **)&m_cu_UMATS_RXYZ_prime,
        // db_cryst.UMATS_RXYZ_prime.size()*sizeof(MAT3)));
        if (db_cryst.UMATS_RXYZ_dbl_prime.size() > 0)
            resize(m_UMATS_RXYZ_dbl_prime, db_cryst.UMATS_RXYZ_dbl_prime.size());
        // gpuErr(cudaMallocManaged((void **)&m_cu_UMATS_RXYZ_dbl_prime,
        // db_cryst.UMATS_RXYZ_dbl_prime.size()*sizeof(MAT3)));

        resize(m_dB_Mats, db_cryst.dB_Mats.size());
        resize(m_dB2_Mats, db_cryst.dB2_Mats.size());
        // gpuErr(cudaMallocManaged((void **)&m_cu_dB_Mats,
        // db_cryst.dB_Mats.size()*sizeof(MAT3)));
        // gpuErr(cudaMallocManaged((void **)&m_cu_dB2_Mats,
        // db_cryst.dB2_Mats.size()*sizeof(MAT3)));

        resize(m_RotMats, db_cryst.RotMats.size());
        resize(m_dRotMats, db_cryst.dRotMats.size());
        resize(m_d2RotMats, db_cryst.d2RotMats.size());
        // gpuErr(cudaMallocManaged((void **)&m_cu_RotMats,
        // db_cryst.RotMats.size()*sizeof(MAT3)));
        // gpuErr(cudaMallocManaged((void **)&m_cu_dRotMats,
        // db_cryst.dRotMats.size()*sizeof(MAT3)));
        // gpuErr(cudaMallocManaged((void **)&m_cu_d2RotMats,
        // db_cryst.d2RotMats.size()*sizeof(MAT3)));

        resize(m_fdet_vectors, db_det.fdet_vectors.size());
        resize(m_sdet_vectors, db_det.sdet_vectors.size());
        resize(m_odet_vectors, db_det.odet_vectors.size());
        resize(m_pix0_vectors, db_det.pix0_vectors.size());
        resize(m_close_distances, db_det.close_distances.size());
        // gpuErr(cudaMallocManaged(&m_cu_fdet_vectors,
        // db_det.fdet_vectors.size()*sizeof(CUDAREAL)));
        // gpuErr(cudaMallocManaged(&m_cu_sdet_vectors,
        // db_det.fdet_vectors.size()*sizeof(CUDAREAL)));
        // gpuErr(cudaMallocManaged(&m_cu_odet_vectors,
        // db_det.fdet_vectors.size()*sizeof(CUDAREAL)));
        // gpuErr(cudaMallocManaged(&m_cu_pix0_vectors,
        // db_det.fdet_vectors.size()*sizeof(CUDAREAL)));
        // gpuErr(cudaMallocManaged(&m_cu_close_distances,
        // db_det.close_distances.size()*sizeof(CUDAREAL)));

        if (db_cryst.fpfdp.size() > 0) {
            resize(m_fpfdp, db_cryst.fpfdp.size());
            resize(m_atom_data, db_cryst.atom_data.size());
            // gpuErr(cudaMallocManaged(&m_cu_fpfdp,
            // db_cryst.fpfdp.size()*sizeof(CUDAREAL)));
            // gpuErr(cudaMallocManaged(&m_cu_atom_data,
            // db_cryst.atom_data.size()*sizeof(CUDAREAL)));
        }
        if (db_cryst.fpfdp_derivs.size() > 0)
            resize(m_fpfdp_derivs, db_cryst.fpfdp_derivs.size());
        // gpuErr(cudaMallocManaged(&m_cu_fpfdp_derivs,
        // db_cryst.fpfdp_derivs.size()*sizeof(CUDAREAL)));

        // already done this in diffBRaggKOKKOS.h
        // gpuErr(cudaMallocManaged(&m_cu_refine_Bmat, 6*sizeof(bool)));
        // gpuErr(cudaMallocManaged(&m_cu_refine_Umat, 3*sizeof(bool)));
        // gpuErr(cudaMallocManaged(&m_cu_refine_Ncells, 3*sizeof(bool)));
        // gpuErr(cudaMallocManaged(&m_cu_refine_panel_origin, 3*sizeof(bool)));
        // gpuErr(cudaMallocManaged(&m_cu_refine_panel_rot, 3*sizeof(bool)));
        // gpuErr(cudaMallocManaged(&m_cu_refine_lambda, 2*sizeof(bool)));

        resize(m_Fhkl, db_cryst.FhklLinear.size());
        // gpuErr(cudaMallocManaged(&m_cu_Fhkl,
        // db_cryst.FhklLinear.size()*sizeof(CUDAREAL)));
        if (db_flags.complex_miller)
            resize(m_Fhkl2, db_cryst.FhklLinear.size());
        // gpuErr(cudaMallocManaged(&m_cu_Fhkl2,
        // db_cryst.FhklLinear.size()*sizeof(CUDAREAL)));

        resize(m_dF_vecs, db_det.dF_vecs.size());
        resize(m_dS_vecs, db_det.dF_vecs.size());
        // gpuErr(cudaMallocManaged((void **)&m_cu_dF_vecs,
        // db_det.dF_vecs.size()*sizeof(VEC3))); gpuErr(cudaMallocManaged((void
        // **)&m_cu_dS_vecs, db_det.dF_vecs.size()*sizeof(VEC3)));

        // gettimeofday(&t3, 0));
        resize(m_floatimage, db_cu_flags.Npix_to_allocate);
        // gpuErr(cudaMallocManaged(&m_cu_floatimage,
        // db_cu_flags.Npix_to_allocate*sizeof(CUDAREAL) ));
        if (db_flags.wavelength_img) {
            resize(m_wavelenimage, db_cu_flags.Npix_to_allocate);
            // gpuErr(cudaMallocManaged(&m_cu_wavelenimage,
            // db_cu_flags.Npix_to_allocate*sizeof(CUDAREAL) ));
        }
        if (db_flags.refine_diffuse) {
            resize(m_d_diffuse_gamma_images, db_cu_flags.Npix_to_allocate * 3);
            resize(m_d_diffuse_sigma_images, db_cu_flags.Npix_to_allocate * 3);
            // gpuErr(cudaMallocManaged(&m_cu_d_diffuse_gamma_images,
            // db_cu_flags.Npix_to_allocate*3*sizeof(CUDAREAL)));
            // gpuErr(cudaMallocManaged(&m_cu_d_diffuse_sigma_images,
            // db_cu_flags.Npix_to_allocate*3*sizeof(CUDAREAL)));
        }
        if (db_flags.refine_fcell) {
            resize(m_d_fcell_images, db_cu_flags.Npix_to_allocate);
            resize(m_d2_fcell_images, db_cu_flags.Npix_to_allocate);
            // gpuErr(cudaMallocManaged(&m_cu_d_fcell_images,
            // db_cu_flags.Npix_to_allocate*1*sizeof(CUDAREAL)));
            // gpuErr(cudaMallocManaged(&m_cu_d2_fcell_images,
            // db_cu_flags.Npix_to_allocate*1*sizeof(CUDAREAL)));
        }
        if (db_flags.refine_eta) {
            resize(m_d_eta_images, db_cu_flags.Npix_to_allocate * 3);
            resize(m_d2_eta_images, db_cu_flags.Npix_to_allocate * 3);
            // gpuErr(cudaMallocManaged(&m_cu_d_eta_images,
            // db_cu_flags.Npix_to_allocate*3*sizeof(CUDAREAL)));
            // gpuErr(cudaMallocManaged(&m_cu_d2_eta_images,
            // db_cu_flags.Npix_to_allocate*3*sizeof(CUDAREAL)));
        }
        if (std::count(db_flags.refine_Umat.begin(), db_flags.refine_Umat.end(), true) > 0) {
            resize(m_d_Umat_images, db_cu_flags.Npix_to_allocate * 3);
            resize(m_d2_Umat_images, db_cu_flags.Npix_to_allocate * 3);
            // gpuErr(cudaMallocManaged(&m_cu_d_Umat_images,
            // db_cu_flags.Npix_to_allocate*3*sizeof(CUDAREAL) ));
            // gpuErr(cudaMallocManaged(&m_cu_d2_Umat_images,
            // db_cu_flags.Npix_to_allocate*3*sizeof(CUDAREAL) ));
        }
        if (std::count(db_flags.refine_Ncells.begin(), db_flags.refine_Ncells.end(), true) > 0 ||
            db_flags.refine_Ncells_def) {
            resize(m_d_Ncells_images, db_cu_flags.Npix_to_allocate * 6);
            resize(m_d2_Ncells_images, db_cu_flags.Npix_to_allocate * 6);
            // gpuErr(cudaMallocManaged(&m_cu_d_Ncells_images,
            // db_cu_flags.Npix_to_allocate*6*sizeof(CUDAREAL)));
            // gpuErr(cudaMallocManaged(&m_cu_d2_Ncells_images,
            // db_cu_flags.Npix_to_allocate*6*sizeof(CUDAREAL)));
        }
        if (std::count(db_flags.refine_panel_rot.begin(), db_flags.refine_panel_rot.end(), true) >
            0)
            resize(m_d_panel_rot_images, db_cu_flags.Npix_to_allocate * 3);
        // gpuErr(cudaMallocManaged(&m_cu_d_panel_rot_images,
        // db_cu_flags.Npix_to_allocate*3*sizeof(CUDAREAL)));
        if (std::count(
                db_flags.refine_panel_origin.begin(), db_flags.refine_panel_origin.end(), true) > 0)
            resize(m_d_panel_orig_images, db_cu_flags.Npix_to_allocate * 3);
        // gpuErr(cudaMallocManaged(&m_cu_d_panel_orig_images,
        // db_cu_flags.Npix_to_allocate*3*sizeof(CUDAREAL)));
        if (std::count(db_flags.refine_lambda.begin(), db_flags.refine_lambda.end(), true) > 0)
            resize(m_d_lambda_images, db_cu_flags.Npix_to_allocate * 2);
        // gpuErr(cudaMallocManaged(&m_cu_d_lambda_images,
        // db_cu_flags.Npix_to_allocate*2*sizeof(CUDAREAL)));
        if (std::count(db_flags.refine_Bmat.begin(), db_flags.refine_Bmat.end(), true) > 0) {
            resize(m_d_Bmat_images, db_cu_flags.Npix_to_allocate * 6);
            resize(m_d2_Bmat_images, db_cu_flags.Npix_to_allocate * 6);
            // gpuErr(cudaMallocManaged(&m_cu_d_Bmat_images,
            // db_cu_flags.Npix_to_allocate*6*sizeof(CUDAREAL)));
            // gpuErr(cudaMallocManaged(&m_cu_d2_Bmat_images,
            // db_cu_flags.Npix_to_allocate*6*sizeof(CUDAREAL)));
        }
        if (db_flags.refine_fp_fdp)
            resize(m_d_fp_fdp_images, db_cu_flags.Npix_to_allocate * 2);
        // gpuErr(cudaMallocManaged(&m_cu_d_fp_fdp_images,
        // db_cu_flags.Npix_to_allocate*2*sizeof(CUDAREAL)));
        if (db_cryst.nominal_hkl.size() > 0)
            resize(m_nominal_hkl, db_cu_flags.Npix_to_allocate * 3);
        // gpuErr(cudaMallocManaged(&m_cu_nominal_hkl,
        // db_cu_flags.Npix_to_allocate*3*sizeof(int)));

        // gettimeofday(&t4, 0);
        // time = (1000000.0*(t4.tv_sec-t3.tv_sec) +
        // t4.tv_usec-t3.tv_usec)/1000.0; printf("TIME SPENT ALLOCATING (IMAGES
        // ONLY):  %3.10f ms \n", time);
        resize(m_panels_fasts_slows, db_cu_flags.Npix_to_allocate * 3);
        // gpuErr(cudaMallocManaged(&m_cu_panels_fasts_slows,
        // db_cu_flags.Npix_to_allocate*3*sizeof(panels_fasts_slows[0])));
        m_npix_allocated = db_cu_flags.Npix_to_allocate;
    }  // END of allocation

    bool ALLOC = !m_device_is_allocated;  // shortcut variable

    gettimeofday(&t2, 0);
    time = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000.0;
    if (TIMERS.recording)
        TIMERS.cuda_alloc += time;
    if (db_flags.verbose > 1)
        printf("TIME SPENT ALLOCATING (TOTAL):  %3.10f ms \n", time);

    // ALLOC = false;
    //  BEGIN COPYING DATA
    gettimeofday(&t1, 0);
    bool FORCE_COPY = true;

    //  END step position

    //  BEGIN sources
    if (db_cu_flags.update_sources || ALLOC || FORCE_COPY) {
        int source_count = db_beam.number_of_sources;
        kokkostbx::transfer_double2kokkos(m_source_X, db_beam.source_X, source_count);
        kokkostbx::transfer_double2kokkos(m_source_Y, db_beam.source_Y, source_count);
        kokkostbx::transfer_double2kokkos(m_source_Z, db_beam.source_Z, source_count);
        kokkostbx::transfer_double2kokkos(m_source_I, db_beam.source_I, source_count);
        kokkostbx::transfer_double2kokkos(m_source_lambda, db_beam.source_lambda, source_count);

        vector_cudareal_t dumb_test = vector_cudareal_t("dumb_test", 3);
        Kokkos::parallel_for(
            "normalize incident vector", source_count, KOKKOS_LAMBDA(const int& i) {
                // VEC3 incident{m_source_X(i), m_source_Y(i), m_source_Z(i)};
                VEC3 incident{3, 5, 2};
                CUDAREAL a = dumb_test(2);
                incident.normalize();
                // m_source_X(i) = incident.x_val();
                // m_source_Y(i) = incident.y_val();
                // m_source_Z(i) = incident.z_val();
            });

        // for (int i=0; i< db_beam.number_of_sources; i++){
        //     VEC3 incident{db_beam.source_X[i], db_beam.source_Y[i],
        //     db_beam.source_Z[i]}; incident /= incident.norm();

        //     m_cu_source_X[i] = incident[0];
        //     m_cu_source_Y[i] = incident[1];
        //     m_cu_source_Z[i] = incident[2];
        //     m_cu_source_I[i] = db_beam.source_I[i];
        //     m_cu_source_lambda[i] = db_beam.source_lambda[i];
        // }
        if (db_flags.verbose > 1)
            printf("H2D sources\n");
    }
    //  END sources

    /*
    //  UMATS
        if (db_cu_flags.update_umats || ALLOC||FORCE_COPY){
            for (int i=0; i< db_cryst.UMATS.size(); i++)
                m_cu_UMATS[i] = db_cryst.UMATS[i];
            for (int i=0; i < db_cryst.UMATS_RXYZ.size(); i++)
                m_cu_UMATS_RXYZ[i] = db_cryst.UMATS_RXYZ[i];
            for (int i=0; i < db_cryst.UMATS_RXYZ_prime.size(); i++)
                m_cu_UMATS_RXYZ_prime[i] = db_cryst.UMATS_RXYZ_prime[i];
            for (int i=0; i < db_cryst.UMATS_RXYZ_dbl_prime.size(); i++)
                m_cu_UMATS_RXYZ_dbl_prime[i] = db_cryst.UMATS_RXYZ_dbl_prime[i];
            if(db_flags.verbose>1)
                printf("H2D Done copying Umats\n") ;
        }
    //  END UMATS


        if (db_cu_flags.update_umats || ALLOC||FORCE_COPY){
        MAT3 Amat_init = db_cryst.eig_U*db_cryst.eig_B*1e10*(db_cryst.eig_O.transpose());
        for(int i_mos =0; i_mos< db_cryst.UMATS_RXYZ.size(); i_mos++){
            m_cu_AMATS[i_mos] = (db_cryst.UMATS_RXYZ[i_mos]*Amat_init).transpose();
                }
            if(db_flags.verbose>1)
                printf("H2D Done copying Amats\n") ;
        }


    //  BMATS
        if(db_cu_flags.update_dB_mats || ALLOC || FORCE_COPY){
            for (int i=0; i< db_cryst.dB_Mats.size(); i++)
                m_cu_dB_Mats[i] = db_cryst.dB_Mats[i];
            for (int i=0; i< db_cryst.dB2_Mats.size(); i++)
                m_cu_dB2_Mats[i] = db_cryst.dB2_Mats[i];
            if(db_flags.verbose>1)
                printf("H2D Done copying dB_Mats\n") ;
        }
    //  END BMATS


    //  ROT MATS
        if(db_cu_flags.update_rotmats || ALLOC || FORCE_COPY){
            for (int i=0; i<db_cryst.RotMats.size(); i++)
                m_cu_RotMats[i] = db_cryst.RotMats[i];
            for (int i=0; i<db_cryst.dRotMats.size(); i++)
                m_cu_dRotMats[i] = db_cryst.dRotMats[i];
            for (int i=0; i<db_cryst.d2RotMats.size(); i++)
                m_cu_d2RotMats[i] = db_cryst.d2RotMats[i];
            if (db_flags.verbose>1)
              printf("H2D Done copying rotmats\n");
        }
    //  END ROT MATS

    //  DETECTOR VECTORS
        if (db_cu_flags.update_detector || ALLOC || FORCE_COPY){
            for (int i=0; i<db_det.fdet_vectors.size(); i++){
                m_cu_fdet_vectors[i] = db_det.fdet_vectors[i];
                m_cu_sdet_vectors[i] = db_det.sdet_vectors[i];
                m_cu_odet_vectors[i] = db_det.odet_vectors[i];
                m_cu_pix0_vectors[i] = db_det.pix0_vectors[i];
            }
            for(int i=0; i < db_det.close_distances.size();i++){
                m_cu_close_distances[i] = db_det.close_distances[i];
            }
            if (db_flags.verbose>1)
              printf("H2D Done copying detector vectors\n");
        }
    //  END  DETECTOR VECTORS

        if ( ALLOC || FORCE_COPY){
          for(int i=0; i< db_cryst.nominal_hkl.size(); i++){
            m_cu_nominal_hkl[i] = db_cryst.nominal_hkl[i];
          }
          for (int i=0; i< db_cryst.atom_data.size(); i++){
            m_cu_atom_data[i] = db_cryst.atom_data[i];
          }
          if (db_flags.verbose>1)
            printf("H2D Done copying atom data\n");
          for(int i=0; i< db_cryst.fpfdp.size(); i++){
            m_cu_fpfdp[i] = db_cryst.fpfdp[i];
          }
          for(int i=0; i< db_cryst.fpfdp_derivs.size(); i++){
            m_cu_fpfdp_derivs[i] = db_cryst.fpfdp_derivs[i];
          }
          if (db_flags.verbose>1)
            printf("H2D Done copying fprime and fdblprime\n");
        }


    //  BEGIN REFINEMENT FLAGS
        if (db_cu_flags.update_refine_flags || ALLOC || FORCE_COPY){
            for (int i=0; i<3; i++){
                m_cu_refine_Umat[i] = db_flags.refine_Umat[i];
                m_cu_refine_Ncells[i] = db_flags.refine_Ncells[i];
                m_cu_refine_panel_origin[i] = db_flags.refine_panel_origin[i];
                m_cu_refine_panel_rot[i] = db_flags.refine_panel_rot[i];
            }
            for(int i=0; i<2; i++)
                m_cu_refine_lambda[i] = db_flags.refine_lambda[i];
            for(int i=0; i<6; i++)
                m_cu_refine_Bmat[i] = db_flags.refine_Bmat[i];
            if (db_flags.verbose>1)
              printf("H2D Done copying refinement flags\n");
        }
    //  END REFINEMENT FLAGS


    //  BEGIN Fhkl
        if (db_cu_flags.update_Fhkl || ALLOC || FORCE_COPY){
            for(int i=0; i < db_cryst.FhklLinear.size(); i++){
              m_cu_Fhkl[i] = db_cryst.FhklLinear[i];
              if (db_flags.complex_miller)
                  m_cu_Fhkl2[i] = db_cryst.Fhkl2Linear[i];
            }
            if (db_flags.verbose>1)
                printf("H2D Done copying step Fhkl\n");
        }
    //  END Fhkl

    //  BEGIN panel derivative vecs
        if(db_cu_flags.update_panel_deriv_vecs || ALLOC || FORCE_COPY){
            for (int i=0; i<db_det.dF_vecs.size(); i++){
                m_cu_dF_vecs[i] = db_det.dF_vecs[i];
                m_cu_dS_vecs[i] = db_det.dS_vecs[i];
            }
            if (db_flags.verbose>1)
                printf("H2D Done copying step panel derivative vectors\n");
        }
    //  END panel derivative vecs

    //  BEGIN panels fasts slows
        if (db_cu_flags.update_panels_fasts_slows || ALLOC || FORCE_COPY){
            for (int i=0; i< panels_fasts_slows.size(); i++)
                m_cu_panels_fasts_slows[i] = panels_fasts_slows[i];
            if (db_flags.verbose>1)
                printf("H2D Done copying panels_fasts_slows\n");
        }
    //  END panels fasts slows


        gettimeofday(&t2, 0);
        time = (1000000.0*(t2.tv_sec-t1.tv_sec) + t2.tv_usec-t1.tv_usec)/1000.0;
        if (TIMERS.recording) TIMERS.cuda_copy_to_dev += time;
        if(db_flags.verbose>1)
            printf("TIME SPENT COPYING DATA HOST->DEV:  %3.10f ms \n", time);

        m_device_is_allocated = true;
        error_msg(cudaGetLastError(), "after copy to device");

        gettimeofday(&t1, 0);

        int Npanels = db_det.fdet_vectors.size()/3;
        int num_atoms = db_cryst.atom_data.size()/5;
        // note cannot use atom data if fpfdp is 0, make this cleaner
        if (db_cryst.fpfdp.size() == 0) {  
            num_atoms = 0;
        }
        // int sm_size = number_of_sources*5*sizeof(CUDAREAL);
        // gpu_sum_over_steps<<<numblocks, blocksize, sm_size >>>(
        bool aniso_eta = db_cryst.UMATS_RXYZ.size() != db_cryst.UMATS_RXYZ_prime.size();
        bool use_nominal_hkl = !db_cryst.nominal_hkl.empty();
        gpu_sum_over_steps<<<numblocks, blocksize>>>(
            Npix_to_model, m_cu_panels_fasts_slows, m_cu_floatimage, m_cu_wavelenimage,
            m_cu_d_Umat_images, m_cu_d2_Umat_images, m_cu_d_Bmat_images, m_cu_d2_Bmat_images,
            m_cu_d_Ncells_images, m_cu_d2_Ncells_images, m_cu_d_fcell_images, m_cu_d2_fcell_images,
            m_cu_d_eta_images, m_cu_d2_eta_images, m_cu_d_lambda_images, m_cu_d2_lambda_images,
            m_cu_d_panel_rot_images, m_cu_d2_panel_rot_images, m_cu_d_panel_orig_images,
            m_cu_d2_panel_orig_images, m_cu_d_fp_fdp_images, db_steps.Nsteps,
            db_flags.printout_fpixel, db_flags.printout_spixel, db_flags.printout,
            db_cryst.default_F, db_det.oversample, db_flags.oversample_omega, db_det.subpixel_size,
            db_det.pixel_size, db_det.detector_thickstep, db_det.detector_thick,
            m_cu_close_distances, db_det.detector_attnlen, db_det.detector_thicksteps,
            db_beam.number_of_sources, db_cryst.phisteps, db_cryst.UMATS.size(),
            db_flags.use_lambda_coefficients, db_beam.lambda0, db_beam.lambda1, db_cryst.eig_U,
            db_cryst.eig_O, db_cryst.eig_B, db_cryst.RXYZ, m_cu_dF_vecs, m_cu_dS_vecs,
            m_cu_UMATS_RXYZ, m_cu_UMATS_RXYZ_prime, m_cu_UMATS_RXYZ_dbl_prime, m_cu_RotMats,
            m_cu_dRotMats, m_cu_d2RotMats, m_cu_UMATS, m_cu_dB_Mats, m_cu_dB2_Mats, m_cu_AMATS,
            m_cu_source_X, m_cu_source_Y, m_cu_source_Z, m_cu_source_lambda, m_cu_source_I,
            db_beam.kahn_factor, db_cryst.Na, db_cryst.Nb, db_cryst.Nc, db_cryst.Nd, db_cryst.Ne,
            db_cryst.Nf, db_cryst.phi0, db_cryst.phistep, db_cryst.spindle_vec,
            db_beam.polarization_axis, db_cryst.h_range, db_cryst.k_range, db_cryst.l_range,
            db_cryst.h_max, db_cryst.h_min, db_cryst.k_max, db_cryst.k_min, db_cryst.l_max,
            db_cryst.l_min, db_cryst.dmin, db_cryst.fudge, db_flags.complex_miller,
            db_flags.verbose, db_flags.only_save_omega_kahn, db_flags.isotropic_ncells,
            db_flags.compute_curvatures, m_cu_Fhkl, m_cu_Fhkl2, m_cu_refine_Bmat,
            m_cu_refine_Ncells, db_flags.refine_Ncells_def, m_cu_refine_panel_origin,
            m_cu_refine_panel_rot, db_flags.refine_fcell, m_cu_refine_lambda, db_flags.refine_eta,
            m_cu_refine_Umat, m_cu_fdet_vectors, m_cu_sdet_vectors, m_cu_odet_vectors,
            m_cu_pix0_vectors, db_flags.nopolar, db_flags.point_pixel, db_beam.fluence,
            db_cryst.r_e_sqr, db_cryst.spot_scale, Npanels, aniso_eta, db_flags.no_Nabc_scale,
            m_cu_fpfdp, m_cu_fpfdp_derivs, m_cu_atom_data, num_atoms, db_flags.refine_fp_fdp,
            m_cu_nominal_hkl, use_nominal_hkl, db_cryst.anisoU, db_cryst.anisoG,
            db_flags.use_diffuse, m_cu_d_diffuse_gamma_images, m_cu_d_diffuse_sigma_images,
            db_flags.refine_diffuse, db_flags.gamma_miller_units, db_flags.refine_Icell,
            db_flags.wavelength_img);

        error_msg(cudaGetLastError(), "after kernel call");

        cudaDeviceSynchronize();
        error_msg(cudaGetLastError(), "after kernel completion");

        if(db_flags.verbose>1)
            printf("KERNEL_COMPLETE gpu_sum_over_steps\n");
        gettimeofday(&t2, 0);
        time = (1000000.0*(t2.tv_sec-t1.tv_sec) + t2.tv_usec-t1.tv_usec)/1000.0;
        if (TIMERS.recording) TIMERS.cuda_kernel += time;
        if(db_flags.verbose>1)
            printf("TIME SPENT(KERNEL):  %3.10f ms \n", time);

        gettimeofday(&t1, 0);
    //  COPY BACK FROM DEVICE
        for (int i=0; i< Npix_to_model; i++){
            floatimage[i] = m_cu_floatimage[i];
        }
        if(db_flags.wavelength_img){
            for (int i=0; i< Npix_to_model; i++){
                d_image.wavelength[i] = m_cu_wavelenimage[i];
            }
        }
        if (db_flags.refine_fcell){
            for (int i=0; i<Npix_to_model; i++){
                d_image.fcell[i] = m_cu_d_fcell_images[i];
                d2_image.fcell[i] = m_cu_d2_fcell_images[i];
            }
        }
        if (std::count(db_flags.refine_Umat.begin(), db_flags.refine_Umat.end(), true) > 0) {
            for (int i = 0; i < 3 * Npix_to_model; i++) {
                d_image.Umat[i] = m_cu_d_Umat_images[i];
                d2_image.Umat[i] = m_cu_d2_Umat_images[i];
            }
        }
        if (std::count(db_flags.refine_panel_rot.begin(), db_flags.refine_panel_rot.end(), true) >
            0) {
            for (int i = 0; i < 3 * Npix_to_model; i++)
                d_image.panel_rot[i] = m_cu_d_panel_rot_images[i];
        }
        if (std::count(
                db_flags.refine_panel_origin.begin(), db_flags.refine_panel_origin.end(), true) >
            0) {
            for (int i = 0; i < 3 * Npix_to_model; i++)
                d_image.panel_orig[i] = m_cu_d_panel_orig_images[i];
        }
        if (db_flags.refine_eta) {
            for (int i = 0; i < 3 * Npix_to_model; i++) {
                d_image.eta[i] = m_cu_d_eta_images[i];
                d2_image.eta[i] = m_cu_d2_eta_images[i];
            }
        }
        if (std::count(db_flags.refine_Ncells.begin(), db_flags.refine_Ncells.end(), true) > 0 ||
            db_flags.refine_Ncells_def) {
            for (int i = 0; i < 6 * Npix_to_model; i++) {
                d_image.Ncells[i] = m_cu_d_Ncells_images[i];
                d2_image.Ncells[i] = m_cu_d2_Ncells_images[i];
            }
        }
        if (db_flags.refine_diffuse) {
            for (int i = 0; i < 3 * Npix_to_model; i++) {
                d_image.diffuse_gamma[i] = m_cu_d_diffuse_gamma_images[i];
                d_image.diffuse_sigma[i] = m_cu_d_diffuse_sigma_images[i];
            }
        }
        if (std::count(db_flags.refine_Bmat.begin(), db_flags.refine_Bmat.end(), true) > 0) {
            for (int i = 0; i < 6 * Npix_to_model; i++) {
                d_image.Bmat[i] = m_cu_d_Bmat_images[i];
                d2_image.Bmat[i] = m_cu_d2_Bmat_images[i];
            }
        }
        if (std::count(db_flags.refine_lambda.begin(), db_flags.refine_lambda.end(), true) > 0) {
            for (int i = 0; i < 2 * Npix_to_model; i++)
                d_image.lambda[i] = m_cu_d_lambda_images[i];
        }

        if (db_flags.refine_fp_fdp) {
            for (int i = 0; i < 2 * Npix_to_model; i++)
                d_image.fp_fdp[i] = m_cu_d_fp_fdp_images[i];
        }

        gettimeofday(&t2, 0);
        time = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000.0;
        if (TIMERS.recording)
            TIMERS.cuda_copy_from_dev += time;
        if (db_flags.verbose > 1)
            printf("TIME SPENT COPYING BACK :  %3.10f ms \n", time);
        error_msg(cudaGetLastError(), "After copy to host");
        */
}
