#include <cstdio>
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

    kokkos_detector local_det(db_det);  
    kokkos_beam local_beam(db_beam);
    kokkos_crystal local_cryst(db_cryst);

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
        // gpuErr(cudaMallocManaged((void **)&m_cu_UMATS_RXYZ,
        // db_cryst.UMATS_RXYZ.size()*sizeof(MAT3))); 
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

        Kokkos::parallel_for(
            "normalize incident vector", source_count, KOKKOS_LAMBDA(const int& i) {
                VEC3 incident{m_source_X(i), m_source_Y(i), m_source_Z(i)};
                incident.normalize();
                m_source_X(i) = incident.x_val();
                m_source_Y(i) = incident.y_val();
                m_source_Z(i) = incident.z_val();
            });

        if (db_flags.verbose > 1)
            printf("H2D sources\n");
    }
    //  END sources

    //  UMATS
    if (db_cu_flags.update_umats || ALLOC || FORCE_COPY) {
        kokkostbx::transfer_vector2kokkos(m_UMATS, db_cryst.UMATS);
        kokkostbx::transfer_vector2kokkos(m_UMATS_RXYZ, db_cryst.UMATS_RXYZ);
        kokkostbx::transfer_vector2kokkos(m_UMATS_RXYZ_prime, db_cryst.UMATS_RXYZ_prime);
        kokkostbx::transfer_vector2kokkos(m_UMATS_RXYZ_dbl_prime, db_cryst.UMATS_RXYZ_dbl_prime);
        if (db_flags.verbose > 1)
            printf("H2D Done copying Umats\n");
    }
    //  END UMATS

    if (db_cu_flags.update_umats || ALLOC || FORCE_COPY) {
        MAT3 Amat_init = db_cryst.eig_U.dot(db_cryst.eig_B);
        Amat_init *= 1e10;
        Amat_init = Amat_init.dot(db_cryst.eig_O.transpose());

        std::vector<MAT3> AMATS(db_cryst.UMATS_RXYZ);
        for (int i = 0; i < AMATS.size(); ++i) {
            AMATS[i] = AMATS[i].dot(Amat_init).transpose();
        }
        kokkostbx::transfer_vector2kokkos(m_AMATS, AMATS);
        if (db_flags.verbose > 1)
            printf("H2D Done copying Amats\n");
    }

    //  BMATS
    if (db_cu_flags.update_dB_mats || ALLOC || FORCE_COPY) {
        kokkostbx::transfer_vector2kokkos(m_dB_Mats, db_cryst.dB_Mats);
        kokkostbx::transfer_vector2kokkos(m_dB2_Mats, db_cryst.dB2_Mats);
        if (db_flags.verbose > 1)
            printf("H2D Done copying dB_Mats\n");
    }
    //  END BMATS

    //  ROT MATS
    if (db_cu_flags.update_rotmats || ALLOC || FORCE_COPY) {
        kokkostbx::transfer_vector2kokkos(m_RotMats, db_cryst.RotMats);
        kokkostbx::transfer_vector2kokkos(m_dRotMats, db_cryst.dRotMats);
        kokkostbx::transfer_vector2kokkos(m_d2RotMats, db_cryst.d2RotMats);
        if (db_flags.verbose > 1)
            printf("H2D Done copying rotmats\n");
    }
    //  END ROT MATS

    //  DETECTOR VECTORS
    if (db_cu_flags.update_detector || ALLOC || FORCE_COPY) {
        kokkostbx::transfer_vector2kokkos(m_fdet_vectors, db_det.fdet_vectors);
        kokkostbx::transfer_vector2kokkos(m_sdet_vectors, db_det.sdet_vectors);
        kokkostbx::transfer_vector2kokkos(m_odet_vectors, db_det.odet_vectors);
        kokkostbx::transfer_vector2kokkos(m_pix0_vectors, db_det.pix0_vectors);
        kokkostbx::transfer_vector2kokkos(m_close_distances, db_det.close_distances);
        if (db_flags.verbose > 1)
            printf("H2D Done copying detector vectors\n");
    }
    //  END  DETECTOR VECTORS

    if (ALLOC || FORCE_COPY) {
        kokkostbx::transfer_vector2kokkos(m_nominal_hkl, db_cryst.nominal_hkl);
        kokkostbx::transfer_vector2kokkos(m_atom_data, db_cryst.atom_data);
        if (db_flags.verbose > 1)
            printf("H2D Done copying atom data\n");

        kokkostbx::transfer_vector2kokkos(m_fpfdp, db_cryst.fpfdp);
        kokkostbx::transfer_vector2kokkos(m_fpfdp_derivs, db_cryst.fpfdp_derivs);
        if (db_flags.verbose > 1)
            printf("H2D Done copying fprime and fdblprime\n");
    }

    //  BEGIN REFINEMENT FLAGS
    if (db_cu_flags.update_refine_flags || ALLOC || FORCE_COPY) {
        kokkostbx::transfer_vector2kokkos(m_refine_Umat, db_flags.refine_Umat);
        kokkostbx::transfer_vector2kokkos(m_refine_Ncells, db_flags.refine_Ncells);
        kokkostbx::transfer_vector2kokkos(m_refine_panel_origin, db_flags.refine_panel_origin);
        kokkostbx::transfer_vector2kokkos(m_refine_panel_rot, db_flags.refine_panel_rot);
        kokkostbx::transfer_vector2kokkos(m_refine_lambda, db_flags.refine_lambda);
        kokkostbx::transfer_vector2kokkos(m_refine_Bmat, db_flags.refine_Bmat);
        if (db_flags.verbose > 1)
            printf("H2D Done copying refinement flags\n");
    }
    //  END REFINEMENT FLAGS

    //  BEGIN Fhkl
    if (db_cu_flags.update_Fhkl || ALLOC || FORCE_COPY) {
        kokkostbx::transfer_vector2kokkos(m_Fhkl, db_cryst.FhklLinear);
        if (db_flags.complex_miller) {
            kokkostbx::transfer_vector2kokkos(m_Fhkl2, db_cryst.Fhkl2Linear);
        }
        if (db_flags.verbose > 1)
            printf("H2D Done copying step Fhkl\n");
    }
    //  END Fhkl

    //  BEGIN panel derivative vecs
    if (db_cu_flags.update_panel_deriv_vecs || ALLOC || FORCE_COPY) {
        kokkostbx::transfer_vector2kokkos(m_dF_vecs, db_det.dF_vecs);
        kokkostbx::transfer_vector2kokkos(m_dS_vecs, db_det.dS_vecs);
        if (db_flags.verbose > 1)
            printf("H2D Done copying step panel derivative vectors\n");
    }
    //  END panel derivative vecs

    //  BEGIN panels fasts slows
    if (db_cu_flags.update_panels_fasts_slows || ALLOC || FORCE_COPY) {
        kokkostbx::transfer_vector2kokkos(m_panels_fasts_slows, panels_fasts_slows);
        if (db_flags.verbose > 1)
            printf("H2D Done copying panels_fasts_slows\n");
    }
    //  END panels fasts slows

    gettimeofday(&t2, 0);
    time = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000.0;
    if (TIMERS.recording)
        TIMERS.cuda_copy_to_dev += time;
    if (db_flags.verbose > 1)
        printf("TIME SPENT COPYING DATA HOST->DEV:  %3.10f ms \n", time);

    m_device_is_allocated = true;
    ::Kokkos::fence("after copy to device");

    gettimeofday(&t1, 0);

    int Npanels = db_det.fdet_vectors.size() / 3;
    int num_atoms = db_cryst.atom_data.size() / 5;
    // note cannot use atom data if fpfdp is 0, make this cleaner
    if (db_cryst.fpfdp.size() == 0) {
        num_atoms = 0;
    }
    // int sm_size = number_of_sources*5*sizeof(CUDAREAL);
    // gpu_sum_over_steps<<<numblocks, blocksize, sm_size >>>(
    bool aniso_eta = db_cryst.UMATS_RXYZ.size() != db_cryst.UMATS_RXYZ_prime.size();
    bool use_nominal_hkl = !db_cryst.nominal_hkl.empty();
/*    kokkos_sum_over_steps<<<numblocks, blocksize>>>(
        Npix_to_model, m_panels_fasts_slows, m_floatimage, m_wavelenimage, m_d_Umat_images,
        m_d2_Umat_images, m_d_Bmat_images, m_d2_Bmat_images, m_d_Ncells_images, m_d2_Ncells_images,
        m_d_fcell_images, m_d2_fcell_images, m_d_eta_images, m_d2_eta_images, m_d_lambda_images,
        m_d2_lambda_images, m_d_panel_rot_images, m_d2_panel_rot_images, m_d_panel_orig_images,
        m_d2_panel_orig_images, m_d_fp_fdp_images, db_steps.Nsteps, db_flags.printout_fpixel,
        db_flags.printout_spixel, db_flags.printout, db_cryst.default_F, db_det.oversample,
        db_flags.oversample_omega, db_det.subpixel_size, db_det.pixel_size,
        db_det.detector_thickstep, db_det.detector_thick, m_close_distances,
        db_det.detector_attnlen, db_det.detector_thicksteps, db_beam.number_of_sources,
        db_cryst.phisteps, db_cryst.UMATS.size(), db_flags.use_lambda_coefficients, db_beam.lambda0,
        db_beam.lambda1, db_cryst.eig_U, db_cryst.eig_O, db_cryst.eig_B, db_cryst.RXYZ, m_dF_vecs,
        m_dS_vecs, m_UMATS_RXYZ, m_UMATS_RXYZ_prime, m_UMATS_RXYZ_dbl_prime, m_RotMats, m_dRotMats,
        m_d2RotMats, m_UMATS, m_dB_Mats, m_dB2_Mats, m_AMATS, m_source_X, m_source_Y, m_source_Z,
        m_source_lambda, m_source_I, db_beam.kahn_factor, db_cryst.Na, db_cryst.Nb, db_cryst.Nc,
        db_cryst.Nd, db_cryst.Ne, db_cryst.Nf, db_cryst.phi0, db_cryst.phistep,
        db_cryst.spindle_vec, db_beam.polarization_axis, db_cryst.h_range, db_cryst.k_range,
        db_cryst.l_range, db_cryst.h_max, db_cryst.h_min, db_cryst.k_max, db_cryst.k_min,
        db_cryst.l_max, db_cryst.l_min, db_cryst.dmin, db_cryst.fudge, db_flags.complex_miller,
        db_flags.verbose, db_flags.only_save_omega_kahn, db_flags.isotropic_ncells,
        db_flags.compute_curvatures, m_Fhkl, m_Fhkl2, m_refine_Bmat, m_refine_Ncells,
        db_flags.refine_Ncells_def, m_refine_panel_origin, m_refine_panel_rot,
        db_flags.refine_fcell, m_refine_lambda, db_flags.refine_eta, m_refine_Umat, m_fdet_vectors,
        m_sdet_vectors, m_odet_vectors, m_pix0_vectors, db_flags.nopolar, db_flags.point_pixel,
        db_beam.fluence, db_cryst.r_e_sqr, db_cryst.spot_scale, Npanels, aniso_eta,
        db_flags.no_Nabc_scale, m_fpfdp, m_fpfdp_derivs, m_atom_data, num_atoms,
        db_flags.refine_fp_fdp, m_nominal_hkl, use_nominal_hkl, db_cryst.anisoU, db_cryst.anisoG,
        db_flags.use_diffuse, m_d_diffuse_gamma_images, m_d_diffuse_sigma_images,
        db_flags.refine_diffuse, db_flags.gamma_miller_units, db_flags.refine_Icell,
        db_flags.wavelength_img);
*/
    ::Kokkos::fence("after kernel call");

    if (db_flags.verbose > 1)
        printf("KERNEL_COMPLETE gpu_sum_over_steps\n");
    gettimeofday(&t2, 0);
    time = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000.0;
    if (TIMERS.recording)
        TIMERS.cuda_kernel += time;
    if (db_flags.verbose > 1)
        printf("TIME SPENT(KERNEL):  %3.10f ms \n", time);

    gettimeofday(&t1, 0);
    //  COPY BACK FROM DEVICE
    kokkostbx::transfer_kokkos2vector(floatimage, m_floatimage);

    if (db_flags.wavelength_img) {
        kokkostbx::transfer_kokkos2vector(d_image.wavelength, m_wavelenimage);
    }
    if (db_flags.refine_fcell) {
        kokkostbx::transfer_kokkos2vector(d_image.fcell, m_d_fcell_images);
        kokkostbx::transfer_kokkos2vector(d2_image.fcell, m_d2_fcell_images);
    }
    if (std::count(db_flags.refine_Umat.begin(), db_flags.refine_Umat.end(), true) > 0) {
        kokkostbx::transfer_kokkos2vector(d_image.Umat, m_d_Umat_images);
        kokkostbx::transfer_kokkos2vector(d2_image.Umat, m_d2_Umat_images);
    }
    if (std::count(db_flags.refine_panel_rot.begin(), db_flags.refine_panel_rot.end(), true) > 0) {
        kokkostbx::transfer_kokkos2vector(d_image.panel_rot, m_d_panel_rot_images);
    }
    if (std::count(db_flags.refine_panel_origin.begin(), db_flags.refine_panel_origin.end(), true) >
        0) {
        kokkostbx::transfer_kokkos2vector(d_image.panel_orig, m_d_panel_orig_images);
    }
    if (db_flags.refine_eta) {
        kokkostbx::transfer_kokkos2vector(d_image.eta, m_d_eta_images);
        kokkostbx::transfer_kokkos2vector(d2_image.eta, m_d2_eta_images);
    }
    if (std::count(db_flags.refine_Ncells.begin(), db_flags.refine_Ncells.end(), true) > 0 ||
        db_flags.refine_Ncells_def) {
        kokkostbx::transfer_kokkos2vector(d_image.Ncells, m_d_Ncells_images);
        kokkostbx::transfer_kokkos2vector(d2_image.Ncells, m_d2_Ncells_images);
    }
    if (db_flags.refine_diffuse) {
        kokkostbx::transfer_kokkos2vector(d_image.diffuse_gamma, m_d_diffuse_gamma_images);
        kokkostbx::transfer_kokkos2vector(d2_image.diffuse_sigma, m_d_diffuse_sigma_images);
    }
    if (std::count(db_flags.refine_Bmat.begin(), db_flags.refine_Bmat.end(), true) > 0) {
        kokkostbx::transfer_kokkos2vector(d_image.Bmat, m_d_Bmat_images);
        kokkostbx::transfer_kokkos2vector(d2_image.Bmat, m_d2_Bmat_images);
    }
    if (std::count(db_flags.refine_lambda.begin(), db_flags.refine_lambda.end(), true) > 0) {
        kokkostbx::transfer_kokkos2vector(d_image.lambda, m_d_lambda_images);
    }
    if (db_flags.refine_fp_fdp) {
        kokkostbx::transfer_kokkos2vector(d_image.fp_fdp, m_d_fp_fdp_images);
    }

    gettimeofday(&t2, 0);
    time = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000.0;
    if (TIMERS.recording)
        TIMERS.cuda_copy_from_dev += time;
    if (db_flags.verbose > 1)
        printf("TIME SPENT COPYING BACK :  %3.10f ms \n", time);
    ::Kokkos::fence("After copy to host");
}
