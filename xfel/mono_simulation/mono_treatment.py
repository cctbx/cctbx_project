from __future__ import absolute_import, division, print_function
from six.moves import range
from cctbx.array_family import flex
import math
from scitbx.matrix import col,sqr
from scitbx.lstbx import normal_eqns
from six.moves import zip
from six import unichr as chr

class refinement_base(object):

  def __init__(OO,self,use_inverse_beam=False):

    OO.parent = self # OO.parent is an instance of the legacy IntegrationMetaProcedure class
    from xfel.mono_simulation import bandpass_gaussian
    from rstbx.bandpass import parameters_bp3

    #take needed parameters from parent
    pxlsz = self.pixel_size # mm/pixel

    detector_origin = col(( -self.inputai.xbeam(),
                            -self.inputai.ybeam(),
                             0.))
    #OO.space_group = self.inputpd["symmetry"].space_group()   #comment this back in as needed for refinement
    indices = flex.miller_index([self.hkllist[pair["pred"]] for pair in self.indexed_pairs])
    OO.reserve_indices = indices
    OO.input_orientation = self.inputai.getOrientation()
    OO.central_wavelength_ang = self.inputai.wavelength
    incident_beam = col((0.,0.,-1.))
    if use_inverse_beam: incident_beam*=-1.

    parameters = parameters_bp3(
       indices=indices, orientation=OO.input_orientation,
       incident_beam=incident_beam,
       packed_tophat=col((1.,1.,0.)),
       detector_normal=col((0.,0.,-1.)),
       detector_fast=col((0.,1.,0.)),detector_slow=col((1.,0.,0.)),
       pixel_size=col((pxlsz,pxlsz,0)),
       pixel_offset=col((0.,0.,0.0)),
       distance=self.inputai.distance(),
       detector_origin=detector_origin
    )
    OO.ucbp3 = bandpass_gaussian(parameters=parameters)

    if "horizons_phil" in OO.parent.__dict__:
      the_tiles = OO.parent.imagefiles.images[OO.parent.image_number
      ].get_tile_manager(OO.parent.horizons_phil
      ).effective_tiling_as_flex_int(
      reapply_peripheral_margin=True,encode_inactive_as_zeroes=True)
      OO.ucbp3.set_active_areas( the_tiles )
    else:
      OO.ucbp3.set_active_areas( [0,0,1700,1700] )
    integration_signal_penetration=0.0 # easier to calculate distance derivatives

    OO.ucbp3.set_sensor_model( thickness_mm = 0.5, mu_rho = 8.36644, # CS_PAD detector at 1.3 Angstrom
      signal_penetration = integration_signal_penetration)

    # test for horizons_phil simply skips the subpixel correction for initial labelit indexing
    if "horizons_phil" in OO.parent.__dict__:
      if OO.parent.horizons_phil.integration.subpixel_joint_model.translations is not None:
        "Subpixel corrections: using joint-refined translation + rotation"
        T = OO.parent.horizons_phil.integration.subpixel_joint_model.translations
        import copy
        resortedT = copy.copy(T)
        for tt in range(0,len(T),2):
          resortedT[tt] = T[tt+1]
          resortedT[tt+1] = T[tt]
        OO.ucbp3.set_subpixel(
            translations = resortedT, rotations_deg = flex.double(
             OO.parent.horizons_phil.integration.subpixel_joint_model.rotations)
          )
    else:
      pass; "Subpixel corrections: none used"

    half_mosaicity_rad = (self.inputai.getMosaicity()/2.) * math.pi/180.
    OO.ucbp3.set_mosaicity(half_mosaicity_rad)
    OO.ucbp3.set_bandpass(OO.central_wavelength_ang - 0.000001, OO.central_wavelength_ang + 0.000001)
    OO.ucbp3.set_domain_size(280. * 17.) # for Holton psI simulation; probably doesn't detract from general case

class refinement(refinement_base):
  grid = range(-20,20)

  def __init__(OO,self,use_inverse_beam=False,mosaic_refinement_target="LSQ",pvr_fix=True):
    OO.mosaic_refinement_target = mosaic_refinement_target # least squares or max-likelihood
    OO.pvr_fix = pvr_fix
    refinement_base.__init__(OO,self,use_inverse_beam)

  def contour_plot_DEPRECATED_DOES_NOT_APPLY_SUBPIXEL_METROLOGY(OO):
      self = OO.parent
      # see if I can reproduce the predicted positions
      pxlsz = self.pixel_size # mm/pixel
      SIGN = -1.

      # get a simple rmsd obs vs predicted
      sumsq = 0.
      nspot = 0
      for pair in self.indexed_pairs:
        deltax = self.spots[pair["spot"]].ctr_mass_x() - self.predicted[pair["pred"]][0]/pxlsz
        deltay = self.spots[pair["spot"]].ctr_mass_y() - self.predicted[pair["pred"]][1]/pxlsz
        sumsq += deltax*deltax + deltay*deltay
        nspot+=1
      print("RMSD obs vs pred in pixels: %7.2f"%(math.sqrt(sumsq/nspot)))

      excursi = flex.double()
      rmsdpos = flex.double()

      for irotx in OO.grid:
        rotx = (0.02 * irotx) * math.pi/180.
        for iroty in OO.grid:
          roty = (0.02 * iroty) * math.pi/180.

          #Not necessary to apply the 3 offset rotations; they have apparently
          #  been applied already.  rotz (0,0,1) is the direct beam
          effective_orientation = OO.input_orientation.rotate_thru((1,0,0),rotx
           ).rotate_thru((0,1,0),roty
           ).rotate_thru((0,0,1),0.0)
          OO.ucbp3.set_orientation(effective_orientation)

          OO.ucbp3.gaussian_fast_slow()
          mean_position = OO.ucbp3.mean_position

          print(mean_position)
          sumsq = 0.
          nspot = 0
          for pair in self.indexed_pairs:
            deltax = mean_position[nspot][1] - self.predicted[pair["pred"]][0]/pxlsz
            deltay = mean_position[nspot][0] - self.predicted[pair["pred"]][1]/pxlsz
            sumsq += deltax*deltax + deltay*deltay
            nspot+=1
          print("RMSD markmodel vs rossmanpred in pixels: %7.2f"%(math.sqrt(sumsq/nspot)))

          #from matplotlib import pyplot as plt
          #plt.plot([mpos[0] for mpos in mean_position],[mpos[1] for mpos in mean_position],"r+")
          #plt.plot([self.predicted[pair["pred"]][1]/pxlsz for pair in indexed_pairs],
          #         [self.predicted[pair["pred"]][0]/pxlsz for pair in indexed_pairs], "b.")
          #plt.show()

          sumsq = 0.
          nspot = 0
          for pair in self.indexed_pairs:
            deltax = self.spots[pair["spot"]].ctr_mass_x() - mean_position[nspot][1]
            deltay = self.spots[pair["spot"]].ctr_mass_y() - mean_position[nspot][0]
            sumsq += deltax*deltax + deltay*deltay
            nspot+=1
          rmsdposition = math.sqrt(sumsq/nspot)
          rmsdpos.append(rmsdposition)
          print("RMSD obs vs markmodel in pixels: %8.4f"%(rmsdposition))

          excursions = flex.double(
            [OO.ucbp3.simple_forward_calculation_spot_position(
            wavelength = OO.central_wavelength_ang,
            observation_no = obsno).rotax_excursion_rad*180./math.pi
            for obsno in range(len(self.indexed_pairs))])

          rmsdexc = math.sqrt(flex.mean(excursions*excursions))
          excursi.append(rmsdexc)
          print("rotx %7.2f roty %7.2f degrees, RMSD excursion %7.3f degrees"%(
          (0.02 * irotx),(0.02 * iroty), rmsdexc))
      return excursi,rmsdpos

  def per_frame_helper_factory(OO):

      class per_frame_helper(normal_eqns.non_linear_ls, normal_eqns.non_linear_ls_mixin):
        def __init__(pfh):
          super(per_frame_helper, pfh).__init__(n_parameters=2)

          pfh.x_0 = flex.double((0.,0.))
          pfh.restart()

        def restart(pfh):
          pfh.x = pfh.x_0.deep_copy()
          pfh.old_x = None

        def step_forward(pfh):
          pfh.old_x = pfh.x.deep_copy()
          pfh.x += pfh.step()

        def step_backward(pfh):
          assert pfh.old_x is not None
          pfh.x, pfh.old_x = pfh.old_x, None

        def parameter_vector_norm(pfh):
          return pfh.x.norm()

        def build_up(pfh, objective_only=False):
          if OO.pvr_fix:
            residuals = pfh.fvec_callable_pvr(pfh.x)
          else:
            residuals = pfh.fvec_callable_NOT_USED_AFTER_BUGFIX(pfh.x)

          pfh.reset()
          if objective_only:
            pfh.add_residuals(residuals, weights=None)
          else:
            grad_r = pfh.jacobian_callable(pfh.x)
            jacobian = flex.double(
              flex.grid(len(OO.parent.indexed_pairs), pfh.n_parameters))
            for j, der_r in enumerate(grad_r):
              jacobian.matrix_paste_column_in_place(der_r,j)
            pfh.add_equations(residuals, jacobian, weights=None)

        def fvec_callable_pvr(pfh,current_values):
          rotx = current_values[0]
          roty = current_values[1]
          effective_orientation = OO.input_orientation.rotate_thru((1,0,0),rotx
           ).rotate_thru((0,1,0),roty
           ).rotate_thru((0,0,1),0.0)
          OO.ucbp3.set_orientation(effective_orientation)
          pfh.last_set_orientation = effective_orientation

          OO.ucbp3.gaussian_fast_slow()

          excursions = flex.double(
            [OO.ucbp3.simple_forward_calculation_spot_position(
            wavelength = OO.central_wavelength_ang,
            observation_no = obsno).rotax_excursion_rad_pvr/(2.*math.pi)
            for obsno in range(len(OO.parent.indexed_pairs))])

          degrees = 360.*excursions
          rmsdexc = math.sqrt(flex.mean(degrees*degrees))
          #print "rotx %7.3f roty %7.3f degrees, -PVR excursion %7.3f degrees"%(
          #(rotx * 180./math.pi),(roty * 180./math.pi), rmsdexc)
          # Note.  Luc Bourhis wants scale to be from 0 to 1. So instead of
          # returning on scale of degrees, use radians/(2*pi)
          # The parameters rotx roty are still expressed in radians
          return excursions

        def jacobian_callable(pfh,current_values):
          rotx = current_values[0]
          roty = current_values[1]
          from scitbx.matrix import sqr
          Ai = sqr(OO.input_orientation.reciprocal_matrix())
          Rx = col((1,0,0)).axis_and_angle_as_r3_rotation_matrix(rotx)
          Ry = col((0,1,0)).axis_and_angle_as_r3_rotation_matrix(roty)
          Rz = col((0,0,1)).axis_and_angle_as_r3_rotation_matrix(0.0)
          dRx_drotx = col((1,0,0)).axis_and_angle_as_r3_derivative_wrt_angle(rotx)
          dRy_droty = col((0,1,0)).axis_and_angle_as_r3_derivative_wrt_angle(roty)
          dA_drotx = Rz * Ry * dRx_drotx * Ai
          dA_droty = Rz * dRy_droty * Rx * Ai

          dexc_drotx = [
            OO.ucbp3.simple_part_excursion_part_rotxy(
            wavelength = OO.central_wavelength_ang,
            observation_no = obsno,
            dA_drotxy = dA_drotx)
            for obsno in range(len(OO.parent.indexed_pairs))]

          dexc_droty = [
            OO.ucbp3.simple_part_excursion_part_rotxy(
            wavelength = OO.central_wavelength_ang,
            observation_no = obsno,
            dA_drotxy = dA_droty)
            for obsno in range(len(OO.parent.indexed_pairs))]
          return flex.double(dexc_drotx)/(2.*math.pi), flex.double(dexc_droty)/(2.*math.pi)

      value = per_frame_helper()
      return value

  def refine_rotx_roty2(OO,enable_rotational_target=True):

      helper = OO.per_frame_helper_factory()
      helper.restart()

      if enable_rotational_target:
        print("Trying least squares minimization of excursions", end=' ')
        from scitbx.lstbx import normal_eqns_solving
        iterations = normal_eqns_solving.naive_iterations(
          non_linear_ls = helper,
          gradient_threshold = 1.E-10)

      results =  helper.x

      print("with %d reflections"%len(OO.parent.indexed_pairs), end=' ')
      print("result %6.2f degrees"%(results[1]*180./math.pi), end=' ')
      print("result %6.2f degrees"%(results[0]*180./math.pi))

      if False: # Excursion histogram
        print("The input mosaicity is %7.3f deg full width"%OO.parent.inputai.getMosaicity())
        # final histogram
        if OO.pvr_fix:
          final = 360.* helper.fvec_callable_pvr(results)
        else:
          final = 360.* helper.fvec_callable_NOT_USED_AFTER_BUGFIX(results)

        rmsdexc = math.sqrt(flex.mean(final*final))
        from matplotlib import pyplot as plt
        nbins = len(final)//20
        n,bins,patches = plt.hist(final,
          nbins, normed=0, facecolor="orange", alpha=0.75)
        plt.xlabel("Rotation on e1 axis, rmsd %7.3f deg"%rmsdexc)
        plt.title("Histogram of cctbx.xfel misorientation")
        plt.axis([-0.5,0.5,0,100])
        plt.plot([rmsdexc],[18],"b|")
        plt.show()

      # Determine optimal mosaicity and domain size model (monochromatic)
      if OO.pvr_fix:
        final = 360.* helper.fvec_callable_pvr(results)
      else:
        final = 360.* helper.fvec_callable_NOT_USED_AFTER_BUGFIX(results)
      #Guard against misindexing -- seen in simulated data, with zone nearly perfectly aligned
      guard_stats = flex.max(final), flex.min(final)
      if False and REMOVETEST_KILLING_LEGITIMATE_EXCURSIONS (guard_stats[0] > 2.0 or guard_stats[1] < -2.0):
        raise Exception("Misindexing diagnosed by meaningless excursion angle (bandpass_gaussian model)");
      print("The mean excursion is %7.3f degrees"%(flex.mean(final)))

      two_thetas = helper.last_set_orientation.unit_cell().two_theta(OO.reserve_indices,OO.central_wavelength_ang,deg=True)
      dspacings = helper.last_set_orientation.unit_cell().d(OO.reserve_indices)
      dspace_sq = dspacings * dspacings
      excursion_rad = final * math.pi/ 180.

      #  First -- try to get a reasonable envelope for the observed excursions.
          ## minimum of three regions; maximum of 50 measurements in each bin
      print("fitting parameters on %d spots"%len(excursion_rad))
      n_bins = min(max(3, len(excursion_rad)//25),50)
      bin_sz = len(excursion_rad)//n_bins
      print("nbins",n_bins,"bin_sz",bin_sz)
      order = flex.sort_permutation(two_thetas)
      two_thetas_env = flex.double()
      dspacings_env = flex.double()
      excursion_rads_env = flex.double()
      for x in range(0,n_bins):
        subset = order[x*bin_sz:(x+1)*bin_sz]
        two_thetas_env.append( flex.mean(two_thetas.select(subset)) )
        dspacings_env.append( flex.mean(dspacings.select(subset)))
        excursion_rads_env.append( flex.max( flex.abs( excursion_rad.select(subset))))

      #  Second -- parameter fit
          ## solve the normal equations
      sum_inv_u_sq = flex.sum(dspacings_env * dspacings_env)
      sum_inv_u    = flex.sum(dspacings_env)
      sum_te_u     = flex.sum(dspacings_env * excursion_rads_env)
      sum_te       = flex.sum(excursion_rads_env)
      Normal_Mat   = sqr((sum_inv_u_sq, sum_inv_u, sum_inv_u, len(dspacings_env)))
      Vector       = col((sum_te_u, sum_te))
      solution     = Normal_Mat.inverse() * Vector
      s_ang = 1./(2*solution[0])
      print("Best LSQ fit Scheerer domain size is %9.2f ang"%(
        s_ang))
      tan_phi_rad = helper.last_set_orientation.unit_cell().d(OO.reserve_indices) / (2. * s_ang)
      tan_phi_deg = tan_phi_rad * 180./math.pi
      k_degrees = solution[1]* 180./math.pi
      print("The LSQ full mosaicity is %8.5f deg; half-mosaicity %9.5f"%(2*k_degrees, k_degrees))
      tan_outer_deg = tan_phi_deg + k_degrees

      if OO.mosaic_refinement_target=="ML":
        from xfel.mono_simulation.max_like import minimizer
        print("input", s_ang,2. * solution[1]*180/math.pi)
        # coerce the estimates to be positive for max-likelihood
        lower_limit_domain_size = math.pow(
         helper.last_set_orientation.unit_cell().volume(),
         1./3.)*20 # 10-unit cell block size minimum reasonable domain

        d_estimate = max(s_ang, lower_limit_domain_size)
        M = minimizer(d_i = dspacings, psi_i = excursion_rad, eta_rad = abs(2. * solution[1]),
                      Deff = d_estimate)
        print("output",1./M.x[0], M.x[1]*180./math.pi)
        tan_phi_rad_ML = helper.last_set_orientation.unit_cell().d(OO.reserve_indices) / (2. / M.x[0])
        tan_phi_deg_ML = tan_phi_rad_ML * 180./math.pi
        # bugfix: Need factor of 0.5 because the plot shows half mosaicity (displacement from the center point defined as zero)
        tan_outer_deg_ML = tan_phi_deg_ML + 0.5*M.x[1]*180./math.pi

      if OO.parent.horizons_phil.integration.mosaic.enable_polychromatic:
        # add code here to perform polychromatic modeling.
        """
        get miller indices DONE
        get model-predicted mono-wavelength centroid S1 vectors
        back-convert S1vec, with mono-wavelength, to detector-plane position, factoring in subpixel correction
        compare with spot centroid measured position
        compare with locus of bodypixels
        """
        print(list(OO.reserve_indices))
        print(len(OO.reserve_indices), len(two_thetas))
        positions = [
              OO.ucbp3.simple_forward_calculation_spot_position(
              wavelength = OO.central_wavelength_ang,
              observation_no = obsno).position
              for obsno in range(len(OO.parent.indexed_pairs))]
        print(len(positions))
        print(positions) # model-predicted positions
        print(len(OO.parent.spots))
        print(OO.parent.indexed_pairs)
        print(OO.parent.spots)
        print(len(OO.parent.spots))
        meas_spots = [OO.parent.spots[pair["spot"]] for pair in OO.parent.indexed_pairs]
  #      for xspot in meas_spots:
  #        xspot.ctr_mass_x(),xspot.ctr_mass_y()
  #        xspot.max_pxl_x()
  #        xspot.bodypixels
  #        xspot.ctr_mass_x()

        # Do some work to calculate an rmsd
        diff_vecs = flex.vec3_double()
        for p,xspot in zip(positions, meas_spots):
          diff_vecs.append((p[0]-xspot.ctr_mass_y(), p[1]-xspot.ctr_mass_x(), 0.0))
        # could use diff_vecs.rms_length()
        diff_vecs_sq = diff_vecs.dot(diff_vecs)
        mean_diff_vec_sq = flex.mean(diff_vecs_sq)
        rmsd = math.sqrt(mean_diff_vec_sq)
        print("mean obs-pred diff vec on %d spots is %6.2f pixels"%(len(positions),rmsd))

        positions_to_fictitious = [
              OO.ucbp3.simple_forward_calculation_spot_position(
              wavelength = OO.central_wavelength_ang,
              observation_no = obsno).position_to_fictitious
              for obsno in range(len(OO.parent.indexed_pairs))]
        # Do some work to calculate an rmsd
        diff_vecs = flex.vec3_double()
        for p,xspot in zip(positions_to_fictitious, meas_spots):
          diff_vecs.append((p[0]-xspot.ctr_mass_y(), p[1]-xspot.ctr_mass_x(), 0.0))
        rmsd = diff_vecs.rms_length()
        print("mean obs-pred_to_fictitious diff vec on %d spots is %6.2f pixels"%(len(positions),rmsd))

        """
        actually, it might be better if the entire set of experimental observations
        is transformed into the ideal detector plane, for the purposes of poly_treatment.


        start here.  Now it would be good to actually implement probability of observing a body pixel given the model.
        We have everything needed right here.
        """
        if OO.parent.horizons_phil.integration.mosaic.enable_AD14F7B:
          # Image plot: obs and predicted positions + bodypixels
          from matplotlib import pyplot as plt
          plt.plot( [p[0] for p in positions_to_fictitious], [p[1] for p in positions_to_fictitious], "r.")
          plt.plot( [xspot.ctr_mass_y() for xspot in meas_spots],
                    [xspot.ctr_mass_x() for xspot in meas_spots], "g.")
          bodypx = []
          for xspot in meas_spots:
            for body in xspot.bodypixels:
              bodypx.append(body)
          plt.plot( [b.y for b in bodypx], [b.x for b in bodypx], "b.")
          plt.axes().set_aspect("equal")
          plt.show()

      print("MEAN excursion",flex.mean(final), end=' ')
      if OO.mosaic_refinement_target=="ML":
        print("mosaicity deg FW=",M.x[1]*180./math.pi)
      else:
        print()
      if OO.parent.horizons_phil.integration.mosaic.enable_AD14F7B: # Excursion vs resolution fit
        AD1TF7B_MAX2T = 30.
        AD1TF7B_MAXDP = 1.
        from matplotlib import pyplot as plt
        fig = plt.figure()
        plt.plot(two_thetas, final, "bo")
        mean = flex.mean(final)
        minplot = flex.min(two_thetas)
        plt.plot([0,minplot],[mean,mean],"k-")
        LR = flex.linear_regression(two_thetas, final)
        #LR.show_summary()
        model_y = LR.slope()*two_thetas + LR.y_intercept()
        plt.plot(two_thetas, model_y, "k-")
        print(helper.last_set_orientation.unit_cell())
        #for sdp,tw in zip (dspacings,two_thetas):
          #print sdp,tw
        if OO.mosaic_refinement_target=="ML":
          plt.title("ML: mosaicity FW=%4.2f deg, Dsize=%5.0fA on %d spots"%(M.x[1]*180./math.pi, 2./M.x[0], len(two_thetas)))
          plt.plot(two_thetas, tan_phi_deg_ML, "r.")
          plt.plot(two_thetas, -tan_phi_deg_ML, "r.")
          plt.plot(two_thetas, tan_outer_deg_ML, "g.")
          plt.plot(two_thetas, -tan_outer_deg_ML, "g.")
        else:
          plt.plot(two_thetas_env, excursion_rads_env *180./math.pi, "r|")
          plt.plot(two_thetas_env, -excursion_rads_env *180./math.pi, "r|")
          plt.plot(two_thetas_env, excursion_rads_env *180./math.pi, "r-")
          plt.plot(two_thetas_env, -excursion_rads_env *180./math.pi, "r-")
          plt.plot(two_thetas, tan_phi_deg, "r.")
          plt.plot(two_thetas, -tan_phi_deg, "r.")
          plt.plot(two_thetas, tan_outer_deg, "g.")
          plt.plot(two_thetas, -tan_outer_deg, "g.")
        plt.xlim([0,AD1TF7B_MAX2T])
        plt.ylim([-AD1TF7B_MAXDP,AD1TF7B_MAXDP])
        OO.parent.show_figure(plt,fig,"psi")
        plt.close()

      from xfel.mono_simulation.util import green_curve_area,ewald_proximal_volume
      if OO.mosaic_refinement_target=="ML":
        OO.parent.green_curve_area = green_curve_area(two_thetas, tan_outer_deg_ML)
        OO.parent.inputai.setMosaicity(M.x[1]*180./math.pi) # full width, degrees
        OO.parent.ML_half_mosaicity_deg = M.x[1]*180./(2.*math.pi)
        OO.parent.ML_domain_size_ang = 1./M.x[0]
        OO.parent.ewald_proximal_volume = ewald_proximal_volume(
            wavelength_ang = OO.central_wavelength_ang,
            resolution_cutoff_ang = OO.parent.horizons_phil.integration.mosaic.ewald_proximal_volume_resolution_cutoff,
            domain_size_ang = 1./M.x[0],
            full_mosaicity_rad = M.x[1])
        return results, helper.last_set_orientation,1./M.x[0] # full width domain size, angstroms
      else:
        assert OO.mosaic_refinement_target=="LSQ"
        OO.parent.green_curve_area = green_curve_area(two_thetas, tan_outer_deg)
        OO.parent.inputai.setMosaicity(2*k_degrees) # full width
        OO.parent.ML_half_mosaicity_deg = k_degrees
        OO.parent.ML_domain_size_ang = s_ang
        OO.parent.ewald_proximal_volume = ewald_proximal_volume(
            wavelength_ang = OO.central_wavelength_ang,
            resolution_cutoff_ang = OO.parent.horizons_phil.integration.mosaic.ewald_proximal_volume_resolution_cutoff,
            domain_size_ang = s_ang,
            full_mosaicity_rad = 2*k_degrees*math.pi/180.)
        return results, helper.last_set_orientation,s_ang # full width domain size, angstroms

  def show_plot(OO,excursi,rmsdpos,minimum):
      excursi.reshape(flex.grid(len(OO.grid), len(OO.grid)))
      rmsdpos.reshape(flex.grid(len(OO.grid), len(OO.grid)))

      from matplotlib import pyplot as plt
      plt.figure()
      CS = plt.contour([i*0.02 for i in OO.grid],[i*0.02 for i in OO.grid], excursi.as_numpy_array())
      plt.clabel(CS, inline=1, fontsize=10, fmt="%6.3f"+chr(176))
      plt.plot([minimum[1]*180./math.pi],[minimum[0]*180./math.pi], "r+")
      plt.title("Rms rotational excursion to reflection condition, degrees")
      plt.axes().set_aspect("equal")
      plt.figure()
      CS = plt.contour([i*0.02 for i in OO.grid],[i*0.02 for i in OO.grid], rmsdpos.as_numpy_array())
      plt.clabel(CS, inline=1, fontsize=10, fmt="%7.4f px")
      plt.title("Rms position shift, obs vs. model, pixels")
      plt.axes().set_aspect("equal")
      plt.show()

class refinement2(refinement_base):
  grid = range(-20,20)

  def __init__(OO,self):

    refinement_base.__init__(OO,self)

  def refine_all(OO):

      class per_frame_helper:
        def __init__(pfh):
          pass
          from rstbx.symmetry.constraints.parameter_reduction import symmetrize_reduce_enlarge
          pfh.convert = symmetrize_reduce_enlarge(space_group=OO.space_group)
          pfh.convert.set_orientation(orientation=OO.input_orientation)

        def fvec_callable(pfh,current_values):
          rotz = current_values[0]
          indep = current_values[1:]

          effective_orientation = OO.input_orientation.rotate_thru((0,0,1),rotz)
          pfh.convert.set_orientation(effective_orientation)
          pfh.convert.forward_independent_parameters()
          effective_orientation = pfh.convert.backward_orientation(independent=indep)

          OO.ucbp3.set_orientation(effective_orientation)
          pfh.last_set_orientation = effective_orientation
          OO.ucbp3.gaussian_fast_slow()

          # note the reversal of x & y with obs vs. predicted
          displacements = flex.double(
            [(
              col(
                  OO.ucbp3.simple_forward_calculation_spot_position(
                  wavelength = OO.central_wavelength_ang,
                  observation_no = obsno).position) -
              col(
                  (OO.parent.spots[OO.parent.indexed_pairs[obsno]["spot"]].ctr_mass_y(),
                   OO.parent.spots[OO.parent.indexed_pairs[obsno]["spot"]].ctr_mass_x(),
                  0.0))
             ).length()
            for obsno in range(len(OO.parent.indexed_pairs))])

          rmsdexc = math.sqrt(flex.mean(displacements*displacements))
          print("rotz %7.3f degrees, RMSD displacement %7.3f pixels"%(
          (rotz * 180./math.pi), rmsdexc))
          return list(displacements)

      helper = per_frame_helper()

      print("Trying least squares minimization of displacements", end=' ')

      results = leastsq(
        func = helper.fvec_callable,
        x0 = [0.] + list(helper.convert.forward_independent_parameters()),
        args = (),
        Dfun = None, #estimate the Jacobian
        full_output = True)

      print("with %d reflections"%len(OO.parent.indexed_pairs), end=' ')
      print("result %6.2f degrees"%(results[0][0]*180./math.pi))
      return results[0], helper.last_set_orientation

def pre_get_predictions(inputai,horizons_phil,raw_image,imageindex,spotfinder,limiting_resolution,domain_size_ang=0):
    from rstbx.apps.slip_helpers import wrapper_of_use_case_bp3
    wrapbp3 = wrapper_of_use_case_bp3( raw_image = raw_image,
      spotfinder = spotfinder, imageindex = imageindex,
      inputai = inputai,
      spot_prediction_limiting_resolution = limiting_resolution,
      phil_params = horizons_phil,
      sub = horizons_phil.integration.use_subpixel_translations)

    ###  XXX  pass this in as a parameter
    bandpass = 1.E-3 # better 1.E-4 than 1.E-3 for fake_psI
    print("THE GET-PREDICTIONS DOMAIN SIZE is %8.3f"%domain_size_ang, end=' ')
    print("FW mosaicity %7.2f deg"%inputai.getMosaicity())
    wrapbp3.set_variables( orientation = inputai.getOrientation(),
                         wave_HI = inputai.wavelength * (1.-(bandpass/2.)),
                         wave_LO = inputai.wavelength * (1.+(bandpass/2.)),
                         half_mosaicity_deg = inputai.getMosaicity()/2.,
                         domain_size = domain_size_ang)
    wrapbp3.ucbp3.picture_fast_slow() # abbreviated
    return wrapbp3

def post_outlier_rejection(parent,image_number,cb_op_to_primitive,horizons_phil,kwargs):
  reserve_wavelength = parent.inputai.wavelength
  verbose = False

  def core_optimization(operational_wavelength):
    parent.inputai.setWavelength(operational_wavelength)
    """parent supplies all the "self" variables referred to above"""
    # first refine rotx and roty
    R = refinement(parent,mosaic_refinement_target=horizons_phil.integration.mosaic.refinement_target,
                   pvr_fix = horizons_phil.integration.mosaic.bugfix2_enable)
    if verbose: excursions,positions = R.contour_plot()
    minimum = R.refine_rotx_roty2(enable_rotational_target =
      horizons_phil.integration.mosaic.enable_rotational_target_highsym)
    if verbose: R.show_plot(excursions,positions,minimum)
    parent.inputai.setOrientation(minimum[1])
    print("RDISTANCE %8.3f X %8.3f Y %8.3f A %8.3f C %8.3f"%(parent.inputai.distance(),
      parent.inputai.xbeam(),parent.inputai.ybeam(),minimum[1].unit_cell().parameters()[0],
      minimum[1].unit_cell().parameters()[2]))

    # now refine rotz, unit cell, distance, beamxy
    refine2 = False
    """As implemented the refine2 seems to skew the excursion vs. resolution plot so as to
       make the domain-size result not meaningful.  Therefore comment this refinement
       out for now.
       P.S. 9/2014 Problem seems to that the code did not implement the subpixel joint model.
    """
    if refine2:
      R2 = refinement2(parent)
      minimum = R2.refine_all()

      parent.inputai.setOrientation(minimum[1])
      print("R2DISTANCE %8.3f X %8.3f Y %8.3f A %8.3f C %8.3f"%(parent.inputai.distance(),
        parent.inputai.xbeam(),parent.inputai.ybeam(),minimum[1].unit_cell().parameters()[0],
        minimum[1].unit_cell().parameters()[2]))

    # last refine rotx and roty
    R = refinement(parent,mosaic_refinement_target=horizons_phil.integration.mosaic.refinement_target,
                   pvr_fix = horizons_phil.integration.mosaic.bugfix2_enable)
    if verbose: excursions,positions = R.contour_plot()
    minimum = R.refine_rotx_roty2(enable_rotational_target =
      horizons_phil.integration.mosaic.enable_rotational_target_highsym)
    return minimum

  from scitbx.simplex import simplex_opt
  class apply_simplex_method(object):
    def __init__(selfOO):
      selfOO.starting_simplex=[]
      selfOO.n = 1
      for ii in range(selfOO.n+1):
        selfOO.starting_simplex.append(flex.random_double(selfOO.n))
      selfOO.optimizer = simplex_opt( dimension=selfOO.n,
                                    matrix  = selfOO.starting_simplex,
                                    evaluator = selfOO,
                                    tolerance=1e-4)
      selfOO.x = selfOO.optimizer.get_solution()

    def target(selfOO, vector):
      selfOO.minimum = core_optimization(
        operational_wavelength = reserve_wavelength + vector[0]*0.001)
      return parent.inputai.getMosaicity()

  if horizons_phil.integration.mosaic.enable_simplex:
    MIN = apply_simplex_method()
    print("MINIMUM=",list(MIN.x))
    minimum = MIN.minimum
  else:
    minimum = core_optimization(reserve_wavelength)

  if verbose: R.show_plot(excursions,positions,minimum)
  parent.inputai.setOrientation(minimum[1])
  print("R3DISTANCE %8.3f X %8.3f Y %8.3f A %8.3f C %8.3f"%(parent.inputai.distance(),
    parent.inputai.xbeam(),parent.inputai.ybeam(),minimum[1].unit_cell().parameters()[0],
    minimum[1].unit_cell().parameters()[2]))


  kwargs["user-reentrant"]=True
  kwargs["domain_size_ang"]=minimum[2]
  parent.integration_concept(image_number,cb_op_to_primitive,False,**kwargs)
