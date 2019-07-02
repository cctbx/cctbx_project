from __future__ import absolute_import, division, print_function
from six.moves import range
import math
from cctbx.array_family import flex
from scitbx import matrix
from libtbx.utils import Sorry

#234567890123456789212345678931234567894123456789512345678961234567897123456789812

from labelit.dptbx.profile_support import show_profile
from rstbx.apps.slip_helpers import slip_callbacks
from rstbx.dials_core.integration_core import integration_core
from six.moves import zip

class IntegrationMetaProcedure(integration_core,slip_callbacks):

  def __init__(self):
    self.block_counter = 0

  def set_up_mask_focus(self,verbose=False):
    self.mask_focus = []
    for frame in self.frame_numbers:
      focus = self.spotfinder.pd['masks'][frame][0:2]
      if len(self.spotfinder.pd['masks'][frame]) < 3 or self.spotfinder.pd['masks'][frame][2] is None:
        self.mask_focus.append( None )
        raise Sorry("No average profile available to set up the integration mask")
        continue; #no average profile; no pred/obs agreement; nothing possible
      average_profile = self.inputpd['masks'][frame][2]
      if verbose:
        box = self.inputpd['masks'][frame][3]
        print(average_profile.focus())
        print(box.focus())
        print("Average Profile:")
        show_profile( average_profile )
        print("Box:")
        show_profile( box )
      self.mask_focus.append( average_profile.focus() )

  def get_predictions_accounting_for_centering(self,cb_op_to_primitive=None,**kwargs):
    # interface requires this function to set current_orientation
    # in the actual setting used for Miller index calculation
    if (self.horizons_phil.known_setting is None or self.horizons_phil.known_setting == self.setting_id ) and \
        self.horizons_phil.integration.model in ["use_case_3_simulated_annealing",
                                                "use_case_3_simulated_annealing_7",
                                                "use_case_3_simulated_annealing_9"]:
      if cb_op_to_primitive==None:
        raise Sorry("Can't use model_3 simulated annealing for non-primitive cells, contact authors.")
      if self.horizons_phil.integration.model=="use_case_3_simulated_annealing":
        self.best_params = dict(zip(("half_mosaicity_deg","wave_HE_ang","wave_LE_ang",
         "reserve_orientation","rotation100_rad","rotation010_rad","rotation001_rad"),
         self.use_case_3_simulated_annealing(self.horizons_phil.integration.use_subpixel_translations))
        )
      elif self.horizons_phil.integration.model=="use_case_3_simulated_annealing_7":
        self.best_params = dict(zip(("half_mosaicity_deg","wave_HE_ang","wave_LE_ang",
         "reserve_orientation","rotation100_rad","rotation010_rad","rotation001_rad",
         "domain_size_ang"),
         self.use_case_3_simulated_annealing_7(self.horizons_phil.integration.use_subpixel_translations))
        )
      elif self.horizons_phil.integration.model=="use_case_3_simulated_annealing_9":
        self.best_params = dict(zip(("half_mosaicity_deg","wave_HE_ang","wave_LE_ang",
         "reserve_orientation","rotation100_rad","rotation010_rad","rotation001_rad",
         "domain_size_ang","ab_factor","c_factor"),
         self.use_case_3_simulated_annealing_9(self.horizons_phil.integration.use_subpixel_translations))
        )
      self.current_orientation = self.best_params["reserve_orientation"]
      self.current_cb_op_to_primitive = cb_op_to_primitive

      BPpredicted = self.bp3_wrapper.ucbp3.selected_predictions_labelit_format()
      BPhkllist = self.bp3_wrapper.ucbp3.selected_hkls()
      self.predicted,self.hkllist = BPpredicted, BPhkllist
      self.partialities = dict(indices=BPhkllist.deep_copy(),
                               data=self.bp3_wrapper.ucbp3.selected_partialities())
      #self.hi = self.bp3_wrapper.ucbp3.selected_hi_predictions()
      #self.lo = self.bp3_wrapper.ucbp3.selected_lo_predictions()
      #print list(self.partialities)
      #for x in range(len(self.hi)):
      #  print "%4d %4d %4d  %7.1f %7.1f %7.1f %7.1f PARTIAL %7.2f"%(
      #  self.hkllist[x][0], self.hkllist[x][1],self.hkllist[x][2],
      #  self.hi[x][0], self.hi[x][1],
      #  self.lo[x][0], self.lo[x][1], self.partialities[x])

      if self.inputai.active_areas != None:
        self.predicted,self.hkllist = self.inputai.active_areas(
                                      self.predicted,self.hkllist,self.pixel_size)
      return

    if self.horizons_phil.integration.model == "user_supplied":

     lower_limit_domain_size = math.pow(
       self.inputai.getOrientation().unit_cell().volume(),
       1./3.)*self.horizons_phil.integration.mosaic.domain_size_lower_limit # default 10-unit cell block size minimum reasonable domain
     actual_used_domain_size = kwargs.get("domain_size_ang",lower_limit_domain_size)
     if cb_op_to_primitive==None:

      from cxi_user import pre_get_predictions
      self.bp3_wrapper = pre_get_predictions(self.inputai, self.horizons_phil,
        raw_image = self.imagefiles.images[self.image_number],
        imageindex = self.frame_numbers[self.image_number],
        spotfinder = self.spotfinder,
        limiting_resolution = self.limiting_resolution,
        domain_size_ang = actual_used_domain_size)
      self.current_orientation = self.inputai.getOrientation()
      self.current_cb_op_to_primitive = cb_op_to_primitive

      BPpredicted = self.bp3_wrapper.ucbp3.selected_predictions_labelit_format()
      BPhkllist = self.bp3_wrapper.ucbp3.selected_hkls()

      self.predicted,self.hkllist = BPpredicted, BPhkllist
      if self.inputai.active_areas != None:
        self.predicted,self.hkllist = self.inputai.active_areas(
                                      self.predicted,self.hkllist,self.pixel_size)
      return

     else:
      self.block_counter+=1
      rot_mat = matrix.sqr(cb_op_to_primitive.c().r().as_double()).transpose()
      centered_orientation = self.inputai.getOrientation()
      self.current_orientation = centered_orientation
      self.current_cb_op_to_primitive = cb_op_to_primitive
      primitive_orientation = centered_orientation.change_basis(rot_mat)
      self.inputai.setOrientation(primitive_orientation)
      from cxi_user import pre_get_predictions
      if self.block_counter < 2:
        KLUDGE = self.horizons_phil.integration.mosaic.kludge1 # bugfix 1 of 2 for protocol 6, equation 2
        self.inputai.setMosaicity(KLUDGE*self.inputai.getMosaicity())
      self.bp3_wrapper = pre_get_predictions(self.inputai, self.horizons_phil,
        raw_image = self.imagefiles.images[self.image_number],
        imageindex = self.frame_numbers[self.image_number],
        spotfinder = self.spotfinder,
        limiting_resolution = self.limiting_resolution,
        domain_size_ang = actual_used_domain_size)

      BPpredicted = self.bp3_wrapper.ucbp3.selected_predictions_labelit_format()
      BPhkllist = self.bp3_wrapper.ucbp3.selected_hkls()

      self.actual = actual_used_domain_size
      self.predicted = BPpredicted
      primitive_hkllist = BPhkllist
      #not sure if matrix needs to be transposed first for outputting HKL's???:
      self.hkllist = cb_op_to_primitive.inverse().apply(primitive_hkllist)
      self.inputai.setOrientation(centered_orientation)
      if self.inputai.active_areas != None:
        self.predicted,self.hkllist = self.inputai.active_areas(
                                      self.predicted,self.hkllist,self.pixel_size)
      if self.block_counter < 2:
         down = self.inputai.getMosaicity()/KLUDGE
         print("Readjusting mosaicity back down to ",down)
         self.inputai.setMosaicity(down)
      return

    if cb_op_to_primitive==None:

      predicted = self.inputai.predict_all(
                  self.image_centers[self.image_number],self.limiting_resolution)
      self.predicted = predicted.vec3() #only good for integrating one frame...
      self.hkllist = predicted.hkl()
      self.current_orientation = self.inputai.getOrientation()
      from cctbx import sgtbx
      self.cb_op_to_primitive = sgtbx.change_of_basis_op()

    else:
      rot_mat = matrix.sqr(cb_op_to_primitive.c().r().as_double()).transpose()
      centered_orientation = self.inputai.getOrientation()
      self.current_orientation = centered_orientation
      self.current_cb_op_to_primitive = cb_op_to_primitive
      primitive_orientation = centered_orientation.change_basis(rot_mat)
      self.inputai.setOrientation(primitive_orientation)
      predicted = self.inputai.predict_all(
                  self.image_centers[self.image_number],self.limiting_resolution)
      self.predicted = predicted.vec3() #only good for integrating one frame...
      primitive_hkllist = predicted.hkl()
      #not sure if matrix needs to be transposed first for outputting HKL's???:
      self.hkllist = cb_op_to_primitive.inverse().apply(primitive_hkllist)
      self.inputai.setOrientation(centered_orientation)
    if self.inputai.active_areas != None:
      self.predicted,self.hkllist = self.inputai.active_areas(
                                    self.predicted,self.hkllist,self.pixel_size)

    if False: #development only; compare the two methods:
      from matplotlib import pyplot as plt
      plt.plot([i[0] for i in BPpredicted],[i[1] for i in BPpredicted],"r.")
      plt.plot([i[0] for i in predicted],[i[1] for i in predicted],"b.")
      plt.show()

  def get_observations_with_outlier_removal(self):
    print("Using spotfinder subset",self.horizons_phil.integration.spotfinder_subset)
    spots = self.spotfinder.images[self.frame_numbers[self.image_number]][self.horizons_phil.integration.spotfinder_subset]
    if getattr(slip_callbacks.slip_callback,"requires_refinement_spots",False):
      from spotfinder.array_family import flex
      self.spotfinder.images[self.frame_numbers[self.image_number]]["refinement_spots"]=flex.distl_spot()
    return spots

  def integration_concept(self,image_number=0,cb_op_to_primitive=None,verbose=False,**kwargs):
    self.image_number = image_number
    NEAR = 10
    pxlsz = self.pixel_size
    self.get_predictions_accounting_for_centering(cb_op_to_primitive,**kwargs)
    FWMOSAICITY = self.inputai.getMosaicity()
    DOMAIN_SZ_ANG = kwargs.get("domain_size_ang",  self.__dict__.get("actual",0)  )
    refineflag = {True:0,False:1}[kwargs.get("domain_size_ang",0)==0]
    self.inputpd["symmetry"].show_summary(prefix="EXCURSION%1d REPORT FWMOS= %6.4f DOMAIN= %6.1f "%(refineflag,FWMOSAICITY,DOMAIN_SZ_ANG))
    from annlib_ext import AnnAdaptor
    self.cell = self.inputai.getOrientation().unit_cell()
    query = flex.double()
    for pred in self.predicted: # predicted spot coord in pixels
      query.append(pred[0]/pxlsz)
      query.append(pred[1]/pxlsz)
    self.reserve_hkllist_for_signal_search = self.hkllist

    reference = flex.double()
    spots = self.get_observations_with_outlier_removal()

    assert len(spots)>NEAR# Can't do spot/pred matching with too few spots
    for spot in spots:
      reference.append(spot.ctr_mass_x())
      reference.append(spot.ctr_mass_y())

    IS_adapt = AnnAdaptor(data=reference,dim=2,k=NEAR)
    IS_adapt.query(query)
    print("Calculate correction vectors for %d observations & %d predictions"%(len(spots),len(self.predicted)))
    indexed_pairs_provisional = []
    correction_vectors_provisional = []
    c_v_p_flex = flex.vec3_double()
    idx_cutoff = float(min(self.mask_focus[image_number]))
    if verbose:
      print("idx_cutoff distance in pixels",idx_cutoff)
    if not self.horizons_phil.integration.enable_one_to_one_safeguard:
     # legacy code, no safeguard against many-to-one predicted-to-observation mapping
     for i in range(len(self.predicted)): # loop over predicteds
      #for n in range(NEAR): # loop over near spotfinder spots
      for n in range(1): # only consider the nearest spotfinder spot
        Match = dict(spot=IS_adapt.nn[i*NEAR+n],pred=i)
        if n==0 and math.sqrt(IS_adapt.distances[i*NEAR+n]) < idx_cutoff:
          indexed_pairs_provisional.append(Match)

          vector = matrix.col(
            [spots[Match["spot"]].ctr_mass_x() - self.predicted[Match["pred"]][0]/pxlsz,
             spots[Match["spot"]].ctr_mass_y() - self.predicted[Match["pred"]][1]/pxlsz])
          correction_vectors_provisional.append(vector)
          c_v_p_flex.append((vector[0],vector[1],0.))
    else:
      one_to_one = {}
      for i in range(len(self.predicted)): # loop over predicteds
        annresultidx = i*NEAR
        obsidx = IS_adapt.nn[annresultidx]
        this_distancesq = IS_adapt.distances[annresultidx]
        if obsidx not in one_to_one or \
           this_distancesq < one_to_one[obsidx]["distancesq"]:
           if math.sqrt(this_distancesq) < idx_cutoff:
             one_to_one[obsidx] = dict(spot=obsidx,pred=i,distancesq=this_distancesq)
      for key,value in one_to_one.items():
        indexed_pairs_provisional.append(value)
        vector = matrix.col(
            [spots[value["spot"]].ctr_mass_x() - self.predicted[value["pred"]][0]/pxlsz,
             spots[value["spot"]].ctr_mass_y() - self.predicted[value["pred"]][1]/pxlsz])
        correction_vectors_provisional.append(vector)
        c_v_p_flex.append((vector[0],vector[1],0.))

    print("... %d provisional matches"%len(correction_vectors_provisional), end=' ')
    print("r.m.s.d. in pixels: %5.2f"%(math.sqrt(flex.mean(c_v_p_flex.dot(c_v_p_flex)))))

    if self.horizons_phil.integration.enable_residual_scatter:
      from matplotlib import pyplot as plt
      fig = plt.figure()
      for cv in correction_vectors_provisional:
        plt.plot([cv[1]],[-cv[0]],"b.")
      plt.title(" %d matches, r.m.s.d. %5.2f pixels"%(len(correction_vectors_provisional),math.sqrt(flex.mean(c_v_p_flex.dot(c_v_p_flex)))))
      plt.axes().set_aspect("equal")
      self.show_figure(plt,fig,"res")
      plt.close()

    if self.horizons_phil.integration.enable_residual_map:
      from matplotlib import pyplot as plt
      fig = plt.figure()
      for match,cv in zip(indexed_pairs_provisional,correction_vectors_provisional):
        plt.plot([spots[match["spot"]].ctr_mass_y()],[-spots[match["spot"]].ctr_mass_x()],"r.")
        plt.plot([self.predicted[match["pred"]][1]/pxlsz],[-self.predicted[match["pred"]][0]/pxlsz],"g.")
        plt.plot([spots[match["spot"]].ctr_mass_y(), spots[match["spot"]].ctr_mass_y() + 10.*cv[1]],
                 [-spots[match["spot"]].ctr_mass_x(), -spots[match["spot"]].ctr_mass_x() - 10.*cv[0]],'b-')
      plt.xlim([0,float(self.inputpd["size2"])])
      plt.ylim([-float(self.inputpd["size1"]),0])
      plt.title(" %d matches, r.m.s.d. %5.2f pixels"%(len(correction_vectors_provisional),math.sqrt(flex.mean(c_v_p_flex.dot(c_v_p_flex)))))
      plt.axes().set_aspect("equal")
      self.show_figure(plt,fig,"map")
      plt.close()
    # insert code here to remove correction length outliers...
    # they are causing terrible
    # problems for finding legitimate correction vectors (print out the list)
    # also remove outliers for the purpose of reporting RMS
    outlier_rejection = True
    cache_refinement_spots = getattr(slip_callbacks.slip_callback,"requires_refinement_spots",False)
    if outlier_rejection:
      correction_lengths = flex.double([v.length() for v in correction_vectors_provisional])
      clorder = flex.sort_permutation(correction_lengths)
      sorted_cl = correction_lengths.select(clorder)

      ACCEPTABLE_LIMIT = 2
      limit = int(0.33 * len(sorted_cl)) # best 1/3 of data are assumed to be correctly modeled.
      if (limit <= ACCEPTABLE_LIMIT):
        raise Sorry("Not enough indexed spots to reject outliers; have %d need >%d" % (limit, ACCEPTABLE_LIMIT))

      y_data = flex.double(len(sorted_cl))
      for i in range(len(y_data)):
        y_data[i] = float(i)/float(len(y_data))

      # ideas are explained in Sauter & Poon (2010) J Appl Cryst 43, 611-616.
      from rstbx.outlier_spots.fit_distribution import fit_cdf,rayleigh
      fitted_rayleigh = fit_cdf(x_data = sorted_cl[0:limit],
                                y_data = y_data[0:limit],
                                distribution=rayleigh)

      inv_cdf = [fitted_rayleigh.distribution.inv_cdf(cdf) for cdf in y_data]

      #print "SORTED LIST OF ",len(sorted_cl), "with sigma",fitted_rayleigh.distribution.sigma
      indexed_pairs = []
      correction_vectors = []
      self.correction_vectors = []
      for icand in range(len(sorted_cl)):
        # somewhat arbitrary sigma = 1.0 cutoff for outliers
        if (sorted_cl[icand]-inv_cdf[icand])/fitted_rayleigh.distribution.sigma > 1.0:
          break
        indexed_pairs.append(indexed_pairs_provisional[clorder[icand]])
        correction_vectors.append(correction_vectors_provisional[clorder[icand]])
        if cache_refinement_spots:
          self.spotfinder.images[self.frame_numbers[self.image_number]]["refinement_spots"].append(
          spots[indexed_pairs[-1]["spot"]])
        if kwargs.get("verbose_cv")==True:
            print("CV OBSCENTER %7.2f %7.2f REFINEDCENTER %7.2f %7.2f"%(
              float(self.inputpd["size1"])/2.,float(self.inputpd["size2"])/2.,
              self.inputai.xbeam()/pxlsz, self.inputai.ybeam()/pxlsz), end=' ')
            print("OBSSPOT %7.2f %7.2f PREDSPOT %7.2f %7.2f"%(
              spots[indexed_pairs[-1]["spot"]].ctr_mass_x(),
              spots[indexed_pairs[-1]["spot"]].ctr_mass_y(),
              self.predicted[indexed_pairs[-1]["pred"]][0]/pxlsz,
              self.predicted[indexed_pairs[-1]["pred"]][1]/pxlsz), end=' ')
            the_hkl = self.hkllist[indexed_pairs[-1]["pred"]]
            print("HKL %4d %4d %4d"%the_hkl,"%2d"%self.setting_id, end=' ')
            radial, azimuthal = spots[indexed_pairs[-1]["spot"]].get_radial_and_azimuthal_size(
              self.inputai.xbeam()/pxlsz, self.inputai.ybeam()/pxlsz)
            print("RADIALpx %5.3f AZIMUTpx %5.3f"%(radial,azimuthal))

        # Store a list of correction vectors in self.
        radial, azimuthal = spots[indexed_pairs[-1]['spot']].get_radial_and_azimuthal_size(
          self.inputai.xbeam()/pxlsz, self.inputai.ybeam()/pxlsz)
        self.correction_vectors.append(
          dict(obscenter=(float(self.inputpd['size1']) / 2,
                          float(self.inputpd['size2']) / 2),
               refinedcenter=(self.inputai.xbeam() / pxlsz,
                              self.inputai.ybeam() / pxlsz),
               obsspot=(spots[indexed_pairs[-1]['spot']].ctr_mass_x(),
                        spots[indexed_pairs[-1]['spot']].ctr_mass_y()),
               predspot=(self.predicted[indexed_pairs[-1]['pred']][0] / pxlsz,
                         self.predicted[indexed_pairs[-1]['pred']][1] / pxlsz),
               hkl=(self.hkllist[indexed_pairs[-1]['pred']][0],
                    self.hkllist[indexed_pairs[-1]['pred']][1],
                    self.hkllist[indexed_pairs[-1]['pred']][2]),
               setting_id=self.setting_id,
               radial=radial,
               azimuthal=azimuthal))

      print("After outlier rejection %d indexed spotfinder spots remain."%len(indexed_pairs))
      if False:
        rayleigh_cdf = [
          fitted_rayleigh.distribution.cdf(x=sorted_cl[c]) for c in range(len(sorted_cl))]
        from matplotlib import pyplot as plt
        plt.plot(sorted_cl,y_data,"r+")
        #plt.plot(sorted_cl,rayleigh_cdf,"g.")
        plt.plot(inv_cdf,y_data,"b.")
        plt.show()
    else:
      indexed_pairs = indexed_pairs_provisional
      correction_vectors = correction_vectors_provisional
    ########### finished with outlier rejection

    self.inputpd["symmetry"].show_summary(prefix="SETTING ")

    is_triclinic = (self.setting_id==1)
    if is_triclinic:
      self.triclinic_pairs = [ dict(pred=self.hkllist[a["pred"]],spot=a["spot"])
        for a in indexed_pairs ]

    if self.horizons_phil.integration.model == "user_supplied":
      if kwargs.get("user-reentrant",None)==None:
        from cxi_user import post_outlier_rejection
        self.indexed_pairs = indexed_pairs
        self.spots = spots
        post_outlier_rejection(self,image_number,cb_op_to_primitive,self.horizons_phil,kwargs)
        return

    ########### finished with user-supplied code

    if self.horizons_phil.integration.spot_shape_verbose:
        from rstbx.new_horizons.spot_shape import spot_shape_verbose
        spot_shape_verbose(rawdata = self.imagefiles.images[self.image_number].linearintdata,
           beam_center_pix = matrix.col((self.inputai.xbeam()/pxlsz, self.inputai.ybeam()/pxlsz)),
           indexed_pairs = indexed_pairs,
           spotfinder_observations = spots,
           distance_mm = self.inputai.distance(),
           mm_per_pixel = pxlsz,
           hkllist = self.hkllist,
           unit_cell = self.cell,
           wavelength_ang = self.inputai.wavelength
        )

    #Other checks to be implemented (future):
    # spot is within active area of detector on a circular detector such as the Mar IP
    # integration masks do not overlap; or deconvolute

    correction_lengths=flex.double([v.length() for v in correction_vectors])
    if verbose:
      print("average correction %5.2f over %d vectors"%(flex.mean(correction_lengths),
      len(correction_lengths)), end=' ')
      print("or %5.2f mm."%(pxlsz*flex.mean(correction_lengths)))
    self.r_residual = pxlsz*flex.mean(correction_lengths)

    #assert len(indexed_pairs)>NEAR # must have enough indexed spots
    if (len(indexed_pairs) <= NEAR):
      raise Sorry("Not enough indexed spots, only found %d, need %d" % (len(indexed_pairs), NEAR))

    reference = flex.double()
    for item in indexed_pairs:
      reference.append(spots[item["spot"]].ctr_mass_x())
      reference.append(spots[item["spot"]].ctr_mass_y())

    PS_adapt = AnnAdaptor(data=reference,dim=2,k=NEAR)
    PS_adapt.query(query)

    self.BSmasks = []
    #self.null_correction_mapping( predicted=self.predicted,
    #                                    correction_vectors = correction_vectors,
    #                                    IS_adapt = IS_adapt,
    #                                    spots = spots)
    self.positional_correction_mapping( predicted=self.predicted,
                                        correction_vectors = correction_vectors,
                                        PS_adapt = PS_adapt,
                                        IS_adapt = IS_adapt,
                                        spots = spots)

    # which spots are close enough to interfere with background?
    MAXOVER=6
    OS_adapt = AnnAdaptor(data=query,dim=2,k=MAXOVER) #six near nbrs
    OS_adapt.query(query)
    if self.mask_focus[image_number] is None:
      raise Sorry("No observed/predicted spot agreement; no Spotfinder masks; skip integration")
    nbr_cutoff = 2.0* max(self.mask_focus[image_number])
    FRAME = int(nbr_cutoff/2)
    #print "The overlap cutoff is %d pixels"%nbr_cutoff
    nbr_cutoff_sq = nbr_cutoff * nbr_cutoff

    #print "Optimized C++ section...",
    self.set_frame(FRAME)
    self.set_background_factor(kwargs["background_factor"])
    self.set_nbr_cutoff_sq(nbr_cutoff_sq)
    self.set_guard_width_sq(self.horizons_phil.integration.guard_width_sq)
    self.set_detector_gain(self.horizons_phil.integration.detector_gain)
    flex_sorted = flex.int()
    for item in self.sorted:
      flex_sorted.append(item[0]);flex_sorted.append(item[1]);

    if self.horizons_phil.integration.mask_pixel_value is not None:
      self.set_mask_pixel_val(self.horizons_phil.integration.mask_pixel_value)

    image_obj = self.imagefiles.imageindex(self.frame_numbers[self.image_number])
    image_obj.read()
    rawdata = image_obj.linearintdata # assume image #1

    if self.inputai.active_areas != None:
      self.detector_xy_draft = self.safe_background( rawdata=rawdata,
                          predicted=self.predicted,
                          OS_adapt=OS_adapt,
                          sorted=flex_sorted,
                          tiles=self.inputai.active_areas.IT,
                          tile_id=self.inputai.active_areas.tile_id);
    else:
      self.detector_xy_draft = self.safe_background( rawdata=rawdata,
                          predicted=self.predicted,
                          OS_adapt=OS_adapt,
                          sorted=flex_sorted);
    for i in range(len(self.predicted)): # loop over predicteds
      B_S_mask = {}
      keys = self.get_bsmask(i)
      for k in range(0,len(keys),2):
        B_S_mask[(keys[k],keys[k+1])]=True
      self.BSmasks.append(B_S_mask)
    #print "Done"
    return

  def show_rejected_spots(self):
    miller = self.get_rejected_miller()
    messag = self.get_rejected_reason()
    for i,j in zip(self.get_rejected_miller(),self.get_rejected_reason()):
      print(i,j)

  def integration_proper(self):
    image_obj = self.imagefiles.imageindex(self.frame_numbers[self.image_number])
    #image_obj.read() #assume image already read
    rawdata = image_obj.linearintdata # assume image #1

    self.integration_proper_fast(rawdata,self.predicted,self.hkllist,self.detector_xy_draft)
    self.integrated_data = self.get_integrated_data()
    self.integrated_sigma= self.get_integrated_sigma()
    self.integrated_miller=self.get_integrated_miller()
    self.detector_xy = self.get_detector_xy()
    self.max_signal = self.get_max_signal()

    for correction_type in self.horizons_phil.integration.absorption_correction:
      if correction_type.apply:
        if correction_type.algorithm=="fuller_kapton":
          print("Absorption correction with %d reflections to correct"%(len(self.detector_xy)))
          from cxi_xdr_xes import absorption
          C = absorption.correction()
          if correction_type.fuller_kapton.smart_sigmas:
            self.fuller_kapton_absorption_correction, self.fuller_kapton_absorption_sigmas = C(
              panel_size_px = (self.inputpd['size1'],self.inputpd['size2']),
              pixel_size_mm = self.pixel_size,
              detector_dist_mm = self.inputai.distance(),
              wavelength_ang = self.inputai.wavelength,
              BSmasks = self.BSmasks,
              get_ISmask_function = self.get_ISmask,
              params = correction_type.fuller_kapton,
              i_no_skip = self.get_integrated_flag(),
              calc_sigmas=True
            )
            # apply corrections and propagate error
            # term1 = (sig(C)/C)^2
            # term2 = (sig(Imeas)/Imeas)^2
            # I' = C*I
            # sig^2(I') = (I')^2*(term1 + term2)
            # sig(I') = sqrt(sig^2(I'))
            term1 = flex.pow(self.fuller_kapton_absorption_sigmas/self.fuller_kapton_absorption_correction, 2)
            term2 = flex.pow(self.integrated_sigma/self.integrated_data, 2)
            self.integrated_data *= self.fuller_kapton_absorption_correction
            integrated_sigma_squared = flex.pow(self.integrated_data, 2) * (term1 + term2)
            self.integrated_sigma = flex.sqrt(integrated_sigma_squared)
            # order is purposeful: the two lines above require that self.integrated_data has already been corrected!
          else:
            self.fuller_kapton_absorption_correction = C(
              panel_size_px = (self.inputpd['size1'],self.inputpd['size2']),
              pixel_size_mm = self.pixel_size,
              detector_dist_mm = self.inputai.distance(),
              wavelength_ang = self.inputai.wavelength,
              BSmasks = self.BSmasks,
              get_ISmask_function = self.get_ISmask,
              params = correction_type.fuller_kapton,
              i_no_skip = self.get_integrated_flag()
            )
            # apply these corrections now
            self.integrated_data *= self.fuller_kapton_absorption_correction
            self.integrated_sigma *= self.fuller_kapton_absorption_correction

    #self.show_rejected_spots()
    return # function has been recoded in C++

  def get_obs(self,space_group_symbol):
    from cctbx.crystal import symmetry
    from cctbx import miller

    xsym = symmetry(unit_cell = self.cell,
                    space_group_symbol=space_group_symbol)

    miller_set = miller.set(crystal_symmetry=xsym,
      indices=self.integrated_miller,anomalous_flag=True)
    miller_array = miller.array(miller_set,self.integrated_data,
      self.integrated_sigma)
    miller_array.set_observation_type_xray_intensity()
    miller_array.set_info("Raw partials from rstbx, not in ASU, no polarization correction")
    return miller_array

  def show_figure(self,plt,fig,tag):
    if self.horizons_phil.integration.graphics_backend=="pdf":

      F=self.imagefiles.frames()
      G=self.imagefiles.imagepath(F[0])
      import os
      fname = os.path.splitext(os.path.basename(G))[0]
      sgi = str(self.inputpd["symmetry"].space_group_info())
      sgi = sgi.replace(" ","")
      sgi = sgi.replace("/","")
      outfile = os.path.join(self.horizons_phil.integration.pdf_output_dir,
        "%s_%02d_%s_%s.pdf"%(fname, self.setting_id, sgi, tag)
      )

      from matplotlib.backends.backend_pdf import PdfPages
      with PdfPages(outfile) as pdf:
        pdf.savefig(fig)
        pdf.savefig()
    else:
        plt.show()
