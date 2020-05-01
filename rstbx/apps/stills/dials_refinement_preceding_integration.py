from __future__ import absolute_import, division, print_function
from six.moves import range
from rstbx.apps.stills.simple_integration import IntegrationMetaProcedure
from rstbx.apps import simple_integration
from scitbx import matrix
import math,copy
from dials.array_family import flex
from six.moves import zip

class integrate_one_frame(IntegrationMetaProcedure):
  def __init__(self, triclinic):
    simple_integration.__init__(self)
    IntegrationMetaProcedure.__init__(self)
    self.triclinic_pairs = triclinic.triclinic_pairs
    self.triclinic = triclinic

  def integration_concept(self, image_number, cb_op_to_primitive, verbose=False, **kwargs):
    if kwargs.get("user-reentrant") != None:
      return self.integration_concept_detail(experiments=kwargs.get("reentrant_experiments"),
                                    reflections=kwargs.get("reentrant_reflections"),
                                    spots=self.triclinic.get_observations_with_outlier_removal(),
                                    image_number=image_number, cb_op_to_primitive=cb_op_to_primitive, **kwargs)

    experiments = self.prepare_dxtbx_models(setting_specific_ai = self.triclinic.inputai, sg="P1")
    reflections = self.prepare_reflection_list(experiments[0].detector)
    self.refine(experiments=experiments, reflections=reflections)

    if self.horizons_phil.isoforms == []: # normal behavior; no isoform constraints
      setting_experiments = self.prepare_dxtbx_models(setting_specific_ai = self.inputai,
      sg = self.inputpd["symmetry"].space_group_info().type().lookup_symbol())
      setting_reflections = copy.deepcopy(reflections)
      setting_reflections["miller_index"] = cb_op_to_primitive.inverse().apply(reflections["miller_index"])
      R = self.refine(experiments=setting_experiments, reflections=setting_reflections)
    else:
      isoform_refineries = []
      setting_reflections = copy.deepcopy(reflections)
      setting_reflections["miller_index"] = cb_op_to_primitive.inverse().apply(reflections["miller_index"])
      look_symbol = self.inputpd["symmetry"].space_group_info().type().lookup_symbol()

      for isoform in self.horizons_phil.isoforms:
        print("Testing isoform %s"%isoform.name,isoform.cell.parameters())
        print("asserting", look_symbol ,"==", isoform.lookup_symbol)

        assert look_symbol == isoform.lookup_symbol
        setting_experiments = self.prepare_dxtbx_models(setting_specific_ai = self.inputai,
          sg = look_symbol, isoform = isoform.cell )
        P = self.refine(experiments=setting_experiments, reflections=setting_reflections, isoform=True)
        print(P.rmsds())
        isoform_refineries.append(P)
      positional_rmsds = [math.sqrt(P.rmsds()[0]**2 + P.rmsds()[1]**2) for P in isoform_refineries]
      print("Positional rmsds for all isoforms:", positional_rmsds)
      minrmsd_mm = min(positional_rmsds)
      minindex = positional_rmsds.index(minrmsd_mm)
      print("The smallest rmsd is %5.1f um from isoform %s"%(1000.*minrmsd_mm,self.horizons_phil.isoforms[minindex].name))
      if self.horizons_phil.isoforms[minindex].rmsd_target_mm is not None:
        print("asserting",minrmsd_mm ,"<", self.horizons_phil.isoforms[minindex].rmsd_target_mm)
        assert minrmsd_mm < self.horizons_phil.isoforms[minindex].rmsd_target_mm
      print("Acceptable rmsd for isoform %s."%(self.horizons_phil.isoforms[minindex].name), end=' ')
      if len (self.horizons_phil.isoforms)==2:
        print("Rmsd gain over the other isoform %5.1f um."%(1000.*abs(positional_rmsds[0] - positional_rmsds[1])))
      else:
        print()
      R = isoform_refineries[minindex]
      # Now one last check to see if direct beam is out of bounds
      if self.horizons_phil.isoforms[minindex].beam_restraint is not None:
        refined_beam = matrix.col(R.get_experiments()[0].detector[0].get_beam_centre(experiments[0].beam.get_s0()))
        known_beam = matrix.col(self.horizons_phil.isoforms[minindex].beam_restraint)
        print("asserting",(refined_beam-known_beam).length(),"<",self.horizons_phil.isoforms[minindex].rmsd_target_mm)
        assert (refined_beam-known_beam).length() < self.horizons_phil.isoforms[minindex].rmsd_target_mm
        # future--circle of confusion could be given as a separate length in mm instead of reusing rmsd_target
      self.identified_isoform = self.horizons_phil.isoforms[minindex].name

    self.integration_concept_detail(experiments=R.get_experiments(), reflections=setting_reflections,
                                    spots=self.triclinic.get_observations_with_outlier_removal(),
                                    image_number=image_number, cb_op_to_primitive=cb_op_to_primitive, **kwargs)

  def integration_concept_detail(self, experiments, reflections, spots,image_number,cb_op_to_primitive,**kwargs):
    detector = experiments[0].detector
    crystal = experiments[0].crystal
    from cctbx.crystal import symmetry
    c_symmetry = symmetry(space_group = crystal.get_space_group(), unit_cell = crystal.get_unit_cell())

    self.image_number = image_number
    NEAR = 10
    pxlsz = detector[0].get_pixel_size()

    Predicted = self.get_predictions_accounting_for_centering(experiments,reflections,cb_op_to_primitive,**kwargs)

    FWMOSAICITY = self.inputai.getMosaicity()
    self.DOMAIN_SZ_ANG = kwargs.get("domain_size_ang",  self.__dict__.get("actual",0)  )
    refineflag = {True:0,False:1}[kwargs.get("domain_size_ang",0)==0]
    c_symmetry.show_summary(prefix="EXCURSION%1d REPORT FWMOS= %6.4f DOMAIN= %6.1f "%(refineflag,FWMOSAICITY,self.DOMAIN_SZ_ANG))
    from annlib_ext import AnnAdaptor
    self.cell = c_symmetry.unit_cell()

    query = flex.double()
    print(len(self.predicted))

    for pred in self.predicted: # predicted spot coord in pixels
      query.append(pred[0]/pxlsz[0])
      query.append(pred[1]/pxlsz[1])

    self.reserve_hkllist_for_signal_search = self.hkllist

    reference = flex.double()

    assert self.length>NEAR# Can't do spot/pred matching with too few spots
    for spot in spots:
      reference.append(spot.ctr_mass_x())
      reference.append(spot.ctr_mass_y())

    IS_adapt = AnnAdaptor(data=reference,dim=2,k=NEAR)
    IS_adapt.query(query)
    idx_cutoff = float(min(self.mask_focus[image_number]))

    from rstbx.apps.slip_helpers import slip_callbacks
    cache_refinement_spots = getattr(slip_callbacks.slip_callback,"requires_refinement_spots",False)

    indexed_pairs_provisional = []
    correction_vectors_provisional = []
    c_v_p_flex = flex.vec3_double()
    this_setting_matched_indices = reflections["miller_index"]
    for j,item in enumerate(this_setting_matched_indices):
      this_setting_index = self.hkllist.first_index(item)
      if this_setting_index:
        Match = dict(spot=j,pred=this_setting_index)
        indexed_pairs_provisional.append(Match)
        vector = matrix.col(
            [reflections["xyzobs.px.value"][j][0] - self.predicted[Match["pred"]][0]/pxlsz[0],
             reflections["xyzobs.px.value"][j][1] - self.predicted[Match["pred"]][1]/pxlsz[1]])
        correction_vectors_provisional.append(vector)
        c_v_p_flex.append((vector[0],vector[1],0.))
    self.N_correction_vectors = len(correction_vectors_provisional)
    self.rmsd_px = math.sqrt(flex.mean(c_v_p_flex.dot(c_v_p_flex)))
    print("... %d provisional matches"%self.N_correction_vectors, end=' ')
    print("r.m.s.d. in pixels: %6.3f"%(self.rmsd_px))

    if self.horizons_phil.integration.enable_residual_scatter:
      from matplotlib import pyplot as plt
      fig = plt.figure()
      for cv in correction_vectors_provisional:
        plt.plot([cv[1]],[-cv[0]],"r.")
      plt.title(" %d matches, r.m.s.d. %5.2f pixels"%(len(correction_vectors_provisional),math.sqrt(flex.mean(c_v_p_flex.dot(c_v_p_flex)))))
      plt.axes().set_aspect("equal")
      self.show_figure(plt,fig,"res")
      plt.close()

    if self.horizons_phil.integration.enable_residual_map:
      from matplotlib import pyplot as plt
      PX = reflections["xyzobs.px.value"]
      fig = plt.figure()
      for match,cv in zip(indexed_pairs_provisional,correction_vectors_provisional):
        plt.plot([PX[match["spot"]][1]],[-PX[match["spot"]][0]],"r.")
        plt.plot([self.predicted[match["pred"]][1]/pxlsz[1]],[-self.predicted[match["pred"]][0]/pxlsz[0]],"g.")
        plt.plot([PX[match["spot"]][1], PX[match["spot"]][1] + 10.*cv[1]],
                 [-PX[match["spot"]][0], -PX[match["spot"]][0] - 10.*cv[0]],'r-')
      if kwargs.get("user-reentrant") != None and self.horizons_phil.integration.spot_prediction == "dials" \
             and self.horizons_phil.integration.enable_residual_map_deltapsi:
        from rstbx.apps.stills.util import residual_map_special_deltapsi_add_on
        residual_map_special_deltapsi_add_on(
          reflections = self.dials_spot_prediction,
          matches = indexed_pairs_provisional, experiments=experiments,
          hkllist = self.hkllist,
          predicted = self.predicted, plot=plt, eta_deg=FWMOSAICITY, deff=self.DOMAIN_SZ_ANG
          )
      plt.xlim([0,detector[0].get_image_size()[1]])
      plt.ylim([-detector[0].get_image_size()[0],0])
      plt.title(" %d matches, r.m.s.d. %5.2f pixels"%(len(correction_vectors_provisional),math.sqrt(flex.mean(c_v_p_flex.dot(c_v_p_flex)))))
      plt.axes().set_aspect("equal")
      self.show_figure(plt,fig,"map")
      plt.close()

    indexed_pairs = indexed_pairs_provisional
    correction_vectors = correction_vectors_provisional
    ########### skip outlier rejection for this derived class

    ### However must retain the ability to write out correction vectiors.
    if True: # at Aaron's request; test later
      correction_lengths = flex.double([v.length() for v in correction_vectors_provisional])
      clorder = flex.sort_permutation(correction_lengths)
      sorted_cl = correction_lengths.select(clorder)
      indexed_pairs = []
      correction_vectors = []
      self.correction_vectors = []
      for icand in range(len(sorted_cl)):
        # somewhat arbitrary sigma = 1.0 cutoff for outliers
        indexed_pairs.append(indexed_pairs_provisional[clorder[icand]])
        correction_vectors.append(correction_vectors_provisional[clorder[icand]])
        if cache_refinement_spots:
          self.spotfinder.images[self.frame_numbers[self.image_number]]["refinement_spots"].append(
          spots[reflections[indexed_pairs[-1]["spot"]]['spotfinder_lookup']])
        if kwargs.get("verbose_cv")==True:
            print("CV OBSCENTER %7.2f %7.2f REFINEDCENTER %7.2f %7.2f"%(
              float(self.inputpd["size1"])/2.,float(self.inputpd["size2"])/2.,
              self.inputai.xbeam()/pxlsz[0], self.inputai.ybeam()/pxlsz[1]), end=' ')
            print("OBSSPOT %7.2f %7.2f PREDSPOT %7.2f %7.2f"%(
              reflections[indexed_pairs[-1]["spot"]]['xyzobs.px.value'][0],
              reflections[indexed_pairs[-1]["spot"]]['xyzobs.px.value'][1],
              self.predicted[indexed_pairs[-1]["pred"]][0]/pxlsz[0],
              self.predicted[indexed_pairs[-1]["pred"]][1]/pxlsz[1]), end=' ')
            the_hkl = self.hkllist[indexed_pairs[-1]["pred"]]
            print("HKL %4d %4d %4d"%the_hkl,"%2d"%self.setting_id, end=' ')
            radial, azimuthal = spots[indexed_pairs[-1]["spot"]].get_radial_and_azimuthal_size(
              self.inputai.xbeam()/pxlsz[0], self.inputai.ybeam()/pxlsz[1])
            print("RADIALpx %5.3f AZIMUTpx %5.3f"%(radial,azimuthal))

        # Store a list of correction vectors in self.
        radial, azimuthal = spots[indexed_pairs[-1]['spot']].get_radial_and_azimuthal_size(
          self.inputai.xbeam()/pxlsz[0], self.inputai.ybeam()/pxlsz[1])
        self.correction_vectors.append(
          dict(obscenter=(float(self.inputpd['size1']) / 2,
                          float(self.inputpd['size2']) / 2),
               refinedcenter=(self.inputai.xbeam() / pxlsz[0],
                              self.inputai.ybeam() / pxlsz[1]),
               obsspot=(reflections[indexed_pairs[-1]['spot']]['xyzobs.px.value'][0],
                        reflections[indexed_pairs[-1]['spot']]['xyzobs.px.value'][1]),
               predspot=(self.predicted[indexed_pairs[-1]['pred']][0] / pxlsz[0],
                         self.predicted[indexed_pairs[-1]['pred']][1] / pxlsz[1]),
               hkl=(self.hkllist[indexed_pairs[-1]['pred']][0],
                    self.hkllist[indexed_pairs[-1]['pred']][1],
                    self.hkllist[indexed_pairs[-1]['pred']][2]),
               setting_id=self.setting_id,
               radial=radial,
               azimuthal=azimuthal))


    self.inputpd["symmetry"] = c_symmetry
    self.inputpd["symmetry"].show_summary(prefix="SETTING ")


    if self.horizons_phil.integration.model == "user_supplied":
      # Not certain of whether the reentrant_* dictionary keys create a memory leak
      if kwargs.get("user-reentrant",None)==None:
        kwargs["reentrant_experiments"] = experiments
        kwargs["reentrant_reflections"] = reflections
        from cxi_user import post_outlier_rejection
        self.indexed_pairs = indexed_pairs
        self.spots = spots
        post_outlier_rejection(self,image_number,cb_op_to_primitive,self.horizons_phil,kwargs)
        return
    ########### finished with user-supplied code


    correction_lengths=flex.double([v.length() for v in correction_vectors])

    self.r_residual = pxlsz[0]*flex.mean(correction_lengths)

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
    # do not use null: self.null_correction_mapping( predicted=self.predicted,
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

  def get_predictions_accounting_for_centering(self,experiments,reflections,cb_op_to_primitive,**kwargs):
    # interface requires this function to set current_orientation
    # in the actual setting used for Miller index calculation
    detector = experiments[0].detector
    crystal = experiments[0].crystal

    if self.horizons_phil.integration.model == "user_supplied":

      lower_limit_domain_size = math.pow(
       crystal.get_unit_cell().volume(),
       1./3.)*self.horizons_phil.integration.mosaic.domain_size_lower_limit # default 10-unit cell block size minimum reasonable domain
      actual_used_domain_size = kwargs.get("domain_size_ang",lower_limit_domain_size)

      self.block_counter+=1
      rot_mat = matrix.sqr(cb_op_to_primitive.c().r().as_double()).transpose()

      from cctbx.crystal_orientation import crystal_orientation, basis_type
      centered_orientation = crystal_orientation(crystal.get_A(),basis_type.reciprocal)
      self.current_orientation = centered_orientation
      self.current_cb_op_to_primitive = cb_op_to_primitive
      primitive_orientation = centered_orientation.change_basis(rot_mat)

      self.inputai.setOrientation(primitive_orientation)
      from cxi_user import pre_get_predictions
      if self.block_counter < 2:
        KLUDGE = self.horizons_phil.integration.mosaic.kludge1 # bugfix 1 of 2 for protocol 6, equation 2
        self.inputai.setMosaicity(KLUDGE*self.inputai.getMosaicity())

      oldbase = self.inputai.getBase()

      #print oldbase.xbeam, oldbase.ybeam
      newbeam = detector[0].get_beam_centre(experiments[0].beam.get_s0())
      newdistance = -detector[0].get_beam_centre_lab(experiments[0].beam.get_s0())[2]

      from labelit.dptbx import Parameters
      base = Parameters(xbeam = newbeam[0], ybeam = newbeam[1], #oldbase.xbeam, ybeam = oldbase.ybeam,
                        distance = newdistance, twotheta = 0.0)
      self.inputai.setBase(base)
      self.inputai.setWavelength(experiments[0].beam.get_wavelength())

      self.bp3_wrapper = pre_get_predictions(self.inputai, self.horizons_phil,
        raw_image = self.imagefiles.images[self.image_number],
        imageindex = self.frame_numbers[self.image_number],
        spotfinder = self.spotfinder,
        limiting_resolution = self.limiting_resolution,
        domain_size_ang = actual_used_domain_size,
        )

      BPpredicted = self.bp3_wrapper.ucbp3.selected_predictions_labelit_format()
      BPhkllist = self.bp3_wrapper.ucbp3.selected_hkls()

      self.actual = actual_used_domain_size
      primitive_hkllist = BPhkllist
      #not sure if matrix needs to be transposed first for outputting HKL's???:
      self.hkllist = cb_op_to_primitive.inverse().apply(primitive_hkllist)

      if self.horizons_phil.integration.spot_prediction == "dials":
        from dials.algorithms.spot_prediction import StillsReflectionPredictor
        predictor = StillsReflectionPredictor(experiments[0])
        Rcalc = flex.reflection_table.empty_standard(len(self.hkllist))
        Rcalc['miller_index'] = self.hkllist
        predictor.for_reflection_table(Rcalc, crystal.get_A())
        self.predicted = Rcalc['xyzcal.mm']
        self.dials_spot_prediction = Rcalc
        self.dials_model = experiments

      elif self.horizons_phil.integration.spot_prediction == "ucbp3":
        self.predicted = BPpredicted

      self.inputai.setOrientation(centered_orientation)
      if self.inputai.active_areas != None:
        self.predicted,self.hkllist = self.inputai.active_areas(
                                      self.predicted,self.hkllist,self.pixel_size)
      if self.block_counter < 2:
         down = self.inputai.getMosaicity()/KLUDGE
         print("Readjusting mosaicity back down to ",down)
         self.inputai.setMosaicity(down)
      return

  def refine(self, experiments, reflections, isoform=None):
    from dials.algorithms.refinement.refiner import phil_scope
    from libtbx.phil import parse

    params = phil_scope.fetch(source=parse('')).extract()
    params.refinement.reflections.weighting_strategy.delpsi_constant=100000.
    params.refinement.reflections.weighting_strategy.override="stills"
    params.refinement.parameterisation.auto_reduction.action="fix"
    #params.refinement.reflections.do_outlier_rejection=True
    #params.refinement.reflections.iqr_multiplier=0.5
    #params.refinement.reflections.minimum_sample_size=50
    #params.refinement.reflections.maximum_sample_size=50
    #params.refinement.reflections.random_seed=1
    if self.horizons_phil.integration.dials_refinement.strategy=="distance":
      params.refinement.parameterisation.beam.fix="all"
      params.refinement.parameterisation.detector.fix_list=["Tau1"] # fix detector rotz, allow distance to refine
    elif self.horizons_phil.integration.dials_refinement.strategy=="wavelength":
      params.refinement.parameterisation.beam.fix="in_spindle_plane,out_spindle_plane"
      params.refinement.parameterisation.detector.fix_list=["Dist","Tau1"] # fix detector rotz and distance
    elif self.horizons_phil.integration.dials_refinement.strategy=="fix":
      params.refinement.parameterisation.beam.fix="all"
      params.refinement.parameterisation.detector.fix_list=["Dist","Tau1"] # fix detector rotz and distance
    if isoform is not None:
      params.refinement.reflections.outlier.algorithm="null"
      params.refinement.parameterisation.crystal.fix="cell"
    from dials.algorithms.refinement.refiner import RefinerFactory
    refiner = RefinerFactory.from_parameters_data_experiments(params,
      reflections, experiments)

    history = refiner.run()
    print(history.keys())
    for item in history["rmsd"]:
      print("%5.2f %5.2f %8.5f"%(item[0],item[1],180.*item[2]/math.pi))

    #for item in history["parameter_vector"]:
    #  print ["%8.5f"%a for a in item]
    print(refiner.selection_used_for_refinement().count(True),"spots used for refinement")
    print(refiner.get_experiments()[0].beam)
    print(refiner.get_experiments()[0].detector)
    print("Distance:", -refiner.get_experiments()[0].detector[0].get_beam_centre_lab(refiner.get_experiments()[0].beam.get_s0())[2])
    print(refiner.get_experiments()[0].crystal)
    return refiner

  def prepare_reflection_list(self,detector):

    spots = self.triclinic.get_observations_with_outlier_removal()
    ordinary_python_list_of_indexed_observations = [
      {
        "id":0,
        "panel":0,
        "miller_index":item["pred"],
        "xyzobs.px.value":(spots[item["spot"]].ctr_mass_x(),spots[item["spot"]].ctr_mass_y(),0.0),
        "xyzobs.px.variance":(0.25,0.25,0.25),
        "spotfinder_lookup":item["spot"]
      }
      for item in self.triclinic_pairs
    ]

    self.length = len(ordinary_python_list_of_indexed_observations)
    R= flex.reflection_table.empty_standard(self.length)

    R['miller_index'] = flex.miller_index([item["miller_index"] for item in ordinary_python_list_of_indexed_observations])
    R['xyzobs.px.value'] = flex.vec3_double([item["xyzobs.px.value"] for item in ordinary_python_list_of_indexed_observations])
    R['xyzobs.px.variance'] = flex.vec3_double([item["xyzobs.px.variance"] for item in ordinary_python_list_of_indexed_observations])
    R['spotfinder_lookup'] = flex.int([item["spotfinder_lookup"] for item in ordinary_python_list_of_indexed_observations])

    R['xyzobs.mm.value'] = flex.vec3_double(self.length)
    R['xyzobs.mm.variance'] = flex.vec3_double(self.length)

    pxlsz = detector[0].get_pixel_size()

    for idx in range(self.length):
      R['xyzobs.mm.value'][idx] = (R['xyzobs.px.value'][idx][0]*pxlsz[0], R['xyzobs.px.value'][idx][1]*pxlsz[1], R['xyzobs.px.value'][idx][2])
      R['xyzobs.mm.variance'][idx] = (R['xyzobs.px.variance'][idx][0]*pxlsz[0], R['xyzobs.px.variance'][idx][1]*pxlsz[1], R['xyzobs.px.variance'][idx][2])

    return R

  def prepare_dxtbx_models(self,setting_specific_ai,sg,isoform=None):

    from dxtbx.model import BeamFactory
    beam = BeamFactory.simple(wavelength = self.inputai.wavelength)

    from dxtbx.model import DetectorFactory
    detector = DetectorFactory.simple(
      sensor = DetectorFactory.sensor("PAD"),
      distance = setting_specific_ai.distance(),
      beam_centre = [setting_specific_ai.xbeam(), setting_specific_ai.ybeam()],
      fast_direction = "+x",
      slow_direction = "+y",
      pixel_size = [self.pixel_size,self.pixel_size],
      image_size = [self.inputpd['size1'],self.inputpd['size1']],
      )

    direct = matrix.sqr(setting_specific_ai.getOrientation().direct_matrix())
    from dxtbx.model import MosaicCrystalKabsch2010
    crystal = MosaicCrystalKabsch2010(
      real_space_a = matrix.row(direct[0:3]),
      real_space_b = matrix.row(direct[3:6]),
      real_space_c = matrix.row(direct[6:9]),
      space_group_symbol = sg,
    )
    crystal.set_mosaicity(setting_specific_ai.getMosaicity())
    if isoform is not None:
      newB = matrix.sqr(isoform.fractionalization_matrix()).transpose()
      crystal.set_B(newB)

    from dxtbx.model import Experiment, ExperimentList
    experiments = ExperimentList()
    experiments.append(Experiment(beam=beam,
                                  detector=detector,
                                  crystal=crystal))

    print(beam)
    print(detector)
    print(crystal)
    return experiments
