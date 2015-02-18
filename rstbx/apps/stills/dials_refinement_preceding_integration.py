from __future__ import division
from rstbx.apps.stills.simple_integration import IntegrationMetaProcedure
from rstbx.apps import simple_integration
from scitbx import matrix
import math,copy
from dials.array_family import flex

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

    setting_experiments = self.prepare_dxtbx_models(setting_specific_ai = self.inputai,
      sg = self.inputpd["symmetry"].space_group_info().type().lookup_symbol())
    setting_reflections = copy.deepcopy(reflections)
    setting_reflections["miller_index"] = cb_op_to_primitive.inverse().apply(reflections["miller_index"])
    R = self.refine(experiments=setting_experiments, reflections=setting_reflections)

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
    DOMAIN_SZ_ANG = kwargs.get("domain_size_ang",  self.__dict__.get("actual",0)  )
    refineflag = {True:0,False:1}[kwargs.get("domain_size_ang",0)==0]
    c_symmetry.show_summary(prefix="EXCURSION%1d REPORT FWMOS= %6.4f DOMAIN= %6.1f "%(refineflag,FWMOSAICITY,DOMAIN_SZ_ANG))
    from annlib_ext import AnnAdaptor
    self.cell = c_symmetry.unit_cell()

    query = flex.double()
    print len(self.predicted)

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
    print "... %d provisional matches"%len(correction_vectors_provisional),
    print "r.m.s.d. in pixels: %5.2f"%(math.sqrt(flex.mean(c_v_p_flex.dot(c_v_p_flex))))

    indexed_pairs = indexed_pairs_provisional
    correction_vectors = correction_vectors_provisional
    ########### skip outlier rejection for this derived class

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

    self.BSmasks = []
    self.null_correction_mapping( predicted=self.predicted,
                                        correction_vectors = correction_vectors,
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
    for i in xrange(len(self.predicted)): # loop over predicteds
      B_S_mask = {}
      keys = self.get_bsmask(i)
      for k in xrange(0,len(keys),2):
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

      from labelit.dptbx import Parameters
      base = Parameters(xbeam = newbeam[0], ybeam = newbeam[1], #oldbase.xbeam, ybeam = oldbase.ybeam,
                        distance = -detector[0].get_distance(), twotheta = 0.0)
      self.inputai.setBase(base)

      #print matrix.col(detector[0].get_slow_axis()).dot(matrix.col((0.,1.,0.)))
      #print matrix.col(detector[0].get_fast_axis()).dot(matrix.col((1.,0.,0.)))

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
         print "Readjusting mosaicity back down to ",down
         self.inputai.setMosaicity(down)
      return

  def refine(self, experiments, reflections):
    from dials.algorithms.refinement.refiner import phil_scope
    from libtbx.phil import parse
    params = phil_scope.fetch(source=parse('')).extract()
    params.refinement.parameterisation.beam.fix="all"
    #params.refinement.parameterisation.detector.fix="orientation"
    params.refinement.parameterisation.detector.fix_list=0,3
    params.refinement.reflections.weighting_strategy.delpsi_constant=100000.
    params.refinement.reflections.weighting_strategy.override="stills"

    from dials.algorithms.refinement.refiner import RefinerFactory
    refiner = RefinerFactory.from_parameters_data_experiments(params,
      reflections, experiments, verbosity=1)

    history = refiner.run()
    print history.keys()
    for item in history["rmsd"]:
      print "%5.2f %5.2f %8.5f"%(item[0],item[1],180.*item[2]/math.pi)

    #for item in history["parameter_vector"]:
    #  print ["%8.5f"%a for a in item]
    print refiner.selection_used_for_refinement().count(True),"spots used for refinement"
    print refiner.get_experiments()[0].detector
    print refiner.get_experiments()[0].crystal
    return refiner

  def prepare_reflection_list(self,detector):

    spots = self.triclinic.get_observations_with_outlier_removal()
    ordinary_python_list_of_indexed_observations = [
      {
        "id":0,
        "panel":0,
        "miller_index":item["pred"],
        "xyzobs.px.value":(spots[item["spot"]].ctr_mass_x(),spots[item["spot"]].ctr_mass_y(),0.0),
        "xyzobs.px.variance":(0.25,0.25,0.25)
      }
      for item in self.triclinic_pairs
    ]

    self.length = len(ordinary_python_list_of_indexed_observations)
    R= flex.reflection_table.empty_standard(self.length)

    R['miller_index'] = flex.miller_index([item["miller_index"] for item in ordinary_python_list_of_indexed_observations])
    R['xyzobs.px.value'] = flex.vec3_double([item["xyzobs.px.value"] for item in ordinary_python_list_of_indexed_observations])
    R['xyzobs.px.variance'] = flex.vec3_double([item["xyzobs.px.variance"] for item in ordinary_python_list_of_indexed_observations])

    R['xyzobs.mm.value'] = flex.vec3_double(self.length)
    R['xyzobs.mm.variance'] = flex.vec3_double(self.length)

    pxlsz = detector[0].get_pixel_size()

    for idx in xrange(self.length):
      R['xyzobs.mm.value'][idx] = (R['xyzobs.px.value'][idx][0]*pxlsz[0], R['xyzobs.px.value'][idx][1]*pxlsz[1], R['xyzobs.px.value'][idx][2])
      R['xyzobs.mm.variance'][idx] = (R['xyzobs.px.variance'][idx][0]*pxlsz[0], R['xyzobs.px.variance'][idx][1]*pxlsz[1], R['xyzobs.px.variance'][idx][2])

    return R

  def prepare_dxtbx_models(self,setting_specific_ai,sg):

    from dxtbx.model.beam import beam_factory
    beam = beam_factory.simple(wavelength = self.inputai.wavelength)

    from dxtbx.model.detector import detector_factory
    detector = detector_factory.simple(
      sensor = detector_factory.sensor("PAD"),
      distance = setting_specific_ai.distance(),
      beam_centre = [setting_specific_ai.xbeam(), setting_specific_ai.ybeam()],
      fast_direction = "+x",
      slow_direction = "+y",
      pixel_size = [self.pixel_size,self.pixel_size],
      image_size = [self.inputpd['size1'],self.inputpd['size1']],
      )

    direct = matrix.sqr(setting_specific_ai.getOrientation().direct_matrix())
    from dxtbx.model.crystal import crystal_model
    crystal = crystal_model(
      real_space_a = matrix.row(direct[0:3]),
      real_space_b = matrix.row(direct[3:6]),
      real_space_c = matrix.row(direct[6:9]),
      space_group_symbol = sg,
      mosaicity = setting_specific_ai.getMosaicity()
    )

    from dxtbx.model.experiment.experiment_list import Experiment, ExperimentList
    experiments = ExperimentList()
    experiments.append(Experiment(beam=beam,
                                  detector=detector,
                                  crystal=crystal))

    print beam
    print detector
    print crystal
    return experiments
