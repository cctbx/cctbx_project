import pickle,math
from labelit.preferences import labelit_commands,labelit_phil
from labelit.command_line.default_param import establish_dict_for_refinement
from labelit.dptbx.autoindex import index_and_refine
from labelit.command_line.stats_index import best_character_to_IndexPrinter
from labelit.command_line.stats_index import AutoIndexOrganizer
from cctbx.array_family import flex
from scitbx import matrix
#234567890123456789212345678931234567894123456789512345678961234567897123456789812

class Empty:pass

class api:
  def __init__(self,*args):
    #application programmer interface for screen, using identical inputs
    # to the command line interface
    args = list(args)
    labelit_phil.merge_command_line(args)
    E = Empty()
    E.argv=['Empty']
    for x in xrange(len(args)):
      E.argv.append(args[x])

    self.establish_dict_for_refinement = establish_dict_for_refinement
    self.index_and_refine = index_and_refine
    self.best_character_to_IndexPrinter = best_character_to_IndexPrinter

    self.Org = AutoIndexOrganizer(
      verbose=labelit_commands.distl.bins.verbose,argument_module=E)
    self.Org.setIndexingDelegate(self.index_and_integrate)
      # legacy:algorithm control could be excercised here

  def __call__(self):
    return self.Org.process()

  def index_and_integrate(self,frames,files,spotfinder_results):
    self.frames = frames
    self.files = files
    self.spotfinder_results = spotfinder_results
    pd = self.establish_dict_for_refinement(frames,spotfinder_results)
    #------------------------------------------------------------
    ai,P = self.index_and_refine(pd,files,spotfinder_results,0)
    self.indexing_ai = ai
    #------------------------------------------------------------
    if labelit_commands.compatibility_allow==False:
      M = self.best_character_to_IndexPrinter(ai,P,pd)
    else:
      from labelit.diffraction.compatibility import \
        best_compatibility_to_IndexPrinter
      M = best_compatibility_to_IndexPrinter(ai,P,pd,files,spotfinder_results)
    #------------------------------------------------------------
    if labelit_commands.__dict__.has_key("writer"):
      labelit_commands.writer.make_image_plots_detail(
        ai=ai,pd=pd,inframes=files,spotfinder_results=spotfinder_results)
    if not labelit_commands.index_only:
      IC = IntegrateCharacters(M,pd)
      IC.write_mosflm_matrices()
      IC.find_best()
      IC.show()
    return pd

  def create_case_only(self,frames,file,spotfinder_results):
    self.pd = self.establish_dict_for_refinement(frames)

  def analyze_one(self,solution):
    inputpd = self.Org.process()

    settings = pickle.load(open("LABELIT_possible","rb"))
    setting = [setting for setting in settings if setting["counter"]==solution][0]

    from labelit.preferences import labelit_commands as param

    pixel_size = float(inputpd['pixel_size'])
    self.pixel_size = pixel_size

    from labelit.dptbx import AutoIndexEngine, Parameters
    ai = AutoIndexEngine(inputpd['endstation'])

    P = Parameters(xbeam=setting["refined x beam"],ybeam=setting["refined y beam"],
             distance=setting["refined distance"],twotheta=float(inputpd["twotheta"]))
    ai.setBase(P)
    ai.setWavelength(float(inputpd['wavelength']))
    ai.setMaxcell(float(inputpd['ref_maxcel']))
    print "Deltaphi is",float(inputpd['deltaphi'])
    ai.setDeltaphi(float(inputpd['deltaphi'])*math.pi/180.)
    ai.setMosaicity(setting["mosaicity"])
    ai.setOrientation(setting["orient"])
    #why aren't hexagonal constraints applied here???
    print inputpd["osc_start"]

    image_centers = [(math.pi/180.)*float(x) for x in inputpd["osc_start"].values()]

    limiting_resolution = param.distl_highres_limit
    print "Limiting resolution",limiting_resolution

    #predict the spots
    pre2m = ai.predict_all(image_centers[0],limiting_resolution)
    self.pre2m = pre2m

    hkllist = ai.hklpredict()
    cell = ai.getOrientation().unit_cell()
    print cell
    for hkl in hkllist:
      #print "%25s %5.2f"%(str(hkl),cell.d(hkl))
      assert cell.d(hkl)>=limiting_resolution
    print "Number of hkls:",(hkllist).size(),
    print "all inside the %4.2f Angstrom limiting sphere."%limiting_resolution
    print "The unit cell is",cell
    self.solution_setting_ai = ai
    self.solution_pd = inputpd
    self.image_centers = image_centers
    return [ai.getOrientation().unit_cell(),hkllist]

  def parameters(self):
    image_centers = [(math.pi/180.)*float(x) for x in self.solution_pd["osc_start"].values()]
    P = dict( image_center_radians = image_centers,
              wavelength = float(self.solution_pd['wavelength']),
              deltaphi = float(self.solution_pd['deltaphi']),
              mosaicity_degrees = self.solution_setting_ai.getMosaicity(),
              orientation = self.solution_setting_ai.getOrientation(),
              )
    return P

def show_observations(obs,out=None):
  if out==None:
    import sys
    out = sys.stdout
  from libtbx.str_utils import format_value

  obs.setup_binner(n_bins = 12)
  result = []
  for i_bin in obs.binner().range_used():
    sel_w = obs.binner().selection(i_bin)
    sel_fo_all = obs.select(sel_w)
    d_max_,d_min_ = sel_fo_all.d_max_min()
    d_range = obs.binner().bin_legend(
      i_bin=i_bin, show_bin_number=False, show_counts=True)
    sel_data = obs.select(sel_w).data()
    sel_sig = obs.select(sel_w).sigmas()

    if(sel_data.size() > 0):
      bin = resolution_bin(
        i_bin        = i_bin,
        d_range      = d_range,
        mean_I       = flex.mean(sel_data),
        n_work       = sel_data.size(),
        mean_I_sigI  = flex.mean(sel_data/sel_sig),
        )
      result.append(bin)
  print >>out, "\n Bin  Resolution Range  Compl.         <I>     <I/sig(I)>"
  for bin in result:
    fmt = " %s %s    %s  %s"
    print >>out,fmt%(
      format_value("%3d",   bin.i_bin),
      format_value("%-17s", bin.d_range),
      format_value("%8.1f", bin.mean_I),
      format_value("%8.1f", bin.mean_I_sigI),
      )
class resolution_bin(object):
  def __init__(self,
               i_bin         = None,
               d_range       = None,
               completeness  = None,
               alpha_work    = None,
               beta_work     = None,
               mean_I        = None,
               mean_I_sigI   = None,
               target_work   = None,
               target_free   = None,
               n_work        = None,
               n_free        = None,
               mean_f_obs    = None,
               fom_work      = None,
               scale_k1_work = None,
               pher_work     = None,
               pher_free     = None,
               sigmaa        = None):
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())

from labelit.dptbx.profile_support import show_profile
from libtbx.development.timers import Timer
from rstbx.apps import simple_integration

class IntegrationMetaProcedure(simple_integration):
  def __init__(self,inputs): # inputs is an instance of class api
    simple_integration.__init__(self)
    self.inputai = inputs.solution_setting_ai #C++ autoindex engine
    self.indexing_ai = inputs.indexing_ai

    #Note...here we may be working with a non-primitive cell, so
    # must work with systematic absences...not implemented yet.
    if self.indexing_ai.getData().size() < 40: return # initial protection

    self.inputpd = inputs.solution_pd #parameter dictionary
    self.inputframes = inputs.frames
    self.imagefiles = inputs.files
    self.spotfinder = inputs.spotfinder_results
    self.image_centers = inputs.image_centers
    print "IMAGE CENTERS",self.image_centers

    # initial resolution from DISTL
    resolution_est = float(self.inputpd['resolution_inspection'])
    print "initial resolution from DISTL",resolution_est

    # resolution limit of the strong spots used for indexing
    resolution_str = self.indexing_ai.high().reciprocal_spacing
    resolution = max(resolution_est,resolution_str)
    print "resolution limit of the strong spots used for indexing",resolution
    self.limiting_resolution = resolution

    #print "resolution: %.2f %.2f"%(resolution_est,resolution_str)
    self.pixel_size = inputs.pixel_size
    self.set_pixel_size(inputs.pixel_size)
    self.set_detector_size(int(self.inputpd["size1"]),int(self.inputpd["size2"]))

    self.pre2m = inputs.pre2m
    self.basic_algorithm()
    self.initialize_increments()
    T = Timer("concept")
    self.integration_concept()
    T = Timer("proper")
    self.integration_proper()

  def basic_algorithm(self,verbose=False):
    Amat = matrix.sqr(self.inputai.getOrientation().direct_matrix())
    self.frames = self.inputpd['osc_start'].keys()
    self.incr_focus = []
    for frame in self.frames:
      focus = self.inputpd['masks'][frame][0:2]
      average_profile = self.inputpd['masks'][frame][2]
      box = self.inputpd['masks'][frame][3]
      if verbose:
        print average_profile.focus()
        print box.focus()
        print "Average Profile:"
        show_profile( average_profile )
        print "Box:"
        show_profile( box )
      self.incr_focus.append( average_profile.focus() )

  """Integration concept.  Focus on each predicted spot position S.
     I(S): the integration mask for S is constructed by superimposing the
           bodypixels of the 10 nearest spotfinder spots; thus getting
           the maximum envelope.
     B(S): the background mask for S is obtained by finding an equal number
           of pixels as are in I(S).  Choose the nearest pixels to I(S)
           that are not within an exclusion distance (guard=3 pixels).
     P(S): A positional correction for I(S)--adjust the position of the mask
           away from the predicted position, before integrating.  This is constructed
           by considering the 10 closest tuples of (spotfinder spot,prediction).
           P(S) is the average positional offset for this set of 10.
     O(S): The set of spots close enough to S that they must be included
           in calculating where B(S) can be sampled.
  """

  def get_predictions_accounting_for_centering(self,cb_op_to_primitive=None):
    if cb_op_to_primitive==None:

      predicted = self.inputai.predict_all(
                  self.image_centers[self.image_number],self.limiting_resolution)
      self.predicted = predicted #only good for integrating one frame...
      self.hkllist = self.inputai.hklpredict()

    else:
      rot_mat = matrix.sqr(cb_op_to_primitive.c().r().as_double()).transpose()
      centered_orientation = self.inputai.getOrientation()
      primitive_orientation = centered_orientation.change_basis(rot_mat)
      self.inputai.setOrientation(primitive_orientation)
      predicted = self.inputai.predict_all(
                  self.image_centers[self.image_number],self.limiting_resolution)
      self.predicted = predicted #only good for integrating one frame...
      primitive_hkllist = self.inputai.hklpredict()
      #not sure if matrix needs to be transposed first for outputting HKL's???:
      self.hkllist = cb_op_to_primitive.inverse().apply(primitive_hkllist)
      self.inputai.setOrientation(centered_orientation)
    if self.inputai.active_areas != None:
      self.predicted,self.hkllist = self.inputai.active_areas(
                                    self.predicted,self.hkllist,self.pixel_size)

  def get_observations_with_outlier_removal(self):
    spots = self.spotfinder.images[self.frames[self.image_number]]["inlier_spots"]
    return spots

  def integration_concept(self,image_number=0,cb_op_to_primitive=None,verbose=False,**kwargs):
    self.image_number = image_number
    NEAR = 10
    pxlsz = self.pixel_size
    self.get_predictions_accounting_for_centering(cb_op_to_primitive)
    from annlib_ext import AnnAdaptor
    self.cell = self.inputai.getOrientation().unit_cell()
    query = flex.double()
    for pred in self.predicted: # predicted spot coord in pixels
      query.append(pred[0]/pxlsz)
      query.append(pred[1]/pxlsz)

    reference = flex.double()
    spots = self.get_observations_with_outlier_removal()

    assert len(spots)>NEAR# Can't do spot/pred matching with too few spots
    for spot in spots:
      reference.append(spot.ctr_mass_x())
      reference.append(spot.ctr_mass_y())

    IS_adapt = AnnAdaptor(data=reference,dim=2,k=NEAR)
    IS_adapt.query(query)

    indexed_pairs_provisional = []
    correction_vectors_provisional = []
    idx_cutoff = float(min(self.inputpd['masks'][self.frames[self.image_number]][0:2]))
    if verbose:
      print "idx_cutoff distance in pixels",idx_cutoff
    for i in xrange(len(self.predicted)): # loop over predicteds
      for n in xrange(NEAR): # loop over near spotfinder spots
        Match = dict(spot=IS_adapt.nn[i*NEAR+n],pred=i)
        if n==0 and math.sqrt(IS_adapt.distances[i*NEAR+n]) < idx_cutoff:
          indexed_pairs_provisional.append(Match)

          vector = matrix.col(
            [spots[Match["spot"]].ctr_mass_x() - self.predicted[Match["pred"]][0]/pxlsz,
             spots[Match["spot"]].ctr_mass_y() - self.predicted[Match["pred"]][1]/pxlsz])
          correction_vectors_provisional.append(vector)

    #insert code here to remove correction length outliers...
    # they are causing terrible
    # problems for finding legitimate correction vectors (print out the list)
    # also remove outliers for the purpose of reporting RMS
    outlier_rejection = True
    if outlier_rejection:
      correction_lengths = flex.double([v.length() for v in correction_vectors_provisional])
      clorder = flex.sort_permutation(correction_lengths)
      sorted_cl = correction_lengths.select(clorder)

      limit = int(0.33 * len(sorted_cl)) # best 1/3 of data are assumed to be correctly modeled.
      y_data = flex.double(len(sorted_cl))
      for i in xrange(len(y_data)):
        y_data[i] = float(i)/float(len(y_data))

      # ideas are explained in Sauter & Poon (2010) J Appl Cryst 43, 611-616.
      from labelit.outlier_spots.fit_distribution import fit_cdf,rayleigh
      fitted_rayleigh = fit_cdf(x_data = sorted_cl[0:limit],
                                y_data = y_data[0:limit],
                                distribution=rayleigh)

      inv_cdf = [fitted_rayleigh.distribution.inv_cdf(cdf) for cdf in y_data]

      #print "SORTED LIST OF ",len(sorted_cl), "with sigma",fitted_rayleigh.distribution.sigma
      indexed_pairs = []
      correction_vectors = []
      for icand in xrange(len(sorted_cl)):
        # somewhat arbitrary sigma = 1.0 cutoff for outliers
        if (sorted_cl[icand]-inv_cdf[icand])/fitted_rayleigh.distribution.sigma > 1.0:
          break
        indexed_pairs.append(indexed_pairs_provisional[clorder[icand]])
        correction_vectors.append(correction_vectors_provisional[clorder[icand]])

        if kwargs.get("verbose_cv")==True:
            print "CV OBSCENTER %7.2f %7.2f REFINEDCENTER %7.2f %7.2f"%(
              float(self.inputpd["size1"])/2.,float(self.inputpd["size2"])/2.,
              self.inputai.xbeam()/pxlsz, self.inputai.ybeam()/pxlsz),
            print "OBSSPOT %7.2f %7.2f PREDSPOT %7.2f %7.2f"%(
              spots[indexed_pairs[-1]["spot"]].ctr_mass_x(),
              spots[indexed_pairs[-1]["spot"]].ctr_mass_y(),
              self.predicted[indexed_pairs[-1]["pred"]][0]/pxlsz,
              self.predicted[indexed_pairs[-1]["pred"]][1]/pxlsz)
      #print "After outlier rejection %d indexed spotfinder spots remain."%len(indexed_pairs)
      if False:
        rayleigh_cdf = [
          fitted_rayleigh.distribution.cdf(x=sorted_cl[c]) for c in xrange(len(sorted_cl))]
        from matplotlib import pyplot as plt
        plt.plot(sorted_cl,y_data,"r+")
        #plt.plot(sorted_cl,rayleigh_cdf,"g.")
        plt.plot(inv_cdf,y_data,"b.")
        plt.show()
    else:
      indexed_pairs = indexed_pairs_provisional
      correction_vectors = correction_vectors_provisional
    ########### finished with outlier rejection

    #Other checks to be implemented (future):
    # spot is within active area of detector on a circular detector such as the Mar IP
    # integration masks do not overlap; or deconvolute

    correction_lengths=flex.double([v.length() for v in correction_vectors])
    if verbose:
      print "average correction %5.2f over %d vectors"%(flex.mean(correction_lengths),
      len(correction_lengths)),
      print "or %5.2f mm."%(pxlsz*flex.mean(correction_lengths))
    self.r_residual = pxlsz*flex.mean(correction_lengths)

    assert len(indexed_pairs)>NEAR # must have enough indexed spots
    reference = flex.double()
    for item in indexed_pairs:
      reference.append(spots[item["spot"]].ctr_mass_x())
      reference.append(spots[item["spot"]].ctr_mass_y())

    PS_adapt = AnnAdaptor(data=reference,dim=2,k=NEAR)
    PS_adapt.query(query)

    self.BSmasks = []
    self.positional_correction_mapping( predicted=self.predicted,
                                        correction_vectors = correction_vectors,
                                        PS_adapt = PS_adapt,
                                        IS_adapt = IS_adapt,
                                        spots = spots)
    """this section is now pushed down to C++
    self.ISmasks = []
    corrections = []
    for i in xrange(len(predicted)): # loop over predicteds
      # calculate the positional correction for this prediction
      # ....average over the 10 nearest positional corrections
      correction = matrix.col([0.0,0.0])
      for n in xrange(NEAR): # loop over near indexed pairs
        correction += correction_vectors[PS_adapt.nn[i*NEAR+n]]
      correction/=NEAR

      I_S_mask = {} #indexed by tuples giving relative xy positions.
      pred = predicted[i]
      predX = pred[0]/pxlsz
      predY = pred[1]/pxlsz
      for n in xrange(NEAR): # loop over near spotfinder spots
        Match = dict(spot=IS_adapt.nn[i*NEAR+n],pred=i)
        #print "       index %6d distance %4.1f"%(
        #      IS_adapt.nn[i*NEAR+n],math.sqrt(IS_adapt.distances[i*NEAR+n]))

        spot = spots[Match["spot"]]
        S_minus_P_vector = matrix.col(
          (spot.ctr_mass_x() - predX,
           spot.ctr_mass_y() - predY))
        for pxl in spot.bodypixels:
          deltaX = pxl.x - spot.ctr_mass_x()
          deltaY = pxl.y - spot.ctr_mass_y()
          I_S_mask[ ( round(predX + deltaX + correction[0]),
                      round(predY + deltaY + correction[1]) )]=True
      self.ISmasks.append(I_S_mask)
      key_pairs = flex.int()
      for key in I_S_mask.keys():
        key_pairs.append(int(key[0])); key_pairs.append(int(key[1]))
      self.append_ISmask(key_pairs)
      corrections.append(correction)
    """

    # which spots are close enough to interfere with background?
    MAXOVER=6
    OS_adapt = AnnAdaptor(data=query,dim=2,k=MAXOVER) #six near nbrs
    OS_adapt.query(query)
    nbr_cutoff = 2.0* max(self.incr_focus[self.image_number])
    FRAME = int(nbr_cutoff/2)
    #print "The overlap cutoff is %d pixels"%nbr_cutoff
    nbr_cutoff_sq = nbr_cutoff * nbr_cutoff

    #print "Optimized C++ section...",
    self.set_frame(FRAME)
    self.set_nbr_cutoff_sq(nbr_cutoff_sq)
    flex_sorted = flex.int()
    for item in self.sorted:
      flex_sorted.append(item[0]);flex_sorted.append(item[1]);

    if self.inputai.active_areas != None:
      self.detector_xy_draft = self.safe_background( predicted=self.predicted,
                          OS_adapt=OS_adapt,
                          sorted=flex_sorted,
                          tiles=self.inputai.active_areas.IT,
                          tile_id=self.inputai.active_areas.tile_id);
    else:
      self.detector_xy_draft = self.safe_background( predicted=self.predicted,
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

    # Never get here...replaced with C++ code
    for i in xrange(len(predicted)): # loop over predicteds
      pred = predicted[i]
      predX = pred[0]/pxlsz
      predY = pred[1]/pxlsz
      correction = corrections[i]
      I_S_mask = self.ISmasks[i]
      # now consider the background
      B_S_mask = {}
      i_bs = 0
      spot_position = matrix.col(( round(predX + correction[0]),
                                  round(predY + correction[1]) ))
      self.detector_xy_draft.append(( round(predX + correction[0]),
                                  round(predY + correction[1]) ))

      #insert a test to make sure spot is within FRAME
      if spot_position[0] > FRAME and spot_position[1] > FRAME and \
         spot_position[0] < int(self.inputpd["size1"]) - FRAME and \
         spot_position[1] < int(self.inputpd["size2"]) - FRAME:

         spot_keys = I_S_mask.keys()
         spot_size = len(spot_keys)

         #Look for potential overlaps
         for n in xrange(MAXOVER):
           distance = OS_adapt.distances[i*MAXOVER+n]
           if distance < nbr_cutoff_sq:
             spot_keys += self.ISmasks[ OS_adapt.nn[i*MAXOVER+n] ].keys()

         for increment in self.sorted:
           candidate_bkgd = spot_position + increment
           b_s_key = (candidate_bkgd[0],candidate_bkgd[1])
           if b_s_key not in spot_keys:
             #eliminate if in guard region
             guard = False
             for key in spot_keys:
               if (b_s_key[0]-key[0])**2 + (b_s_key[1]-key[1])**2  < 10:
                 guard = True
                 break
             if guard: continue
             i_bs += 1
             B_S_mask[b_s_key] = True
           if i_bs == spot_size: break
      self.BSmasks.append(B_S_mask)
      # private interface.  If B_S_mask is the empty dictionary, the spot
      # was out of boundary and it is not possible to integrate

  def integration_proper(self):
    rawdata = self.imagefiles.images[self.image_number].linearintdata # assume image #1
    self.integration_proper_fast(rawdata,self.predicted,self.hkllist,self.detector_xy_draft)
    self.integrated_data = self.get_integrated_data()
    self.integrated_sigma= self.get_integrated_sigma()
    self.integrated_miller=self.get_integrated_miller()
    self.detector_xy = self.get_detector_xy()
    return # function has been recoded in C++

    self.integrated_data = flex.double()
    self.integrated_sigma= flex.double()
    self.integrated_miller=flex.miller_index()
    self.detector_xy = flex.vec2_double()
    from rstbx.diffraction import corrected_backplane
    for i in xrange(len(self.predicted)):
      signal = flex.double()
      bkgrnd = flex.double()
      if self.BSmasks[i].keys()==[]:continue # out-of-boundary spots
      #print "Integrating ",self.hkllist[i],
      keys = self.get_ISmask(i)
      smask = []
      for ks in xrange(0,len(keys),2):
        smask.append((keys[ks],keys[ks+1]))
      bmask = self.BSmasks[i]
      for spixel in smask:
        ispixel = int(spixel[0]),int(spixel[1])
        signal.append(rawdata[ispixel])
      BP = corrected_backplane(0,0)
      for bpixel in bmask:
        ibpixel = int(bpixel[0]),int(bpixel[1])
        bkgrnd.append(rawdata[ibpixel])
        BP.accumulate(int(bpixel[0]),int(bpixel[1]),rawdata[ibpixel])
      BP.finish()
      #print "<background>=%7.2f"%flex.mean(bkgrnd),
      #print "<signal>=%7.2f"%flex.mean(signal),

      corr_signal = flex.double()
      corr_bkgrnd = flex.double()
      for spixel in smask:
        ispixel = int(spixel[0]),int(spixel[1])
        corr_signal.append(rawdata[ispixel] -
                           BP.localmean(ispixel[0],ispixel[1]))
      for bpixel in bmask:
        ibpixel = int(bpixel[0]),int(bpixel[1])
        corr_bkgrnd.append(rawdata[ibpixel] - BP.localmean(ibpixel[0],ibpixel[1]))
      #print "<corr. background>=%7.2f"%flex.mean(corr_bkgrnd),
      #print "<corr. signal>=%7.2f"%flex.mean(corr_signal)

      summation_intensity = flex.sum(corr_signal)
      uncorrected_signal = flex.sum(signal)
      mcount = len(signal)
      ncount = len(bkgrnd)
      sum_background = flex.sum(bkgrnd)
      # Use gain == 1
      # variance formula from Andrew Leslie, Int Tables article
      variance = uncorrected_signal + sum_background*(float(mcount)/ncount)**2
      sigma = math.sqrt(variance)
      self.integrated_data.append(summation_intensity)
      self.integrated_sigma.append(sigma)
      self.integrated_miller.append(self.hkllist[i])
      self.detector_xy.append(self.detector_xy_draft[i])

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

  def integration_masks_as_xy_tuples(self):
    values = []
    for imsk in xrange(len(self.BSmasks)):
      smask_keys = self.get_ISmask(imsk)
      for ks in xrange(0,len(smask_keys),2):
        values.append((smask_keys[ks],smask_keys[ks+1]))
    return values

  def background_masks_as_xy_tuples(self):
    values = []
    for imsk in xrange(len(self.BSmasks)):
      bmask = self.BSmasks[imsk]
      for key in bmask.keys():
        values.append((key[0],key[1]))
    return values

  def user_callback(self,dc,wxpanel,wx):
    # arguments are a wx Device Context, an Xray Frame, and the wx Module itself
    # BLUE: predictions
    for ix,pred in enumerate(self.predicted):
        if self.BSmasks[ix].keys()==[]:continue
        x,y = wxpanel._img.image_coords_as_screen_coords(
          pred[1]/self.pixel_size,
          pred[0]/self.pixel_size)
        dc.SetPen(wx.Pen('blue'))
        dc.SetBrush(wx.BLUE_BRUSH)
        dc.DrawCircle(x,y,1)

    for imsk in xrange(len(self.BSmasks)):
      smask_keys = self.get_ISmask(imsk)
      bmask = self.BSmasks[imsk]
      if len(bmask.keys())==0: continue

      # CYAN: integration mask
      for ks in xrange(0,len(smask_keys),2):
        x,y = wxpanel._img.image_coords_as_screen_coords(smask_keys[ks+1],
                                                         smask_keys[ks])
        dc.SetPen(wx.Pen('cyan'))
        dc.SetBrush(wx.CYAN_BRUSH)
        dc.DrawCircle(x,y,1)

      # YELLOW: background mask
      for key in bmask.keys():
        x,y = wxpanel._img.image_coords_as_screen_coords(key[1],key[0])
        dc.SetPen(wx.Pen('yellow'))
        dc.SetBrush(wx.CYAN_BRUSH)
        dc.DrawCircle(x,y,1)

    for spot in self.spotfinder.images[self.frames[self.image_number]]["inlier_spots"]:
      # RED: spotfinder spot pixels
      for pxl in spot.bodypixels:
        x,y = wxpanel._img.image_coords_as_screen_coords(
          pxl.y,
          pxl.x)
        dc.SetPen(wx.Pen('red'))
        dc.SetBrush(wx.RED_BRUSH)
        dc.DrawCircle(x,y,1)

      # GREEN: spotfinder centers of mass
      x,y = wxpanel._img.image_coords_as_screen_coords(
        spot.ctr_mass_y(),
        spot.ctr_mass_x())
      dc.SetPen(wx.Pen('green'))
      dc.SetBrush(wx.GREEN_BRUSH)
      dc.DrawCircle(x,y,1)

  def user_callback1(self,dc,wxpanel,wx):
    x,y = wxpanel._img.image_coords_as_screen_coords(100,100)
    dc.SetPen(wx.Pen('green'))
    dc.SetBrush(wx.GREEN_BRUSH)
    dc.DrawCircle(x,y,10)

  def initialize_increments(self,image_number=0):
    #initialize a data structure that contains possible vectors
    # background_pixel - spot_center
    # consider a large box 4x as large as the presumptive mask.
    from scitbx.array_family import flex
    Incr = []
    Distsq = flex.double()
    for i in xrange(-self.incr_focus[image_number][0],1+self.incr_focus[image_number][0]):
      for j in xrange(-self.incr_focus[image_number][1],1+self.incr_focus[image_number][1]):
        Incr.append(matrix.col((i,j)))
        Distsq.append(i*i+j*j)
    order = flex.sort_permutation(Distsq)
    self.sorted = [] # a generic list of points close in distance to a central point
    for i in xrange(len(order)):
      #print i,order[i],Distsq[order[i]],Incr[order[i]]
      self.sorted.append(Incr[order[i]])
