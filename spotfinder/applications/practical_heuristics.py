from __future__ import division
from spotfinder.array_family import flex
import types,math
from spotfinder.exception import SpotfinderError
from spotfinder.core_toolbox import Distl,SpotFilterAgent,SingleMask
from spotfinder.math_support import stats_profile
from spotfinder.math_support import scitbx_stats
from spotfinder.applications.heuristic_tbx.spotreporter import spotreporter
from libtbx.development.timers import Timer, Profiler
from libtbx.utils import Sorry

TALLY2=0
OVERLAY = 1
DEBUG=0
VERBOSE_COUNT=0

def spot_filter(spotlist,function):
  return [x for x in spotlist if function(x)]

def resolution_comparisons(a,b):
  #gives -1,0,1 depending on a>b, a==b, a<b
  if a.resolution<b.resolution: return 1
  if a.resolution==b.resolution: return 0
  return -1

def resol_to_radius(resolution,pd):
  distance = float(pd['distance'])
  wavelength = float(pd['wavelength'])
  theta = math.asin(wavelength/2.0/resolution)
  return distance * math.tan(2.0*theta)

from spotfinder.diffraction.geometry import radius_to_resol

class pSpot:
  def __init__(self,x,y,intensity):
    self.mx = x
    self.my = y
    self.mi = intensity
  def x(self):
    return self.mx
  def y(self):
    return self.my
  def intensity(self):
    return self.mi

class SaturationMeasure(object):
  def forward_message(self):
    return "calculate saturation on %d of %d"%(self.n_resolution_spots,
            self.n_goodspots)
  def saturation(self):
    return self.p_saturation
  def message(self):
    return "%%Saturation, Top %d Peaks"%self.n_sample
  def format(self):
    return "%.2f"%(100*self.p_saturation)
  def message2(self):
    return "In-Resolution Ovrld Spots"
  def format2(self):
    return "%d"%(self.OverCount)

class ListNode:
  def __init__(self,nodelist,nodekey,parent=0):
    self.data = nodelist
    self.descriptor = nodekey
    self.parent = parent

class ListManager(dict):
  """From an initial spotlist, pick qualified subsets"""
  def __init__(self,masterlist,masterkey,spotfilter):
    self.master = masterlist
    self.nodes = {1:ListNode(masterlist,masterkey)}
    self.key = {masterkey:1}
    self.spotfilter = spotfilter

  def alias(self,oldkey,newkey):
    self.key[newkey]=self.key[oldkey]

  def add_child(self,oldkey,newkey,childselection):
    parent = self.key[oldkey];    nextkey = max(self.nodes.keys())+1
    self.nodes[nextkey]=ListNode(childselection,newkey,parent)
    self.key[newkey]=nextkey

  def spot_filter(self,oldkey,newkey,apfunction):
    #parent
    parent = self.key[oldkey]
    if parent == 1:
      parentselection = xrange(self.master.size())
    else:
      parentselection = self.nodes[parent].data
    nextkey = max(self.nodes.keys())+1
    childselection = flex.int()
    for x in parentselection:
      if apfunction(self.master[x]):
        childselection.append(x)
    self.nodes[nextkey]=ListNode(childselection,newkey,parent)
    self.key[newkey]=nextkey

  def c_spot_filter(self,oldkey,newkey,apfunction,arguments=[]):
    parent = self.key[oldkey]
    if parent == 1:
      parentselection = xrange(self.master.size())
    else:
      parentselection = self.nodes[parent].data
    nextkey = max(self.nodes.keys())+1
    self.spotfilter.set_arguments(flex.double(arguments))
    childselection = self.spotfilter.filter(
      self.master,parentselection,apfunction)
    self.nodes[nextkey]=ListNode(childselection,newkey,parent)
    self.key[newkey]=nextkey

  def precompute_resolution(self,twotheta,rotation_axis,camera_convention):
    self.spotfilter.precompute_resolution(self.master,twotheta,
      flex.double(rotation_axis),camera_convention)

  def single_mask(self,oldkey,phil_params):
    parent = self.key[oldkey]
    if parent == 1:
      parentselection = xrange(self.master.size())
    else:
      parentselection = self.nodes[parent].data
    Mask = SingleMask(master=self.master,selection=parentselection,
           minimum_spot_count =
           phil_params.codecamp.minimum_spot_count or 25)
    return Mask

  def resolution_sort(self,oldkey):
    parentselection = self.get_parentselection(oldkey)
    sort_order = self.spotfilter.resolution_sort(self.master,parentselection)
    sorted_resolutions=self.spotfilter.get_resolution(self.master,sort_order)
    return sort_order,sorted_resolutions

  def resolution_sort_nztt(self,oldkey):
    parentselection = self.get_parentselection(oldkey)
    sort_order = self.spotfilter.resolution_sort_nztt(self.master,parentselection)
    sorted_resolutions=self.spotfilter.get_resolution(self.master,sort_order)
    return sort_order,sorted_resolutions

  def order_by_criterion(self,oldkey,criterion):
    parentselection = self.get_parentselection(oldkey)
    sort_order = self.spotfilter.order_by(self.master,parentselection,criterion)
    spotlist = flex.distl_spot()
    for number in sort_order:
        spotlist.append(self.master[number])
    return spotlist

  def most_recent_child(self):
    return self.nodes[max(self.nodes.keys())].descriptor

  def get_indices(self,key):
    if self.key[key]==1: return xrange(len(self.master))
    return self.nodes[self.key[key]].data

  def get_indices_without(self,key,subtractive_key):
    initial_set = list(self.get_indices(key))
    subtractive_set = list(self.get_indices(subtractive_key))
    initial_set.sort()
    subtractive_set.sort()
    #use an algorithm to subtract list B from list A
    C = []
    ptrA = 0
    ptrB = 0
    lenA = len(initial_set)
    lenB = len(subtractive_set)
    while ptrA < lenA and ptrB <= lenB:
      if ptrB == lenB:
        C.append(initial_set[ptrA]); ptrA+=1; continue
      if initial_set[ptrA]==subtractive_set[ptrB]:
        ptrA+=1; ptrB+=1; continue
      if initial_set[ptrA]<subtractive_set[ptrB]:
        C.append(initial_set[ptrA]); ptrA+=1; continue
      ptrB+=1;
    return flex.int(C) # returned indices will no longer be sorted by resolution

  def get_property(self,key,property):
    parentselection = self.get_parentselection(key)
    return self.spotfilter.get_property(self.master,parentselection,property)

  def get_parentselection(self,key):
    parent = self.key[key]
    if parent == 1:
      parentselection = xrange(self.master.size())
    else:
      parentselection = self.nodes[parent].data
    return parentselection

  def __getitem__(self,key):
    if key.find("N_")==0:
      newkey=key[2:]
      return self.nodes[self.key[newkey]].data.size()
    if key=='spotoutput':
      return self
    if key=='resolution_spots':#for testing only
      spotlist = []
      for number in self.nodes[self.key[key]].data:
        spotlist.append(self.master[number])
      return spotlist
    if key in self.key.keys():
      if self.key[key]==1: return self.master
      spotlist = flex.distl_spot()
      for number in self.nodes[self.key[key]].data:
        spotlist.append(self.master[number])
      return spotlist
    return dict.__getitem__(self,key)

  def __delitem__(self,key):
    if key in self.key.keys(): pass
    else:
      dict.__delitem__(self,key)

  def keys(self):
    builtin = dict.keys(self)
    for key in self.key.keys():
      builtin.append(key); builtin.append("N_"+key)
    return builtin

  def has_key(self,key): return key in self.keys()

def pickle_safe_spotcenter(spot,algorithm):
  #Image overlays and autoindexing maintain backward compatibility
  #  with spot information pickled under the old maximum_pixel method.
  if (algorithm == 'maximum_pixel' or not spot.com_valid()):
        return 'maximum_pixel'
  elif algorithm == 'center_of_mass':
        return 'center_of_mass'

class heuristics_base(object):

  def __init__(self,pd,phil_params):
    self.pd = pd
    self.phil_params = phil_params
    self.pd['resolution_inspection']='100.0'
    self.pd['ref_maxcel']='10.0' # default maximum cell dimension (always override this)
    self.NspotMin = self.phil_params.distl_minimum_number_spots_for_indexing
    self.NspotMax = self.phil_params.distl_maximum_number_spots_for_indexing
    self.BinMin = 25
    self.errormessage = None
    self.images = {}
    self.protocol = 'tnear2'
    self.overlapping = False #flag indicates whether the special procedure was used
    self.force_detail = False #flag indicates whether percent_overlap > force_detail cutoff

  def register_frames(self,frameinfo,imagefilesinstance):
    if type(frameinfo) in [types.IntType,types.LongType]:
      frames = [int(frameinfo),]
    elif type(frameinfo) in [types.TupleType,types.ListType]:
      frames = frameinfo
    for x in xrange(len(frames)):
      self.images[frames[x]] = self.oneImage(frames[x],self.pd,
                                       imagefilesinstance.imageindex(frames[x]))
      self.determine_maxcell(frames[x],self.pd)
      self.images[frames[x]]['spotoutput']['relpath']=imagefilesinstance.imagepath(frames[x])

  def oneImage(self,framenumber,pd,image):
    # The only way to get pixel size & pixel dimensions (currently) is from header
    pimage = image
    pimage.read()
    #print "Detector type",type(pimage)
    if not pd.has_key('endstation'):
      from iotbx.detectors.context import endstation
      pd['endstation']=endstation.EndStation_from_ImageObject(pimage,self.phil_params)
      for key in pd['endstation'].mosflm():
        pd[key]=pd['endstation'].mosflm()[key]
    pd['vendortype'] = pimage.vendortype
    pd['binning']      = "%d"%pimage.bin
    pd['pixel_size'] = "%f"%pimage.pixel_size
    self.pixel_size = float(pd['pixel_size'])

    pd['size1']      = "%d"%pimage.size1
    self.size1       = float(pd['size1'])
    pd['size2']      = "%d"%pimage.size2
    self.size2       = float(pd['size2'])
    if not pd.has_key('osc_start'): pd['osc_start'] = {}
    pd['osc_start'][framenumber] = "%f"%pimage.osc_start
    if not pd.has_key('file'): pd['file'] = {}
    pd['file'][framenumber] = pimage.filename
    self.two_theta_degrees = float(pd['twotheta'])

    if not pd.has_key('xbeam') or not pd.has_key('ybeam'):
      raise SpotfinderError("Deprecation warning: inputs had no beam position",pd)
    self.complex_nominal_center = complex(float(pd["xbeam"]),float(pd["ybeam"]))

    arguments = "" # distl_aggressive no longer supported, eg "-s2 5 -s3 8"

    #allow for self.phil_params.distl_highres_limit
    if self.phil_params.distl_highres_limit != None:
      arguments = arguments + " -ro %.3f"%self.phil_params.distl_highres_limit

    #  test implementation of parallel processing for spotfinder
    #from labelit.webice_support import parallel_distl
    #pimage,pd = parallel_distl.split_image(pimage,pd)

    try: sf = Distl(arguments,pimage,pd,
            report_overloads=self.phil_params.distl_report_overloads,
            params=self.phil_params)
    except Sorry, e:
      raise e
    except Exception, e:
      raise SpotfinderError("Spotfinder cannot analyze image %s"%frame,pd)

    #To support sublattice detection, make pixel-wise Z-scores persistent
    if self.phil_params.distl_keep_Zdata:
      pimage.linear_Z_data = sf.Z_data() #potentially uses a lot of memory

    #************************************************************
    #
    #  Very important.  For the mar image plate (and allother circular detectors,
    #  must specify that these are embedded circular detectors and so adjust the
    #  percentage underloads.  Or else libdistl must be changed so as not to
    # search in these regions:  yes, this would be better because I propose to
    # change the search loop anyway.  However, the getUnderload function must
    # also be modified!!!
    #
    #  ice ring search must also be modified to take this into account, but this
    #  part must be more precise than the above.
    #
    #************************************************************

    if self.two_theta_degrees==0.0:
      mm_minimum_radius = resol_to_radius(
      self.phil_params.distl_lowres_limit,pd)

    sfa = SpotFilterAgent(pixel_size = pimage.pixel_size,
                          xbeam = float(pd["xbeam"]),
                          ybeam = float(pd["ybeam"]),
                          distance = float(pd['distance']),
                          wavelength = float(pd['wavelength']),
                          icerings = sf.icerings,)
    #from libtbx import easy_pickle
    #easy_pickle.dump("file.dmp",sfa)
    #easy_pickle.load("file.dmp")

    fstats = ListManager(
      masterlist=sf.spots,masterkey='spots_total',spotfilter=sfa)

    fstats['distl_resolution'] = sf.imgresol()

    # 1. Get all spots

    fstats.alias(oldkey = 'spots_total',newkey = 'goodspots')
    if VERBOSE_COUNT: print "total DISTL spots",fstats['N_spots_total']

    fstats.c_spot_filter('goodspots','spots_non-ice','ice_ring_test')
    fstats['ice-ring_impact'] = sf.nicerings()
    fstats['ice-ring_bounds'] = [(sf.icerings[i].lowerresol,sf.icerings[i].upperresol)
      for i in xrange(fstats['ice-ring_impact']) ]

    #**********************************************************************
    #  Known parts of code that are inefficient: use 35000-spot HK97 example
    #  1. 3.8 seconds: sorting the resolution spots (item #4 in tnear2) (Corrected)
    #  2. 9.7 seconds: spreadsheet bookkeeping; method2_resolution.py, xrow loop (Partly corrected 5/09)
    #  3. 4.4 seconds: ice2::RingFinder::filtered()  (Corrected)

    # 3. omit spots too close to the beamstop
    if self.two_theta_degrees==0.0:
      fstats.c_spot_filter('spots_non-ice','hi_pass_resolution_spots',
                                         'resolution_test',
                                         arguments=[mm_minimum_radius,])
    else:
      fstats.precompute_resolution(
        self.two_theta_degrees * math.pi / 180., #two theta radians
        pd['endstation'].rotation_axis(),
        pd['endstation'].camera_convention())

      fstats.c_spot_filter('spots_non-ice','hi_pass_resolution_spots',
                                           'resolution_test_nztt',
                          arguments=[self.phil_params.distl_lowres_limit])


    if VERBOSE_COUNT:
      print "after lowres filter",fstats["N_hi_pass_resolution_spots"]
#start here.
#In the end, make sure these work:
#interface with mosflm: "TWOTHETA" keyword fails; "TILT" fix works with fudge factor
#James Holton's spots index (make unit test); anything with resolution_mm
#diffimage display!:  ring_markup() method of webice_support/__init__
#report the two theta value in stats index
    # 4. Calculate resolution cutoff
    if self.two_theta_degrees==0.0:
      sorted_order,sorted_resolutions = fstats.resolution_sort(
      "hi_pass_resolution_spots")
    else:
      sorted_order,sorted_resolutions = fstats.resolution_sort_nztt(
      "hi_pass_resolution_spots")
    # targetBinNumber: first try number of candidate Bragg spots per bin
    targetBinNumber = max(self.BinMin, len(sorted_order)//20)

    # cutoff threshhold: expected number spots in highest resolution bin
    cutoff_ratio = self.phil_params.tnear_resolution_divisor
    fstats['resolution_divisor']=cutoff_ratio
    lowerCutoffBinNumber = targetBinNumber / cutoff_ratio

    if len(sorted_order) < targetBinNumber:
      # So few spots that there is only one bin
      # no resolution determination possible
      fstats['resolution'] = None
      fstats['resolution_detail'] = len(sorted_order)

    elif False and sorted_resolutions[self.BinMin-1] < 4.0:
      '''This filter (that essentially says you must have low
      resolution spots) hindered the ability to index the case
      ana/procrun0000084148/sphN1_*.mar2300.  Further investigation
      showed that not a single one of the 94 reference cases in the
      regression database was affected by this test, so the test is
      being provisionally removed.  It is not known if the test has a
      beneficial effect on cases with no protein diffraction, i.e., just
      ice rings.'''
      # This test indicates that there are so few spots at low
      #  resolution that the first bin extends past 4.0 Angstroms,
      #  into the region where ice rings might be found
      #  since there are too few low resolution data
      #  conclude that these are false Bragg spots
      fstats['resolution'] = None
      sorted_resolutions=[]
      # deprecated; would need to be recoded:
      fstats['N_spots_resolution'] = 0

    else:
      from spotfinder.diffraction.geometry import Geom2d
      frac_calc = Geom2d(pd)
      #spot-based ice-ring filtering, added Aug 2004
      '''The case ana/procrun0000084148/sphN1_*.mar2300 (image 090) shows the
      limits of this filter as presently implemented.  In that particular
      case, this ice-ring filter eliminates spots from a large area in the
      2-3 Angstrom resolution range.  However, to be believed, the purported
      ice-spots should define a circle or ellipse centered at the beam position
      with a resolution spread that is very narrow.  Code could be
      written much more effectively to search and elimate these rings'''
      from spotfinder.applications.heuristic_tbx.ice2 import RingFinder
      from spotfinder.applications.heuristic_tbx.ice_nztt import RingFinder_nztt

      if (abs(float(pd['twotheta'])) > 0.0):
        # Very inefficient--8 seconds per call; will need to be optimized
        #PP = Profiler("nztt")
        Ring = RingFinder_nztt(sorted_resolutions,targetBinNumber,frac_calc,
                               image)
        #del PP
      else:
        Ring = RingFinder(sorted_resolutions,targetBinNumber,frac_calc)

      fstats.add_child("hi_pass_resolution_spots",
                       "ice_free_resolution_spots",
                       Ring.filtered(sorted_order))
      fstats['ice-ring_impact'] += Ring.ice_ring_impact()
      fstats['ice-ring_bounds'] += Ring.ice_ring_bounds()

      #*************************the filtering of existing spots

      if VERBOSE_COUNT:
       print "after spot_based ice-ring filter",fstats[
             'N_ice_free_resolution_spots']

      if fstats['N_ice_free_resolution_spots'] < targetBinNumber:
        # So few spots that there is only one bin
        # no resolution determination possible
        fstats['resolution'] = None
        fstats['resolution_detail'] = fstats['N_ice_free_resolution_spots']
      else:
        from spotfinder.applications.heuristic_tbx.method2_resolution\
          import ResolutionShells
        sorted_resolutions = fstats.spotfilter.get_resolution(
                             fstats.master,fstats.get_indices('ice_free_resolution_spots'))

        Shell = ResolutionShells(sorted_resolutions,
          sorted_resolutions[targetBinNumber-1],frac_calc)
        #Shell.show()

        if self.phil_params.distl.bins.verbose:
          ShellR = spotreporter(sorted_resolutions,self.phil_params,
            sorted_resolutions[targetBinNumber-1],
            wavelength = float(pd['wavelength']),
            max_total_rows=self.phil_params.distl.bins.N,
            fractionCalculator=frac_calc,
            use_binning_of=self)
          ShellR.total_signal(fstats,fstats.get_indices('ice_free_resolution_spots'))
          ShellR.background(sf)
          ShellR.sigma_analysis()
          print
          ShellR.show(message="Analysis of spots after ice removal, but prior to resolution cutoff, for image \n%s"%pimage.filename)

        # slight adjustment so that we don't focus on the lowest
        #  resolution data shell (spends too much time in Ewald sphere)
        #  But don't make this adjustment if we are limited by low spot count
        if Shell.rows() > 2 and \
            Shell.Population[1] > 2*lowerCutoffBinNumber and \
            Shell.Population[1] > self.BinMin:
          lowerCutoffBinNumber = Shell.Population[1] / cutoff_ratio

        # first determination of cutoff ignoring corner effect
        for x in xrange(Shell.rows()):
          idx = Shell.rows()-x-1
          if Shell.Population[idx] > lowerCutoffBinNumber:
            lastshell = idx
            break

        # eliminate pathological case where there are a lot of low resolution
        # spots and then some ice rings at high resolution.  If there are ten
        # empty shells in a row (empty defined as having fewer than
        # VetterCut spots), redetermine lastshell.
        #
        # Originally, VetterCut := lowerCutoffBinNumber
        # Modify the heuristic to cover two extremes:
        # Case 1.  Based on the HK97 virus work; there is a class of diffraction
        #    cases showing a distinct dip in diffraction intensity in the
        #    5-Angstrom regime.  If lowerCutoffBinNumber is used (usually
        #    20% of the spot count in the 2nd bin), the algorithm can
        #    misbehave and reject all the high-resolution Bragg spots.
        #    For one case in particular a tiny change in beam position
        #    flips the resolution cutoff from a 6.5- to 4.4-Angstrom.
        # Case 2.  There are many cases where the diffraction clearly
        #    drops off, and at higher resolutions there are random-signal
        #    spots plus ice rings.  The spot count in these high-res bins
        #    is typically small, <20.  Use "20" as a heuristic cutoff
        #    between Case 1 & Case 2.
        # In future, it may be productive to focus on smarter algorithms to
        # execute ellipse recognition in the "Ringfinder" section above

        if lowerCutoffBinNumber <= 20:
          VetterCut = lowerCutoffBinNumber
        else:
          VetterCut = 20 + lowerCutoffBinNumber // 5
        lastshell = Shell.vetter(VetterCut,lastshell)

        #option for overriding the resolution analysis based on falloff
        # of spot count.  Force spots at least this far out
        if self.phil_params.force_method2_resolution_limit is not None:
          for x in xrange(Shell.rows()):
            if Shell.Limit[x]<self.phil_params.force_method2_resolution_limit and lastshell<x:
              lastshell = x
              break

        # extend resolution cutoff taking into account corner cutoff
        while Shell.Fract[lastshell] < 1.0 and lastshell+1 < Shell.rows():
          nextshell = lastshell+1
          if Shell.adjustPop[nextshell] > lowerCutoffBinNumber:
            lastshell+=1
          else: break

        fstats['resolution'] = Shell.Limit[lastshell]
        if self.phil_params.distl_highres_limit!=None:
           fstats['resolution']=max(fstats['resolution'],
             self.phil_params.distl_highres_limit)

        for x in xrange(lastshell+1):
          if self.phil_params.spotfinder_verbose: print "(%.2f,%d)"%(Shell.Limit[x],Shell.Population[x])

        if self.two_theta_degrees==0.0:
          fstats.c_spot_filter(          #hi-resolution radius
          'ice_free_resolution_spots',
          'lo_pass_resolution_spots',
          'lo_pass_resolution_test',
          arguments=[resol_to_radius(fstats['resolution'],self.pd),])
        else:
          fstats.c_spot_filter(          #hi-resolution radius
          'ice_free_resolution_spots',
          'lo_pass_resolution_spots',
          'lo_pass_resolution_test_nztt',
          arguments=[fstats['resolution'],])

        if VERBOSE_COUNT:
          print "ice_free_resolution_spots ",fstats['N_ice_free_resolution_spots']

        fstats['shells']=Shell

        if VERBOSE_COUNT:
          print "after resolution-shell cutoff",fstats['N_lo_pass_resolution_spots']
        fstats['resolution_mm']=resol_to_radius(fstats['resolution'],pd)

    fstats['saturation'] = self.calculate_saturation(fstats,image)

    if self.phil_params.distl.bins.verbose:
     try:
      subset = fstats.spotfilter.get_resolution(
               fstats.master,fstats.get_indices('lo_pass_resolution_spots'))
      ShellR = spotreporter(subset,self.phil_params,use_binning_of=self)
      ShellR.total_signal(fstats,fstats.get_indices('lo_pass_resolution_spots'))
      ShellR.background(sf)
      ShellR.sigma_analysis()
      print
      ShellR.show(message="Analysis of spots after resolution filtering, but prior to spot quality heuristics, for image \n%s"%pimage.filename)
     except Exception: # in case the low-pass filter step was skipped; e.g., blank image.
      pass

    # 5. eliminate multi-modal spots & punctate spots
    fstats.c_spot_filter(
      fstats.most_recent_child(),
      'spots_unimodal',
      'modal_test',
      arguments=[self.phil_params.distl_profile_bumpiness,
                 self.phil_params.distl.minimum_spot_area])
      # parameters are bumpiness(max number of local maxima in peak),minimum pixels

    if VERBOSE_COUNT: print "not bumpy & not punctate",fstats['N_spots_unimodal']

    # 6. Compute distributions for outlier rejection
    #Y = Timer("Inliers")#try to get this down from 232 seconds to 8 seconds
                         #(For HK97 sz=4)

    inlier_idx_raw = flex.int()
    inlier_neigh=[]
    unimodal_spots = fstats['spots_unimodal']

    if fstats['N_spots_unimodal']>2: # avoids divide-by-zero error in stats calls, below
      intensities = flex.double([s.intensity() for s in unimodal_spots])
      areas       = flex.double([s.bodypixels.size() for s in unimodal_spots])
      # Deprecate shapes because code is unstable; in rare cases spots can have
      # body pixels but no border pixels ( RAW/APS_07_2005/pEI4/pEI4_8_1.0001 )
      #shapes      = flex.double([s.shape() for s in unimodal_spots])
      eccentricity= flex.double([s.model_eccentricity() for s in unimodal_spots])
      skewness    = fstats.get_property(fstats.most_recent_child(),'skewness')

      fstats['intensity'] = (i_ave,i_std) = scitbx_stats(intensities)
      fstats['area']      = (a_ave,a_std) = scitbx_stats(areas)
      #fstats['shape']     = (s_ave,s_std) = scitbx_stats(shapes)
      fstats['eccen']     = (e_ave,e_std) = scitbx_stats(eccentricity)
      fstats['skewness']  = (k_ave,k_std) = scitbx_stats(skewness)

      #Filtering on skewness is important when using center-of-mass positions
      # for autoindexing.  The goal is to eliminate pairs of unresolved Bragg
      # spots.  Condition is flagged when the max_pixel to center-of-mass
      # vector is more than 45% of the semi-major axis
      fstats.c_spot_filter(
          fstats.most_recent_child(),
          'spots_low_skew',
          'low_skew',
          arguments=[min(0.45, 2.0 * k_std + k_ave),])

      #Filter on the intensity distribution.
      fstats.c_spot_filter(
          fstats.most_recent_child(),
          'spots_good_intensity',
          'intensity_inlier',
          arguments=[i_ave - 5.0*i_std, i_ave + 5.0*i_std, image.saturation,])

      #special code for expecting a high MOSFLM RESID based on the
      # presence of very high signal/noise ratio (essentially lysozyme strength).
      # This is all just a supposition based on one dataset, offset_1 from Ana.
      # The assumption may be wrong; e.g. the high resid may be due to very low
      # background instead of high s/n; or due to some geometrical distortion.
      #This breaks encapsulation and will have to be re-organized in the future.
      if not pd.has_key('special_resid'): pd['special_resid']=''
      if fstats['intensity'][0]>160.: pd['special_resid']='RESID 10.0 #High s/n'
      #The most unusual thing about the procrun0000077831/TMC114_WT2_run77831_1_001
      #dataset is its large differential between spot areas of the largest spots
      # and the smallest spots.  Hypothesize that this also leads to high
      # weighted residuals.
      if stats_profile(areas)>10.: pd['special_resid']='RESID 10.0 #High s/n'

      from annlib_ext import AnnAdaptor
      data = fstats.get_property('goodspots',
             self.phil_params.distl_spotcenter_algorithm)
      query = fstats.get_property(fstats.most_recent_child(),
              self.phil_params.distl_spotcenter_algorithm)

      A = AnnAdaptor(data,2)       # construct k-d tree for reference set
      A.query(query)               # find nearest neighbors of query points

      neighbors = (self.pixel_size) * flex.sqrt(A.distances)

      """Explanation: Distance to the nearest neighbor is math.sqrt(A.distances[i]), in
         units of pixels.  The vector to the nearest neighbor, in units of pixels, is
         ( fstats.master[A.nn[i]].x()-query[2*i], fstats.master[A.nn[i]].y()-query[2*i+1])
      """

      fstats['neighbor']  = (n_ave,n_std) = scitbx_stats(neighbors)

      sep_input_spots = fstats[fstats.most_recent_child()]
      sep_input_indices = fstats.get_indices(fstats.most_recent_child())

      overlapping_count = 0
      for idx in xrange(len(sep_input_spots)):
        try:
          if float(pd['pixel_size'])*sep_input_spots[idx].majoraxis() * \
             self.phil_params.overlapping_spot_criterion > neighbors[idx]:
            overlapping_count+=1
        except Exception:
          pass
      percent_overlap = 100*overlapping_count/len(sep_input_spots)
      #print "overlap %2.0f%% vs. cutoff %2.0f%%"%(percent_overlap,self.phil_params.percent_overlap_forcing_detail)

      from spotfinder.core_toolbox.close_spots_detail import NearNeighborVectors
      Afull = AnnAdaptor(data,2)
      Afull.query(data)
      NV = NearNeighborVectors(ave_area = a_ave, query_centers = data,
                               fstats = fstats, ann_adaptor = Afull)
      #NV.show_vector_map()
      #NV.show_maxima()

      DEVELOP_ASSERT = True
      self.force_detail = percent_overlap > self.phil_params.percent_overlap_forcing_detail
      if DEVELOP_ASSERT or self.force_detail:
        self.overlapping = True
        #extraordinary procedure, when many spots are close.  Include closely
        #  spaced spots in autoindexing, if it appears that they truly
        #  reflect lattice spacing.
        if self.phil_params.spotfinder_verbose:
          print len(sep_input_spots),"spot count before close neighbor analysis;",
          print overlapping_count,"(%2.0f%%) rejected on neighbors;"%(percent_overlap)

        pmax = NV.vectors()
        #filter out spots that are potentially large enough to contain
        #two spots, now that we know a candidate projected unit-cell vector
        # expect big time savings if this section is pushed down to C++
        compact_idx = []

        if 0<len(pmax)<=3:
          sq_vectors = [v[0]*v[0] + v[1]*v[1] for v in pmax]
          for idx in xrange(len(sep_input_spots)):
            spot_compact = True
            for iv,vector in enumerate(pmax):
              thisspot = sep_input_spots[idx]
              #more efficient code--disallowed spots (based on nearest neighbor
              # vector) are now flagged based on the width of the ellipse model
              # rather than a more time-consuming body-pixel match:

              #calculate angle theta between major axis and neighbor vector
              dot = thisspot.eigenvector(1)[0]*vector[0] + thisspot.eigenvector(1)[1]*vector[1]
              costheta = dot / math.sqrt(sq_vectors[iv])
              sintheta = math.sqrt(1 - costheta*costheta)
              spota = thisspot.a()
              spotb = thisspot.b()
              a_sin_theta = spota * sintheta
              b_cos_theta = spotb * costheta
              #evaluate the ellipse full width along neighbor vector direction
              width_sq = 4.*(spota*spota*spotb*spotb/
                (a_sin_theta*a_sin_theta + b_cos_theta*b_cos_theta) )
              if width_sq >= sq_vectors[iv]:
                spot_compact = False
                break

            if spot_compact: compact_idx.append(idx)
            else: pass #print "eliminate the spot at",thisspot.x(),thisspot.y()
        else:
          compact_idx = xrange(len(sep_input_spots))
        if self.phil_params.spotfinder_verbose: print len(sep_input_spots)-len(compact_idx),"large spots rejected"

        # finally, allow certain close spots (but not all of them) to be included
        for idx in compact_idx:
          try:
            S =  ( int( fstats.master[A.nn[idx]].max_pxl_x()-query[2*idx] ),
                   int( fstats.master[A.nn[idx]].max_pxl_y()-query[2*idx+1]) )
            #accept spot by default
            reject_spot = False
            #if close to nearest neighbor reject
            if sep_input_spots[idx].majoraxis() * \
               self.phil_params.overlapping_spot_criterion > math.sqrt(S[0]*S[0]+S[1]*S[1]):
               # with normal procedure, spot would be rejected at this point
               reject_spot = True
               if len(pmax)<=3:
                 for vector in pmax:
                   if abs(vector[0]-S[0])<=1 and abs(vector[1]-S[1])<=1:
                     #but allow if the proximity is to a candidate cell near neighbor
                     reject_spot=False
                     break
            if reject_spot: continue

            inlier_idx_raw.append(sep_input_indices[idx])
            inlier_neigh.append(neighbors[idx])
          except Exception:pass

        if self.phil_params.spotfinder_verbose:
          print len(compact_idx)-len(inlier_idx_raw),"close spots rejected"
          print len(sep_input_spots),"input for spot separation analysis;",
          jj = len(sep_input_spots)-len(inlier_idx_raw)
          print jj,"(%2.0f%%) rejected on special criteria;"%(100.*jj/len(sep_input_spots))

      else: #normal procedure, assuming not too many close spots
        for idx in xrange(len(sep_input_spots)):
          try:
            if math.fabs(intensities[idx]-i_ave) <= 5.0*i_std and \
               intensities[idx]<image.saturation and \
               float(pd['pixel_size'])*sep_input_spots[idx].majoraxis() * \
               self.phil_params.overlapping_spot_criterion<=neighbors[idx]:
               inlier_idx_raw.append(sep_input_indices[idx])
               inlier_neigh.append(neighbors[idx])
          except Exception:
            #print "REJECT spot on exception"
            pass #sometimes throw an error when majoraxis is requested (edge spots)

        if len(inlier_idx_raw)<self.NspotMin:
          for idx in xrange(len(sep_input_spots)):
            try:
              proximal_radius = 2.0 * self.phil_params.overlapping_spot_criterion * float(pd['pixel_size'])*sep_input_spots[idx].majoraxis()
              if math.fabs(intensities[idx]-i_ave) <= 5.0*i_std and \
                 intensities[idx]<image.saturation and \
                 float(pd['pixel_size'])*sep_input_spots[idx].majoraxis() *self.phil_params.overlapping_spot_criterion > neighbors[idx] and \
                 sf.isIsolated(sep_input_spots[idx],proximal_radius):
                 #print "Very few Bragg spots; forced to accept spot with neighbor distance",neighbors[idx]/(float(pd['pixel_size'])*sep_input_spots[idx].majoraxis())
                 inlier_idx_raw.append(sep_input_indices[idx])
                 inlier_neigh.append(neighbors[idx])
            except Exception:
              print "REJECT spot on exception"
              pass #sometimes throw an error when majoraxis is requested (edge spots)

    if fstats.has_key('lo_pass_resolution_spots'):
      fstats.alias('lo_pass_resolution_spots','spots_resolution')
    elif fstats.has_key('ice_free_resolution_spots'):
      fstats.alias('ice_free_resolution_spots','spots_resolution')
    else:
      fstats.alias('hi_pass_resolution_spots','spots_resolution')

    fstats.add_child('spots_unimodal','spots_separated',
                     inlier_idx_raw)

    fstats.alias(fstats.most_recent_child(),'spots_inlier')
    fstats.alias('spots_inlier','inlier_spots')#for backward compatibility

    if self.phil_params.distl.bins.verbose:
     try:
      subset = fstats.spotfilter.get_resolution(
               fstats.master,fstats.get_indices('inlier_spots'))
      ShellR = spotreporter(subset,self.phil_params,use_binning_of=self)
      ShellR.total_signal(fstats,fstats.get_indices('inlier_spots'))
      ShellR.background(sf)
      ShellR.sigma_analysis()
      print
      ShellR.show(message="Analysis of good Bragg spots after quality heuristics, for image \n%s"%pimage.filename)
     except Exception:
      pass #if there aren't enough spots, just skip the tabular printout

    if not pd.has_key('masks'):  pd['masks']={}
    pd['masks'][framenumber] = None
    if fstats['N_spots_inlier'] > ( #guard against C++ hard-coded minimum
      self.phil_params.codecamp.minimum_spot_count or 25):
      Msk = fstats.single_mask('inlier_spots',self.phil_params)
      pd['masks'][framenumber] = [Msk.x,Msk.y]

    if VERBOSE_COUNT: print "inlier spots",fstats["N_inlier_spots"]

    #need this for later calculation of maxcell
    fstats['neighbors'] = flex.double(inlier_neigh)

    return fstats

  def calculate_saturation(self,fstats,image):
    S = SaturationMeasure()
    resolution_spots = fstats[fstats.most_recent_child()]
    S.n_resolution_spots = len(resolution_spots)
    S.n_goodspots = fstats['N_goodspots']
    ispots = flex.int()
    for i in xrange(S.n_resolution_spots):
      spt = resolution_spots[i]
      ispots.append( image.linearintdata[(spt.max_pxl_x(),spt.max_pxl_y())] )
    S.OverCount = (ispots >= int(image.saturation)).count(True)
    perm = flex.sort_permutation(ispots,True)
    spot_peaks = flex.int([ispots[p] for p in perm[0:min(50,len(perm))]])
    S.n_sample = len(spot_peaks)
    if S.n_sample < 2: #no peaks on the image
      S.p_saturation = 0.0; return S
    S.average = float(reduce(lambda x,y:x+y, spot_peaks))/len(spot_peaks)
    S.p_saturation = S.average/image.saturation
    # for Raxis-II and Mar Image Plates, it has been noticed that image.linearintdata
    # can contain pixel values greater than image.saturation.  No doubt this is
    # a shortcoming in how the pixel values are decoded for overloaded values (not fixed)
    return S

  def show(self):
    for frame in self.images.keys():
      for key in ['N_spots_total','N_spots_non-ice','N_spots_resolution','N_spots_unimodal',
                  'N_spots_inlier','intensity','area',
                  'neighbor','maxcel','resolution',
                  'distl_resolution','ice-ring_impact']:
        if self.images[frame].has_key(key):
          print "\t"+key,self.images[frame][key]
      for key in ['eccen']:
        if self.images[frame].has_key('eccen'):
          if self.images[frame].has_key(key):
            print "\t"+key,self.images[frame][key]

  def determine_maxcell(self,frame,pd):
    if self.phil_params.codecamp.maxcell != None:
      self.images[frame]['maxcel']=self.phil_params.codecamp.maxcell
      return
    if self.images[frame]['N_spots_inlier']>2:
      neighbors   = self.images[frame]['neighbors']
      n_ave,n_std = scitbx_stats(neighbors)

      average_nearest_neighbor = n_ave

      NNBIN = self.NspotMin//2 # recommended bin size for nearest neighbor histogram
      peak_of_interest = n_ave

      for pss in xrange(1):  #make two passes thru histo

        fineness = max( [self.pixel_size / 2.,
          peak_of_interest / (1.0+float(len(neighbors))/float(NNBIN))] )
        min_neighbors = min(neighbors)

        def histo_index_to_mm(index):
          return min_neighbors+index*fineness

        #neighbor_histogram
        histogram = [0]*int(peak_of_interest*4/fineness)
        for y in neighbors:
          ibin = int( (y-min_neighbors)/fineness )
          if ibin < len(histogram):
            histogram[ibin]+=1

        if TALLY2:
          for row in xrange(len(histogram)):
            low = histo_index_to_mm(row)
            hi  = histo_index_to_mm(row+1)
            print "%.2f to %.2f, %5.1f Ang"%(low,hi,radius_to_resol(histo_index_to_mm(row+.5),pd)),
            print "*"*histogram[row],
            print

        most_probable_neighbor = histo_index_to_mm(0.5 + histogram.index(max(histogram)))

        #Compute yet another measure of unit cell--first peak in histogram
        peak1 = 0; peak1i = 0
        for peakpt in xrange(len(histogram)):
          if histogram[peakpt]>peak1i:
            peak1=peakpt; peak1i=histogram[peakpt]
          if histogram[peakpt]<0.5*peak1i and peak1i> 0.1*(max(histogram)):
            break

        first_peak_if_any = histo_index_to_mm(0.5 + peak1)
        peak_of_interest = min(first_peak_if_any,most_probable_neighbor)

        # Another trial measure: 5th percentile neighbor
        sort_perm = flex.sort_permutation(neighbors)
        cutoff_mm = neighbors[sort_perm[ int(0.05*len(sort_perm)) ]]

      #Caveat: this rough calculation becomes less accurate for non-zero two-theta
      self.images[frame]['neighboring_spot_separation']=min(most_probable_neighbor , cutoff_mm)
      MAXTOL = 1.5 # Margin of error for max unit cell estimate
                   # 11/19/02 old value 1.4; new value 2.0 to accomodate a
                   #   pathological case where zone is down the large unit cell.
                   #   Larger value is now tolerable because of post-mosflm
                   #   check for systematic absences.
                   # 9/28/03 change back to 1.5 with new autoindexer because
                   #   value of 2.0 can (infrequently) lead to misindexing;
                   #   based on experience preparing figure 4.
      maxcell = max(MAXTOL * radius_to_resol(most_probable_neighbor,pd),
                    0.0    * radius_to_resol(min(neighbors),pd),
                    0.0    * radius_to_resol(first_peak_if_any,pd),
                    MAXTOL * radius_to_resol(cutoff_mm,pd))
      self.images[frame]['maxcel']=maxcell

  def setError(self,newmessage):
    if self.errormessage == None: self.errormessage=newmessage

  def get_resolution_inspection(self):
    all_frames = self.images.keys()
    all_resolutions=flex.double([self.images[f]['resolution'] for f in all_frames])
    ave_resolution=flex.mean(all_resolutions)
    self.pd['resolution_inspection']='%f'%(ave_resolution)
    return self.pd['resolution_inspection']
