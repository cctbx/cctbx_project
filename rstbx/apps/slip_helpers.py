import math
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from scitbx import lbfgs

class minimizer(object):
  def __init__(self,data,stuff,frame,max_iterations=50):

    self.n = len(data)
    self.values = data
    self.initial = data.deep_copy()


    #mosaicity,hi energy wave, lo energy wave
    #want mosaicity between 0.1 and 1.0 of input value.
    self.lower_bound = flex.double([0.1*data[0],0.98*data[1],(data[1]+data[2])/2.])
    upper_bound = flex.double([10.*data[0],(data[1]+data[2])/2.,1.02*data[2]])
    mean_value = (upper_bound - self.lower_bound)/2.
    self.full_range = upper_bound - self.lower_bound
    starting_params = flex.tan(math.pi*(((data - self.lower_bound)/self.full_range)-0.5))
    self.x = starting_params.deep_copy()
    print "staerting params",list(self.x)


    self.stuff = stuff
    self.frame = frame
    # optimize parameters
    self.minimizer = lbfgs.run(target_evaluator=self,
      termination_params = lbfgs.termination_parameters(max_calls=20))

  def compute_functional_and_gradients(self):
    print "trying",list(self.x)

    # try to constrain rather than restrain.  Use arctan function to get this right
    # minimizer can choose any value from -inf to inf but we'll keep values carefully in range
    # mos

    #if self.x[0]/self.initial[0] >2.: self.x[0]=2.*self.initial[0]
    #assert self.x[1]<self.x[2]
    #if self.x[1]/self.initial[1] <0.99: self.x[1]=0.99*self.initial[1]
    #if self.x[2]/self.initial[2] >1.01: self.x[2]=1.01*self.initial[2]
    #print "adjust",list(self.x)

    def get_score(x):
      unpacked = self.lower_bound + self.full_range*(
        0.5 + (flex.atan(x)/math.pi)
      )
      #print "          trying",list(unpacked)
      return self.stuff.use_case_3_score_only(
      self.frame,unpacked[0],unpacked[1],unpacked[2])

    #f = self.stuff.use_case_3_score_only(
    #  self.frame,self.values[0],self.values[1],self.values[2])
    f = get_score(self.x)

    gradients = flex.double(len(self.x))
    for i in xrange(self.n):
      #factors = [1.000]*self.n
     # factors[i] *= 1.001
     # D_i = 0.001 * self.x[i]
      #f_i_plus_D_i = self.stuff.use_case_3_score_only(
     #   self.frame,self.x[0]*factors[0],self.x[1]*factors[1],self.x[2]*factors[2])
      trial_x = self.x.deep_copy()
      trial_x[i]+=1.
      f_i_plus_D_i = get_score(trial_x)
      df_di = (f_i_plus_D_i - f)/1.
      gradients[i]=df_di
    return f,gradients

  def unpacked(self,x):
    return self.lower_bound + self.full_range*(
        0.5 + (flex.atan(x)/math.pi)
      )

class wrapper_of_use_case_bp3(object):
  def __init__(self, raw_image, spotfinder, imageindex, inputai, limiting_resolution, phil_params, sub=None):
    """MODEL:  polychromatic beam with top hat bandpass profile.
               isotropic mosaicity with top hat half-width; spots are brought into reflecting condition
               by a finite rotation about the axis that is longitudinal to the projection of the q-vector
               onto the detector plane.
    """
    from rstbx.bandpass import use_case_bp3,parameters_bp3
    from scitbx.matrix import col
    from math import pi
    from cctbx.crystal import symmetry

    self.detector_origin = col((-inputai.getBase().xbeam, -inputai.getBase().ybeam, 0.))
    crystal = symmetry(unit_cell=inputai.getOrientation().unit_cell(),space_group = "P1")
    indices = crystal.build_miller_set(anomalous_flag=True, d_min = limiting_resolution)
    parameters = parameters_bp3(
       indices=indices.indices(), orientation=inputai.getOrientation(),
       incident_beam=col((0.,0.,1.)),
       packed_tophat=col((1.,1.,0.)),
       detector_normal=col((0.,0.,-1.)), detector_fast=col((0.,1.,0.)),detector_slow=col((1.,0.,0.)),
       pixel_size=col((raw_image.pixel_size,raw_image.pixel_size,0)),
       pixel_offset=col((0.5,0.5,0.0)), distance=inputai.getBase().distance,
       detector_origin=self.detector_origin
    )

    self.ucbp3 = use_case_bp3(parameters=parameters)
    self.ucbp3.set_active_areas(
      raw_image.get_tile_manager(phil_params).effective_tiling_as_flex_int(
      reapply_peripheral_margin=True))
    if sub != None:  self.ucbp3.set_subpixel( flex.double(sub) )
    # Reduce Miller indices to a manageable set.  NOT VALID if the crystal rotates significantly
    self.ucbp3.prescreen_indices(inputai.wavelength)
    # done with Miller set reduction
    from annlib_ext import AnnAdaptorSelfInclude as AnnAdaptor
    body_pixel_reference = flex.double()
    for spot in spotfinder.images[imageindex]["goodspots"]:
      for pxl in spot.bodypixels:
        body_pixel_reference.append(pxl.y + 0.5)
        body_pixel_reference.append(pxl.x + 0.5)

    self.N_bodypix = body_pixel_reference.size()//2
    self.ucbp3.set_adaptor(body_pixel_reference)

  def set_variables(self, orientation, wave_HI, wave_LO, half_mosaicity_deg):
    half_mosaicity_rad = half_mosaicity_deg * math.pi/180.
    self.ucbp3.set_mosaicity(half_mosaicity_rad)
    self.ucbp3.set_bandpass(wave_HI,wave_LO)
    self.ucbp3.set_orientation(orientation)

  def score_only(self):
    self.ucbp3.picture_fast_slow()
    # not sure why x and y origin shifts are swapped here, but this seemed to work
    swapped_origin = (-self.detector_origin[1],-self.detector_origin[0],0.)
    self.ucbp3.spot_rectangles(swapped_origin)
    self.ucbp3.spot_rectregions(swapped_origin,1.0)
    self.ucbp3.enclosed_pixels_and_margin_pixels()

    return self.ucbp3.score_only_detail(weight=50.)


class slip_callbacks:
  def slip_callback(self,frame):
    #best_params=self.use_case_3_simulated_annealing()
    #best_params = self.use_case_3_grid_refine(frame)
    #self.inputai.setOrientation(best_params[3])
    #self.use_case_3_refactor(frame,best_params[0],best_params[1], best_params[2])

    normal = True
    # BLUE: predictions
    blue_data = []
    for ix,pred in enumerate(self.predicted):
        if self.BSmasks[ix].keys()==[]:continue
        x,y = frame.pyslip.tiles.picture_fast_slow_to_map_relative(
          (pred[1]/self.pixel_size) +0.5,
          (pred[0]/self.pixel_size) +0.5)
        blue_data.append((x,y))
    if normal: self.blue_layer = frame.pyslip.AddPointLayer(
          blue_data, color="blue", name="<blue_layer>",
          radius=2,
          renderer = frame.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])

    yellow_data = []; cyan_data = []
    for imsk in xrange(len(self.BSmasks)):
      smask_keys = self.get_ISmask(imsk)
      bmask = self.BSmasks[imsk]
      if len(bmask.keys())==0: continue

      # CYAN: integration mask
      for ks in xrange(0,len(smask_keys),2):
        cyan_data.append(
          frame.pyslip.tiles.picture_fast_slow_to_map_relative(
           smask_keys[ks+1] + 0.5,smask_keys[ks] + 0.5))

      # YELLOW: background mask
      for key in bmask.keys():
        yellow_data.append(
          frame.pyslip.tiles.picture_fast_slow_to_map_relative(
            key[1] + 0.5 ,key[0] + 0.5))
    if normal: self.cyan_layer = frame.pyslip.AddPointLayer(
          cyan_data, color="cyan", name="<cyan_layer>",
          radius=1.5,
          renderer = frame.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
    if normal: self.yellow_layer = frame.pyslip.AddPointLayer(
          yellow_data, color="yellow", name="<yellow_layer>",
          radius=1.5,
          renderer = frame.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])

    red_data = []; green_data = []

    for spot in self.spotfinder.images[self.frames[self.image_number]]["goodspots"]:
      # RED: spotfinder spot pixels
      for pxl in spot.bodypixels:
        red_data.append(
          frame.pyslip.tiles.picture_fast_slow_to_map_relative(
            pxl.y + 0.5, pxl.x + 0.5))

      # GREEN: spotfinder centers of mass
      green_data.append(
          frame.pyslip.tiles.picture_fast_slow_to_map_relative(
            spot.ctr_mass_y() + 0.5, spot.ctr_mass_x() + 0.5))

    self.red_layer = frame.pyslip.AddPointLayer(
          red_data, color="red", name="<red_layer>",
          radius=1.5,
          renderer = frame.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
    if normal: self.green_layer = frame.pyslip.AddPointLayer(
          green_data, color="green", name="<green_layer>",
          radius=1.5,
          renderer = frame.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])

  def use_case_1(self,frame):

    # A rehash of the spot-prediction algorithm for pedagogical use.
    # Use an "ewald proximity" filter as suggested by Ralf.
    # aside from a few extra spots due to ewald proximity, this is exactly the
    # same spot model developed initially for the Sept/Dec 2011 CXI runs.

    orange_data = []
    from scitbx.matrix import col,sqr
    print "wavelength",self.inputai.wavelength
    print "orientation",self.inputai.getOrientation()
    A = sqr(self.inputai.getOrientation().reciprocal_matrix())
    print "base",self.inputai.getBase()
    print "pixel size",self.pixel_size
    detector_origin = col((-self.inputai.getBase().xbeam, -self.inputai.getBase().ybeam, 0.))
    detector_fast = col((0.,1.,0.))
    detector_slow = col((1.,0.,0.))
    distance = self.inputai.getBase().distance

    s0 = col((0.,0.,1/self.inputai.wavelength))
    s0_length = s0.length()
    detector_normal = col((0.,0.,-1.))

    from cctbx.crystal import symmetry
    crystal = symmetry(unit_cell=self.inputai.getOrientation().unit_cell(),space_group = "P1")
    indices = crystal.build_miller_set(anomalous_flag=True, d_min = self.limiting_resolution)
    for H in indices.indices():
      s = A * H
      q = (s + s0)
      #print q.length(), s0_length
      if abs(q.length() - s0_length) > 0.001: continue
      q_unit = q.normalize()

      # check if diffracted ray parallel to detector face

      q_dot_n = q_unit.dot(detector_normal)

      if q_dot_n >= 0: continue

      r = (q_unit * distance / q_dot_n) - detector_origin

      x = r.dot(detector_fast)
      y = r.dot(detector_slow)
      print x,y

      orange_data.append( frame.pyslip.tiles.picture_fast_slow_to_map_relative(
          (x/self.pixel_size) +0.5,
          (y/self.pixel_size) +0.5))

    self.orange_layer = frame.pyslip.AddPointLayer(
          orange_data, color="orange", name="<orange_layer>",
          radius=3.0,
          renderer = frame.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])


  def use_case_2(self,frame):

    # Extend the model.  Assume monochromatic beam but finite radial mosaicity.  Dispense
    # with the "ewald proximity" mechanism; now spots are brought into reflecting condition
    # by a finite rotation about the axis that is longitudinal to the projection of the q-vector
    # onto the detector plane.

    orange_data = []
    from scitbx.matrix import col,sqr
    from math import pi
    print "Moasicity degrees, half",0.1
    mosaicity_rad = 0.1 * pi/180.  #half-width top-hat mosaicity
    A = sqr(self.inputai.getOrientation().reciprocal_matrix())
    detector_origin = col((-self.inputai.getBase().xbeam, -self.inputai.getBase().ybeam, 0.))
    detector_fast = col((0.,1.,0.))
    detector_slow = col((1.,0.,0.))
    distance = self.inputai.getBase().distance

    #s0:  parallel to the direction of incident radiation
    s0 = col((0.,0.,1/self.inputai.wavelength))
    s0_length = s0.length()
    s0_unit = s0.normalize()
    detector_normal = col((0.,0.,-1.))
    #  Cn, the circular section through the Ewald sphere.
    Cncenter = -s0
    Cnradius_squared = s0.length_sq()
    #  Taking a page from mathworld.wolfram.com, calculate the distance d
    #  between the centers of Co and Cn,
    d = s0_length

    from cctbx.crystal import symmetry
    crystal = symmetry(unit_cell=self.inputai.getOrientation().unit_cell(),space_group = "P1")
    indices = crystal.build_miller_set(anomalous_flag=True, d_min = self.limiting_resolution)
    for H in indices.indices():
      s = A * H
      rotax = s.normalize().cross(s0_unit) #The axis that most directly brings the Bragg spot onto Ewald sphere
      s_rad = s.length()
      s_rad_sq = s.length_sq()

      # take a page from ewald_sphere.cpp, determine intersection of two coplanar circles
      #  Co, the circle centered on reciprocal origin and containing the point s,
      #  Cn, the circle centered on -s0 (ewald sphere center) of radius (1/lambda) with normal rotax.
      # Consider the intersection of two circles:
      #  Co, the circle of rotation of H.
      # Cocenter = 0; so it falls out of the equations

      #   The chord of intersection between Co and Cn lies a
      #   distance x along the (Cocenter - Cncenter) vector
      chord_direction =      (rotax.cross( - Cncenter)).normalize();

      a = s.length_sq()/(2.*s0_length)
      b = math.sqrt(s.length_sq() - (a*a))      #  Calculate half-length of the chord of intersection
      #  Two intersection points
      intersections_0p = -a * s0_unit+ b*chord_direction
      intersections_1p = -a * s0_unit- b*chord_direction
      iangle_0= math.acos (intersections_0p.dot(s) / (s_rad_sq))
      iangle_1= math.acos (intersections_1p.dot(s) / (s_rad_sq))

      assert approx_equal((intersections_0p+s0).length()-s0_length,0. )

      if iangle_0 < mosaicity_rad:
        intersection = intersections_0p
      elif iangle_1 < mosaicity_rad:
        intersection = intersections_1p
      else: continue

      q = (intersection + s0)
      q_unit = q.normalize()

      # check if diffracted ray parallel to detector face

      q_dot_n = q_unit.dot(detector_normal)

      if q_dot_n >= 0: continue
      print "IANGLES",iangle_0 * 180./pi, iangle_1 * 180./pi

      r = (q_unit * distance / q_dot_n) - detector_origin

      x = r.dot(detector_fast)
      y = r.dot(detector_slow)

      orange_data.append( frame.pyslip.tiles.picture_fast_slow_to_map_relative(
          (x/self.pixel_size) +0.5,
          (y/self.pixel_size) +0.5))

    self.orange_layer = frame.pyslip.AddPointLayer(
          orange_data, color="orange", name="<orange_layer>",
          radius=3.0,
          renderer = frame.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])

  def use_case_2cpp(self,frame):
    from rstbx.bandpass import use_case_bp2_picture_fast_slow
    # Extend the model.  Assume monochromatic beam but finite radial mosaicity.  Dispense
    # with the "ewald proximity" mechanism; now spots are brought into reflecting condition
    # by a finite rotation about the axis that is longitudinal to the projection of the q-vector
    # onto the detector plane.

    from scitbx.matrix import col
    from math import pi

    detector_origin = col((-self.inputai.getBase().xbeam, -self.inputai.getBase().ybeam, 0.))

    from cctbx.crystal import symmetry
    crystal = symmetry(unit_cell=self.inputai.getOrientation().unit_cell(),space_group = "P1")
    indices = crystal.build_miller_set(anomalous_flag=True, d_min = self.limiting_resolution)

    picture_fast_slow = use_case_bp2_picture_fast_slow(
       indices=indices.indices(), orientation=self.inputai.getOrientation(),
       incident_beam=col((0.,0.,1.)), wavelength=self.inputai.wavelength,
       detector_normal=col((0.,0.,-1.)), detector_fast=col((0.,1.,0.)),detector_slow=col((1.,0.,0.)),
       pixel_size=col((self.pixel_size,self.pixel_size,0)),
       pixel_offset=col((0.5,0.5,0.0)), distance=self.inputai.getBase().distance,
       detector_origin=detector_origin,
       half_mosaicity_rad=0.1 * pi/180.
    )

    map_relative = frame.pyslip.tiles.vec_picture_fast_slow_to_map_relative(picture_fast_slow)

    self.orange_layer = frame.pyslip.AddPointLayer(
          map_relative, color="orange", name="<orange_layer>",
          radius=3.0,
          renderer = frame.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])

  def use_case_3cpp(self,frame):
    from rstbx.bandpass import use_case_bp3_picture_fast_slow
    # Extend the model.  Assume polychromatic beam with top hat profile.  Assume finite radial mosaicity.  Dispense
    # with the "ewald proximity" mechanism; now spots are brought into reflecting condition
    # by a finite rotation about the axis that is longitudinal to the projection of the q-vector
    # onto the detector plane.

    from scitbx.matrix import col
    from math import pi

    detector_origin = col((-self.inputai.getBase().xbeam, -self.inputai.getBase().ybeam, 0.))

    from cctbx.crystal import symmetry
    crystal = symmetry(unit_cell=self.inputai.getOrientation().unit_cell(),space_group = "P1")
    indices = crystal.build_miller_set(anomalous_flag=True, d_min = self.limiting_resolution)

    cpp_results = use_case_bp3_picture_fast_slow(
       indices=indices.indices(), orientation=self.inputai.getOrientation(),
       incident_beam=col((0.,0.,1.)),
       #tophat=col((self.inputai.wavelength,self.inputai.wavelength+0.00001,0.1*pi/180.)),
       tophat=col((self.inputai.wavelength*0.9975,self.inputai.wavelength*1.0025,0.1*pi/180.)),
       detector_normal=col((0.,0.,-1.)), detector_fast=col((0.,1.,0.)),detector_slow=col((1.,0.,0.)),
       pixel_size=col((self.pixel_size,self.pixel_size,0)),
       pixel_offset=col((0.5,0.5,0.0)), distance=self.inputai.getBase().distance,
       detector_origin=detector_origin
    )
    picture_fast_slow = cpp_results[0].select(cpp_results[2])
    map_relative = frame.pyslip.tiles.vec_picture_fast_slow_to_map_relative(picture_fast_slow)

    self.yellow_layer = frame.pyslip.AddPointLayer(
          map_relative, color="yellow", name="<yellow_layer>",
          radius=3.0,
          renderer = frame.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
    picture_fast_slow = cpp_results[1].select(cpp_results[2])
    map_relative = frame.pyslip.tiles.vec_picture_fast_slow_to_map_relative(picture_fast_slow)

    self.red_layer = frame.pyslip.AddPointLayer(
          map_relative, color="red", name="<red_layer>",
          radius=3.0,
          renderer = frame.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])

  def use_case_3box(self,frame):
    from rstbx.bandpass import use_case_bp3_picture_fast_slow
    # Extend the model.  Assume polychromatic beam with top hat profile.  Assume finite radial mosaicity.  Dispense
    # with the "ewald proximity" mechanism; now spots are brought into reflecting condition
    # by a finite rotation about the axis that is longitudinal to the projection of the q-vector
    # onto the detector plane.

    from scitbx.matrix import col
    from math import pi

    detector_origin = col((-self.inputai.getBase().xbeam, -self.inputai.getBase().ybeam, 0.))

    from cctbx.crystal import symmetry
    crystal = symmetry(unit_cell=self.inputai.getOrientation().unit_cell(),space_group = "P1")
    indices = crystal.build_miller_set(anomalous_flag=True, d_min = self.limiting_resolution)
    half_mosaicity_rad = 0.1*pi/180.
    cpp_results = use_case_bp3_picture_fast_slow(
       indices=indices.indices(), orientation=self.inputai.getOrientation(),
       incident_beam=col((0.,0.,1.)),
       tophat=col((self.inputai.wavelength*0.9975,self.inputai.wavelength*1.0025,half_mosaicity_rad)),
       detector_normal=col((0.,0.,-1.)), detector_fast=col((0.,1.,0.)),detector_slow=col((1.,0.,0.)),
       pixel_size=col((self.pixel_size,self.pixel_size,0)),
       pixel_offset=col((0.5,0.5,0.0)), distance=self.inputai.getBase().distance,
       detector_origin=detector_origin
    )
    picture_fast_slow = cpp_results[0].select(cpp_results[2])
    map_relative_hi = frame.pyslip.tiles.vec_picture_fast_slow_to_map_relative(picture_fast_slow)
    picture_fast_slow = cpp_results[1].select(cpp_results[2])
    map_relative_lo = frame.pyslip.tiles.vec_picture_fast_slow_to_map_relative(picture_fast_slow)

    # not sure if I've swapped x/y correctly
    beam_coor = frame.pyslip.tiles.vec_picture_fast_slow_to_map_relative(
      [(0.5 + self.inputai.getBase().ybeam/self.pixel_size,
        0.5 + self.inputai.getBase().xbeam/self.pixel_size)])

    polydata = []
    beam_pos = col(beam_coor[0])
    for idx in xrange(len(map_relative_hi)):
      hi_pos = col(map_relative_hi[idx])
      lo_pos = col(map_relative_lo[idx])
      radial_vector = (hi_pos-beam_pos)
      radial_unit_vec = radial_vector.normalize()
      radius = radial_vector.length()
      tangential_unit_vec = col((-radial_unit_vec[1],radial_unit_vec[0])) # 90-degree rotation
      tangential_excursion = tangential_unit_vec * radius * half_mosaicity_rad
      polydata.append( ([(hi_pos + tangential_excursion).elems,
                         (hi_pos - tangential_excursion).elems,
                         (lo_pos - tangential_excursion).elems,
                         (lo_pos + tangential_excursion).elems,
                         (hi_pos + tangential_excursion).elems
                     ],{}) )

    self.red_layer = frame.pyslip.AddPolygonLayer( # needs to be changed for Linux (antialiasing removed)
          polydata, color="red", name="<red_layer>",
          width=1.0,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])

  def use_case_3_refactor(self,frame,half_deg,wave_HI, wave_LO):
    from rstbx.bandpass import use_case_bp3,parameters_bp3
    # Extend the model.  Assume polychromatic beam with top hat profile.  Assume finite radial mosaicity.  Dispense
    # with the "ewald proximity" mechanism; now spots are brought into reflecting condition
    # by a finite rotation about the axis that is longitudinal to the projection of the q-vector
    # onto the detector plane.

    from scitbx.matrix import col
    from math import pi

    detector_origin = col((-self.inputai.getBase().xbeam, -self.inputai.getBase().ybeam, 0.))

    from cctbx.crystal import symmetry
    crystal = symmetry(unit_cell=self.inputai.getOrientation().unit_cell(),space_group = "P1")
    indices = crystal.build_miller_set(anomalous_flag=True, d_min = self.limiting_resolution)
    half_mosaicity_rad = half_deg*pi/180.
    parameters = parameters_bp3(
       indices=indices.indices(), orientation=self.inputai.getOrientation(),
       incident_beam=col((0.,0.,1.)),
       packed_tophat=col((wave_HI,wave_LO,half_mosaicity_rad)),
       detector_normal=col((0.,0.,-1.)), detector_fast=col((0.,1.,0.)),detector_slow=col((1.,0.,0.)),
       pixel_size=col((self.pixel_size,self.pixel_size,0)),
       pixel_offset=col((0.5,0.5,0.0)), distance=self.inputai.getBase().distance,
       detector_origin=detector_origin
    )

    cpp_results = use_case_bp3(parameters=parameters)
    cpp_results.set_active_areas(
      frame.pyslip.tiles.raw_image.get_tile_manager(frame.inherited_params).effective_tiling_as_flex_int(
      reapply_peripheral_margin=True))
    cpp_results.picture_fast_slow()
    picture_fast_slow = cpp_results.hi_E_limit.select(cpp_results.observed_flag)
    map_relative_hi = frame.pyslip.tiles.vec_picture_fast_slow_to_map_relative(picture_fast_slow)
    picture_fast_slow = cpp_results.lo_E_limit.select(cpp_results.observed_flag)
    map_relative_lo = frame.pyslip.tiles.vec_picture_fast_slow_to_map_relative(picture_fast_slow)

    poly = cpp_results.spot_rectangles((self.inputai.getBase().ybeam,self.inputai.getBase().xbeam,0.))
    map_relative_poly = frame.pyslip.tiles.vec_picture_fast_slow_to_map_relative(poly)
    cpp_polydata = []
    for idx in xrange(0,len(map_relative_poly),5):
      cpp_polydata.append( ([ map_relative_poly[idx+0],
                              map_relative_poly[idx+1],
                              map_relative_poly[idx+2],
                              map_relative_poly[idx+3],
                              map_relative_poly[idx+4]
                     ],{}) )

    self.red_layer = frame.pyslip.AddPolygonLayer( # needs to be changed for Linx (antialiasing removed)
          cpp_polydata, color="red", name="<red_layer>",
          width=1.0,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])

    poly = cpp_results.spot_rectregions((self.inputai.getBase().ybeam,self.inputai.getBase().xbeam,0.),1.0)
    map_relative_poly = frame.pyslip.tiles.vec_picture_fast_slow_to_map_relative(poly)
    cpp_polydata = []
    for idx in xrange(0,len(map_relative_poly),5):
      cpp_polydata.append( ([ map_relative_poly[idx+0],
                              map_relative_poly[idx+1],
                              map_relative_poly[idx+2],
                              map_relative_poly[idx+3],
                              map_relative_poly[idx+4]
                     ],{}) )

    self.pink_layer = frame.pyslip.AddPolygonLayer( # needs to be changed for Linx (antialiasing removed)
          cpp_polydata, color="pink", name="<pink_layer>",
          width=1.0,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
    print "Entering C++ pixels"
    cpp_results.enclosed_pixels_and_margin_pixels()
    print "Done with C++ pixels"
    internal = cpp_results.enclosed_px
    map_relative_pixels = frame.pyslip.tiles.vec_picture_fast_slow_to_map_relative(internal)
    self.yellow_layer = frame.pyslip.AddPointLayer(
          map_relative_pixels, color="yellow", name="<yellow_layer>",
          radius=1.0,
          renderer = frame.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
    internal = cpp_results.margin_px
    map_relative_pixels = frame.pyslip.tiles.vec_picture_fast_slow_to_map_relative(internal)
    self.green_layer = frame.pyslip.AddPointLayer(
          map_relative_pixels, color="green", name="<green_layer>",
          radius=1.0,
          renderer = frame.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])

    map_relative_pixels = frame.pyslip.tiles.vec_picture_fast_slow_to_map_relative(
                          cpp_results.selected_predictions())

    self.blue_layer = frame.pyslip.AddPointLayer(
          map_relative_pixels, color="blue", name="<blue_layer>",
          radius=2.0,
          renderer = frame.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])

    # Now figure out how to implement the scoring function. Do this in picture fast low coordinates
    # FIRST Set:  spot bodypixels:
    from annlib_ext import AnnAdaptorSelfInclude as AnnAdaptor
    body_pixel_reference = flex.double()
    for spot in self.spotfinder.images[self.frames[self.image_number]]["goodspots"]:
      for pxl in spot.bodypixels:
        body_pixel_reference.append(pxl.y + 0.5)
        body_pixel_reference.append(pxl.x + 0.5)
    self.adapt = AnnAdaptor(data=body_pixel_reference,dim=2,k=1)

    N_bodypix = body_pixel_reference.size()//2

    # second set:  predict box
    enclosed_pixels = cpp_results.enclosed_px
    N_enclosed = enclosed_pixels.size()
    N_enclosed_body_pixels = 0
    query = flex.double()
    for pixel in enclosed_pixels:
      query.append(pixel[0]); query.append(pixel[1])
    self.adapt.query(query)
    import math
    for p in xrange(N_enclosed):
      if math.sqrt(self.adapt.distances[p]) < 0.1:
        N_enclosed_body_pixels += 1

    # third set:  marginal
    marginal_pixels = cpp_results.margin_px
    margin_distances = cpp_results.margin_distances

    N_marginal = marginal_pixels.size()
    N_marginal_body_pixels = 0
    marginal_body = 0
    marginal_nonbody = 0
    query = flex.double()
    for pixel in marginal_pixels:
      query.append(pixel[0]); query.append(pixel[1])
    self.adapt.query(query)
    for p in xrange(N_marginal):
      if math.sqrt(self.adapt.distances[p]) < 0.1:
        N_marginal_body_pixels += 1
        marginal_body += 0.5 + 0.5 * math.cos (-math.pi * margin_distances[p]) #taking MARGIN==1
      else:
        marginal_nonbody += 0.5 + 0.5 * math.cos (math.pi * margin_distances[p])
    print "marginal body/nonbody",marginal_body, marginal_nonbody

    print "There are %d body pixels of which %d are enclosed and %d are marginal leaving %d remote"%(
      N_bodypix,N_enclosed_body_pixels,N_marginal_body_pixels,
      N_bodypix-N_enclosed_body_pixels-N_marginal_body_pixels)
    print "There are %d enclosed pixels of which %d are body, %d are nonbody"%(
      N_enclosed,N_enclosed_body_pixels,N_enclosed-N_enclosed_body_pixels)
    print "There are %d marginal pixels of which %d are body, %d are nonbody"%(
      N_marginal,N_marginal_body_pixels,N_marginal-N_marginal_body_pixels)
    Score = 0
    # the scoring function to account for these spots:
    #    pink -- spot body pixels inside predict box = 0
    #    red  -- spot body pixels > 2 pxl away from predict box = 1
    Score += N_bodypix-N_enclosed_body_pixels-N_marginal_body_pixels
    #    gradation -- body pixels within the marginal zone
    Score += marginal_body + marginal_nonbody
    #    blank -- nonbody pixels outside of 2 pixel margin = 0
    #    yellow -- nonbody pixels inside predict box = 1
    Score += N_enclosed-N_enclosed_body_pixels
    #    gradation -- in between zone, within margin
    print "The score is",Score

  def use_case_3_score_only(self,frame,half_deg,wave_HI, wave_LO):
    from rstbx.bandpass import use_case_bp3,parameters_bp3
    # Extend the model.  Assume polychromatic beam with top hat profile.  Assume finite radial mosaicity.  Dispense
    # with the "ewald proximity" mechanism; now spots are brought into reflecting condition
    # by a finite rotation about the axis that is longitudinal to the projection of the q-vector
    # onto the detector plane.

    from scitbx.matrix import col
    from math import pi

    detector_origin = col((-self.inputai.getBase().xbeam, -self.inputai.getBase().ybeam, 0.))

    from cctbx.crystal import symmetry
    crystal = symmetry(unit_cell=self.inputai.getOrientation().unit_cell(),space_group = "P1")
    indices = crystal.build_miller_set(anomalous_flag=True, d_min = self.limiting_resolution)
    half_mosaicity_rad = half_deg*pi/180.
    parameters = parameters_bp3(
       indices=indices.indices(), orientation=self.inputai.getOrientation(),
       incident_beam=col((0.,0.,1.)),
       packed_tophat=col((wave_HI,wave_LO,half_mosaicity_rad)),
       detector_normal=col((0.,0.,-1.)), detector_fast=col((0.,1.,0.)),detector_slow=col((1.,0.,0.)),
       pixel_size=col((self.pixel_size,self.pixel_size,0)),
       pixel_offset=col((0.5,0.5,0.0)), distance=self.inputai.getBase().distance,
       detector_origin=detector_origin
    )

    cpp_results = use_case_bp3(parameters=parameters)
    cpp_results.set_active_areas(
      frame.pyslip.tiles.raw_image.get_tile_manager(frame.inherited_params).effective_tiling_as_flex_int(
      reapply_peripheral_margin=True))
    cpp_results.picture_fast_slow()
    poly = cpp_results.spot_rectangles((self.inputai.getBase().ybeam,self.inputai.getBase().xbeam,0.))

    poly = cpp_results.spot_rectregions((self.inputai.getBase().ybeam,self.inputai.getBase().xbeam,0.),1.0)
    cpp_results.enclosed_pixels_and_margin_pixels()

    # Now figure out how to implement the scoring function. Do this in picture fast low coordinates
    # FIRST Set:  spot bodypixels:
    from annlib_ext import AnnAdaptorSelfInclude as AnnAdaptor
    body_pixel_reference = flex.double()
    for spot in self.spotfinder.images[self.frames[self.image_number]]["goodspots"]:
      for pxl in spot.bodypixels:
        body_pixel_reference.append(pxl.y + 0.5)
        body_pixel_reference.append(pxl.x + 0.5)
    self.adapt = AnnAdaptor(data=body_pixel_reference,dim=2,k=1)

    N_bodypix = body_pixel_reference.size()//2

    # second set:  predict box
    enclosed_pixels = cpp_results.enclosed_px
    N_enclosed = enclosed_pixels.size()
    N_enclosed_body_pixels = 0
    query = flex.double()
    for pixel in enclosed_pixels:
      query.append(pixel[0]); query.append(pixel[1])
    self.adapt.query(query)
    import math
    for p in xrange(N_enclosed):
      if math.sqrt(self.adapt.distances[p]) < 0.1:
        N_enclosed_body_pixels += 1

    # third set:  marginal
    marginal_pixels = cpp_results.margin_px
    margin_distances = cpp_results.margin_distances
    WGT = 50.
    N_marginal = marginal_pixels.size()
    N_marginal_body_pixels = 0
    marginal_body = 0
    marginal_nonbody = 0
    query = flex.double()
    for pixel in marginal_pixels:
      query.append(pixel[0]); query.append(pixel[1])
    self.adapt.query(query)
    for p in xrange(N_marginal):
      if math.sqrt(self.adapt.distances[p]) < 0.1:
        N_marginal_body_pixels += 1
        marginal_body += 0.5 + 0.5 * math.cos (-math.pi * margin_distances[p]) #taking MARGIN==1
      else:
        marginal_nonbody += 0.5 + 0.5 * math.cos (math.pi * margin_distances[p])
    marginal_body *=WGT
    if False:
      print "marginal body/nonbody",marginal_body, marginal_nonbody
      print "There are %d body pixels of which %d are enclosed and %d are marginal leaving %d remote"%(
        N_bodypix,N_enclosed_body_pixels,N_marginal_body_pixels,
        N_bodypix-N_enclosed_body_pixels-N_marginal_body_pixels)
      print "There are %d enclosed pixels of which %d are body, %d are nonbody"%(
        N_enclosed,N_enclosed_body_pixels,N_enclosed-N_enclosed_body_pixels)
      print "There are %d marginal pixels of which %d are body, %d are nonbody"%(
        N_marginal,N_marginal_body_pixels,N_marginal-N_marginal_body_pixels)
    Score = 0
    # the scoring function to account for these spots:
    #    pink -- spot body pixels inside predict box = 0
    #    red  -- spot body pixels > 2 pxl away from predict box = 1
    Score += WGT*(N_bodypix-N_enclosed_body_pixels-N_marginal_body_pixels)
    #    gradation -- body pixels within the marginal zone
    Score += marginal_body + marginal_nonbody
    #    blank -- nonbody pixels outside of 2 pixel margin = 0
    #    yellow -- nonbody pixels inside predict box = 1
    Score += N_enclosed-N_enclosed_body_pixels
    #    gradation -- in between zone, within margin
    return Score

  def use_case_3_grid_refine(self,frame):
    reserve_orientation = self.inputai.getOrientation()

    wrapbp3 = wrapper_of_use_case_bp3( raw_image = frame.pyslip.tiles.raw_image,
      spotfinder = self.spotfinder, imageindex = self.frames[self.image_number],
      inputai = self.inputai,
      limiting_resolution = self.limiting_resolution,
      phil_params = frame.inherited_params)
    wrapbp3.set_variables( orientation = self.inputai.getOrientation(),
                           wave_HI = self.inputai.wavelength*0.9975,
                           wave_LO = self.inputai.wavelength*1.0025,
                           half_mosaicity_deg = 0.1)
    #print "score...",wrapbp3.score_only()

    wave_HI = self.inputai.wavelength*0.9975
    wave_LO = self.inputai.wavelength*1.0025
    low_score = None
    for half_deg in [0.06, 0.08, 0.10, 0.12, 0.14]:
      for bandpass in [0.004, 0.005, 0.006, 0.007, 0.008]:
        for mean_multiplier in [0.9990, 1.0000, 1.0010, 1.0020]:
#               A1=0.;A2=0.;A3=0.
          for A1 in (math.pi/180.)*flex.double([-0.1,0.0,0.1]):
            for A2 in (math.pi/180.)*flex.double([-0.1,0.0,0.1]):
              for A3 in (math.pi/180.)*flex.double([-0.1,0.0,0.1]):
                ori = reserve_orientation.rotate_thru((1,0,0),A1).rotate_thru((0,1,0),A2).rotate_thru((0,0,1),A3)
                self.inputai.setOrientation(ori)
                HI = self.inputai.wavelength*(mean_multiplier-(bandpass/2.))
                LO = self.inputai.wavelength*(mean_multiplier+(bandpass/2.))
                #score = self.use_case_3_score_only(
                #       frame,half_deg,HI,LO)
                wrapbp3.set_variables( orientation = self.inputai.getOrientation(),
                           wave_HI = HI,wave_LO = LO,half_mosaicity_deg = half_deg)

                score = wrapbp3.score_only()
                if low_score == None or score < low_score:
                  low_score = score
                  best_params = (half_deg,HI,LO,ori,A1,A2,A3)
                  print "wave %7.4f - %7.4f bandpass %.2f half %7.4f score %7.1f"%(HI,LO,100.*(LO-HI)/LO,half_deg,score)
    print "Rendering image with wave %7.4f - %7.4f bandpass %.2f half %7.4f score %7.1f"%(
      best_params[1],best_params[2],100.*(best_params[2]-best_params[1])/best_params[1],best_params[0],low_score)
    print "rotation angles",best_params[4],best_params[5],best_params[6]
    return best_params

  def use_case_3_simulated_annealing(self,subpixel=None):
    reserve_orientation = self.inputai.getOrientation()

    wrapbp3 = wrapper_of_use_case_bp3( raw_image = self.imagefiles.images[self.image_number],
      spotfinder = self.spotfinder, imageindex = self.frames[self.image_number],
      inputai = self.inputai,
      limiting_resolution = self.limiting_resolution,
      phil_params = self.horizons_phil,
      sub = subpixel)

    from rstbx.bandpass.simulated_annealing import SALight

    SA = SALight()

    # Half mosaicity in degrees
    # Mid-wavelength adjustment factor
    # Bandpass fractional full width
    # adjustment angle in degrees
    # adjustment angle in degrees
    # adjustment angle in degrees

    # starting values; likely expected values
    SA.x = flex.double([0.1,1.00,0.006,0.0,0.0,0.0])
    SA.initial = SA.x.deep_copy()

    # reasonable length scale (expected interval, half width)
    SA.L = flex.double([0.02,0.001,0.001,0.05,0.05,0.05])

    SA.format = "Mosaicity %6.3f Wave mean %7.4f bandpass %7.4f Angles %8.5f %8.5f %8.5f"

    def set_variables_from_sa_x(x):
      ori = reserve_orientation.rotate_thru((1,0,0),(math.pi/180.)*x[3]
                              ).rotate_thru((0,1,0),(math.pi/180.)*x[4]
                              ).rotate_thru((0,0,1),(math.pi/180.)*x[5])
      self.inputai.setOrientation(ori)
      mean_multiplier = x[1]
      bandpass = x[2]
      HI = self.inputai.wavelength*(mean_multiplier-(bandpass/2.))
      LO = self.inputai.wavelength*(mean_multiplier+(bandpass/2.))
      wrapbp3.set_variables( orientation = self.inputai.getOrientation(),
                           wave_HI = HI,wave_LO = LO,half_mosaicity_deg = x[0])
      #pack into format for calling function
      these_params = (x[0],HI,LO,ori,(math.pi/180.)*x[3],(math.pi/180.)*x[4],(math.pi/180.)*x[5])
      return these_params

    set_variables_from_sa_x(SA.x)
    last_score = wrapbp3.score_only()
    low_score = last_score + 0 # makes a copy
    Tstart = 600
    for T in xrange(Tstart, 1, -1):
      decreasing_increment = (T/Tstart)*SA.random_increment()
      last_x = SA.x.deep_copy()
      test_params = SA.x + decreasing_increment
      if test_params[2]<=0: continue # can't have negative bandpass; unphysical!
      if test_params[0]<=0: continue # can't have negative mosaicity; unphysical!
      SA.x += decreasing_increment
      print T, SA.format%tuple(SA.x),
      set_variables_from_sa_x(SA.x)
      new_score = wrapbp3.score_only()
      print "Score %8.1f"%new_score,

      if new_score < low_score:
        low_score = 1.0*new_score
      if new_score < last_score:
        probability_of_acceptance=1.0
      else:
        probability_of_acceptance = math.exp(-(new_score-last_score)/(2.5*T))
      if flex.random_double(1)[0] < probability_of_acceptance:
        #new position accepted
        last_score = 1.0*new_score
        print "accepted"
      else:
        SA.x = last_x.deep_copy()
        print "rejected"

    print "Final"
    print T, SA.format%tuple(SA.x),"Score %8.1f"%last_score,"final"

    #these three lines set the bp3 wrapper so it can be used from the calling class (simple_integration.py)
    best_params = set_variables_from_sa_x(SA.x)
    wrapbp3.score_only()
    self.bp3_wrapper = wrapbp3

    print "Rendering image with wave %7.4f - %7.4f bandpass %.2f half %7.4f score %7.1f"%(
      best_params[1],best_params[2],100.*(best_params[2]-best_params[1])/best_params[1],best_params[0],last_score)
    print "rotation angles",best_params[4],best_params[5],best_params[6]
    return best_params
