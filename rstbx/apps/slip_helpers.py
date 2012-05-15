import math
from libtbx.test_utils import approx_equal

class slip_callbacks:
  def slip_callback(self,frame):

    self.use_case_3box(frame)

    # BLUE: predictions
    blue_data = []
    for ix,pred in enumerate(self.predicted):
        if self.BSmasks[ix].keys()==[]:continue
        x,y = frame.pyslip.tiles.picture_fast_slow_to_map_relative(
          (pred[1]/self.pixel_size) +0.5,
          (pred[0]/self.pixel_size) +0.5)
        blue_data.append((x,y))
    self.blue_layer = frame.pyslip.AddPointLayer(
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
    self.cyan_layer = frame.pyslip.AddPointLayer(
          cyan_data, color="cyan", name="<cyan_layer>",
          radius=1.5,
          renderer = frame.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
    self.yellow_layer = frame.pyslip.AddPointLayer(
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
    self.green_layer = frame.pyslip.AddPointLayer(
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


# are there more spots than inliers?  all good spots?
# can we refine the 3-parameters to get a better model based on spot positions?
#figure out a good correlation coefficient matching predictions & spotfinder bodypixels.
#can the cc be continuous-valued instead of discrete.
