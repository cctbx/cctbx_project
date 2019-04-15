
# TODO:
#  - cached scenes

from __future__ import division
from libtbx.utils import Sorry
import libtbx.phil
from libtbx import object_oriented_patterns as oop
from math import sqrt
import math

def generate_systematic_absences(array,
                                  expand_to_p1=False,
                                  generate_bijvoet_mates=False):
  from cctbx import crystal
  sg = array.space_group()
  base_symm = crystal.symmetry(
    unit_cell=array.unit_cell(),
    space_group=sg.build_derived_reflection_intensity_group(False))
  base_array = array.customized_copy(crystal_symmetry=base_symm)
  complete_set = base_array.complete_set().customized_copy(
    crystal_symmetry=array.crystal_symmetry())
  absence_array = complete_set.sys_absent_flags()
  if (generate_bijvoet_mates):
    absence_array = absence_array.generate_bijvoet_mates()
  if (expand_to_p1):
    absence_array = absence_array.expand_to_p1().customized_copy(
      crystal_symmetry=array.crystal_symmetry())
      #niggli_cell().expand_to_p1()
  return absence_array


class scene(object):
  """
  Data for visualizing a Miller array graphically, either as a 3D view or
  a 2D zone of reciprocal space.  Currently used in phenix.data_viewer, but
  easily extensible to any other graphics platform (e.g. a PNG embedded in
  a web page).
  """
  def __init__(self, miller_array, settings, merge=None, foms=None):
    self.miller_array = miller_array
    self.renderscale = 100.0
    self.foms = foms
    if self.miller_array.is_complex_array():
      # want to display map coefficient as circular colours but weighted with FOMS
      # so copy any provided foms in the empty list of sigmas
      if self.foms:
        assert ( self.miller_array.size() == self.foms.size() )
        #self.miller_array._sigmas = foms.deep_copy()
      #else:
      #  self.miller_array._sigmas = None
    self.settings = settings
    self.merge_equivalents = merge
    from cctbx import miller
    from cctbx import crystal
    from cctbx.array_family import flex
    self.multiplicities = None
    self.process_input_array()
    array = self.work_array
    uc = array.unit_cell()
    self.unit_cell = uc
    self.slice_selection = None
    self.axis_index = None
    if (settings.slice_mode):
      self.axis_index = ["h","k","l"].index(self.settings.slice_axis)
      self.slice_selection = miller.simple_slice(
        indices=array.indices(),
        slice_axis=self.axis_index,
        slice_index=settings.slice_index)
      #if (self.slice_selection.count(True) == 0):
        #raise ValueError("No data selected!")
    index_span = array.index_span()
    self.d_min = array.d_min()
    self.hkl_range = index_span.abs_range()
    self.axes = [ uc.reciprocal_space_vector((self.hkl_range[0],0,0)),
                  uc.reciprocal_space_vector((0,self.hkl_range[1],0)),
                  uc.reciprocal_space_vector((0,0,self.hkl_range[2])) ]
    self.generate_view_data()
    if (self.slice_selection is not None):
      self.indices = self.work_array.indices().select(self.slice_selection)
      self.data = self.data.select(self.slice_selection)
      self.phases = self.phases.select(self.slice_selection)
      self.radians = self.radians.select(self.slice_selection)
      self.ampl = self.ampl.select(self.slice_selection)
      self.foms = self.foms.select(self.slice_selection)
    else :
      self.indices = array.indices()
    self.points = uc.reciprocal_space_vector(self.indices) * self.renderscale
    self.missing_flags = flex.bool(self.radii.size(), False)
    self.sys_absent_flags = flex.bool(self.radii.size(), False)
    if (settings.show_missing):
      self.generate_missing_reflections()
    if (settings.show_systematic_absences) and (not settings.show_only_missing):
      self.generate_systematic_absences()
    # XXX hack for process_pick_points
    self.visible_points = flex.bool(self.points.size(), True)
    n_points = self.points.size()
    assert (self.colors.size() == n_points)
    assert (self.indices.size() == n_points)
    assert (self.radii.size() == n_points)
    assert (self.missing_flags.size() == n_points)
    assert (self.sys_absent_flags.size() == n_points)
    assert (self.data.size() == n_points)
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    assert (self.phases.size() == n_points)
    assert (self.radians.size() == n_points)
    assert (self.ampl.size() == n_points)
    if foms:
      assert (self.foms.size() == n_points)
    self.clear_labels()


  def process_input_array(self):
    from cctbx.array_family import flex
    array = self.miller_array.deep_copy()
    multiplicities = None
    if self.merge_equivalents :
      if self.settings.show_anomalous_pairs:
        from cctbx import miller
        merge = array.merge_equivalents()
        multiplicities = merge.redundancies()
        asu, matches = multiplicities.match_bijvoet_mates()
        mult_plus, mult_minus = multiplicities.hemispheres_acentrics()
        anom_mult = flex.int(
          min(p, m) for (p, m) in zip(mult_plus.data(), mult_minus.data()))
        #flex.min_max_mean_double(anom_mult.as_double()).show()
        anomalous_multiplicities = miller.array(
          miller.set(asu.crystal_symmetry(),
                     mult_plus.indices(),
                     anomalous_flag=False), anom_mult)
        anomalous_multiplicities = anomalous_multiplicities.select(
          anomalous_multiplicities.data() > 0)

        array = anomalous_multiplicities
        multiplicities = anomalous_multiplicities

      else:
        merge = array.merge_equivalents()
        array = merge.array()
        multiplicities = merge.redundancies()
    settings = self.settings
    data = array.data()
    self.missing_set = oop.null()
    #if (array.is_xray_intensity_array()):
    #  data.set_selected(data < 0, flex.double(data.size(), 0.))
    if (array.is_unique_set_under_symmetry()) and (settings.map_to_asu):
      array = array.map_to_asu()
      if (multiplicities is not None):
        multiplicities = multiplicities.map_to_asu()
    if (settings.d_min is not None):
      array = array.resolution_filter(d_min=settings.d_min)
      if (multiplicities is not None):
        multiplicities = multiplicities.resolution_filter(
          d_min=settings.d_min)
    self.filtered_array = array.deep_copy()
    if (settings.expand_anomalous):
      array = array.generate_bijvoet_mates()
      if (multiplicities is not None):
        multiplicities = multiplicities.generate_bijvoet_mates()
    if (self.settings.show_missing):
      self.missing_set = array.complete_set().lone_set(array)
      if self.settings.show_anomalous_pairs:
        self.missing_set = self.missing_set.select(
          self.missing_set.centric_flags().data(), negate=True)
    if (settings.expand_to_p1):
      original_symmetry = array.crystal_symmetry()
      if array.is_complex_array() and self.foms:
        tmp_millarr = miller.array( miller.set(original_symmetry, array.indices() ),
                                data=self.foms)
        self.foms = self.tmp_millarr.expand_to_p1().customized_copy(
          crystal_symmetry=original_symmetry).data().deep_copy()
      array = array.expand_to_p1().customized_copy(
        crystal_symmetry=original_symmetry)
      #array = array.niggli_cell().expand_to_p1()
      #self.missing_set = self.missing_set.niggli_cell().expand_to_p1()
      self.missing_set = self.missing_set.expand_to_p1().customized_copy(
        crystal_symmetry=original_symmetry)
      if (multiplicities is not None):
        multiplicities = multiplicities.expand_to_p1().customized_copy(
            crystal_symmetry=original_symmetry)
    data = array.data()
    self.r_free_mode = False
    self.phases = flex.double(data.size(), float('nan'))
    self.radians = flex.double(data.size(), float('nan'))
    self.ampl = flex.double(data.size(), float('nan'))
    if not self.foms:
      self.foms = flex.double(data.size(), float('nan'))
    self.sigmas = None
    #self.sigmas = flex.double(data.size(), float('nan'))
    if isinstance(data, flex.bool):
      self.r_free_mode = True
      data_as_float = flex.double(data.size(), 0.0)
      data_as_float.set_selected(data==True, flex.double(data.size(), 1.0))
      data = data_as_float
      self.data = data.deep_copy()
    else :
      if isinstance(data, flex.double):
        self.data = data.deep_copy()
      elif isinstance(data, flex.complex_double):
        self.data = data.deep_copy()
        self.ampl = flex.abs(data)
        self.phases = flex.arg(data) * 180.0/math.pi
        # purge nan values from array to avoid crash in fmod_positive()
        b = flex.bool([bool(math.isnan(e)) for e in self.phases])
        # replace the nan values with an arbitrary float value
        self.phases = self.phases.set_selected(b, 42.4242)
        # indicate the coresponding phase/radian is completely undetermnined
        self.foms = self.foms.set_selected(b, 0.0)
        # Now cast negative degrees to equivalent positive degrees
        self.phases = flex.fmod_positive(self.phases, 360.0)
        self.radians = flex.arg(data)
        # replace the nan values with an arbitrary float value
        self.radians = self.radians.set_selected(b, 0.424242)
      elif hasattr(array.data(), "as_double"):
        self.data = array.data().as_double()
      else:
        raise RuntimeError("Unexpected data type: %r" % data)
      if (settings.show_data_over_sigma):
        if (array.sigmas() is None):
          raise Sorry("sigmas not defined.")
        sigmas = array.sigmas()
        non_zero_sel = sigmas != 0
        array = array.select(non_zero_sel)
        array = array.customized_copy(data=array.data()/array.sigmas())
        self.data = array.data()
        if (multiplicities is not None):
          multiplicities = multiplicities.select(non_zero_sel)
      if array.sigmas() is not None:
        self.sigmas = array.sigmas()
      else:
        self.sigmas = None
    self.work_array = array
    self.work_array.set_info(self.miller_array.info() )
    self.multiplicities = multiplicities


  def generate_view_data(self):
    from scitbx.array_family import flex
    from scitbx import graphics_utils
    settings = self.settings
    data_for_colors = data_for_radii = None
    data = self.data #self.work_array.data()
    sigmas = self.sigmas
    if (isinstance(data, flex.double) and data.all_eq(0)):
      data = flex.double(data.size(), 1)
    if ((self.multiplicities is not None) and
        (settings.scale_colors_multiplicity)):
      data_for_colors = self.multiplicities.data().as_double()
      assert data_for_colors.size() == data.size()
    elif (settings.sqrt_scale_colors) and (isinstance(data, flex.double)):
      data_for_colors = flex.sqrt(flex.abs(data))
    elif isinstance(data, flex.complex_double):
      data_for_colors = self.radians
      # when using map coefficients the sigmas are filled with foms if provided
      #if self.sigmas:
      #  self.foms = self.sigmas
      foms_for_colours = self.foms
    elif (settings.sigma_color) and sigmas is not None:
      data_for_colors = sigmas.as_double()
    else :
      data_for_colors = flex.abs(data.deep_copy())
    if ((self.multiplicities is not None) and
        (settings.scale_radii_multiplicity)):
      #data_for_radii = data.deep_copy()
      data_for_radii = self.multiplicities.data().as_double()

      if (settings.sigma_radius) and sigmas is not None:
        data_for_radii = sigmas * self.multiplicities.as_double()
        #print "sigmas: " + self.miller_array.info().label_string()
      assert data_for_radii.size() == data.size()
    #elif (settings.sqrt_scale_radii) and (isinstance(data, flex.double)):
    #  data_for_radii = flex.sqrt(flex.abs(data))
    elif (settings.sigma_radius) and sigmas is not None:
      data_for_radii = sigmas.as_double()
      #print "sigmas: " + self.miller_array.info().label_string()
    else :
      #data_for_radii = flex.abs(data.deep_copy())
      data_for_radii = flex.pow(flex.abs(data.deep_copy()), settings.nth_root_scale_radii)

    if (settings.slice_mode):
      data = data.select(self.slice_selection)
      if (not settings.keep_constant_scale):
        data_for_radii = data_for_radii.select(self.slice_selection)
        data_for_colors = data_for_colors.select(self.slice_selection)
        foms_for_colours = foms_for_colours.select(self.slice_selection)
    if isinstance(data, flex.complex_double):
      colors = graphics_utils.colour_by_phi_FOM(data_for_colors, foms_for_colours)
    elif (settings.color_scheme in ["rainbow", "heatmap", "redblue"]):
      colors = graphics_utils.color_by_property(
        properties=data_for_colors,
        selection=flex.bool(data_for_colors.size(), True),
        color_all=False,
        gradient_type=settings.color_scheme)
    elif (settings.color_scheme == "grayscale"):
      colors = graphics_utils.grayscale_by_property(
        properties=data_for_colors,
        selection=flex.bool(data_for_colors.size(), True),
        shade_all=False,
        invert=settings.black_background)
    else :
      if (settings.black_background):
        base_color = (1.0,1.0,1.0)
      else :
        base_color = (0.0,0.0,0.0)
      colors = flex.vec3_double(data_for_colors.size(), base_color)
    if (settings.slice_mode) and (settings.keep_constant_scale):
      colors = colors.select(self.slice_selection)
      data_for_radii = data_for_radii.select(self.slice_selection)
    uc = self.work_array.unit_cell()
    #abc = uc.parameters()[0:3]
    min_dist = min(uc.reciprocal_space_vector((1,1,1)))
    min_radius = 0.20 * min_dist
    max_radius = 40 * min_dist

    #if (settings.sqrt_scale_radii) and (not settings.scale_radii_multiplicity):
    #  data_for_radii = flex.sqrt(flex.abs(data_for_radii))
    if len(data_for_radii):
      max_value = flex.max(data_for_radii)
      scale = max_radius / max_value
      radii = data_for_radii * (scale * self.settings.scale)
      too_small = radii < min_radius
      if (too_small.count(True) > 0):
        radii.set_selected(too_small, flex.double(radii.size(), min_radius))
      assert radii.size() == colors.size()
    else:
      radii = flex.double()
      max_radius = 0
    self.radii = radii
    self.max_radius = max_radius
    self.colors = colors
    if isinstance(data, flex.complex_double):
      self.foms = foms_for_colours


  def generate_missing_reflections(self):
    from cctbx import miller
    from cctbx.array_family import flex
    settings = self.settings
    array = self.work_array
    uc = array.unit_cell()
    if (settings.show_only_missing):
      self.colors = flex.vec3_double()
      self.points = flex.vec3_double()
      self.radii = flex.double()
      self.missing_flags = flex.bool()
      self.indices = flex.miller_index()
      self.data = flex.double()
      self.sys_absent_flags = flex.bool()
    if (settings.slice_mode):
      slice_selection = miller.simple_slice(
        indices=self.missing_set.indices(),
        slice_axis=self.axis_index,
        slice_index=settings.slice_index)
      missing = self.missing_set.select(slice_selection).indices()
    else :
      missing = self.missing_set.indices()
    n_missing = missing.size()
    if (n_missing > 0):
      points_missing = uc.reciprocal_space_vector(missing) * 100.
      self.points.extend(points_missing)
      if (settings.color_scheme == "heatmap"):
        self.colors.extend(flex.vec3_double(n_missing, (0.,1.,0.)))
      elif (not settings.color_scheme in ["rainbow","redblue"]):
        self.colors.extend(flex.vec3_double(n_missing, (1.,0,0)))
      else :
        self.colors.extend(flex.vec3_double(n_missing, (1.,1.,1.)))
      self.radii.extend(flex.double(n_missing, self.max_radius / 2))
      self.missing_flags.extend(flex.bool(n_missing, True))
      self.indices.extend(missing)
      self.data.extend(flex.double(n_missing, -1.))
      self.sys_absent_flags.extend(flex.bool(n_missing, False))


  def generate_systematic_absences(self):
    from cctbx import miller
    from cctbx import crystal
    from cctbx.array_family import flex
    settings = self.settings
    array = self.filtered_array # XXX use original array here!
    absence_array = generate_systematic_absences(
      array=self.filtered_array,
      expand_to_p1=settings.expand_to_p1,
      generate_bijvoet_mates=settings.expand_anomalous)
    if (settings.slice_mode):
      slice_selection = miller.simple_slice(
        indices=absence_array.indices(),
        slice_axis=self.axis_index,
        slice_index=settings.slice_index)
      absence_array = absence_array.select(slice_selection)
    absence_flags = absence_array.data()
    if (absence_flags.count(True) == 0):
      print "No systematic absences found!"
    else :
      new_indices = flex.miller_index()
      for i_seq in absence_flags.iselection():
        hkl = absence_array.indices()[i_seq]
        #if (hkl in self.indices):
        #  j_seq = self.indices.index(hkl)
        #  self.colors[j_seq] = (1.0, 0.5, 1.0)
        #  self.radii[j_seq] = self.max_radius
        #  self.missing_flags[j_seq] = False
        #  self.sys_absent_flags[j_seq] = True
        #else :
        new_indices.append(hkl)
      if (new_indices.size() > 0):
        uc = self.work_array.unit_cell()
        points = uc.reciprocal_space_vector(new_indices) * 100
        self.points.extend(points)
        n_sys_absent = new_indices.size()
        self.radii.extend(flex.double(new_indices.size(), self.max_radius))
        self.indices.extend(new_indices)
        self.missing_flags.extend(flex.bool(new_indices.size(), False))
        self.sys_absent_flags.extend(flex.bool(new_indices.size(), True))
        self.data.extend(flex.double(n_sys_absent, -1.))
        if (settings.color_scheme == "redblue"):
          self.colors.extend(flex.vec3_double(new_indices.size(), (1.,1.0,0.)))
        else :
          self.colors.extend(flex.vec3_double(new_indices.size(), (1.,0.5,1.)))


  def clear_labels(self):
    self.label_points = set([])


  def get_resolution_at_point(self, k):
    hkl = self.indices[k]
    return self.unit_cell.d(hkl)


  def get_reflection_info(self, k):
    hkl = self.indices[k]
    d_min = self.unit_cell.d(hkl)
    if (self.missing_flags[k]) or (self.sys_absent_flags[k]):
      value = None
    else :
      value = self.data[k]
    return (hkl, d_min, value)


class render_2d(object):
  def __init__(self, scene, settings):
    self.scene = scene
    self.settings = settings
    self.setup_colors()


  def GetSize(self):
    raise NotImplementedError()


  def get_center_and_radius(self):
    w, h = self.GetSize()
    r = (min(w,h) // 2) - 20
    center_x = max(w // 2, r + 20)
    center_y = max(h // 2, r + 20)
    return center_x, center_y, r


  def setup_colors(self):
    if (self.settings.black_background):
      self._background = (0.,0.,0.)
      self._foreground = (1.,1.,1.)
      if (self.settings.color_scheme == "heatmap"):
        self._missing = (0.,1.,0.)
      elif (not self.settings.color_scheme in ["rainbow", "redblue"]):
        self._missing = (1.,0.,0.)
      else :
        self._missing = (1.,1.,1.)
    else :
      self._background = (1.,1.,1.)
      self._foreground = (0.,0.,0.)
      if (self.settings.color_scheme == "heatmap"):
        self._missing = (0.,1.,0.)
      elif (not self.settings.color_scheme in ["rainbow", "redblue"]):
        self._missing = (1.,0.,0.)
      else :
        self._missing = (0.,0.,0.)


  def get_scale_factor(self):
    return 100.


  def render(self, canvas):
    self._points_2d = []
    self._radii_2d = []
    assert (self.settings.slice_mode)
    if (self.settings.slice_axis == "h"):
      i_x, i_y = 1, 2
      axes = ("k", "l")
    elif (self.settings.slice_axis == "k"):
      i_x, i_y = 0, 2
      axes = ("h", "l")
    else :
      i_x, i_y = 0, 1
      axes = ("h", "k")
    center_x, center_y, r = self.get_center_and_radius()
    x_max = self.scene.axes[i_x][i_x] * 100.
    y_max = self.scene.axes[i_y][i_y] * 100.
    if (self.settings.show_axes):
      # FIXME dimensions not right?
      x_end = self.scene.axes[i_x][i_x], self.scene.axes[i_x][i_y]
      y_end = self.scene.axes[i_y][i_x], self.scene.axes[i_y][i_y]
      x_len = sqrt(x_end[0]**2 + x_end[1]**2)
      y_len = sqrt(y_end[0]**2 + y_end[1]**2)
      x_scale = (r+10) / x_len
      y_scale = (r+10) / y_len
      x_end = (x_end[0] * x_scale, x_end[1] * x_scale)
      y_end = (y_end[0] * y_scale, y_end[1] * y_scale)
      self.draw_line(canvas, center_x, center_y, center_x+x_end[0],
        center_y-x_end[1])
      self.draw_line(canvas, center_x, center_y, center_x+y_end[0],
        center_y-y_end[1])
      self.draw_text(canvas, axes[0], center_x + x_end[0] - 6,
        center_y - x_end[1] - 20)
      self.draw_text(canvas, axes[1], center_x + y_end[0] + 6,
        center_y - y_end[1])
    max_radius = self.scene.max_radius * r / max(x_max, y_max)
    max_radius *= self.settings.scale
    r_scale = ( 1/ self.scene.d_min) * self.get_scale_factor() # FIXME
    for k, hkl in enumerate(self.scene.points):
      x_, y_ = hkl[i_x], hkl[i_y]
      x = center_x + r * x_ / r_scale
      y = center_y - r * y_ / r_scale
      r_point = self.scene.radii[k] * r / max(x_max, y_max)
      if (self.settings.uniform_size):
        r_point = max_radius
      else :
        r_point = max(0.5, r_point)
      r_point *= self.settings.scale
      self._points_2d.append((x,y))
      self._radii_2d.append(r_point)
      if (self.scene.missing_flags[k]):
        self.draw_open_circle(canvas, x, y, r_point)
      elif (self.scene.sys_absent_flags[k]):
        self.draw_open_circle(canvas, x, y, r_point, self.scene.colors[k])
      else :
        self.draw_filled_circle(canvas, x, y, r_point, self.scene.colors[k])


  def draw_line(self, canvas, x1, y1, x2, y2):
    raise NotImplementedError()


  def draw_text(self, canvas, text, x, y):
    raise NotImplementedError()


  def draw_open_circle(self, canvas, x, y, radius, color=None):
    raise NotImplementedError()


  def draw_filled_circle(self, canvas, x, y, radius, color):
    raise NotImplementedError()


philstr = """
  data = None
    .type = path
    .optional = False
  labels = None
    .type = str
  symmetry_file = None
    .type = str
  black_background = True
    .type = bool
  show_axes = True
    .type = bool
  show_data_over_sigma = False
    .type = bool
  nth_root_scale_radii = 1.0
    .type = int
  sqrt_scale_colors = False
    .type = bool
  phase_color = False
    .type = bool
  sigma_color = False
    .type = bool
  sigma_radius = False
    .type = bool
  expand_to_p1 = False
    .type = bool
  expand_anomalous = False
    .type = bool
  spheres = True
    .type = bool
  show_missing = False
    .type = bool
  show_only_missing = False
    .type = bool
  show_systematic_absences = False
    .type = bool
  sphere_detail = 20
    .type = int
  pad_radii = False
    .type = bool
  slice_mode = False
    .type = bool
  keep_constant_scale = True
    .type = bool
  slice_axis = *h k l
    .type = choice
  slice_index = 0
    .type = int
  color_scheme = *rainbow heatmap redblue grayscale mono
    .type = choice
  show_labels = True
    .type = bool
  uniform_size = False
    .type = bool
  d_min = None
    .type = float
  scale = 1
    .type = int
  map_to_asu = False
    .type = bool
  scale_radii_multiplicity = False
    .type = bool
  scale_colors_multiplicity = False
    .type = bool
  show_anomalous_pairs = False
    .type = bool
"""

master_phil = libtbx.phil.parse( philstr )
params = master_phil.fetch().extract()

def reset_settings():
  """
  Reset all display settings to their default values as specified in the phil definition string
  """
  global master_phil
  global philstr
  global params
  params = master_phil.fetch(source = libtbx.phil.parse( philstr) ).extract()


def settings():
  """
  Get a global phil parameters object containing the display settings
  """
  global params
  return params


