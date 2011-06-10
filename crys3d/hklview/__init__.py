
# TODO:
#  - cached scenes

from libtbx.utils import Sorry
from math import sqrt

class scene (object) :
  def __init__ (self, miller_array, settings) :
    self.miller_array = miller_array
    from cctbx import miller
    from scitbx import graphics_utils
    from cctbx.array_family import flex
    array = self.miller_array.deep_copy()
    if (array.is_unique_set_under_symmetry()) :
      array = array.map_to_asu()
    if (array.is_xray_intensity_array()) :
      data = array.data()
      data.set_selected(data < 0, flex.double(data.size(), 0.))
      array = array.customized_copy(data=data)
    if (settings.show_data_over_sigma) :
      if (array.sigmas() is None) :
        raise Sorry("sigmas not defined.")
      sigmas = array.sigmas()
      data = array.data()
      array = array.select(sigmas != 0).customized_copy(data=data/sigmas)
    if (settings.d_min is not None) :
      array = array.resolution_filter(d_min=settings.d_min)
    uc = array.unit_cell()
    self.unit_cell = uc
    if (settings.expand_to_p1) :
      array = array.expand_to_p1()
    if (settings.expand_anomalous) :
      array = array.generate_bijvoet_mates()
    slice_selection = None
    if (settings.slice_mode) :
      slice_selection = miller.simple_slice(
        indices=array.indices(),
        slice_axis=settings.slice_axis,
        slice_index=settings.slice_index)
    index_span = array.index_span()
    self.hkl_range = index_span.abs_range()
    self.axes = [ uc.reciprocal_space_vector((self.hkl_range[0],0,0)),
                  uc.reciprocal_space_vector((0,self.hkl_range[1],0)),
                  uc.reciprocal_space_vector((0,0,self.hkl_range[2])) ]
    indices = array.indices()
    data = array.data()
    assert isinstance(data, flex.double) or isinstance(data, flex.bool)
    if isinstance(data, flex.bool) :
      data = flex.double(data.size(), 1.0)
    if (settings.sqrt_scale_colors) :
      data_for_colors = flex.sqrt(data)
    else :
      data_for_colors = data
    if (settings.sqrt_scale_radii) :
      data_for_radii = flex.sqrt(data)
    else :
      data_for_radii = data
    if (settings.color_scheme == 0) :
      colors = graphics_utils.color_by_property(
        properties=data_for_colors,
        selection=flex.bool(data.size(), True),
        color_all=False,
        use_rb_color_gradient=False)
    elif (settings.color_scheme == 1) :
      colors = graphics_utils.grayscale_by_property(
        properties=data_for_colors,
        selection=flex.bool(data.size(), True),
        shade_all=False,
        invert=settings.black_background)
    else :
      if (settings.black_background) :
        base_color = (1.0,1.0,1.0)
      else :
        base_color = (0.0,0.0,0.0)
      colors = flex.vec3_double(data.size(), base_color)
    if (slice_selection is not None) :
      data = data.select(slice_selection)
      if (data.size() == 0) :
        raise ValueError("No data selected!")
      indices = indices.select(slice_selection)
      if (settings.keep_constant_scale) :
        colors = colors.select(slice_selection)
      elif (settings.color_scheme == 0) :
        colors = graphics_utils.color_by_property(
          properties=data_for_colors.select(slice_selection),
          selection=flex.bool(data.size(), True),
          color_all=False,
          use_rb_color_gradient=False,
          grayscale=(settings.color_scheme == 1))
      elif (settings.color_scheme == 1) :
        colors = graphics_utils.grayscale_by_property(
          properties=data_for_colors.select(slice_selection),
          selection=flex.bool(data.size(), True),
          shade_all=False,
          invert=settings.black_background)
      elif (settings.color_scheme == 2) :
        colors = flex.vec3_double(data.size(), (1.0,1.0,1.0))
    self.colors = colors
    self.points = uc.reciprocal_space_vector(indices) * 100.
    self.indices = indices
    abc = uc.parameters()[0:3]
    min_radius = 0.20 / max(abc)
    max_radius = 50 / max(abc)
    scale = max_radius / flex.max(data_for_radii)
    if (settings.sqrt_scale_radii) :
      data = flex.sqrt(data)
    radii = data * scale
    too_small = radii < min_radius
    radii.set_selected(too_small, flex.double(radii.size(), min_radius))
    self.radii = radii
    self.max_radius = flex.max(radii)
    self.missing = flex.bool(self.radii.size(), False)
    if (settings.show_missing_reflections) :
      if (settings.show_only_missing) :
        self.colors = flex.vec3_double()
        self.points = flex.vec3_double()
        self.radii = flex.double()
        self.missing = flex.bool()
        self.indices = flex.miller_index()
      complete_set = array.complete_set()
      if (settings.slice_mode) :
        slice_selection = miller.simple_slice(
          indices=complete_set.indices(),
          slice_axis=settings.slice_axis,
          slice_index=settings.slice_index)
        missing = complete_set.select(slice_selection).lone_set(array).indices()
      else :
        missing = complete_set.lone_set(array).indices()
      n_missing = missing.size()
      if (n_missing > 0) :
        points_missing = uc.reciprocal_space_vector(missing) * 100.
        self.points.extend(points_missing)
        if (settings.color_scheme != 0) :
          self.colors.extend(flex.vec3_double(n_missing, (1.,0,0)))
        else :
          self.colors.extend(flex.vec3_double(n_missing, (1.,1.,1.)))
        self.radii.extend(flex.double(n_missing, max_radius / 2))
        self.missing.extend(flex.bool(n_missing, True))
        self.indices.extend(missing)
    # XXX hack for process_pick_points
    self.visible_points = flex.bool(self.points.size(), True)
    assert (self.colors.size() == self.points.size() == self.indices.size() ==
            self.radii.size() == self.missing.size())
    self.clear_labels()

  def clear_labels (self) :
    self.label_points = set([])

  def get_resolution_at_point (self, k) :
    hkl = self.indices[k]
    return self.unit_cell.d(hkl)

class settings (object) :
  def __init__ (self) :
    self.black_background = True
    self.show_axes = True
    self.show_data_over_sigma = False
    self.sqrt_scale_radii = True
    self.sqrt_scale_colors = False
    self.expand_to_p1 = False
    self.expand_anomalous = False
    self.display_as_spheres = True
    self.show_missing_reflections = False
    self.show_only_missing = False
    self.sphere_detail = 20
    self.pad_radii = False
    self.slice_mode = False
    self.keep_constant_scale = True
    self.slice_axis = 0
    self.slice_index = 0
    self.color_scheme = 0
    self.show_labels = True
    self.uniform_size = False
    self.d_min = None
