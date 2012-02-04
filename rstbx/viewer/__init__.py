
import rstbx.utils
from libtbx.math_utils import ifloor, iceil
from libtbx.str_utils import format_value
import math

pi_over_180 = math.pi / 180

class screen_params (object) :
  """
  Manager for all display parameters: this is independent of the actual image
  data, although it stores various attributes such as detector dimensions.
  The primary function is to convert between different coordinate systems and
  determine which part of the image to display.
  """
  def __init__ (self, img_w=None, img_h=None) :
    self.img_w = img_w
    self.img_h = img_h
    self.screen_w = None
    self.screen_h = None
    self.thumb_w = None
    self.thumb_h = None
    self.thumb_ratio = None
    self.img_x_offset = 0
    self.img_y_offset = 0
    self.screen_x_start = 0
    self.screen_y_start = 0
    self.detector_pixel_size = 0
    self.zoom = 0
    self.last_thumb_x = 0 # NKS: hooks for keeping pan position while
    self.last_thumb_y = 0 #  rendering the Prev or Next image

  def inherit_params (self, params) :
    """
    Adopts the current screen parameters (to preserve zoom level, offset,
    etc.), but only if the image dimensions are the same.
    """
    if (self.img_w != params.img_w) or (self.img_h != params.img_h) :
      return
    self.img_x_offset = params.img_x_offset
    self.img_y_offset = params.img_y_offset
    self.screen_x_start = params.screen_x_start
    self.screen_y_start = params.screen_y_start
    self.last_thumb_x = params.last_thumb_x
    self.last_thumb_y = params.last_thumb_y
    self.zoom = params.zoom

  def set_zoom (self, zoom) :
    assert (zoom >= 0)
    # XXX don't do anything fancy if image is uninitialized (for zoom view)
    if (None in [self.screen_w, self.screen_h, self.img_w, self.img_h]) :
      self.zoom = zoom
      return
    # XXX adjust offsets to preserve current center
    x0, y0, w0, h0 = self.get_bitmap_params()
    increase_zoom = (zoom > self.zoom)
    decrease_zoom = (zoom < self.zoom)
    self.zoom = zoom
    center_x, center_y = int(x0 + w0/2), int(y0 + h0/2)
    x, y, w, h = self.get_bitmap_params()
    if (increase_zoom) :
      if ((x + w) < self.img_w) :
        self.img_x_offset = max(min(int(center_x - w/2), self.img_w - w), 0)
      if ((y + h) < self.img_h) :
        self.img_y_offset = max(min(int(center_y - h/2), self.img_h - h), 0)
    elif (decrease_zoom) :
      self.img_x_offset = max(int(center_x - w/2), 0)
      self.img_y_offset = max(int(center_y - h/2), 0)

  def set_screen_size (self, w, h) :
    self.screen_w = w
    self.screen_h = h
    scale = self.get_scale()
    w_img_pixels = int(w / scale)
    h_img_pixels = int(h / scale)
    if (w_img_pixels > (self.img_w - self.img_x_offset)) :
      self.img_x_offset = max(0, self.img_w - w_img_pixels)
    if (h_img_pixels > (self.img_h - self.img_y_offset)) :
      self.img_y_offset = max(0, self.img_h - w_img_pixels)

  def set_image_size (self, w, h) :
    self.img_w = w
    self.img_h = h

  def set_thumbnail_size (self, w, h, ratio) :
    self.thumb_w = w
    self.thumb_h = h
    self.thumb_ratio = ratio

  def get_image_size (self) :
    assert (not None in [self.img_w, self.img_h])
    return self.img_w, self.img_h

  def get_thumbnail_size (self) :
    assert (not None in [self.thumb_w, self.thumb_h])
    return (self.thumb_w, self.thumb_h)

  def set_detector_resolution (self, p) :
    self.detector_pixel_size = p

  def get_scale (self) :
    assert (not None in [self.screen_w, self.screen_h, self.img_w, self.img_h])
    if (self.zoom == 0) :
      return min(self.screen_w / self.img_w, self.screen_h / self.img_h)
    else :
      return self.zoom

  def get_detector_dimensions (self) :
    assert (not None in [self.img_w, self.img_h, self.detector_pixel_size])
    return (self.img_w * self.detector_pixel_size,
            self.img_h * self.detector_pixel_size)

  def adjust_screen_coordinates (self, x, y) :
    xi, yi, w, h = self.get_bitmap_params()
    scale = self.get_scale()
    x_ = x + max(0, (self.screen_w - (w*scale)) / 2)
    y_ = y + max(0, (self.screen_h - (h*scale)) / 2)
    return (int(x_), int(y_))

  def get_bitmap_params (self) :
    scale = self.get_scale()
    x0 = max(self.img_x_offset, 0)
    y0 = max(self.img_y_offset, 0)
    w_scaled = min(self.img_w, (self.screen_w / scale)) #- x0
    h_scaled = min(self.img_h, (self.screen_h / scale))
    return (x0, y0, w_scaled, h_scaled)

  def get_thumbnail_box (self) :
    x, y, w, h = self.get_bitmap_params()
    tr = self.thumb_ratio
    return int(x / tr), int(y / tr), int(w / tr), int(h / tr)

  def get_zoom_box (self, x, y, boxsize=400, mag=16) :
    #assert ((boxsize % mag) == 0)
    n_pixels = iceil(boxsize / mag)
    x0 = min(self.img_w - n_pixels, ifloor(x - (n_pixels / 2)))
    y0 = min(self.img_h - n_pixels, ifloor(y - (n_pixels / 2)))
    return (x0, y0, n_pixels, n_pixels)

  def translate_image (self, delta_x, delta_y) :
    """
    Translate the viewport to a different area of the image.  Arguments are
    in pixels.
    """
    scale = self.get_scale()
    x_new = max(0, ifloor(self.img_x_offset - (delta_x / scale)))
    y_new = max(0, ifloor(self.img_y_offset - (delta_y / scale)))
    max_x = ifloor(self.img_w - (self.screen_w / scale))
    max_y = ifloor(self.img_h - (self.screen_h / scale))
    self.img_x_offset = min(x_new, max_x)
    self.img_y_offset = min(y_new, max_y)

  def center_view_from_thumbnail (self, x, y) :
    """
    Translate the viewport to center on the X,Y coordinates equivalent to the
    point clicked in the thumbnail view.  Arguments are in screen coordinates
    relative to the upper left corner of the thumbnail (which is assumed to be
    displayed in its entirety).
    """
    if (self.zoom == 0) : return
    self.last_thumb_x = x
    self.last_thumb_y = y
    x0, y0, w, h = self.get_bitmap_params()
    img_x = max(0, ifloor((x * self.thumb_ratio) - (w / 2)))
    img_y = max(0, ifloor((y * self.thumb_ratio) - (h / 2)))
    scale = self.get_scale()
    max_x = ifloor(self.img_w - (self.screen_w / scale))
    max_y = ifloor(self.img_h - (self.screen_h / scale))
    self.img_x_offset = min(img_x, max_x)
    self.img_y_offset = min(img_y, max_y)

  def image_coords_as_screen_coords (self, x, y) :
    """
    Convert image pixel coordinates to viewport pixel coordinates.
    """
    scale = self.get_scale()
    x1 = self.screen_x_start + ((x+0.5) - self.img_x_offset) * scale
    y1 = self.screen_y_start + ((y+0.5) - self.img_y_offset) * scale
    xi, yi, w, h = self.get_bitmap_params()
    x2 = x1 + max(0, (self.screen_w - (w*scale)) / 2)
    y2 = y1 + max(0, (self.screen_h - (h*scale)) / 2)
    return ((x2), (y2))

  def screen_coords_as_image_coords (self, x, y) :
    """
    Convert pixel coordinates in the viewport to pixel coordinates in the
    raw image.
    """
    scale = self.get_scale()
    xi, yi, w, h = self.get_bitmap_params()
    x1 = x - max(0, (self.screen_w - (w*scale)) / 2)
    y1 = y - max(0, (self.screen_h - (h*scale)) / 2)
    x2 = self.img_x_offset + (x1 / scale)
    y2 = self.img_y_offset + (y1 / scale)
    return (ifloor(x2) + 1, ifloor(y2) + 1)

  def image_coords_as_array_coords (self, x, y) :
    """
    Convert image pixel coordinates to indexing values in the FlexImage
    object.
    """
    return y-1, x-1

  def array_coords_as_detector_coords (self, x, y) :
    """
    Convert array indices to points on the detector surface.  Used in the
    calculation of approximate lattice dimensions based on peaks in a
    user-drawn line in the viewer.
    """
    x_, y_ = y+1, x+1
    return self._raw.image_coords_as_detector_coords(x_, y_)

  def distance_between_points (self, x1, y1, x2, y2) :
    """
    Given a pair of image pixel coordinates, calculate the distance between
    them on the detector face in mm.
    """
    x1_mm, y1_mm = self._raw.image_coords_as_detector_coords(x1, y1)
    x2_mm, y2_mm = self._raw.image_coords_as_detector_coords(x2, y2)
    return math.sqrt((x1_mm - x2_mm)**2 + (y1_mm - y2_mm)**2)

class image (screen_params) :
  def __init__ (self, file_name) :
    screen_params.__init__(self)
    self.file_name = file_name
    from iotbx.detectors import ImageFactory
    img = ImageFactory(file_name)
    img.read()
    self._raw = img
    print img.show_header()
    self.set_image_size(
      w=self._raw.size2,
      h=self._raw.size1)
    self.set_detector_resolution(self._raw.pixel_size)
    from spotfinder.command_line.signal_strength import master_params
    params = master_params.extract()
    self._raw.initialize_viewer_properties(params)
    self._beam_center = None
    self._integration = None
    self._spots = None
    #self.update_image()
    #self.convert_to_bitmap()

  def set_integration_results (self, integration) :
    self._integration = integration
    mp = integration['mapped_predictions']
    print "%d spot predictions loaded" % mp.size()
    print "max. resolution is %g" % integration['resolution']

  def set_spots (self, spots) :
    self._spots = spots

  def set_beam_center (self, xbeam, ybeam) :
    self._beam_center = (xbeam, ybeam)

  def create_flex_image (self,
                         brightness=100,
                         color_scheme=0,
                         binning=1) :
    fi = self._raw.get_flex_image(
      binning=binning,
      brightness=brightness / 100.,
    )
    #from scitbx.array_family import flex
    #print flex.max(self._raw.linearintdata), flex.min(self._raw.linearintdata)
    fi.setWindow(0.0, 0.0, 1)
    fi.adjust(color_scheme=color_scheme)
    fi.prep_string()
    return fi

  def update_settings (self, **kwds) :
    self._wx_img = None
    self.update_image(**kwds)

  def update_image (self, brightness=100, color_scheme=0) :
    """
    Re-process the image to adjust brightness and colors, and generate a new
    wx.Image object and corresponding thumbnail image.
    """
    import wx
    self._img = self.create_flex_image(
      brightness=brightness,
      color_scheme=color_scheme)
    w = self._img.ex_size2()
    h = self._img.ex_size1()
    self.set_image_size(w, h)
    wx_image = wx.EmptyImage(w, h)
    wx_image.SetData(self._img.export_string)
    self._wx_img = wx_image
    binning = 8
    if (w > 2560) :
      binning = 16
    fi_thumb = self.create_flex_image(brightness=brightness,
      color_scheme=color_scheme,
      binning=binning)
    w = fi_thumb.ex_size2()
    h = fi_thumb.ex_size1()
    wx_thumb = wx.EmptyImage(w, h)
    wx_thumb.SetData(fi_thumb.export_string)
    self.set_thumbnail_size(w, h, binning)
    self._wx_thumb = wx_thumb
    self._wx_thumb_bmp = wx_thumb.ConvertToBitmap()

  def get_bitmap (self) :
    """
    Returns the primary wx.Image scaled and clipped to the current screen
    parameters for display in the main canvas.
    """
    import wx
    x, y, w, h = self.get_bitmap_params()
    scale = self.get_scale()
    img = self._wx_img.GetSubImage((x, y, w, h))
    img = img.Scale(w * scale, h * scale, wx.IMAGE_QUALITY_NORMAL)
    return img.ConvertToBitmap()

  def get_thumbnail_bitmap (self) :
    """
    Returns the thumbnail image (without any further processing).
    """
    return self._wx_thumb_bmp #.ConvertToBitmap()

  def get_zoomed_bitmap (self, x, y, boxsize=400, mag=16) :
    """
    Returns a zoomed-in view of the image, centered around the clicked
    position.
    """
    import wx
    x0, y0, w, h = self.get_zoom_box(x, y, boxsize, mag)
    assert (w == h)
    img = self._wx_img.GetSubImage((x0, y0, w, h))
    return img.Scale(boxsize, boxsize, wx.IMAGE_QUALITY_NORMAL)

  # XXX does this need to be in C++?
  def get_intensities_in_box (self, x, y, boxsize=400, mag=16) :
    x0, y0, w, h = self.get_zoom_box(x, y, boxsize, mag)
    i, j = self.image_coords_as_array_coords(x0, y0)
    d = self._raw.linearintdata
    format = " ".join([ "%d" for n in range(w) ])
    values = []
    for u in range(1, h+1) : # XXX why can't I start at 0?
      values.append([])
      for v in range(1, w+1) :
        intensity = d[i+u, j+v]
        values[u-1].append(intensity)
    #for row in values :
      #print format % tuple(row)
    return values

  def get_drawable_spots (self) :
    """
    Given an array of spotfinder results (generated separately), determine
    which of these are within the current bounds of the viewport.
    """
    if (self._spots is None) : return []
    x, y, w, h = self.get_bitmap_params()
    all_spots = []
    for spot in self._spots :
      all_spots.append(( spot.ctr_mass_x(), spot.ctr_mass_y() ))
    spots_out = self._get_drawable_points(all_spots)
    return spots_out

  def _get_drawable_points (self, points) :
    points_out = []
    x, y, w, h = self.get_bitmap_params()
    for ym,xm in points :
      if ((x+w) >= xm >= x) and ((y+h) >= ym >= y) :
        xm_, ym_ = self.image_coords_as_screen_coords(xm, ym)
        points_out.append((xm_, ym_))
    return points_out

  def get_drawable_background_mask (self) :
    if (self._integration is None) : return []
    points_out = self._get_drawable_points(
      self._integration['background_masks_xy'])
    return points_out

  def get_drawable_predictions (self) :
    if (self._integration is None) : return []
    points_out = self._get_drawable_points(
      self._integration['mapped_predictions'])
    return points_out

  def get_drawable_integration_mask (self) :
    if (self._integration is None) : return []
    points_out = self._get_drawable_points(
      self._integration['integration_masks_xy'])
    return points_out

  def set_beam_center_from_screen_coords (self, x, y) :
    """
    Reposition the beam center for the current image - this is not saved, but
    it will override the beam center in the image header.  Arguments are
    screen pixel coordinates in the main viewport.
    """
    x_im, y_im = self.screen_coords_as_image_coords(x, y)
    if ((x_im <= 0) or (y_im <= 0) or
        (x_im > self.img_w) or (y_im > self.img_h)) :
      raise Sorry("Coordinates are out of image!")
    x_point, y_point = self._raw.image_coords_as_detector_coords(x_im, y_im)
    old_x, old_y = self.get_beam_center_mm()
    self._beam_center = (x_point, y_point)
    return (old_x, old_y, x_point, y_point)

  def reset_beam_center (self) :
    self._beam_center = None

  def get_beam_center_mm (self) :
    if (self._beam_center is not None) :
      center_x, center_y = self._beam_center
    else:
      center_x, center_y = self._raw.get_beam_center_mm()
    return center_x, center_y

  def get_beam_center (self) :
    center_x, center_y = self.get_beam_center_mm()
    return self._raw.detector_coords_as_image_coords(center_x, center_y)

  def get_detector_distance (self) :
    dist = self._raw.distance
    twotheta = self.get_detector_2theta()
    if (twotheta == 0.0) :
      return dist
    else :
      return dist / math.cos(twotheta)

  def get_detector_2theta (self) :
    two_theta = self._raw.twotheta
    return two_theta * pi_over_180

  def get_wavelength (self) :
    return self._raw.wavelength

  def get_point_info (self, x, y) :
    """
    Determine the intensity, resolution, and array indices of a pixel.
    Arguments are in image pixel coordinates (starting from 1,1).
    """
    from spotfinder import core_toolbox
    x_point, y_point = self._raw.image_coords_as_detector_coords(x, y)
    x0, y0 = self._raw.detector_coords_as_image_coords(x_point, y_point)
    center_x, center_y = self.get_beam_center_mm()
    dist = self.get_detector_distance()
    two_theta = self.get_detector_2theta()
    wavelength = self.get_wavelength()
    """
    calc = core_toolbox.resolution_on_image(
                xbeam_mm = center_x,
                ybeam_mm = center_y,
                distance_mm = self._raw.distance,
                wavelength_ang = wavelength,
                twotheta_rad = two_theta)
    d_min = calc.resolution_at_point(x_point, y_point)
    """ # future generalization

    if (dist > 0) :
      scattering_angle = rstbx.utils.get_scattering_angle(
        x=x_point,
        y=y_point,
        center_x=center_x,
        center_y=center_y,
        distance=dist,
        detector_two_theta=two_theta,
        distance_is_corrected=True)
      if (scattering_angle == 0.0) :
        d_min = None
      else :
        d_min = wavelength / (2 * math.sin(scattering_angle / 2))
    else:
      d_min = None
    slow, fast = self.image_coords_as_array_coords(x, y)
    try :
      intensity = self._raw.linearintdata[slow, fast]
    except IndexError :
      return None
    else :
      return point_info(slow, fast, intensity, d_min)

  def line_between_points (self, x1, y1, x2, y2, n_values=100) :
    """
    Given two points on the image, sample intensities along a line connecting
    them (using linear interpolation).  This also calculates the coordinates
    of each sample point, which is used for lattice dimension calculations
    once peaks have been identified.  Arguments are in image pixel coordinates
    (starting at 1,1).
    """
    x1_, y1_ = self.image_coords_as_array_coords(x1, y1)
    x2_, y2_ = self.image_coords_as_array_coords(x2, y2)
    n_values = ifloor(math.sqrt((x2_-x1_)**2 + (y2_-y1_)**2))
    delta_x = (x2_ - x1_) / (n_values - 1)
    delta_y = (y2_ - y1_) / (n_values - 1)
    vals = []
    img_coords = []
    d = self._raw.linearintdata
    # TODO remarkably, this is reasonably fast in Python, but it would
    # probably be more at home in scitbx.math
    for n in range(n_values) :
      x = x1_ + (n * delta_x)
      y = y1_ + (n * delta_y)
      xd, yd = self.array_coords_as_detector_coords(x, y)
      img_coords.append((xd,yd))
      x_1 = ifloor(x)
      x_2 = iceil(x)
      y_1 = ifloor(y)
      y_2 = iceil(y)
      v11 = d[(x_1, y_1)]
      v12 = d[(x_1, y_2)]
      v21 = d[(x_2, y_1)]
      v22 = d[(x_2, y_2)]
      if (x_2 == x_1) :
        if (y_2 == y_1) :
          vxy = v11
        else :
          vxy = ((v12 * (y - y_1)) + (v11 * (y_2 - y))) / (y_2 - y_1)
      elif (y_2 == y_1) :
        vxy =  ((v21 * (x - x_1)) + (v11 * (x_2 - x))) / (x_2 - x_1)
      else :
        dxdy = (y_2 - y_1) * (x_2 - x_1)
        vxy = ((v11 / dxdy) * (x_2 - x) * (y_2 - y)) + \
              ((v21 / dxdy) * (x - x_1) * (y_2 - y)) + \
              ((v12 / dxdy) * (x_2 - x) * (y - y_1)) + \
              ((v22 / dxdy) * (x - x_1) * (y - y_1))
      vals.append(vxy)
    lattice_length = None
    if (len(vals) > 5) :
      # first find peaks in the profile
      peaks = []
      avg = sum(vals) / len(vals)
      filtered_vals = []
      for x in vals :
        if (x <= avg*3) :
          filtered_vals.append(x)
      background = sum(filtered_vals) / len(filtered_vals)
      i = 2
      while (i < len(vals) - 2) :
        x = vals[i]
        if (x <= background) :
          pass
        elif ((x > vals[i-1]) and (x > vals[i-2]) and
              (x > vals[i+1]) and (x > vals[i+2])) :
          peaks.append(i)
        i += 1
      if (len(peaks) > 0) :
        # calculate the average lattice length
        center_x, center_y = self.get_beam_center_mm()
        distances = []
        i = 1
        while (i < len(peaks)) :
          x1,y1 = img_coords[peaks[i-1]]
          x2,y2 = img_coords[peaks[i]]
          rs_distance = rstbx.utils.reciprocal_space_distance(x1, y1, x2, y2,
            wavelength=self.get_wavelength(),
            center_x=center_x,
            center_y=center_y,
            distance=self.get_detector_distance(),
            detector_two_theta=self.get_detector_2theta(),
            distance_is_corrected=True)
          assert (rs_distance > 0)
          distances.append(1 / rs_distance)
          i += 1
        lattice_length = sum(distances) / len(distances)
    distance = self.distance_between_points(x1, y1, x2, y2)
    return line_profile(vals, distance, lattice_length)

class point_info (object) :
  """
  Container for storing attributes of a pixel, for display by the main
  frame (currently on the statusbar).
  """
  def __init__ (self, slow, fast, intensity, d_min) :
    self.slow = slow
    self.fast = fast
    self.intensity = intensity
    self.d_min = d_min

  def format (self) :
    return "resolution = %s  intensity = %.2f  slow=%d  fast=%d" % (
      format_value("%.2f A", self.d_min), self.intensity, self.slow, self.fast)

class line_profile (object) :
  def __init__ (self, values, distance, lattice_length) :
    self.values = values
    self.distance = distance
    self.lattice_length = lattice_length

# TODO replace this with libtbx.phil
class settings (object) :
  def __init__ (self) :
    self.zoom_level = 0
    self.brightness = 100
    self.show_beam_center = True
    self.invert_beam_center_axes = False
    self.show_spotfinder_spots = True
    self.show_integration = True
    self.enable_collect_values = True
    self.color_scheme = 0
