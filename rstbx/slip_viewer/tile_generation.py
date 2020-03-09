# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$

from __future__ import absolute_import, division, print_function
from six.moves import range

import math
import wx

######
# Base class for a tile object - handles access to tiles.
######
'''Jan '12 ... make it so the rstbx widget actually zooms in DONE
           ... draggable area large enough for whole image DONE
           ... remove the dependency on map tile library geological DONE
           ... provide cache  DONE
           ... put in negative zoom levels  DONE
           ... tighten up the flex image  -- refactor  DONE
           ... connect the open filename command line option to the tile generator DONE
           ... check in the code so far.  Put most of it in cxi xdr xes DONE
           ... provide cache invalidation DONE
           ... give mechanism to change brightness; color scheme
'''

def _get_flex_image(
    data, vendortype, binning=1, brightness=1.0, saturation=65535.0,
    show_untrusted=False, color_scheme=0):
  # This is a combination of the get_data_type() and get_flex_image()
  # functions from iotbx.detectors.detectorbase.  XXX This may turn
  # out to be generally useful (see
  # e.g. rstbx.viewer.create_flex_image()), but where to place it?
  # dxtbx Format class?
  typehash = str(data.__class__)
  if typehash.find('int') >= 0:
    from iotbx.detectors import FlexImage
  elif typehash.find('double') >= 0:
    from iotbx.detectors import FlexImage_d as FlexImage

  return FlexImage(
    binning=binning,
    brightness=brightness,
    rawdata=data,
    saturation=int(round(saturation)),
    vendortype=vendortype,
    show_untrusted=show_untrusted,
    color_scheme=color_scheme
  )


def _get_flex_image_multipanel(panels, raw_data, beam, brightness=1.0,
                               binning=1, show_untrusted=False, color_scheme=0):
  # From xfel.cftbx.cspad_detector.readHeader() and
  # xfel.cftbx.cspad_detector.get_flex_image().  XXX Is it possible to
  # merge this with _get_flex_image() above?  XXX Move to dxtbx Format
  # class (or a superclass for multipanel images)?

  from math import ceil

  from iotbx.detectors import generic_flex_image
  from libtbx.test_utils import approx_equal
  from scitbx.array_family import flex
  from scitbx.matrix import col, rec, sqr
  from xfel.cftbx.detector.metrology import get_projection_matrix

  assert len(panels) == len(raw_data), (len(panels), len(raw_data))

  # Determine next multiple of eight of the largest panel size.
  data_max_focus = None
  for data in raw_data:
      if data_max_focus is None:
          data_max_focus = data.focus()
      else:
          data_max_focus = (
              max(data_max_focus[0], data.focus()[0]),
              max(data_max_focus[1], data.focus()[1]),
          )
  data_padded = (
      8 * int(ceil(data_max_focus[0] / 8)),
      8 * int(ceil(data_max_focus[1] / 8)),
  )

  # Assert that all saturated values are equal and not None.  While
  # dxtbx records a separated trusted_range for each panel,
  # generic_flex_image supports only accepts a single common value for
  # the saturation.
  saturation = None
  for panel in panels:
      if saturation is None:
          saturation = panel.get_trusted_range()[1]
      else:
          assert approx_equal(saturation, panel.get_trusted_range()[1])
  assert saturation is not None

  # Create rawdata and my_flex_image before populating it.
  rawdata = flex.double(
    flex.grid(len(panels) * data_padded[0], data_padded[1]))
  my_flex_image = generic_flex_image(
    rawdata=rawdata,
    binning=binning,
    size1_readout=data_max_focus[0],
    size2_readout=data_max_focus[1],
    brightness=brightness,
    saturation=saturation,
    show_untrusted=show_untrusted,
    color_scheme=color_scheme
  )

  # Calculate the average beam center across all panels, in meters
  # not sure this makes sense for detector which is not on a plane?
  beam_center = col((0, 0, 0))
  npanels = 0
  for panel in panels:
    try:
      beam_center += col(panel.get_beam_centre_lab(beam.get_s0()))
      npanels += 1
    except RuntimeError as e: # catch DXTBX_ASSERT for no intersection
      pass
  beam_center /= (npanels / 1e-3)

  # XXX If a point is contained in two panels simultaneously, it will
  # be assigned to the panel defined first.  XXX Use a Z-buffer
  # instead?
  for i in range(len(panels)):
    # Determine the pixel size for the panel (in meters), as pixel
    # sizes need not be identical.
    data = raw_data[i]
    panel = panels[i]
    pixel_size = (panel.get_pixel_size()[0] * 1e-3,
                  panel.get_pixel_size()[1] * 1e-3)

    if len(panels) == 24 and panels[0].get_image_size() == (2463,195):
      rawdata.matrix_paste_block_in_place(
        block=data.as_double(),
        i_row=i * data_padded[0],
        i_column=0)
      # XXX hardcoded panel height and row gap
      my_flex_image.add_transformation_and_translation((1,0,0,1),
                                                       (-i*(195+17),0))

      continue

    elif len(panels) == 120 and panels[0].get_image_size() == (487,195):
      i_row = i // 5
      i_col = i % 5
      rawdata.matrix_paste_block_in_place(
        block=data.as_double(),
        i_row=i * data_padded[0],
        i_column=0)
      # XXX hardcoded panel height and row gap
      my_flex_image.add_transformation_and_translation(
        (1,0,0,1), (-i_row*(195+17),-i_col*(487+7)))

      continue

    # Get unit vectors in the fast and slow directions, as well as the
    # the locations of the origin and the center of the panel, in
    # meters. The origin is taken w.r.t. to average beam center of all
    # panels. This avoids excessive translations that can result from
    # rotations around the laboratory origin. Related to beam centre above
    # and dials#380 not sure this is right for detectors which are not
    # coplanar since system derived from first panel...
    fast = col(panel.get_fast_axis())
    slow = col(panel.get_slow_axis())
    origin = col(panel.get_origin()) * 1e-3  - beam_center

    center = origin \
             + (data.focus()[0] - 1) / 2 * pixel_size[1] * slow \
             + (data.focus()[1] - 1) / 2 * pixel_size[0] * fast
    normal = slow.cross(fast).normalize()

    # Determine rotational and translational components of the
    # homogeneous transformation that maps the readout indices to the
    # three-dimensional laboratory frame.
    Rf = sqr((  fast(0, 0),   fast(1, 0),   fast(2, 0),
               -slow(0, 0),  -slow(1, 0),  -slow(2, 0),
              normal(0, 0), normal(1, 0), normal(2, 0)))
    tf = -Rf * center
    Tf = sqr((Rf(0, 0), Rf(0, 1), Rf(0, 2), tf(0, 0),
              Rf(1, 0), Rf(1, 1), Rf(1, 2), tf(1, 0),
              Rf(2, 0), Rf(2, 1), Rf(2, 2), tf(2, 0),
              0,        0,        0,        1))

    # E maps picture coordinates onto metric Cartesian coordinates,
    # i.e. [row, column, 1 ] -> [x, y, z, 1].  Both frames share the
    # same origin, but the first coordinate of the screen coordinate
    # system increases downwards, while the second increases towards
    # the right.  XXX Is this orthographic projection the only one
    # that makes any sense?
    E = rec(elems=[0, +pixel_size[1], 0,
                   -pixel_size[0], 0, 0,
                   0, 0, 0,
                   0, 0, 1],
            n=[4, 3])

    # P: [x, y, z, 1] -> [row, column, 1].  Note that data.focus()
    # needs to be flipped to give (horizontal, vertical) size,
    # i.e. (width, height).
    Pf = get_projection_matrix(
      pixel_size, (data.focus()[1], data.focus()[0]))[0]

    rawdata.matrix_paste_block_in_place(
      block=data.as_double(),
      i_row=i * data_padded[0],
      i_column=0)

    # Last row of T is always [0, 0, 0, 1].
    T = Pf * Tf * E
    R = sqr((T(0, 0), T(0, 1),
             T(1, 0), T(1, 1)))
    t = col((T(0, 2), T(1, 2)))
    my_flex_image.add_transformation_and_translation(R, t)
  my_flex_image.followup_brightness_scale()
  return my_flex_image


class _Tiles(object):
    # maximum number of tiles held in each level cache
    MaxTileList = 512

    def __init__(self, filename):
        (self.tile_size_x, self.tile_size_y) = (256,256)
        self.levels = [-3,-2,-1,0,1,2,3,4,5]

        # set min and max tile levels
        self.min_level = -3
        self.max_level = 5
        self.extent = (-180.0, 180., -166.66 , 166.66) #longitude & latitude limits
        self.set_image(filename)
        self.current_brightness = 1.0
        self.current_color_scheme = 0
        self.user_requests_antialiasing = False

        self.show_untrusted = False

    def set_image(self, file_name_or_data, metrology_matrices=None, get_raw_data=None):

        self.reset_the_cache()
        if file_name_or_data is None:
          self.raw_image = None
          return
        if type(file_name_or_data) is type(""):
          from iotbx.detectors import ImageFactory
          self.raw_image = ImageFactory(file_name_or_data)
          self.raw_image.read()
        else:
          try:
            self.raw_image = file_name_or_data._raw
          except AttributeError:
            self.raw_image = file_name_or_data
        #print "SETTING NEW IMAGE",self.raw_image.filename

        # XXX Since there doesn't seem to be a good way to refresh the
        # image (yet), the metrology has to be applied here, and not
        # in frame.py.

        detector = self.raw_image.get_detector()

        if len(detector) > 1 and metrology_matrices is not None:
          self.raw_image.apply_metrology_from_matrices(metrology_matrices)

        if get_raw_data is not None:
          self.raw_image.set_raw_data(get_raw_data(self.raw_image))
        raw_data = self.raw_image.get_raw_data()
        if not isinstance(raw_data, tuple):
          raw_data = (raw_data,)

        if len(detector) > 1:
          self.flex_image = _get_flex_image_multipanel(
            brightness=self.current_brightness / 100,
            panels=detector,
            show_untrusted=self.show_untrusted,
            raw_data=raw_data,
            beam = self.raw_image.get_beam(),
            color_scheme=self.current_color_scheme)
        else:
          self.flex_image = _get_flex_image(
            brightness=self.current_brightness / 100,
            data=raw_data[0],
            saturation=self.raw_image.get_detector()[0].get_trusted_range()[1],
            vendortype=self.raw_image.get_vendortype(),
            show_untrusted=self.show_untrusted,
            color_scheme=self.current_color_scheme
          )

        if self.zoom_level >= 0:
          self.flex_image.adjust(color_scheme=self.current_color_scheme)

    def set_image_data(self, raw_image_data):
      self.reset_the_cache()
      # XXX Since there doesn't seem to be a good way to refresh the
      # image (yet), the metrology has to be applied here, and not
      # in frame.py.

      detector = self.raw_image.get_detector()
      self.raw_image.set_raw_data(raw_image_data)
      if len(detector) == 1 and len(raw_image_data) == 1:
        raw_image_data = raw_image_data[0]

      if len(detector) > 1:
        self.flex_image = _get_flex_image_multipanel(
          brightness=self.current_brightness / 100,
          panels=detector,
          raw_data=raw_image_data,
          beam = self.raw_image.get_beam())
      else:
        self.flex_image = _get_flex_image(
          brightness=self.current_brightness / 100,
          data=raw_image_data,
          saturation=self.raw_image.get_detector()[0].get_trusted_range()[1],
          vendortype=self.raw_image.get_vendortype(),
          show_untrusted=self.show_untrusted
        )

      self.flex_image.adjust(color_scheme=self.current_color_scheme)

    def update_brightness(self,b,color_scheme=0):
        raw_data = self.raw_image.get_raw_data()
        if not isinstance(raw_data, tuple):
          raw_data = (raw_data,)

        if len(self.raw_image.get_detector()) > 1:
          # XXX Special-case read of new-style images until multitile
          # images are fully supported in dxtbx.
          self.flex_image = _get_flex_image_multipanel(
            brightness=b / 100,
            panels=self.raw_image.get_detector(),
            show_untrusted=self.show_untrusted,
            raw_data=raw_data,
            beam=self.raw_image.get_beam(),
            color_scheme=color_scheme)
        else:
          self.flex_image = _get_flex_image(
            brightness=b / 100,
            data=raw_data[0],
            saturation=self.raw_image.get_detector()[0].get_trusted_range()[1],
            vendortype=self.raw_image.get_vendortype(),
            show_untrusted=self.show_untrusted,
            color_scheme=color_scheme
          )

        self.reset_the_cache()
        self.UseLevel(self.zoom_level)
        self.current_color_scheme = color_scheme
        self.current_brightness = b
        self.flex_image.adjust(color_scheme)

    def update_color_scheme(self,color_scheme=0):
        self.flex_image.adjust(color_scheme)
        self.reset_the_cache()
        self.UseLevel(self.zoom_level)
        self.current_color_scheme = color_scheme

    def reset_the_cache(self):

        # setup the tile caches and Least Recently Used lists
        self.cache = {}
        self.lru = {}
        for l in self.levels:
            self.cache[l] = {}
            self.lru[l] = []

    def flex_image_get_tile(self,x,y):
      # The supports_rotated_tiles_antialiasing_recommended flag in
      # the C++ FlexImage class indicates whether the underlying image
      # instance supports tilted readouts.  Anti-aliasing only makes
      # sense if it does.
      if self.raw_image is not None and \
        self.zoom_level >=2 and \
        self.flex_image.supports_rotated_tiles_antialiasing_recommended and \
        self.user_requests_antialiasing:
        # much more computationally intensive to prepare nice-looking pictures of tilted readout
        self.flex_image.setZoom(self.zoom_level+1)
        fraction = 512./self.flex_image.size1()/(2**(self.zoom_level+1))
        self.flex_image.setWindowCart(  y, x, fraction )
        self.flex_image.prep_string()
        w,h = self.flex_image.ex_size2(), self.flex_image.ex_size1()
        assert w==512
        assert h==512
        wx_image = wx.EmptyImage(w/2,h/2)
        try:
          import PIL.Image as Image
        except ImportError:
          import Image
        try:
          I = Image.fromstring("RGB",(512,512),self.flex_image.export_string)
        except NotImplementedError:
          I = Image.frombytes("RGB",(512,512),self.flex_image.export_string)
        J = I.resize((256,256),Image.ANTIALIAS)
        wx_image.SetData(J.tostring())
        return wx_image.ConvertToBitmap()
      elif self.raw_image is not None:
        self.flex_image.setZoom(self.zoom_level)
        fraction = 256./self.flex_image.size1()/(2**self.zoom_level)
        self.flex_image.setWindowCart(  y, x, fraction )
        self.flex_image.prep_string()
        w,h = self.flex_image.ex_size2(), self.flex_image.ex_size1()
        assert w==256
        assert h==256
        wx_image = wx.EmptyImage(w,h)
        wx_image.SetData(self.flex_image.export_string)
        return wx_image.ConvertToBitmap()
      else:
        wx_image = wx.EmptyImage(256,256)
        return wx_image.ConvertToBitmap()

    def get_binning(self):
      if self.zoom_level>=0: return 1.
      return 2.**-self.zoom_level

    def UseLevel(self, n):
        """Prepare to serve tiles from the required level.

        n  The required level

        Returns a tuple (map_width, map_height, ppd_x, ppd_y) if successful,
        else None.  The width/height values are pixels.  The ppd_? values are
        pixels-per-degree values for X and Y direction.
        """
        # try to get cache for this level, no cache means no level
        #print "IN USE LEVEL",n
        try:
            self.tile_cache = self.cache[n]
            self.tile_list = self.lru[n]
        except KeyError:
            return None
        self.zoom_level = n
        if self.raw_image is None: #dummy values if there is no image
          self.center_x_lon = self.center_y_lat = 500.
          return (1024,1024,1.,1.)
        self.num_tiles_x = int(math.ceil((self.flex_image.size1()*(2**self.zoom_level))/256.))
        self.num_tiles_y = int(math.ceil((self.flex_image.size2()*(2**self.zoom_level))/256.))
        self.ppd_x = 2.**self.zoom_level
        self.ppd_y = 2.**self.zoom_level
        #print "USELEVEL %d # tiles: %d %d"%(n,self.num_tiles_x,self.num_tiles_y)
        #print "USELEVEL %d returning"%n,(self.tile_size_x * self.num_tiles_x,
        #        self.tile_size_y * self.num_tiles_y,
        #        self.ppd_x, self.ppd_y)
        # The longitude & latitude coordinates at the image center:
        self.center_x_lon = self.extent[0] + (1./self.ppd_x) * (0 +
          self.flex_image.size2() * (2**self.zoom_level) / 2.
          )
        self.center_y_lat = self.extent[3] - (1./self.ppd_y) * (0 +
          self.flex_image.size1() * (2**self.zoom_level) / 2.
          )
        # The 2+num_tiles is just a trick to get PySlip to think the map is
        # slightly larger, allowing zoom level -3 to be properly framed:
        # ....for larger display sizes it is necessary to increase this...
        # ....can tile_generator get the display size & figure it out?
        return (self.tile_size_x * (2+self.num_tiles_x),
                self.tile_size_y * (2+self.num_tiles_y),
                self.ppd_x, self.ppd_y)

    def get_initial_instrument_centering_within_picture_as_lon_lat(self):
      import sys
      detector = self.raw_image.get_detector()
      if sys.platform.lower().find("linux") >= 0:
        if len(detector) > 1:
          return 0.,0.
        else:
          return self.center_x_lon-self.extent[0], self.center_y_lat-self.extent[3]
      else:
        if len(detector) > 1:
          return self.extent[0], self.extent[3]
        else:
          return self.center_x_lon, self.center_y_lat

    def GetTile(self, x, y):
        #from libtbx.development.timers import Timer
        #T = Timer("get tile")
        """Get bitmap for tile at tile coords (x, y).

        x  X coord of tile required (tile coordinates)
        y  Y coord of tile required (tile coordinates)

        Returns bitmap object containing the tile image.
        Tile coordinates are measured from map top-left.
        """
        try:
            # if tile in cache, return it from there
            pic = self.tile_cache[(x, y)]
            index = self.tile_list.index((x, y))
            del self.tile_list[index]
        except KeyError:
            pic = self.flex_image_get_tile(x,y)
            self.tile_cache[(x, y)] = pic

        self.tile_list.insert(0, (x, y))
        if len(self.tile_cache)>=self.MaxTileList:
          del self.tile_cache[self.tile_list[-1]]
          del self.tile_list[-1]
        return pic

    def get_flex_pixel_coordinates(self, lon, lat):
      fast_picture_coord_pixel_scale, slow_picture_coord_pixel_scale = \
        self.lon_lat_to_picture_fast_slow(lon,lat)
      if self.flex_image.supports_rotated_tiles_antialiasing_recommended: # for generic_flex_image
        tilted = self.flex_image.picture_to_readout(
          slow_picture_coord_pixel_scale,fast_picture_coord_pixel_scale)
        return tilted
      else: # standard flex_image
        return slow_picture_coord_pixel_scale,fast_picture_coord_pixel_scale

    def lon_lat_to_picture_fast_slow(self,longitude,latitude):
      # input latitude and longitude in degrees (conceptually)
      # output fast and slow picture coordinates in units of detector pixels
      # slow is pointing down (x).  fast is pointing right (y).

      detector = self.raw_image.get_detector()
      if len(detector) == 1:
        (size2, size1) = detector[0].get_image_size()
      else:
        # XXX Special-case until multitile detectors fully supported.
        (size1, size2) = (self.flex_image.size1(), self.flex_image.size2())

      return \
        (size2/2.) - (self.center_x_lon - longitude),  \
        (size1/2.) - (latitude - self.center_y_lat)

    def picture_fast_slow_to_lon_lat(self,pic_fast_pixel,pic_slow_pixel):
      # inverse of the preceding function

      detector = self.raw_image.get_detector()
      if detector.num_panels() == 1:
        (size1, size2) = detector.get_image_size()
      else:
        # XXX Special-case until multitile detectors fully supported.
        (size1, size2) = (self.flex_image.size1(), self.flex_image.size2())

      return \
        (size2/2.) - self.center_x_lon - pic_fast_pixel, \
        (size1/2.) + self.center_y_lat - pic_slow_pixel

    def picture_fast_slow_to_map_relative(self,pic_fast_pixel,pic_slow_pixel):
      #return up/down, left/right map relative coords for pyslip layers
      return pic_fast_pixel+self.extent[0],-pic_slow_pixel+self.extent[3]

    def map_relative_to_picture_fast_slow(self, map_rel_vert, map_rel_horiz):
      # return fast, slow picture coords
      return map_rel_vert-self.extent[0],-map_rel_horiz+self.extent[3]

    def vec_picture_fast_slow_to_map_relative(self,vector):
      value = []
      for vec in vector:
        value.append(self.picture_fast_slow_to_map_relative(vec[0],vec[1]))
      return value

    def get_spotfinder_data(self, params):

      pointdata = []
      test_pattern = False
      if test_pattern is True and self.raw_image.__class__.__name__.find("CSPadDetector") >= 0:
        key_count = -1
        for key, asic in self.raw_image._tiles.iteritems():
          key_count += 1
          focus = asic.focus()
          for slow in range(0,focus[0],20):
            for fast in range(0,focus[1],20):
              slowpic,fastpic = self.flex_image.tile_readout_to_picture(key_count,slow,fast)
              mr1,mr2 = self.picture_fast_slow_to_map_relative(fastpic,slowpic)
              pointdata.append((mr1,mr2,{"data":key}))

      elif (self.raw_image.__class__.__name__.find("CSPadDetector") >= 0):
        from cxi_xdr_xes.cftbx.spotfinder.speckfinder import spotfind_readout

        key_count = -1
        for key, asic in self.raw_image._tiles.iteritems():
          key_count += 1
          indexing = spotfind_readout(
              readout=asic,
              peripheral_margin=params.spotfinder.peripheral_margin)

          for spot in indexing:
            slow = int(round(spot[0]))
            fast = int(round(spot[1]))

            slowpic,fastpic = self.flex_image.tile_readout_to_picture(key_count,slow,fast)
            mr1,mr2 = self.picture_fast_slow_to_map_relative(fastpic,slowpic)
            pointdata.append((mr1,mr2,{"data":key}))

      else:
        from spotfinder.command_line.signal_strength import master_params
        working_params = master_params.fetch(sources = []) #placeholder for runtime mods
        working_params.show(expert_level=1)
        distl_params = working_params.extract()

        spotfinder,frameno = self.raw_image.get_spotfinder(distl_params)
        spots = spotfinder.images[frameno]["spots_total"]
        for spot in spots:
          mr = self.picture_fast_slow_to_map_relative(
              spot.max_pxl_y() + 0.5, spot.max_pxl_x() + 0.5)
#             spot.ctr_mass_y() + 0.5, spot.ctr_mass_x() + 0.5)
          pointdata.append(mr)
      return pointdata

    def get_effective_tiling_data(self, params):
      box_data = []
      text_data = []
      if hasattr(self.raw_image, 'get_tile_manager'):
        IT = self.raw_image.get_tile_manager(params).effective_tiling_as_flex_int()
        for i in range(len(IT) // 4):
          tile = IT[4*i:4*i+4]
          attributes = {'color': '#0000FFA0', 'width': 1, 'closed': False}
          box_data.append(
              ((self.picture_fast_slow_to_map_relative(tile[1], tile[0]),
                self.picture_fast_slow_to_map_relative(tile[1], tile[2])),
               attributes))
          box_data.append(
              ((self.picture_fast_slow_to_map_relative(tile[1], tile[0]),
                self.picture_fast_slow_to_map_relative(tile[3], tile[0])),
               attributes))
          box_data.append(
              ((self.picture_fast_slow_to_map_relative(tile[1], tile[2]),
                self.picture_fast_slow_to_map_relative(tile[3], tile[2])),
               attributes))
          box_data.append(
              ((self.picture_fast_slow_to_map_relative(tile[3], tile[0]),
                self.picture_fast_slow_to_map_relative(tile[3], tile[2])),
               attributes))
          txt_x, txt_y = self.picture_fast_slow_to_map_relative(
              (tile[1]+tile[3])//2, (tile[0]+tile[2])//2)
          text_data.append((txt_x, txt_y, "%i" %i))
      return box_data, text_data

    def get_resolution(self, x, y, readout=None):
        """
        Determine the resolution of a pixel.
        Arguments are in image pixel coordinates (starting from 1,1).
        """

        d_min = None
        detector = self.raw_image.get_detector()
        beam = self.raw_image.get_beam()
        if detector is None or beam is None:
          return None
        beam = beam.get_s0()

        if len(detector) > 1:
          if readout is None:
            return None

          panel = detector[readout]

        else:
          panel = detector[0]

        if abs(panel.get_distance()) > 0:
          return panel.get_resolution_at_pixel(beam, (x, y))
        else:
          return None

    def get_detector_distance(self):
        detector = self.raw_image.get_detector()
        if len(detector) == 1:
          dist = abs(detector[0].get_distance())
        else:
          # XXX Special-case until multitile detectors fully
          # supported.
          dist = self.raw_image.distance
        twotheta = self.get_detector_2theta()
        if (twotheta == 0.0):
            return dist
        else :
            return dist / math.cos(twotheta)

    def get_detector_2theta(self):
        from scitbx.matrix import col

        detector = self.raw_image.get_detector()
        if len(detector) == 1:
          n = col(detector[0].get_normal())
          s0 = col(self.raw_image.get_beam().get_unit_s0())
          two_theta = s0.angle(n, deg=False)
        else:
          # XXX Special-case until multitile detectors fully
          # supported.
          try:
            two_theta = self.raw_image.twotheta * math.pi / 180
          except AttributeError:
            two_theta = 0

        return two_theta
