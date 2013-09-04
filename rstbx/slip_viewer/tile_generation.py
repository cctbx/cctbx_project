# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$

from __future__ import division

import math
import wx
import rstbx.utils

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

    def set_image(self, file_name_or_data, metrology_matrices=None):

        self.reset_the_cache()
        if file_name_or_data is None:
          self.raw_image = None
          return
        if type(file_name_or_data) is type(""):
          from iotbx.detectors import ImageFactory
          self.raw_image = ImageFactory(file_name_or_data)
        else:
          self.raw_image = file_name_or_data._raw
        self.raw_image.read()
        #print "SETTING NEW IMAGE",self.raw_image.filename

        # XXX Since there doesn't seem to be a good way to refresh the
        # image (yet), the metrology has to be applied here, and not
        # in frame.py.

        if (self.raw_image.implements_metrology_matrices() and
            metrology_matrices is not None):
            self.raw_image.apply_metrology_from_matrices(metrology_matrices)

        self.flex_image = self.raw_image.get_flex_image(
            brightness=self.current_brightness / 100)
        self.flex_image.adjust(color_scheme=self.current_color_scheme)

    def update_brightness(self,b,color_scheme=0):
        self.flex_image= self.raw_image.get_flex_image(brightness=b/100.)
        self.update_color_scheme(color_scheme)
        self.current_brightness = b

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
        import Image
        I = Image.fromstring("RGB",(512,512),self.flex_image.export_string)
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

        Returns a tuple (map_width, map_height, ppd_x, ppd_y) if succesful,
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
      if sys.platform.lower().find("linux") >= 0:
        if self.raw_image.__class__.__name__.find('FormatPYmultitile') >= 0:
          return 0.,0.
        else:
          return self.center_x_lon-self.extent[0], self.center_y_lat-self.extent[3]
      else:
        if self.raw_image.__class__.__name__.find('FormatPYmultitile') >= 0:
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
      return \
        (float(self.raw_image.size2)/2.) - \
        (self.center_x_lon - float(longitude)),  \
        (float(self.raw_image.size1)/2.) - \
        (float(latitude) - self.center_y_lat)

    def picture_fast_slow_to_lon_lat(self,pic_fast_pixel,pic_slow_pixel):
      # inverse of the preceding function
      return \
        (float(self.raw_image.size2)/2.) - self.center_x_lon - pic_fast_pixel, \
        (float(self.raw_image.size1)/2.) + self.center_y_lat - pic_slow_pixel

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
          for slow in xrange(0,focus[0],20):
            for fast in xrange(0,focus[1],20):
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
        for i in xrange(len(IT) // 4):
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

    def get_resolution (self, x, y, readout=None) :
        """
        Determine the resolution of a pixel.
        Arguments are in image pixel coordinates (starting from 1,1).
        """
        x_point, y_point = self.raw_image.image_coords_as_detector_coords(x, y,readout)
        x0, y0 = self.raw_image.detector_coords_as_image_coords(x_point, y_point)
        center_x, center_y = self.raw_image.get_beam_center_mm()
        dist = self.get_detector_distance()
        two_theta = self.get_detector_2theta()
        wavelength = self.raw_image.wavelength

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

        return d_min

    def get_detector_distance (self) :
        dist = self.raw_image.distance
        twotheta = self.get_detector_2theta()
        if (twotheta == 0.0) :
            return dist
        else :
            return dist / math.cos(twotheta)

    def get_detector_2theta (self) :
        try:
            two_theta = self.raw_image.twotheta
            return two_theta * (math.pi / 180)
        except AttributeError:
            return 0
