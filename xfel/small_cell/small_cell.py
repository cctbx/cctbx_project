from __future__ import division
#-*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

""" Series of functions used by cctbx.small_cell """

import numpy as np
from scipy.optimize import minimize
import math
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import itertools
from scitbx.matrix import col, sqr
from cctbx.crystal import symmetry
import cctbx.miller
from cctbx.miller import flex
from cctbx.uctbx import unit_cell
from cctbx import sgtbx
import operator

import wx
app = wx.App(0)
wx.SystemOptions.SetOptionInt("osx.openfiledialog.always-show-types", 1)

""" Calculates the euclidian distance between two 2D points """
def measure_distance (a,b): return math.sqrt((math.pow(b[0]-a[0],2)+math.pow(b[1]-a[1],2)))

def d_in_pixels (d_spacing, wavelength, distance, pixel_size):
  '''Calculate distance in pixels from beam center for a given d spacing using image's
     center, wavelength and detector distance'''
  return distance * math.tan(2 * math.asin(wavelength/(2*d_spacing))) / pixel_size

class small_cell_spot:
  """ Bucket for tracking parameters for a specific reflection """
  def __init__(self, spot_dict, ID):
    '''
    @p spot_dict: dictionary from a dials reflection table
    @p ID: identifier for the spot
    '''
    self.ID = ID
    self.spot_dict = spot_dict
    self.hkls = [] # used during ambiguous HKL resolution
    self.hkl = None # only set when an hkl combination is locked in


    self.x, self.y, self.z = self.spot_dict['xyzrecip']

    self.xyz = col([self.x,self.y,self.z])

    self.pred = None # used in calculating final RMSD

    self.peak_pixels = flex.vec3_int()
    l, r, t, b, z0, z1 = self.spot_dict['bbox']
    z = z1-z0
    for my, y in zip(xrange(b-t),xrange(t,b)):
      for mx, x in zip(xrange(r-l),xrange(l,r)):
        self.peak_pixels.append((x,y,z))


class small_cell_hkl(object):
  """ Class for storing an asymmetric unit hkl and an
  original hkl """
  def __init__(self, ahkl, ohkl):
    """
    @param ahkl Asymmetric unit hkl, aka the hkl with no symops applied
    @param ohkl Original hkl, aka the hkl with symops applied
    """
    self._ahkl = ahkl
    self._ohkl = ohkl
    self.flipped = False
    self.connections = []

  def __eq__(self, other):
    if hasattr(other, "ohkl") and self.ohkl == other.ohkl:
      return True
    return False

  def __getattr__(self, name):
    if name == "ohkl":
      if self.flipped:
        return -self._ohkl
      else:
        return self._ohkl
    elif name == "ahkl":
      if self.flipped:
        return -self._ahkl
      else:
        return self._ahkl
    else:
      raise AttributeError()

  def __setattr__(self, name, value):
    if name == "ohkl":
      self._ohkl = value
    elif name == "ahkl":
      self._ohkl = value
    else:
      object.__setattr__(self, name, value)

  def get_ohkl_str(self):
    """ Get a nicely formatted string for this hkl's original hkl """
    if self.ohkl is None:
      return None
    return "[% 2d,% 2d,% 2d]"%(self.ohkl.elems[0],self.ohkl.elems[1],self.ohkl.elems[2])

  def get_ahkl_str(self):
    """ Get a nicely formatted string for this hkl's asymmetric unit hkl """
    if self.ahkl is None:
      return None
    return "[% 2d,% 2d,% 2d]"%(self.ahkl.elems[0],self.ahkl.elems[1],self.ahkl.elems[2])

class small_cell_connection(object):
  """ Represents a connection between two reflections """
  def __init__(self, hkl1, hkl2, spot1, spot2, dobs, dcalc):
    """
    @param hkl1 small_cell_hkl object 1
    @param hkl2 small_cell_hkl object 2
    @param spot1 small_cell_spot 1
    @param spot2 small_cell_spot 2
    @param dobs observed reciprocal space distance between spot1 and spot2
    @param dcalc calculated reciprocal space distance between hkl1 and hkl2
    """
    self.hkl1 = hkl1
    self.hkl2 = hkl2
    self.spot1 = spot1
    self.spot2 = spot2
    self.dobs = dobs
    self.dcalc = dcalc

def test_spot_connection(hklA,hklB,xyzA,xyzB,metrical_matrix,phil):
  """ Given two hkls and two xyzs, determine if the reflections are 'connected', meaning
  the observed distance between the refelctions is similar to the calculated distance
  to the refelctions.
  See Brewster et. al. (2015), equation 4
  @param hklA small_cell_hkl object A
  @param hklB small_cell_hkl object B
  @param xyzA col object 1, observed xyz in reciprocal space
  @param xyzB col object 2, observed xyz in reciprocal space
  @param metrical_matrix metrical matrix of the unit cell (see equation 3)
  @param phil parsed small_cell phil parameters
  """
  dH = hklA.ohkl[0] - hklB.ohkl[0]
  dK = hklA.ohkl[1] - hklB.ohkl[1]
  dL = hklA.ohkl[2] - hklB.ohkl[2]
  column = col([dH,dK,dL])
  delta_calc = math.sqrt((column.transpose() * metrical_matrix * column)[0])
  delta_obv  = (xyzA-xyzB).length()

  return approx_equal(delta_calc, delta_obv, out=None, eps=phil.small_cell.spot_connection_epsilon), delta_obv, delta_calc

from rstbx.slip_viewer.frame import XrayFrame as SlipXrayFrame
class SmallCellXrayFrame(SlipXrayFrame):
  """ Extention of cctbx's image viewer for displaying small cell results """
  # inline method for drawing a circle.  From Johan's ring_frame.py in the slip viewer
  def _draw_ring_layer(self, dc, data, map_rel):
    """Draw a points layer.

    dc       the device context to draw on
    data     an iterable of point tuples:
             (x, y, place, radius, colour, x_off, y_off, pdata)
    map_rel  points relative to map if True, MUST BE TRUE for lightweight
    Assumes all points are the same colour, saving 100's of ms.
    """

    assert map_rel is True
    if len(data)==0:
      return
    (lon, lat, place, radius, colour, x_off, y_off, pdata) = data[0]

    scale = 2**self.pyslip.tiles.zoom_level

    # Draw points on map/view, using transparency if implemented.
    try:
      dc = wx.GCDC(dc)
    except NotImplementedError:
      pass
    dc.SetPen(wx.Pen(colour))
    dc.SetBrush(wx.Brush(colour, wx.TRANSPARENT))
    for (lon, lat, place, radius, colour, x_off, y_off, pdata) in data:
      (x, y) = self.pyslip.ConvertGeo2View((lon, lat))
      dc.DrawCircle(x, y, radius * scale)

class small_cell_callbacks:
  """ Worker class used by the wx app """
  def slip_callback(self,frame):
    """ Draw the predictions.  Code taken from slip_helpers.py """

    show_spotfinder_spots = True
    show_predictions = True
    show_connected_spots = True
    show_results = True
    show_rings = True
    show_text = True

    #BLUE: predictions
    blue_data = []; blue_text = []; blue_squares = []
    for pred, idx in zip(self.predicted, self.indicies):
      #if self.BSmasks[ix].keys()==[]:continue
      #x,y = frame.pyslip.tiles.picture_fast_slow_to_map_relative(
        #(pred[1]/self.pixel_size) +0.5,
        #(pred[0]/self.pixel_size) +0.5)
      #blue_data.append((x,y))
      x, y = frame.pyslip.tiles.picture_fast_slow_to_map_relative(pred[0],pred[1])
      blue_data.append((x,y))
      s1 = (x-5,y-5)
      s2 = (x-5,y+5)
      s3 = (x+5,y+5)
      s4 = (x+5,y-5)
      square = (s1,s2,s3,s4)
      blue_squares.append((square,{'closed':True}))
      blue_text.append((s2[0],s2[1],"%d %d %d"%(idx[0],idx[1],idx[2])))
    if show_predictions:
      if show_text:
        self.blue_text = frame.pyslip.AddTextLayer(blue_text, name="<blue_text>")
      self.blue_layer = frame.pyslip.AddPointLayer(
          blue_data, color="blue", name="<blue_layer>",
          radius=2,
          renderer = frame.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
      self.blue_squares = frame.pyslip.AddPolygonLayer(
          blue_squares, name = "<blue_square>",
          colour = "blue")

    if show_rings:
      #draw the d'rings
      sym = symmetry(unit_cell=self.ori.unit_cell(),space_group=self.horiz_phil.small_cell.spacegroup)
      hkl_list = cctbx.miller.build_set(sym, False, d_min=self.horiz_phil.small_cell.high_res_limit)

      #comment out the next four lines to see all the rings regardles of diffraction condition
      miller_array = hkl_list.indices().__class__()
      for idx in self.indicies:
        miller_array.append(idx)
      hkl_list = cctbx.miller.set(sym,miller_array)

      twotheta = hkl_list.two_theta(wavelength = frame._img._raw.wavelength)
      #L_mm = []
      L_pixels = []
      ring_indices = []
      bc = col((self.beam_x, self.beam_y))
      d = self.img.get_detector()[0]
      s0 = self.img.get_beam().get_s0()
      for tt in twotheta:
        #L_mm.append(self.distance * math.tan(tt[1]))
        lmm = self.distance * math.tan(tt[1])
        lpx = d.millimeter_to_pixel(col((lmm,0)) + bc)
        L_pixels.append(measure_distance(lpx,d.get_beam_centre_px(s0)))
        ring_indices.append(tt[0])
      #for lmm in L_mm:
        #L_pixels.append(lmm/frame._img._raw.pixel_size)

      bcx, bcy = frame.pyslip.tiles.picture_fast_slow_to_map_relative(
        (self.beam_x / frame._img._raw.pixel_size), (self.beam_y / frame._img._raw.pixel_size))
      #bcx, bcy = frame._img._raw.detector_coords_as_image_coords_float(self.beam_y, self.beam_x)
      #bcx, bcy = frame.pyslip.tiles.picture_fast_slow_to_map_relative(bcx + 0.5, bcy + 0.5)
      ring_text = []
      for pxl, i, idx in zip(L_pixels, xrange(len(L_pixels)), ring_indices):
        if pxl > 100: continue

        ring_data = [(bcx, bcy,
          {"color": "red", "radius": pxl})]

        frame.pyslip.AddPointLayer(
              ring_data,
              map_rel=True,
              visible=True,
              show_levels=[-3, -2, -1, 0, 1, 2, 3, 4, 5],
              renderer=frame._draw_ring_layer,
              name="<ring_layer%d>"%i)

        ring_text.append((bcx,bcy+pxl,"%d %d %d"%(idx[0],idx[1],idx[2])))

      if show_text:
        self.ring_text = frame.pyslip.AddTextLayer(ring_text, name="<ring_text>", colour = "red")

    red_data = []; springgreen_data = []; yellow_data = []#; pink_data = []
    green_data = []; green_text = []; green_squares = []
    purple_data = []; purple_squares = []; purple_text = []

    for spot in self.spots:
      # RED: spotfinder spot pixels
      for pxl in spot.peak_pixels:
        red_data.append(
          frame.pyslip.tiles.picture_fast_slow_to_map_relative(
            pxl[0] + 0.5, pxl[1] + 0.5))

      # GREEN: spotfinder centers of mass
      x, y = frame.pyslip.tiles.picture_fast_slow_to_map_relative(
            spot.spot_dict['xyzobs.px.value'][0], spot.spot_dict['xyzobs.px.value'][1])
      green_data.append((x,y))
      s1 = (x-5,y-5)
      s2 = (x-5,y+5)
      s3 = (x+5,y+5)
      s4 = (x+5,y-5)
      square = (s1,s2,s3,s4)
      green_squares.append((square,{'closed':True}))
      green_text.append((s3[0],s3[1],"%d"%(len(green_data)-1)))


    # Yellow: buffer strip around spot finder pixels
    for spot in self.buffers:
      for pixel in spot:
        yellow_data.append(
            frame.pyslip.tiles.picture_fast_slow_to_map_relative(pixel[0],pixel[1]))

    # Spring green: background pixels
    for spot in self.backgrounds:
      for pixel in spot:
        springgreen_data.append(
            frame.pyslip.tiles.picture_fast_slow_to_map_relative(pixel[0],pixel[1]))

    for spot in self.connected_spots:
      # PURPLE: spotfinder centers of the connected spots
      x, y = frame.pyslip.tiles.picture_fast_slow_to_map_relative(
        spot.spot_dict['xyzobs.px.value'][0], spot.spot_dict['xyzobs.px.value'][1])
      purple_data.append((x,y))

      s1 = (x-5,y-5)
      s2 = (x-5,y+5)
      s3 = (x+5,y+5)
      s4 = (x+5,y-5)
      square = (s1,s2,s3,s4)
      purple_squares.append((square,{'closed':True}))

      purple_text.append((s2[0],s2[1],"%d %d %d"% \
                          (spot.hkl.ohkl[0],spot.hkl.ohkl[1],spot.hkl.ohkl[2])))

    if show_spotfinder_spots:
      self.red_layer = frame.pyslip.AddPointLayer(
            red_data, color="red", name="<red_layer>",
            radius=1.5,
            renderer = frame.pyslip.LightweightDrawPointLayer,
            show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
      self.green_squares = frame.pyslip.AddPolygonLayer(
                green_squares, name = "<green_square>",
                colour = "green")
      if show_text:
        self.green_text = frame.pyslip.AddTextLayer(green_text, name="<green_text>")
      self.green_layer = frame.pyslip.AddPointLayer(
            green_data, color="green", name="<green_layer>",
            radius=1.5,
            renderer = frame.pyslip.LightweightDrawPointLayer,
            show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])

    if show_connected_spots:
      if show_text:
        self.purple_text = frame.pyslip.AddTextLayer(purple_text, name="<purple_text>")
      self.purple_squares = frame.pyslip.AddPolygonLayer(
            purple_squares, name = "<purple_square>",
            colour = "purple")
      self.purple_layer = frame.pyslip.AddPointLayer(
            purple_data, color="purple", name="<purple_layer>",
            radius=2,
            renderer = frame.pyslip.LightweightDrawPointLayer,
            show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
    if show_results:
      self.yellow_layer = frame.pyslip.AddPointLayer(
            yellow_data, color="yellow", name="<yellow_layer>",
            radius=2,
            renderer = frame.pyslip.LightweightDrawPointLayer,
            show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
      self.springgreen_layer = frame.pyslip.AddPointLayer(
            springgreen_data, color="spring green", name="<springgreen_layer>",
            radius=2,
            renderer = frame.pyslip.LightweightDrawPointLayer,
            show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])

def filter_indicies(ori,wavelength,resolution,phil,img):
  """ Given a unit cell, determine reflections in the diffracting condition, assuming the mosiaicity
  passed in the target phil file. Include their locations in reciprocal space given a crystal
  orientaiton.
  @param ori crystal orientation
  @param wavelength incident beam wavelength (angstroms)
  @param resolution limiting resolution to determine miller indices
  @param phil parsed small cell phil parameters
  @param img dxtbx format object
  @return list of original indices, list of asymmetric indices
  """
  beam = img.get_beam()
  sym = symmetry(unit_cell=ori.unit_cell(),space_group=phil.small_cell.spacegroup)
  ops = []
  for op in sym.space_group().expand_inv(sgtbx.tr_vec((0,0,0))).all_ops(): # this gets the spots related by inversion, aka Bijvoet mates
    r = op.r().as_hkl()
    subops = r.split(',')
    tmpop = [1,1,1]
    if '-' in subops[0]: tmpop[0] = -1
    if '-' in subops[1]: tmpop[1] = -1
    if '-' in subops[2]: tmpop[2] = -1

    if tmpop not in ops:
      ops.append(tmpop)

  asu_indices = sym.build_miller_set(anomalous_flag=False, d_min = resolution)
  #original_indicies = sym.build_miller_set(anomalous_flag=True, d_min = resolution)
  asu_indices_with_dups = []
  original_indicies = []
  for idx in asu_indices.indices():
    for op in ops:
      orig_idx = (idx[0]*op[0],idx[1]*op[1],idx[2]*op[2])
      if orig_idx not in original_indicies:
        original_indicies.append(orig_idx)
        asu_indices_with_dups.append(idx)
#  asu_indices = original_indicies.map_to_asu().indices()

  A = sqr(ori.reciprocal_matrix())
  s0 = col(beam.get_s0())

  ret_orig = []
  ret_asu = []
#  for index_o, index_a in zip(original_indicies.indices(),asu_indices):
  for index_o, index_a in zip(original_indicies,asu_indices_with_dups):
    s = A * col(index_o)
    q = s + s0
    ratio = q.length() / s0.length()
    if ratio > 1.0 - phil.small_cell.faked_mosaicity and ratio < 1.0 + phil.small_cell.faked_mosaicity:
      ret_orig.append(index_o)
      ret_asu.append(index_a)
  return (ret_orig, ret_asu)

def write_cell (ori,wavelength,max_clique,phil,img):
  """ Dump a series of useful debugging files viewable by gnuplot
  @param ori crystal orientation
  @param wavelength incident beam wavelength (angstroms)
  @param max_clique final maximum clique (list of small_cell_spot objects)
  @param phil parsed small cell phil parameters
  @param img dxtbx format object
   """
  A = sqr(ori.reciprocal_matrix())
  abasis = A * col((1,0,0))
  bbasis = A * col((0,1,0))
  cbasis = A * col((0,0,1))

  f = open("spots.dat",'w')
  for index in filter_indicies(ori,wavelength,phil.small_cell.high_res_limit,phil,img)[0]:
    v = A * col(index)
    f.write(" % 6.3f % 6.3f % 6.3f\n"%(v[0],v[1],v[2]))
  f.close()

  f = open("cell.dat",'w')
  dots = 10
  for h in xrange(-dots,dots):
    for k in xrange(-4,4):
      for l in xrange(-dots,dots):
        v = A * col([h,k,l])
        f.write(" % 6.3f % 6.3f % 6.3f\n"%(v[0],v[1],v[2]))
  f.close()

  f = open("clique_hkls.dat",'w')
  for spot in max_clique:
    v = A * col(spot.hkl.ohkl)
    f.write(" % 6.3f % 6.3f % 6.3f\n"%(v[0],v[1],v[2]))
  f.close()

  f = open("clique_xyzs.dat",'w')
  for spot in max_clique:
    v = spot.xyz
    f.write(" % 6.3f % 6.3f % 6.3f\n"%(v[0],v[1],v[2]))
  f.close()

  f = open("clique_hats.dat",'w')
  for spot in max_clique:
    v = (abasis*spot.hkl.ohkl[0])+(bbasis*spot.hkl.ohkl[1])+(cbasis*spot.hkl.ohkl[2])
    f.write(" % 6.3f % 6.3f % 6.3f\n"%(v[0],v[1],v[2]))
  f.close()

  f = open("arrows.p",'w')
  f.write("set arrow to % 6.3f, % 6.3f, % 6.3f lc rgb 'red'  \n"%(abasis[0],abasis[1],abasis[2]))
  f.write("set arrow to % 6.3f, % 6.3f, % 6.3f lc rgb 'green'\n"%(bbasis[0],bbasis[1],bbasis[2]))
  f.write("set arrow to % 6.3f, % 6.3f, % 6.3f lc rgb 'blue' \n"%(cbasis[0],cbasis[1],cbasis[2]))
  f.close()

  f = open("gp.p", 'w')
  f.write("set parametric\n")
  f.write("set ticslevel 0\n")
  f.write("""splot [-pi:pi][-pi/2:pi/2] (cos(u)*cos(v)/%f), (sin(u)*cos(v)/%f), (1/%f)-(sin(v)/%f), "cell.dat", "spots.dat", "clique_hkls.dat", "clique_xyzs.dat", "clique_hats.dat"\n"""%(wavelength,wavelength,wavelength,wavelength))
  f.write("""load "arrows.p"\n""")
  f.close()

"""
Not used
def hkl_to_xy_new (ori,wavelength,hkl,distance,bc,pixel_size):
  rv = ori.unit_cell().reciprocal_space_vector(hkl)
  rvre = sqr(ori.crystal_rotation_matrix()) * rv
  rvre = (rv[0],rv[1],rv[2] - (1/wavelength)) # direct beam anti-parallel (0,0,1)

  #dpx = 1765
  #dpy = 1765

  #dsx = dpx * pixel_size
  #dsy = dpy * pixel_size

  d = -distance / rvre[2]
  dx = rvre[0] * d
  dy = rvre[1] * d

  pxf = dx/pixel_size + 0.5 + bc[0]/pixel_size
  pyf = dy/pixel_size + 0.5 + bc[1]/pixel_size



  #pxf = (dx/dsx + 0.5) * dpx
  #pyf = (dy/dsy + 0.5) * dpy

  #pxf = pxf - (dpx/2 - bc[0]/pixel_size)
  #pyf = pyf - (dpy/2 - bc[1]/pixel_size)

  return(pxf,pyf)
  #return(dx,dy)
"""

def hkl_to_xy (ori,wavelength,hkl,distance,bc,pixel_size,img):
  """ Given an hkl, crystal orientation, and sufficient experimental parameters, compute
  the refelction's predicted xy position on a given image
  @param ori crystal orientation
  @param wavelength incident beam wavelength (angstroms)
  @param hkl small_cell_hkl object
  @param distance detector distance (mm)
  @param bc beam center (tuple, mm)
  @param pixel_size pixel size in mm
  @param img dxtbx format object
  """
  detector = img.get_detector()[0]
  beam = img.get_beam()

  detector_normal = col(detector.get_normal())
  detector_fast   = col(detector.get_fast_axis())
  detector_slow   = col(detector.get_slow_axis())
  detector_origin = col([bc[0], bc[1], 0.])
  pixel_size      = col([pixel_size,pixel_size,0])

  A = sqr(ori.reciprocal_matrix())

  #s0:  parallel to the direction of incident radiation
  s0 = col(beam.get_s0())
  s0_length = s0.length()
  s0_unit = s0.normalize();
  assert s0_length > 0.

  s = (A * hkl) #s, the reciprocal space coordinates, lab frame, of the oriented Miller index
  s_rad_sq = s.length_sq()
  assert s_rad_sq > 0.
  rotax = s.normalize().cross(s0_unit) #The axis that most directly brings the Bragg spot onto Ewald sphere
  chord_direction = (rotax.cross(s0)).normalize()

  a = s.length_sq()/(2.*s0_length) # see diagram
  b = math.sqrt(s.length_sq() - (a*a))#  Calculate half-length of the chord of intersection

  intersection = (-a * s0_unit) - (b * chord_direction)

  #acos_argument = intersection.dot(s) / s_rad_sq
  #iangle_1 = math.acos ( min(1.0,acos_argument) ) #avoid math domain error
  #assert approx_equal((intersection+s0).length()-s0_length,0. )

  #if (iangle_1 < effective_half_mosaicity_rad) { assume that it is

  q = intersection + s0
  q_unit = q.normalize()

  #check if diffracted ray parallel to detector face
  q_dot_n = q_unit.dot(detector_normal)
  assert q_dot_n != 0.

  return detector.get_ray_intersection_px(q)

def small_cell_index(path, horiz_phil):
  """ Index an image with a few spots and a known, small unit cell,
  with unknown basis vectors """

  import os,math

  # Load the dials and small cell parameters
  from dxtbx.datablock import DataBlockFactory
  from dials.algorithms.peak_finding.spotfinder_factory import SpotFinderFactory
  from dials.model.serialize.dump import reflections as reflections_dump

  print "Loading %s"%path

  # this image object isn't used until later when the predictions are made
  from dxtbx.format.Registry import Registry
  format_class = Registry.find(path)
  img = format_class(path)
  detector = img.get_detector()[0]
  beam = img.get_beam()

  # load the image
  try:
    # if it's a pickle file, read it directly
    import pickle
    DATA = pickle.load(open(path, 'rb'))
  except KeyError:
    # otherwise, we have to build the DATA dictionary manually
    from xfel.cxi.cspad_ana import cspad_tbx
    beam_x, beam_y = detector.get_beam_centre(beam.get_s0())
    DATA = cspad_tbx.dpack(active_areas=None,
                           address=None,
                           beam_center_x=beam_x,
                           beam_center_y=beam_x,
                           ccd_image_saturation=detector.get_trusted_range()[1],
                           data=img.get_raw_data(),
                           distance=detector.get_distance(),
                           pixel_size=detector.get_pixel_size()[0],
                           saturated_value=detector.get_trusted_range()[1],
                           timestamp=None,
                           wavelength=beam.get_wavelength())

  if not DATA.has_key('DETECTOR_ADDRESS') or DATA['DETECTOR_ADDRESS'] is None:
    DATA['DETECTOR_ADDRESS'] = 'CxiDs1-0|Cspad-0'

  if horiz_phil.small_cell.override_beam_x is not None:
    DATA['BEAM_CENTER_X'] = horiz_phil.small_cell.override_beam_x

  if horiz_phil.small_cell.override_beam_y is not None:
    DATA['BEAM_CENTER_Y'] = horiz_phil.small_cell.override_beam_y

  if horiz_phil.small_cell.override_distance is not None:
    DATA['DISTANCE'] = horiz_phil.small_cell.override_distance

  if horiz_phil.small_cell.override_wavelength is not None:
    DATA['WAVELENGTH'] = horiz_phil.small_cell.override_wavelength

  beam_x = DATA['BEAM_CENTER_X']
  beam_y = DATA['BEAM_CENTER_Y']
  beam_center = col((beam_x, beam_y))
  pixel_size = DATA['PIXEL_SIZE']
  wavelength = DATA['WAVELENGTH']
  distance = DATA['DISTANCE']
  raw_data = img.get_raw_data()

  print "Using beam center %s, %s, distance %s, and wavelength %s"%(beam_x, beam_y, distance, wavelength)

  # create the spot finder
  find_spots = SpotFinderFactory.from_parameters(horiz_phil)

  # spotfind
  datablock = DataBlockFactory.from_in_memory([img])
  reflections = find_spots(datablock)

  # filter the reflections for those near asic boundries
  from dials.algorithms.shoebox import MaskCode
  mask_peak = MaskCode.Valid|MaskCode.Foreground
  print "Filtering %s reflections by proximity to asic boundries..."%len(reflections),

  sel = flex.bool()
  for sb in reflections['shoebox']:
    focus = sb.mask.focus()
    l, r, t, b, z0, z1 = sb.bbox
    coords = []
    for my, y in zip(xrange(b-t),xrange(t,b)):
      for mx, x in zip(xrange(r-l),xrange(l,r)):
        if sb.mask[0,my,mx] == mask_peak:
          coords.append((x,y))
    test = flex.bool([is_bad_pixel(raw_data,c) for c in coords])
    sel.append(test.count(True) == 0)
  reflections = reflections.select(sel)

  reflections_dump(reflections, "spotfinder.pickle")
  print "saved %d"%len(reflections)

  recip_coords = flex.vec3_double()
  radial_sizes = flex.double()
  azimuthal_sizes = flex.double()
  for ref in reflections:
    # calculate reciprical space coordinates
    x, y, z = ref['xyzobs.px.value']
    xyz = col(detector.get_pixel_lab_coord((x,y)))
    xyz /= wavelength * xyz.length()
    xyz -= col(beam.get_s0()) # translate to origin of reciprocal space
    recip_coords.append(xyz)


    # Calculate unit-length radial and azimuthal (tangential) direction
    # vectors, r and a, respectively.  The azimuthal direction vector
    # is the radial vector rotated by 90 degrees counter-clockwise.

    radial = (col((x,y)) - beam_center).normalize()
    azimuthal = col((-radial[1],+radial[0]))

    # Determine the extent of the spot along the radial and azimuthal
    # directions from its center.

    a_max = float('-inf')
    a_min = float('+inf')
    r_max = float('-inf')
    r_min = float('+inf')
    l, r, t, b, z0, z1 = ref['shoebox'].bbox
    coords = []

    for my, y in zip(xrange(b-t),xrange(t,b)):
      for mx, x in zip(xrange(r-l),xrange(l,r)):
        if ref['shoebox'].mask[0,my,mx] == mask_peak:
          p = col((x,y)) - beam_center
          pa = p.dot(azimuthal)
          pr = p.dot(radial)

          if pa > a_max:
            a_max = pa
          if pa < a_min:
            a_min = pa
          if pr > r_max:
            r_max = pr
          if pr < r_min:
            r_min = pr

    radial_sizes.append(r_max - r_min)
    azimuthal_sizes.append(a_max - a_min)

  reflections['xyzrecip'] = recip_coords
  reflections['radial_size'] = radial_sizes
  reflections['azimuthal_size'] = azimuthal_sizes

  all_spots = []
  for i, ref in enumerate(reflections):
    all_spots.append(small_cell_spot(ref, i))

  # Unit cell calculated from indexed virtual powder diffraction
  sym = symmetry(unit_cell=horiz_phil.small_cell.powdercell,
                space_group_symbol=horiz_phil.small_cell.spacegroup)
  hkl_list = cctbx.miller.build_set(sym, False, d_min=horiz_phil.small_cell.high_res_limit)

  spacings = hkl_list.d_spacings()

  rcparams = sym.unit_cell().reciprocal().parameters()
  a = rcparams[0]
  b = rcparams[1]
  c = rcparams[2]
  alpha = rcparams[3] * math.pi / 180
  beta  = rcparams[4] * math.pi / 180
  gamma = rcparams[5] * math.pi / 180

  mm = sym.unit_cell().reciprocal().metrical_matrix()
  mm = sqr([mm[0],mm[3],mm[4],
            mm[3],mm[1],mm[5],
            mm[4],mm[5],mm[2]])

  # for every combination of spots examined, test possible translation and inversions
  # based on the symmetry of cell in question
  ops = []
  for op in sym.space_group().expand_inv(sgtbx.tr_vec((0,0,0))).all_ops(): # this gets the spots related by inversion, aka Bijvoet mates
    r = op.r().as_hkl()
    subops = r.split(',')
    tmpop = [1,1,1]
    if '-' in subops[0]: tmpop[0] = -1
    if '-' in subops[1]: tmpop[1] = -1
    if '-' in subops[2]: tmpop[2] = -1

    if tmpop not in ops:
      ops.append(tmpop)

  # make a list of the spots and the d-spacings they fall on
  spots_on_drings = []

  for spot in all_spots:
    dist = measure_distance((spot.spot_dict['xyzobs.px.value'][0],spot.spot_dict['xyzobs.px.value'][1]),(beam_x/pixel_size, beam_y/pixel_size))
    inner = dist - (spot.spot_dict['radial_size']/2)
    outer = dist + (spot.spot_dict['radial_size']/2)


    found_one = False
    for d in spacings:
      dpix = d_in_pixels(d[1], wavelength, distance, pixel_size)
      if dpix <= outer and dpix >= inner:
        # we will only examine asymmetric unit HKLs first.  Later we will try and determine original HKLs
        spot.hkls.append(small_cell_hkl(col(d[0]),col(d[0])))
        found_one = True

    if found_one:
      spots_on_drings.append(spot)

  overlap_limit = 5; overlap_count = 0
  for spot in spots_on_drings:
    if len(spot.hkls) > overlap_limit:
      spots_on_drings.remove(spot)
      overlap_count += 1
  print "Removed %d spots that overlaped more than %d rings."%(overlap_count,overlap_limit)

  print "Spots on d-rings:  %d"%len(spots_on_drings)
  print "Total sf spots:    %d"%len(reflections)

  max_clique_len = 0
  integrated_count = 0

  print "Finding spot connections..."

  # test every pair of spots to see if any of thier possible HKLs are connected
  spots_count = 2
  if len(spots_on_drings) >= spots_count:
    count = 0
    for i in itertools.permutations(range(len(spots_on_drings)),spots_count):
      count += 1
      spotA = spots_on_drings[i[0]]
      spotB = spots_on_drings[i[1]]

      for hklA in spotA.hkls:
        # don't test the same hklb twice.  This can happen if there is a zero in the index.
        tested_B = []
        for hklB_a in spotB.hkls:
          for op in ops:
            hklB = small_cell_hkl(hklB_a.ahkl, col([hklB_a.ahkl[0]*op[0],
                                                    hklB_a.ahkl[1]*op[1],
                                                    hklB_a.ahkl[2]*op[2]]))

            if hklA == hklB or hklB in tested_B:
              continue
            tested_B.append(hklB)

            approx_eq, delta_obv, delta_calc = test_spot_connection(hklA,hklB,spotA.xyz,spotB.xyz,mm,horiz_phil)

            if approx_eq:
              hklA.connections.append(small_cell_connection(hklA,hklB,spotA,spotB,delta_obv,delta_calc))
              hklB.connections.append(small_cell_connection(hklB,hklA,spotB,spotA,delta_obv,delta_calc))
              #print "EQUAL: spot %2d [% 2d,% 2d,% 2d] - spot %2d [% 2d,% 2d,% 2d] = %6.6f (obv), %6.6f (calc)"% \
                    #(spotA.ID, hklA.ohkl[0], hklA.ohkl[1], hklA.ohkl[2],
                     #spotB.ID, hklB.ohkl[0], hklB.ohkl[1], hklB.ohkl[2],
                     #delta_obv, delta_calc)

    # if I want to print out the full graph, i would do it here using spots_on_drings and test the connections attribute of each spot
    for spot in spots_on_drings:
      for hkl in spot.hkls:
        for con in hkl.connections:
          print "Spot ID", spot.ID, "OHKL", con.hkl1.get_ohkl_str(), "is connected to spot ID", con.spot2.ID, "OHKL", con.hkl2.get_ohkl_str()

    # Now, figure out which spot/hkl combo is the most connected spot/hkl
    most_connected_spot = None
    most_connected_hkl = None
    tie = []
    for spot in spots_on_drings:
      for hkl in spot.hkls:
        if len(hkl.connections) <= 0:
          continue
        if most_connected_spot is None or len(hkl.connections) > len(most_connected_hkl.connections):
          most_connected_spot = spot
          most_connected_hkl = hkl
          tie = []
        elif len(hkl.connections) == len(most_connected_hkl.connections):
          if len(tie) == 0:
            tie.append((most_connected_spot,most_connected_hkl))
          tie.append((spot,hkl))

    if most_connected_spot is None or most_connected_hkl is None:
      print "No spots can be connected to each other to resolve HKL ambiguities."
      print "IMAGE STATS %s: spots %5d, max clique: %5d, integrated %5d spots"%(path,len(reflections),max_clique_len,integrated_count)
      return

    if len(tie) > 0:
      print "TIE!  Picking the most connected spot with the smallest connection differences"
      most_connected_spot = None
      most_connected_hkl = None
      best_dist = float("inf")
      for spot, hkl in tie:
        dists = flex.double()
        for conn in hkl.connections:
          print spot.ID, conn.hkl1.ohkl.elems,conn.hkl2.ohkl.elems,
          dist = abs(conn.dobs - conn.dcalc)
          if not dist in dists:
            dists.append(dist)
            print "%32.32f"%dist
          else:
            print "DUP" #assume the same dist wouldn't occur twice, unless it's a symmetry operation of the same asu hkl
        avg_dist = sum(dists) / len(dists)
        print spot.ID, "avgdist", avg_dist
        #assert avg_dist != best_dist  # i mean, i guess this could happen, but not really prepared to deal with it
        if avg_dist < best_dist:
          best_dist = avg_dist
          most_connected_spot = spot
          most_connected_hkl = hkl
      assert most_connected_spot is not None and most_connected_hkl is not None

    print "Most connected spot: %3d %s with %d connections"%(most_connected_spot.ID, most_connected_hkl.ohkl.elems, len(most_connected_hkl.connections))

    most_connected_spot.hkl = most_connected_hkl # we get one for free

    print "Building clique graph..."

    # first zero out all the spot hkls arrays as we are going to re-assign them based on the most connected spot
    for spot in spots_on_drings:
      spot.hkls = []

    mapping = [] # 2ples.  (index into sub_clique array, index into spot's hkl array)
    sub_clique = []
    for conn in most_connected_hkl.connections:
      print "SPOT %3d (% 3d, % 3d, % 3d) <--> SPOT %3d (% 3d, % 3d, % 3d) dObs: %6.6f, dCalc: %6.6f, diff: %6.6f"% \
            (conn.spot1.ID, conn.hkl1.ohkl[0],conn.hkl1.ohkl[1],conn.hkl1.ohkl[2],
             conn.spot2.ID, conn.hkl2.ohkl[0],conn.hkl2.ohkl[1],conn.hkl2.ohkl[2],
             conn.dobs, conn.dcalc, abs(conn.dobs - conn.dcalc))

      conn.hkl2.connections = [] # zero these out as well so we can re-form them
      conn.spot2.hkls.append(conn.hkl2)
      if conn.spot2 in sub_clique:
        mapping.append((sub_clique.index(conn.spot2),len(conn.spot2.hkls)-1))
      else:
        sub_clique.append(conn.spot2)
        mapping.append((len(sub_clique)-1,len(conn.spot2.hkls)-1))

    # re-calculate the connections
    global degrees
    degrees = []
    for e1 in mapping:
      spot1 = sub_clique[e1[0]]
      hkl1  = spot1.hkls[e1[1]]

      approx_eq, delta_obv, delta_calc = test_spot_connection(hkl1,most_connected_hkl,
                                                              spot1.xyz,most_connected_spot.xyz,mm,horiz_phil)
      hkl1.connections.append(small_cell_connection(hkl1,most_connected_hkl,
                                                    spot1,most_connected_spot,delta_obv,delta_calc))

      for e2 in mapping:
        if e1 == e2:
          continue

        spot2 = sub_clique[e2[0]]
        hkl2  = spot2.hkls[e2[1]]

        approx_eq, delta_obv, delta_calc = test_spot_connection(hkl1,hkl2,spot1.xyz,spot2.xyz,mm,horiz_phil)
        if approx_eq:
          hkl1.connections.append(small_cell_connection(hkl1,hkl2,spot1,spot2,delta_obv,delta_calc))

      degrees.append(len(hkl1.connections))

    # sort the mapping based on degeneracy. this should speed clique finding.
    mapping = sorted(mapping, cmp = lambda x, y: cmp(len(sub_clique[x[0]].hkls[x[1]].connections),
                                                     len(sub_clique[y[0]].hkls[y[1]].connections)))
    degrees = flex.size_t(sorted(degrees))


    #build the clique graph
    graph = []
    for e1 in mapping:
      row = []
      spot1 = sub_clique[e1[0]]
      hkl1  = spot1.hkls[e1[1]]

      for e2 in mapping:
        if e1 == e2:
          row.append(0)
          continue

        spot2 = sub_clique[e2[0]]
        hkl2  = spot2.hkls[e2[1]]

        row.append(0)
        for conn in hkl1.connections:
          if conn.spot2 is spot2 and conn.hkl2 == hkl2:
            row[-1] = 1
            break

      graph.append(row)

    print mapping
    for row in graph:
      print row

    graph_lines = []

    #reported = []

    for j in range(len(mapping)):
      conn_count = 0
      spotA = sub_clique[mapping[j][0]]
      hklA = spotA.hkls[mapping[j][1]]
      line = "%d(%d,%d,%d) typeA "%(spotA.ID,hklA.ohkl[0],hklA.ohkl[1],hklA.ohkl[2])
      print "Spot %d %s is connected to "%(spotA.ID,hklA.ohkl.elems),
      #for i in range(j, len(mapping)):
      for i in range(j+1): #range(len(mapping)):
        if graph[j][i]:
          #a2b = (spotA.ID,hklA.ohkl.elems,spotB.ID,spotB.ohkl.elems)
          #b2a = (spotA.ID,hklA.ohkl.elems,spotB.ID,spotB.ohkl.elems)

          conn_count = conn_count + 1
          spotB = sub_clique[mapping[i][0]]
          hklB = spotB.hkls[mapping[i][1]]
          print "[%d, %s]"%(spotB.ID, hklB.ohkl.elems),
          line += "%d(%d,%d,%d) "%(spotB.ID,hklB.ohkl[0],hklB.ohkl[1],hklB.ohkl[2])
      print "Conn count:", conn_count
      graph_lines.append(line + "\n")

    print "converting to flex"
    graph_flex = flex.bool(flex.grid(len(graph),len(graph)))
    for j in xrange(len(graph)):
      for i in xrange(len(graph)):
        graph_flex[i,j] = bool(graph[i][j])

    #calcuate maximum size cliques using the Bron-Kerbosch algorithm:
    #http://en.wikipedia.org/wiki/Bron-Kerbosch_algorithm
    #code re-written from here (originally by Andy Hayden at #organizationName):
    #http://stackoverflow.com/questions/13904636/implementing-bronkerbosch-algorithm-in-python
    #choose the pivot to be the node with highest degree in the union of P and X, based on this paper:
    #http://www.sciencedirect.com/science/article/pii/S0304397508003903

    print "starting to find max clique of ", path

    _range = flex.size_t(range(len(mapping)))

    unmapped_cliques = []

    global total_calls
    total_calls = 0
    def bronk2(R, P, X, g):
        global degrees, total_calls
        total_calls = total_calls + 1
        if not any((P, X)):
            unmapped_cliques.append(R)
            return

        assert list(P.intersection(X)) == []

        u = P.concatenate(X)
        max_index, max_value = max(enumerate(degrees.select(u)), key=operator.itemgetter(1))
        pivot = u[max_index]

        n = N(pivot,g)
        b = flex.bool(len(P),True)
        for v in n:
          b = b & (P != v)

        subset =  P.select(b)
        for v in subset:
            R_v = R.concatenate(flex.size_t([v]))
            P_v = P.intersection(N(v,g))
            X_v = X.intersection(N(v,g))
            bronk2(R_v, P_v, X_v, g)
            P = P.select(P != v)
            X.append(v)
            X = flex.sorted(X)
    # N for neighbors
    def N(v, g):
        row = g[v:v+1,0:g.focus()[0]].as_1d()
        return _range.select(row)

    bronk2(flex.size_t(),flex.size_t(range(len(mapping))),flex.size_t(),graph_flex)

    print "Total calls to bronk: ", total_calls

    # map the cliques to the spots
    cliques = []
    for row in unmapped_cliques:
      new_row = []
      for column in row:
        new_row.append(mapping[column])
      cliques.append(new_row)
    print "cliques",
    print list(cliques)

    #find the biggest clique
    biggest = -1
    max_clique = []
    for clique in cliques:
      #print len(clique)
      if len(clique) > biggest:
        max_clique = clique
        biggest = len(clique)
    print "max clique:", max_clique

    max_clique_spots = []
    max_clique_spots.append(most_connected_spot)

    for entry in max_clique:
      spot = sub_clique[entry[0]]
      if spot in max_clique_spots:
        print "Duplicate spot in the max_clique, can't continue."
        print "IMAGE STATS %s: spots %5d, max clique: %5d, integrated %5d spots"%(path,len(reflections),max_clique_len,integrated_count)
        return

      assert spot.hkl is None
      spot.hkl = spot.hkls[entry[1]]
      max_clique_spots.append(spot)
      #print spot.ID, spot.hkl.ohkl.elems


    # resolve the ambiguity where two spots can have the same index, by looking at the distances observed and calculated and finding the best spot
    # build a dictionary where the keys are the original indices and the values are lists of spots.
    matched = {}
    for spot in max_clique_spots:
      key = spot.hkl.ohkl.elems
      if not key in matched.keys():
        matched[key] = []
      matched[key].append(spot)

    # an ambiguous entry will have multiple spots for the same index.  likely the spots are very close in reciprocal space.
    ambig_keys = []
    for key in matched:
      assert len(matched[key]) > 0
      if len(matched[key]) > 1:
        ambig_keys.append(key)

    for key in ambig_keys:
      print "Resolving ambiguity in", key, ": ", len(matched[key]), "spots with the same hkl"
      best_spot = None
      best_dist = float("inf")
      for spot in matched[key]:
        avg_dist = 0
        num_conns = 0
        for conn in spot.hkl.connections:
          if conn.spot2.hkl is not None and conn.spot2.hkl.ohkl.elems not in ambig_keys:
            print spot.ID, conn.hkl1.ohkl.elems,conn.hkl2.ohkl.elems,
            dist = abs(conn.dobs - conn.dcalc)
            print dist
            avg_dist = avg_dist + dist
            num_conns = num_conns + 1
        avg_dist = avg_dist / num_conns
        print spot.ID, "avgdist", avg_dist
        assert avg_dist != best_dist  # i mean, i guess this could happen, but not really prepared to deal with it
        if avg_dist < best_dist:
          best_dist = avg_dist
          best_spot = spot
      assert best_spot is not None
      for spot in matched[key]:
        if spot is not best_spot:
          max_clique_spots.remove(spot)

    #for spot in max_clique_spots:
      #if spot.ID == 11 and spot.hkl.ohkl.elems == (0,1,-1):
        #print "BREAK"

    if len(max_clique_spots) > 4:
      print "############################"
    print "Final resolved clique spots:"
    for spot in max_clique_spots:
      print spot.ID, spot.hkl.ohkl.elems
    if len(max_clique_spots) > 4:
      print "############################"
      gfile = open(os.path.basename(path).rstrip(".pickle") + ".sif", 'w')
      for line in graph_lines:
        gfile.write(line)
      gfile.close()

    working_set = []
    for spot in max_clique_spots:
      working_set.append(spot)

    max_clique_len = len(max_clique_spots)
    ok_to_integrate = False

    # loop, adding new spots to the clique and re-refining the unit cell paramters until no new spots can be added
    loop_count = 0
    while True:

      #calculate the basis vectors
      loop_count = loop_count + 1
      result = get_crystal_orientation(working_set, sym, img, beam_x, beam_y, distance, True, loop_count)
      if result is None:
        print "Couldn't get basis vectors for max clique"
        break
      ok_to_integrate = True

      # here I should test the angles too
      ori, refined_bcx, refined_bcy, refined_wavelength, refined_distance = result

      if approx_equal(ori.unit_cell().reciprocal().parameters()[0], a, out=None, eps=1.e-2) and \
         approx_equal(ori.unit_cell().reciprocal().parameters()[1], b, out=None, eps=1.e-2) and \
         approx_equal(ori.unit_cell().reciprocal().parameters()[2], c, out=None, eps=1.e-2):

        print "cell parameters approx. equal"
      else:
        print "cell parameters NOT APPROX EQUAL"

      sym.unit_cell().show_parameters()
      ori.unit_cell().show_parameters()

      # get the indicies that satisfy the difraction condition
      indicies_orig, indicies_asu = filter_indicies(ori,refined_wavelength,horiz_phil.small_cell.high_res_limit,horiz_phil,img)

      # caluculate predicted detector locations for these indicies
      f = open("preds.txt", "w")
      predicted = []
      for index in indicies_orig:
        xy = hkl_to_xy(ori,refined_wavelength,col(index),refined_distance,(refined_bcx,refined_bcy),pixel_size,img)
        if xy[0] > 0 and xy[0] <= detector.get_image_size()[0] and xy[1] > 0 and xy[1] <= detector.get_image_size()[0]:
          f.write("%d %d %d %f %f\n"%(index[0],index[1],index[2],xy[0],xy[1]))
        predicted.append(xy)
      f.close()

      # find spotfinder spots within a certain distance of predictions
      possibles = []
      for spot in working_set:
        possibles.append(spot)
        spot.pred = hkl_to_xy(ori,refined_wavelength,col(spot.hkl.ohkl.elems),refined_distance,(refined_bcx,refined_bcy),pixel_size,img)

      for spot in all_spots:
        found_it = False
        for ws_spot in working_set:
          if spot.ID == ws_spot.ID:
            found_it = True
            break
        if found_it:
          continue

        s_x = spot.spot_dict['xyzobs.px.value'][0]
        s_y = spot.spot_dict['xyzobs.px.value'][1]
        for index_o, index_a, pred in zip(indicies_orig, indicies_asu, predicted):
          if measure_distance((s_x, s_y),(pred[1],pred[0])) <= horiz_phil.small_cell.interspot_distance:
            # don't add the working_set spots here.  They were added above
            found_it = False
            for ws_spot in working_set:
              if ws_spot.hkl.ohkl.elems == index_o:
                found_it = True

            if not found_it:
              spot.hkl = small_cell_hkl(col(index_a),col(index_o))
              possibles.append(spot)

              spot.pred = pred

      # throw out any spots which share the same original index.  This can happen if a prediction is within 10 pixels of two spots.
      indexed = []
      for possible in possibles:
        count = 0
        for test in possibles:
          if test.hkl == possible.hkl:
            count += 1
        if count == 1:
          indexed.append(possible)

      indexed_rmsd = spots_rmsd(indexed)
      working_rmsd = spots_rmsd(working_set)

      print "Working set: %d spots, RMSD: %f"%(len(working_set),working_rmsd)
      print "Indexed set: %d spots, RMSD: %f"%(len(indexed),indexed_rmsd)
      print "Working set: ",
      for s in working_set: print s.ID,
      print
      print "Possible set: ",
      for s in possibles: print s.ID,
      print
      print "Indexed set: ",
      for s in indexed: print s.ID,
      print

      if len(working_set) < len(indexed):# and working_rmsd * 1.1 < indexed_rmsd: # allow a small increase in RMSD
        working_set = indexed
        print "Doing another round of unit cell refinement"
      else:
        print "Done refining unit cell.  No new spots to add."
        break
      # end finding preds and crystal orientation matrix refinement loop

    if result is not None:
      write_cell(ori,refined_wavelength,indexed,horiz_phil,img)

    indexed_hkls = flex.vec2_double()
    indexed_intensities = flex.double()
    indexed_sigmas = flex.double()

    if ok_to_integrate:
      results = []
      buffers = []
      backgrounds = []
      indexed_hkls = flex.miller_index()
      indexed_intensities = flex.double()
      indexed_sigmas = flex.double()
      mapped_predictions = flex.vec2_double()
      max_signal = flex.double()

      rmsd = 0
      rmsd_n = 0
      for spot in indexed:
        peakpix = []
        peakvals = []
        tmp = []
        is_bad = False
        for p in spot.peak_pixels:
          if is_bad_pixel(raw_data,p):
            is_bad = True
            break
          p = (p[0]+.5,p[1]+.5)
          peakpix.append(p)
          tmp.append(p)
          peakvals.append(raw_data[int(p[1]),int(p[0])])
        if is_bad: continue

        buffers.append(grow_by(peakpix,1))

        tmp.extend(buffers[-1])
        backgrounds.append(grow_by(tmp,1))
        tmp.extend(backgrounds[-1])
        backgrounds[-1].extend(grow_by(tmp,1))

        background = []
        bg_vals = []
        raw_bg_sum = 0
        for p in backgrounds[-1]:
          i = raw_data[int(p[1]),int(p[0])]
          if i is not None and i > 0:
            background.append(p)
            bg_vals.append(i)
            raw_bg_sum += i

        ret = reject_background_outliers(background, bg_vals)
        if ret is None:
          print "Not enough background pixels to integrate spot %d"%spot.ID
          continue
        background, bg_vals = ret
        backgrounds[-1] = background

        bp_a,bp_b,bp_c = get_background_plane_parameters(bg_vals, background)

        intensity = 0
        bg_peak = 0
        for v,p in zip(peakvals,peakpix):
          intensity += v - (bp_a*p[0] + bp_b*p[1] + bp_c)
          bg_peak += bp_a*p[0] + bp_b*p[1] + bp_c

        gain = 7.5
        sigma = math.sqrt(gain * (intensity + bg_peak + ((len(peakvals)/len(bg_vals))**2) * raw_bg_sum))

        print "ID: %3d, ohkl: %s, ahkl: %s, I: %9.1f, sigI: %9.1f, RDiff: %9.6f"%( \
          spot.ID, spot.hkl.get_ohkl_str(), spot.hkl.get_ahkl_str(), intensity, sigma,
          (sqr(ori.reciprocal_matrix())*spot.hkl.ohkl - spot.xyz).length()) #list(sqr(ori.reciprocal_matrix()).inverse() * col(spot.rs_spot))
        #print list(sqr(ori.reciprocal_matrix())*spot.hkl.ohkl)
        #print list(spot.xyz), "Diff", list(sqr(ori.reciprocal_matrix())*spot.hkl.ohkl - spot.xyz)

        max_sig = raw_data[int(spot.spot_dict['xyzobs.px.value'][1]),int(spot.spot_dict['xyzobs.px.value'][0])]

        s = "Orig HKL: % 4d % 4d % 4d "%(spot.hkl.ohkl.elems)
        s = s + "Asu HKL: % 4d % 4d % 4d "%(spot.hkl.ahkl.elems)
        s = s + "I: % 10.1f sigI: % 8.1f I/sigI: % 8.1f "%(intensity, sigma, intensity/sigma)
        s = s + "Size (pix): %3d Max pix val: %6d\n"%(len(spot.peak_pixels),max_sig)
        results.append(s)

        if spot.pred is None:
          #mapped_predictions.append((spot.spot_dict['xyzobs.px.value'][1] + 0.5, spot.spot_dict['xyzobs.px.value'][0] + 0.5))
          mapped_predictions.append((spot.spot_dict['xyzobs.px.value'][0], spot.spot_dict['xyzobs.px.value'][1]))
        else:
          mapped_predictions.append((spot.pred[0],spot.pred[1]))

        indexed_hkls.append(spot.hkl.ohkl.elems)
        indexed_intensities.append(intensity)
        indexed_sigmas.append(sigma)
        max_signal.append(max_sig)

        if spot.pred is not None:#and spot.ID != 5: #and spot.ID != 0 and spot.ID != 15:
          rmsd_n += 1
          rmsd += measure_distance(col((spot.spot_dict['xyzobs.px.value'][0],spot.spot_dict['xyzobs.px.value'][1])),col(spot.pred))**2

      if len(results) >= horiz_phil.small_cell.min_spots_to_integrate:
        f = open(path + ".int","w")
        for line in results:
          f.write(line)
        f.close()

        info = dict(
          xbeam = refined_bcx,
          ybeam = refined_bcy,
          distance = distance,
          wavelength = wavelength,
          #residual = local["r_residual"],
          #mosaicity = local["r_mosaicity"],
          pointgroup = horiz_phil.small_cell.spacegroup,
          observations = [cctbx.miller.set(sym,indexed_hkls).array(indexed_intensities,indexed_sigmas)],
          mapped_predictions = [mapped_predictions],
          model_partialities = [None],
          sa_parameters = [None],
          max_signal = [max_signal],
          current_orientation = [ori],
          current_cb_op_to_primitive = [sgtbx.change_of_basis_op()], #identity.  only support primitive lattices.
        )
        G = open(os.path.join(path.rstrip(os.path.basename(path)),"int-" + os.path.basename(path).strip()),"wb")
        pickle.dump(info,G,pickle.HIGHEST_PROTOCOL)

        print "cctbx.small_cell: integrated %d spots."%len(results),
        integrated_count = len(results)
      else:
        print "cctbx.small_cell: not enough spots to integrate (%d)."%len(results),

      if rmsd_n > 0:
        print " RMSD: %f"%math.sqrt((1/rmsd_n)*rmsd)
      else:
        print " Cannot calculate RMSD.  Not enough integrated spots or not enough clique spots near predictions."

      if True:# and len(results) >= horiz_phil.small_cell.min_spots_to_integrate:# and len(max_clique_spots) > 4:
        from rstbx.command_line.slip_viewer import master_str as slip_params
        from libtbx import phil
        from spotfinder import phil_str
        from spotfinder.command_line.signal_strength import additional_spotfinder_phil_defs

        work_phil = phil.process_command_line("",master_string=slip_params + phil_str + additional_spotfinder_phil_defs)
        work_params = work_phil.work.extract()

        frame = SmallCellXrayFrame(None, -1, "X-ray image display", size=(800,720))

        # Update initial settings with default values.  Needs
        # to be done before image is loaded (but after the frame is
        # instantiated).
        frame.params = work_params
        frame.pyslip.tiles.user_requests_antialiasing = work_params.anti_aliasing
        frame.settings_frame.panel.center_ctrl.SetValue(True)
        frame.settings_frame.panel.integ_ctrl.SetValue(True)
        frame.settings_frame.panel.spots_ctrl.SetValue(False)
        frame.settings.show_effective_tiling = work_params.show_effective_tiling
        frame.settings_frame.panel.collect_values()
        #paths = work_phil.remaining_args

        #from rstbx.apps.slip_helpers import slip_callbacks
        worker = small_cell_callbacks()
        worker.horiz_phil = horiz_phil
        worker.predicted = predicted
        worker.indicies = indicies_orig
        #worker.predicted = rectangle_data
        #worker.hkllist = hkllist
        worker.pixel_size = pixel_size
        worker.spots = all_spots
        #worker.spots = spots_on_drings
        #worker.connected_spots = connected_spots
        worker.connected_spots = max_clique_spots
        #worker.spot_pixels = spot_pixels
        worker.buffers = buffers
        worker.backgrounds = backgrounds

        worker.ori = ori
        worker.beam_x = refined_bcx
        worker.beam_y = refined_bcy
        worker.distance = distance
        worker.img = img

        frame.user_callback = worker.slip_callback
        frame.load_image(format_class(path))
        frame.Show()
        app.MainLoop()

    print "IMAGE STATS %s: spots %5d, max clique: %5d, integrated %5d spots"%(path,len(all_spots),max_clique_len,integrated_count)
    return max_clique_len

def spots_rmsd(spots):
  """ Calculate the rmsd for a series of small_cell_spot objects
  @param list of small_cell_spot objects
  @param RMSD (pixels) of each spot
  """
  rmsd = 0
  for spot in spots:
    rmsd += measure_distance(col((spot.spot_dict['xyzobs.px.value'][0],spot.spot_dict['xyzobs.px.value'][1])),col(spot.pred))**2
  return math.sqrt(rmsd/len(spots))

def hkl_to_xyz(hkl,abasis,bbasis,cbasis):
  """ Compute reciprocal space coordinates of an hkl given a set of basis vectors
  @param hkl miller index (tuple)
  @param abasis vector a
  @param bbasis vector b
  @param cbasis vector c
  @return reciprocal space coordinates of the hkl
  """
  return (hkl[0]*abasis) + (hkl[1]*bbasis) + (hkl[2]*cbasis)

def get_crystal_orientation(spots, sym, img, beam_x, beam_y, distance, use_minimizer=True, loop_count = 0):
  """ given a set of refelctions and input geometry, determine a crystal orientation using a set
  of linear equations, then refine it.
  @param spots list of small cell spot objects
  @param sym cctbx symmetry object
  @param img dxtbx format object
  @param beam_x beam center x position (mm)
  @param beam_y beam center y position (mm)
  @param distance detector distance (mm)
  @param use_minimizer if true, refine final orientation using a miminizer
  @param loop_count output during printout (debug use only)
  @return a tuple containing the orientation and refined beam x, beam y, wavelength, and distance
  """

  # determine initial orientation matrix from the set of reflections
  miller_indices = flex.vec3_double([i.hkl.ohkl for i in spots])
  u_vectors = flex.vec3_double([i.xyz for i in spots])
  from xfel.small_cell.solve_orientation import small_cell_orientation
  solver = small_cell_orientation(miller_indices, u_vectors, sym)
  ori = solver.unrestrained_setting()

  from cctbx import crystal_orientation
  F = sqr(sym.unit_cell().fractionalization_matrix()).transpose()

  try:
    Amat_start = sqr(ori.crystal_rotation_matrix()) * F
    ori_start = crystal_orientation.crystal_orientation(Amat_start, crystal_orientation.basis_type.reciprocal)

    print "powder  cell and residuals round %3d from %3d spots "%(loop_count,len(spots)),"[%.7f, %.7f, %.7f]"%(0,0,0),         ;sym.unit_cell().show_parameters()
    print "derived cell and residuals round %3d from %3d spots "%(loop_count,len(spots)),"[%.7f, %.7f, %.7f]"%tuple(solver.residuals),;ori.unit_cell().show_parameters()

    det = sqr(ori_start.crystal_rotation_matrix()).determinant()
    if not approx_equal(det, 1.0, out=None, eps=1.e-3):
      # some attempts to deal with matrices with non-one determinates is preserved here, commented out
      if approx_equal(det, -1.0, out=None, eps=1.e-3):
        #print "Reflected rotation matrix.  Fixing..."
        print "Note, reflected rotation matrix."

        #ori_start = ori_start.make_positive() # reflect the basis vectors

        #det = sqr(ori_start.crystal_rotation_matrix()).determinant()
        #if not approx_equal(det, 1.0, out=None, eps=1.e-3):
          #print "Couldn't fix the rotation matrix.  Cannot use this image."
          #return None

        #for spot in spots: # fix the signs of the hkls in the clique using this new basis
          #spot.hkl.flipped = True

      else:
        print "Invalid rotation matrix.  Determinant: ", det
        return None
  except Exception,e: # can fail here w/ a corrupt metrical matrix or math domain error
    print e.message
    return None

  if not use_minimizer:
    return (ori_start, beam_x, beam_y, img.get_beam().get_wavelength(), distance)

  #try:
    #o = optimise_basis(spots,ori_start) # call to old minimizer not used anymore
  #except Exception, e:
    #print e.message
    #return None

  from scitbx.math import r3_rotation_unit_quaternion_as_matrix,r3_rotation_axis_and_angle_from_matrix,euler_angles_as_matrix
  #r3 = sqr(r3_rotation_unit_quaternion_as_matrix(col(o.x[0:4])))
  #aa = r3_rotation_axis_and_angle_from_matrix(r3)
  #ori_rot = ori_start.rotate_thru(col(aa.axis),aa.angle())
  ori_rot = ori_start

  """
  Several minimizers are listed here as possible canidates for minimizing the orientation matrix. Only one is actually
  used (see the below call to minimize)
  FIXME: remove scipy dependency
  """

  # minimize using scipy: http://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html
  #params_start = sym.unit_cell().parameters() #these three lines commented out because not optimizing unit cell anymore
  #constraints = list(params_start[0:3])
  #constraints.append(params_start[4])

  #constraints.append(beam_x)
  #constraints.append(beam_y)
  #constraints.append(img.wavelength)
  #constraints.append(img.distance)

  def simplex_uc(x):
    for i in xrange(4):
      if x[i] < constraints[i] * 0.95:
        x[i] = constraints[i] * 0.95
      elif x[i] > constraints[i] * 1.05:
        x[i] = constraints[i] * 1.05

    ##x[6] = img.wavelength
    #if x[6] < constraints[6] * 0.99:
      #x[6] = constraints[6] * 0.99
    #elif x[6] > constraints[6] * 1.01:
      #x[6] = constraints[6] * 1.01

    a = x[0]
    b = x[1]
    c = x[2]
    alpha = 90
    beta  = x[3]
    gamma = 90
    #bcx = x[4]
    #bcy = x[5]
    #wavelength = float(x[6])
    #distance = x[7]
    e1 = x[4]
    e2 = x[5]
    e3 = x[6]

    rotation = euler_angles_as_matrix((e1,e2,e3),deg=True)

    sym_working = symmetry(unit_cell="%f, %f, %f, %f, %f, %f"%(a,b,c,alpha,beta,gamma),
                           space_group_symbol=SPACE_GROUP)
    F_working = sqr(sym_working.unit_cell().fractionalization_matrix()).transpose()
    Amat_working = sqr(ori_rot.crystal_rotation_matrix()) * rotation * F_working

    ori_working = crystal_orientation.crystal_orientation(Amat_working, crystal_orientation.basis_type.reciprocal)

    M = sqr(ori_working.reciprocal_matrix())

    f = 0
    for spot in spots:
      #spot_obs = col((spot.sf_spot.ctr_mass_y(),spot.sf_spot.ctr_mass_x()))
      #spot_calc = col(hkl_to_xy(ori_working,wavelength,spot.hkl.ohkl,distance,(bcx,bcy),img.pixel_size))

      #print "%.1f"%(spot_obs - spot_calc).length(),
      #f += ((spot_obs - spot_calc).length())**2
      diff = spot.xyz - (M * spot.hkl.ohkl)
      f += diff.dot(diff)

    #print "f: %f, unit_cell = %f, %f, %f, %f, %f, %f, bc: (%f,%f), wavelength: %f, dist: %f, euler: %f, %f, %f"%\
              #(f,a,b,c,alpha,beta,gamma,bcx,bcy,wavelength,distance,e1,e2,e3)

    return f

  def simplex_uc_ori_only(x, sym):

    e1 = x[0]
    e2 = x[1]
    e3 = x[2]

    rotation = euler_angles_as_matrix((e1,e2,e3),deg=True)

    sym_working = sym
    F_working = sqr(sym_working.unit_cell().fractionalization_matrix()).transpose()
    Amat_working = rotation * sqr(ori_rot.crystal_rotation_matrix()) * F_working

    ori_working = crystal_orientation.crystal_orientation(Amat_working, crystal_orientation.basis_type.reciprocal)

    M = sqr(ori_working.reciprocal_matrix())

    f = 0
    for spot in spots:
      diff = spot.xyz - (M * spot.hkl.ohkl)
      f += diff.dot(diff)

    #print "f: %f, unit_cell = %f, %f, %f, %f, %f, %f, bc: (%f,%f), wavelength: %f, dist: %f, euler: %f, %f, %f"%\
              #(f,a,b,c,alpha,beta,gamma,bcx,bcy,wavelength,distance,e1,e2,e3)

    #print "f: %f, euler: %f, %f, %f"%(f,e1,e2,e3)

    return f

  #p = sym.unit_cell().parameters() # DON'T REFINE UNIT CELL ANYMORE
  #x0 = np.array([p[0],p[1],p[2],p[4],beam_x,beam_y,img.wavelength,img.distance,0,0,0]) # DON'T REFINE UNIT CELL ANYMORE
  #res = minimize(simplex_uc, x0, method='nelder-mead',    # DON'T REFINE UNIT CELL ANYMORE
  #               options={'xtol': 1e-8, 'disp': True})    # DON'T REFINE UNIT CELL ANYMORE
  #res = minimize(simplex_uc, x0, method='BFGS', jac=None,
                   #options={'disp': True})


  def simplex_uc_rotz_only(x):

      rotz = x[0]

      #rotation = euler_angles_as_matrix((e1,e2,e3),deg=True)

      #sym_working = symmetry(unit_cell=SYM_STR, space_group_symbol=SPACE_GROUP)
      #F_working = sqr(sym_working.unit_cell().fractionalization_matrix()).transpose()
      #Amat_working = rotation * sqr(ori_rot.crystal_rotation_matrix()) * F_working

      #ori_working = crystal_orientation.crystal_orientation(Amat_working, crystal_orientation.basis_type.reciprocal)

      #M = sqr(ori_working.reciprocal_matrix())
      M = sqr(ori_rot.reciprocal_matrix())

      Z = sqr([math.cos(rotz),-math.sin(rotz),0,
               math.sin(rotz), math.cos(rotz),0,
               0,0,1])

      f = 0
      for spot in spots:
        diff = spot.xyz - (Z * M * spot.hkl.ohkl)
        f += diff.dot(diff)

      #print "f: %f, unit_cell = %f, %f, %f, %f, %f, %f, bc: (%f,%f), wavelength: %f, dist: %f, euler: %f, %f, %f"%\
                #(f,a,b,c,alpha,beta,gamma,bcx,bcy,wavelength,distance,e1,e2,e3)

      #print "f: %f, euler: %f, %f, %f"%(f,e1,e2,e3)

      return f

  x0 = np.array([0,0,0])
  res = minimize(simplex_uc_ori_only, x0, args=(sym,), method='nelder-mead',
                 options={'xtol': 1e-8, 'disp': True})

  #rotation_final = euler_angles_as_matrix((res.x[4],res.x[5],res.x[6]),deg=True) # DON'T REFINE UNIT CELL ANYMORE
  rotation_final = euler_angles_as_matrix((res.x[0],res.x[1],res.x[2]),deg=True)
  #print "Refined beam x: %f, beam y: %f, wavelength: %f, distance: %f"%(res.x[4],res.x[5],res.x[6],res.x[7])

  #sym_final = symmetry(unit_cell="%f, %f, %f, %f, %f, %f"%(res.x[0],res.x[1],res.x[2],90,res.x[3],90),# DON'T REFINE UNIT CELL ANYMORE
                       #space_group_symbol=SPACE_GROUP)# DON'T REFINE UNIT CELL ANYMORE
  sym_final = sym
  F_final = sqr(sym_final.unit_cell().fractionalization_matrix()).transpose()
  #Z = sqr([math.cos(res.x[0]),-math.sin(res.x[0]),0,
           #math.sin(res.x[0]), math.cos(res.x[0]),0,
           #0,0,1])
  #Amat_final = Z * sqr(ori_rot.crystal_rotation_matrix()) * F_final
  Amat_final = rotation_final * sqr(ori_rot.crystal_rotation_matrix()) * F_final
  #Amat_final = rotation_final * F_final

  ori_final = crystal_orientation.crystal_orientation(Amat_final, crystal_orientation.basis_type.reciprocal)

  m = ori_final.crystal_rotation_matrix()
  print "Final crystal rotation matrix: ", list(m)

  # not sure this is the correct conversion to euler angles...
  m = m.as_list_of_lists()
  print "Final euler angles: ", math.atan2(m[2][0],m[2][1])*180/math.pi,math.acos(m[2][2])*180/math.pi,math.atan2(m[0][2],m[1][2])*180/math.pi

  #return (ori_final, float(res.x[4]),float(res.x[5]),float(res.x[6]),float(res.x[7]))
  return (ori_final, beam_x, beam_y, img.get_beam().get_wavelength(), distance)


from scitbx import lbfgs, matrix
class optimise_basis(object):
  """Given a set of points in three-space and a unit cell's initial position,
  find the best orientation of the unit cell to match the spots"""
  # this minimizer no longer used

  def __init__(self, spots, ori_start):
    """
    @param spots
    @param ori_start
    """

    #self.x = flex.double((.5,.5,.5,.5,1))
    #self.x = flex.double((1,0,0,0,1))
    self.x = flex.double((1,0,0,0))
    #p = .95 # very close to identity rotation
    #a =  math.sqrt(p)
    #bcd = math.sqrt((1-p)/3)
    #self.x = flex.double((a,bcd,bcd,bcd,1))

    self.spots = spots
    for spot in spots:
      spot.start_xyz = sqr(ori_start.reciprocal_matrix()) * spot.hkl.ohkl

    print "------"
    self.minimizer = lbfgs.run(target_evaluator=self)


  def compute_functional_and_gradients(self):
    print "Pre-norm:     constraint: %6.3f, "%(self.x[0]**2+self.x[1]**2+self.x[2]**2+self.x[3]**2),  \
      "x: % 8.8f % 8.8f % 8.8f % 8.8f"%tuple(self.x)

    norm = col((self.x[0],self.x[1],self.x[2],self.x[3]))
    norm = norm.normalize()
    self.x[0] = norm[0]
    self.x[1] = norm[1]
    self.x[2] = norm[2]
    self.x[3] = norm[3]

    a = self.x[0]
    b = self.x[1]
    c = self.x[2]
    d = self.x[3]
    w = col((self.x[1],self.x[2],self.x[3]))
    #lagrange = self.x[4]

    f = 0
    dfda = 0
    dfdb = 0
    dfdc = 0
    dfdd = 0

    # finite differences used to verify gradients
    fd_delta = 1e-7
    fda  = a+fd_delta
    fdwb = col((b+fd_delta,c,d))
    fdwc = col((b,c+fd_delta,d))
    fdwd = col((b,c,d+fd_delta))
    #fdl  = lagrange+fd_delta

    fdga = fdgb = fdgc = fdgd = 0#fdgl = 0

    for spot in self.spots:
      f += optimise_basis.f(spot.start_xyz,a,w,spot.xyz)

      fdga += optimise_basis.f(spot.start_xyz,fda,w ,spot.xyz)
      fdgb += optimise_basis.f(spot.start_xyz,a,fdwb,spot.xyz)
      fdgc += optimise_basis.f(spot.start_xyz,a,fdwc,spot.xyz)
      fdgd += optimise_basis.f(spot.start_xyz,a,fdwd,spot.xyz)

      x = spot.start_xyz[0]
      y = spot.start_xyz[1]
      z = spot.start_xyz[2]

      #o: observed
      xo = spot.xyz[0]
      yo = spot.xyz[1]
      zo = spot.xyz[2]

      p = xo-x-2*a*z*c+2*a*y*d-2*c*y*b+2*x*c**2+2*x*d**2-2*z*b*d
      q = yo-y-2*a*x*d+2*a*z*b-2*d*z*c+2*y*d**2+2*y*b**2-2*b*x*c
      r = zo-z-2*a*y*b+2*a*x*c-2*b*x*d+2*z*b**2+2*z*c**2-2*c*y*d

      dfda += 2*p*(-2*z*c+2*y*d)+2*q*(-2*x*d+2*z*b)+2*r*(-2*y*b+2*x*c)
      dfdb += 2*p*(-2*c*y-2*z*d)+2*q*(2*a*z+4*y*b-2*x*c)+2*r*(-2*a*y-2*x*d+4*z*b)
      dfdc += 2*p*(-2*a*z-2*y*b +4*x*c)+2*q*(-2*d*z-2*b*x)+2*r*(2*a*x+4*z*c-2*y*d)
      dfdd += 2*p*(2*a*y+4*x*d-2*z*b)+2*q*(-2*a*x-2*z*c+4*y*d)+2*r*(-2*b*x-2*c*y)

    #l_step = lagrange * (a**2 + b**2 + c**2 + d**2 - 1)
    #f += l_step
    #dfda += 2*lagrange*a
    #dfdb += 2*lagrange*b
    #dfdc += 2*lagrange*c
    #dfdb += 2*lagrange*d

    #dfdl = a**2+b**2+c**2+d**2-1

    #fdga += l_step
    #fdgb += l_step
    #fdgc += l_step
    #fdgd += l_step

    fdga = (fdga - f)/fd_delta
    fdgb = (fdgb - f)/fd_delta
    fdgc = (fdgc - f)/fd_delta
    fdgd = (fdgd - f)/fd_delta

    #fdgl = ((f-l_step+(fdl * (a**2 + b**2 + c**2 + d**2 - 1))) - f)/fd_delta


    g  = flex.double((dfda,dfdb,dfdc,dfdd))#,dfdl))
    fd = flex.double((fdga,fdgb,fdgc,fdgd))#,fdgl))

    print "F: %8.7f, constraint: %6.3f, "%(f,a**2+b**2+c**2+d**2), "x: % 8.8f % 8.8f % 8.8f % 8.8f"%tuple(self.x)
    print "Gradients (analytically derived)    : % 8.8f % 8.8f % 8.8f % 8.8f"%tuple(g)
    print "Gradients (frac. diffs. estimates)  : % 8.8f % 8.8f % 8.8f % 8.8f"%tuple(fd)
    print "------"
    return (f, g)

  @staticmethod
  def f(start_xyz,a,w,obs_xyz):
    delta = (2*a*w.cross(start_xyz)) + (2*w.cross(w.cross(start_xyz)))
    diff = obs_xyz - (start_xyz + delta)
    return diff.dot(diff)

def grow_by(pixels, amt):
  """
  Given a list of pixels, grow it contiguously by the given
  number of pixels
  """
  ret = []
  tested = []

  def recurse(pixel, depth):
    if depth > amt or pixel in ret or pixel in tested:
      return
    if not pixel in pixels:
      ret.append(pixel)
    if depth == 0:
      tested.append(pixel)

    recurse((pixel[0]-1,pixel[1]  ),depth+1)
    recurse((pixel[0]+1,pixel[1]  ),depth+1)
    recurse((pixel[0]  ,pixel[1]-1),depth+1)
    recurse((pixel[0]  ,pixel[1]+1),depth+1)

  for pix in pixels:
    recurse(pix,0)

  return ret

def get_background_value(img, center):
  square_radius = 5
  data = flex.double()
  for j in xrange(2*square_radius+1):
    pass # THIS FUNCTION NOT USED YET

def reject_background_outliers(bg_pixels, bg_vals):
  """
  Given a set of background pixel coordinates and values, determine the 80% that
  are most similar, taking account a background plane gradient. Recursive.
  @param bg_pixels list of 2D pixel values
  @param bg_values corresponding pixel values
  @return the culled list of background pixels
  """
  assert len(bg_vals) == len(bg_pixels)

  pairs = sorted(zip(bg_pixels,bg_vals),key=lambda pair:pair[1])
  bg_pixels = []
  bg_vals = []
  for bgc, bgv in pairs:
    bg_pixels.append(bgc)
    bg_vals.append(bgv)

  eighty_percent = int(0.8*len(pairs))
  if len(pairs[0:eighty_percent]) <= 0:
    return None

  bp_a,bp_b,bp_c = get_background_plane_parameters(bg_vals[0:eighty_percent], bg_pixels[0:eighty_percent])

  fit_vals = []

  for p, q, v in [(a[0][1],a[0][0],a[1]) for a in pairs]:
    fit_vals.append(v-(p*bp_a + q*bp_b + bp_c))
  stddev = flex.double(fit_vals).standard_deviation_of_the_sample()

  culled_bg = []
  culled_bg_vals = []
  for i in xrange(len(fit_vals)):
    if abs(fit_vals[i]) < 3 * stddev:
      culled_bg.append(bg_pixels[i])
      culled_bg_vals.append(bg_vals[i])


  if len(culled_bg) == len(bg_pixels):
    return bg_pixels, bg_vals
  else:
    #print "Rejected %d background pixels"%(len(bg_pixels)-len(culled_bg))
    return reject_background_outliers(culled_bg,culled_bg_vals)

def get_background_plane_parameters(bgvals,bgpixels):
  """ Given a set of pixels, determine the best fit plane assuming they are background
  see http://journals.iucr.org/d/issues/1999/10/00/ba0027/index.html
  @param bg_pixels list of 2D pixel values
  @param bg_values corresponding pixel values
  @return the background plane parameters
  """
  assert len(bgvals) == len(bgpixels)

  M1 = [0.,0.,0.,0.,0.,0.,0.,0.,0.]
  M2 = [0.,0.,0.]
  n = len(bgvals)

  for v, pix in zip(bgvals,bgpixels):
    p = pix[1]
    q = pix[0]

    M1[0] += p**2
    M1[1] += p*q
    M1[2] += p
    M1[3] += p*q
    M1[4] += q**2
    M1[5] += q
    M1[6] += p
    M1[7] += q
    M1[8] += n

    M2[0] += p*v
    M2[1] += q*v
    M2[2] += v

  return sqr(M1).inverse() * col(M2)


def is_bad_pixel(raw_data, pix):
  """ if a pixel is in the below diagram and is white, then it is too close to the edge of the tile:
  --X--
  -XXX-
  XX*XX
  -XXX-
  --X--

  *: location of @p pix
  @param pixel aray
  @param pixel coordinate
  @return bool indicating whether the pixel is bad
  """
  pixels = [(-2, 0),
            (-1, 1),
            (-1, 0),
            (-1,-1),
            ( 0,-2),
            ( 0,-1),
            ( 0, 0),
            ( 0, 1),
            ( 0, 2),
            ( 1,-1),
            ( 1, 0),
            ( 1, 1),
            ( 2, 0)]
  for p in pixels:
    i = raw_data[pix[1]+p[1],pix[0]+p[0]]
    if i == None or i <= 0:
      return True

  return False
