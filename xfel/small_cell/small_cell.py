from __future__ import division
from __future__ import print_function
from six.moves import range
#-*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

""" Series of functions used by cctbx.small_cell """

import numpy as np
from scipy.optimize import minimize
import math
from dials.array_family import flex
from libtbx.test_utils import approx_equal
import itertools
from scitbx.matrix import col, sqr
from cctbx.crystal import symmetry
import cctbx.miller
from cctbx.uctbx import unit_cell
from cctbx import sgtbx
import operator
from dxtbx.model.experiment_list import ExperimentListFactory, ExperimentListDumper

from dials.algorithms.shoebox import MaskCode
mask_peak = MaskCode.Valid|MaskCode.Foreground

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
    for my, y in zip(range(b-t),range(t,b)):
      for mx, x in zip(range(r-l),range(l,r)):
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

def filter_indicies(ori,beam,resolution,phil):
  """ Given a unit cell, determine reflections in the diffracting condition, assuming the mosiaicity
  passed in the target phil file. Include their locations in reciprocal space given a crystal
  orientation.
  @param ori crystal orientation
  @param dxtbx beam object
  @param resolution limiting resolution to determine miller indices
  @param phil parsed small cell phil parameters
  @return list of original indices, list of asymmetric indices
  """
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
  asu_indices_with_dups = []
  original_indicies = []
  for idx in asu_indices.indices():
    for op in ops:
      orig_idx = (idx[0]*op[0],idx[1]*op[1],idx[2]*op[2])
      if orig_idx not in original_indicies:
        original_indicies.append(orig_idx)
        asu_indices_with_dups.append(idx)

  A = sqr(ori.reciprocal_matrix())
  s0 = col(beam.get_s0())

  ret_orig = []
  ret_asu = []
  for index_o, index_a in zip(original_indicies,asu_indices_with_dups):
    s = A * col(index_o)
    q = s + s0
    ratio = q.length() / s0.length()
    if ratio > 1.0 - phil.small_cell.faked_mosaicity and ratio < 1.0 + phil.small_cell.faked_mosaicity:
      ret_orig.append(index_o)
      ret_asu.append(index_a)
  return (ret_orig, ret_asu)

def write_cell (ori,beam,max_clique,phil):
  """ Dump a series of useful debugging files viewable by gnuplot
  @param ori crystal orientation
  @param beam dxtbx beam object
  @param max_clique final maximum clique (list of small_cell_spot objects)
  @param phil parsed small cell phil parameters
  @param img dxtbx format object
   """
  wavelength = beam.get_wavelength()
  A = sqr(ori.reciprocal_matrix())
  abasis = A * col((1,0,0))
  bbasis = A * col((0,1,0))
  cbasis = A * col((0,0,1))

  f = open("spots.dat",'w')
  for index in filter_indicies(ori,beam,phil.small_cell.high_res_limit,phil)[0]:
    v = A * col(index)
    f.write(" % 6.3f % 6.3f % 6.3f\n"%(v[0],v[1],v[2]))
  f.close()

  f = open("cell.dat",'w')
  dots = 10
  for h in range(-dots,dots):
    for k in range(-4,4):
      for l in range(-dots,dots):
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

def hkl_to_xy (ori,hkl,detector,beam):
  """ Given an hkl, crystal orientation, and sufficient experimental parameters, compute
  the refelction's predicted xy position on a given image
  @param ori crystal orientation
  @param detector dxtbx detector object
  @param beam dxtbx beam object
  """

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
  q = intersection + s0

  try:
    panel_id, xy = detector.get_ray_intersection(q)
  except RuntimeError:
    return None, None
  xy = detector[panel_id].millimeter_to_pixel(xy)
  return panel_id, xy

def ori_to_crystal(ori, spacegroup):
  from dxtbx.model import MosaicCrystalSauter2014
  direct_matrix = ori.direct_matrix()
  real_a = direct_matrix[0:3]
  real_b = direct_matrix[3:6]
  real_c = direct_matrix[6:9]
  crystal = MosaicCrystalSauter2014(real_a, real_b, real_c, spacegroup)
  crystal.set_domain_size_ang(100) # hardcoded here, but could be refiend using nave_parameters
  crystal.set_half_mosaicity_deg(0.05) # hardcoded here, but could be refiend using nave_parameters
  return crystal

def small_cell_index(path, horiz_phil):
  """ Index an image with a few spots and a known, small unit cell,
  with unknown basis vectors """

  # Load the dials and small cell parameters
  from dxtbx.model.experiment_list import ExperimentListFactory
  from dials.algorithms.spot_finding.factory import SpotFinderFactory
  from dials.model.serialize.dump import reflections as reflections_dump

  print("Loading %s"%path)

  # load the image
  from dxtbx.format.Registry import Registry
  format_class = Registry.find(path)
  img = format_class(path)
  imageset = img.get_imageset([path])
  raw_data = img.get_raw_data()

  # create the spot finder
  find_spots = SpotFinderFactory.from_parameters(horiz_phil)

  # spotfind
  experiments = ExperimentListFactory.from_imageset_and_crystal(imageset, None)[0]
  reflections = find_spots(experiments)

  # filter the reflections for those near asic boundries
  print("Filtering %s reflections by proximity to asic boundries..."%len(reflections), end=' ')

  sel = flex.bool()
  for sb in reflections['shoebox']:
    focus = sb.mask.focus()
    l, r, t, b, z0, z1 = sb.bbox
    coords = []
    for my, y in zip(range(b-t),range(t,b)):
      for mx, x in zip(range(r-l),range(l,r)):
        if sb.mask[0,my,mx] == mask_peak:
          coords.append((x,y))
    test = flex.bool([is_bad_pixel(raw_data,c) for c in coords])
    sel.append(test.count(True) == 0)
  reflections = reflections.select(sel)

  reflections_dump(reflections, "spotfinder.pickle")
  print("saved %d"%len(reflections))

  max_clique_len, experiments, refls = small_cell_index_detail(experiments, reflections, horiz_phil)
  return max_clique_len

def small_cell_index_detail(experiments, reflections, horiz_phil, write_output = True):
  """ Index an image with a few spots and a known, small unit cell,
  with unknown basis vectors """
  import os,math

  imagesets = experiments.imagesets()
  assert len(imagesets) == 1
  imageset = imagesets[0]
  path = imageset.paths()[0]

  detector = imageset.get_detector()
  beam = imageset.get_beam()

  if horiz_phil.small_cell.override_wavelength is not None:
    beam.set_wavelength(horiz_phil.small_cell.override_wavelength)
  wavelength = beam.get_wavelength()
  s0 = col(beam.get_s0())
  s0u = s0.normalize()

  raw_data = imageset[0]
  if not isinstance(raw_data, tuple):
    raw_data = (raw_data,)

  recip_coords = flex.vec3_double()
  radial_labs = flex.vec3_double()
  radial_sizes = flex.double()
  azimuthal_sizes = flex.double()
  s0_projs = flex.vec3_double()
  for ref in reflections:
    # calculate reciprical space coordinates
    x, y, z = ref['xyzobs.px.value']
    panel = detector[ref['panel']]
    xyz_lab = col(panel.get_pixel_lab_coord((x,y)))
    xyz = xyz_lab / (wavelength * xyz_lab.length())
    xyz -= col(beam.get_s0()) # translate to origin of reciprocal space
    recip_coords.append(xyz)

    # Calculate unit-length radial and azimuthal (tangential) direction
    # vectors, r and a, respectively.  The azimuthal direction vector
    # is the radial vector rotated by 90 degrees counter-clockwise.

    s0_proj = xyz_lab.length()*math.cos(xyz_lab.angle(s0)) * s0u
    radial_lab = xyz_lab - s0_proj
    radial = radial_lab.normalize()
    azimuthal = radial.cross(s0u)

    # Determine the extent of the spot along the radial and azimuthal
    # directions from its center.

    a_max = float('-inf')
    a_min = float('+inf')
    r_max = float('-inf')
    r_min = float('+inf')
    l, r, t, b, z0, z1 = ref['shoebox'].bbox
    coords = []

    for my, y in zip(range(b-t),range(t,b)):
      for mx, x in zip(range(r-l),range(l,r)):
        if ref['shoebox'].mask[0,my,mx] == mask_peak:
          p = col(panel.get_pixel_lab_coord((x,y)))
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

    radial_labs.append(radial_lab)
    radial_sizes.append(r_max - r_min)
    azimuthal_sizes.append(a_max - a_min)
    s0_projs.append(s0_proj)

  reflections['xyzrecip'] = recip_coords
  reflections['radial_lab'] = radial_labs
  reflections['radial_size'] = radial_sizes
  reflections['azimuthal_size'] = azimuthal_sizes
  reflections['s0_proj'] = s0_projs

  from dials.algorithms.indexing import index_reflections
  reflections['imageset_id'] = reflections['id']
  reflections.centroid_px_to_mm(detector)
  reflections.map_centroids_to_reciprocal_space(detector, beam)

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
    dist = col(spot.spot_dict['radial_lab']).length()
    inner = dist - (spot.spot_dict['radial_size']/2)
    outer = dist + (spot.spot_dict['radial_size']/2)

    #L = 2dsinT
    inner_angle = math.atan2(inner, col(spot.spot_dict['s0_proj']).length())
    outer_angle = math.atan2(outer, col(spot.spot_dict['s0_proj']).length())
    outer_d = wavelength/2/math.sin(inner_angle/2) # inner becomes outer
    inner_d = wavelength/2/math.sin(outer_angle/2) # outer becomes inner

    found_one = False
    for d in spacings:
      if d[1] <= outer_d and d[1] >= inner_d:
        # we will only examine asymmetric unit HKLs first.  Later we will try and determine original HKLs
        spot.hkls.append(small_cell_hkl(col(d[0]),col(d[0])))
        found_one = True

    if found_one:
      spots_on_drings.append(spot)

  overlap_limit = horiz_phil.small_cell.d_ring_overlap_limit; overlap_count = 0
  if overlap_limit is None:
    print("Accepting all spots on d-rings")
  else:
    for spot in spots_on_drings:
      if len(spot.hkls) > overlap_limit:
        spots_on_drings.remove(spot)
        overlap_count += 1
    print("Removed %d spots that overlaped more than %d rings."%(overlap_count,overlap_limit))

  print("Spots on d-rings:  %d"%len(spots_on_drings))
  print("Total sf spots:    %d"%len(reflections))

  max_clique_len = 0
  integrated_count = 0

  print("Finding spot connections...")

  # test every pair of spots to see if any of their possible HKLs are connected
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
          print("Spot ID", spot.ID, "OHKL", con.hkl1.get_ohkl_str(), "is connected to spot ID", con.spot2.ID, "OHKL", con.hkl2.get_ohkl_str())

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
      print("No spots can be connected to each other to resolve HKL ambiguities.")
      print("IMAGE STATS %s: spots %5d, max clique: %5d, integrated %5d spots"%(path,len(reflections),max_clique_len,integrated_count))
      return

    if len(tie) > 0:
      print("TIE!  Picking the most connected spot with the smallest connection differences")
      most_connected_spot = None
      most_connected_hkl = None
      best_dist = float("inf")
      for spot, hkl in tie:
        dists = flex.double()
        for conn in hkl.connections:
          print(spot.ID, conn.hkl1.ohkl.elems,conn.hkl2.ohkl.elems, end=' ')
          dist = abs(conn.dobs - conn.dcalc)
          if not dist in dists:
            dists.append(dist)
            print("%32.32f"%dist)
          else:
            print("DUP") #assume the same dist wouldn't occur twice, unless it's a symmetry operation of the same asu hkl
        avg_dist = sum(dists) / len(dists)
        print(spot.ID, "avgdist", avg_dist)
        #assert avg_dist != best_dist  # i mean, i guess this could happen, but not really prepared to deal with it
        if avg_dist < best_dist:
          best_dist = avg_dist
          most_connected_spot = spot
          most_connected_hkl = hkl
      assert most_connected_spot is not None and most_connected_hkl is not None

    print("Most connected spot: %3d %s with %d connections"%(most_connected_spot.ID, most_connected_hkl.ohkl.elems, len(most_connected_hkl.connections)))

    most_connected_spot.hkl = most_connected_hkl # we get one for free

    print("Building clique graph...")

    # first zero out all the spot hkls arrays as we are going to re-assign them based on the most connected spot
    for spot in spots_on_drings:
      spot.hkls = []

    mapping = [] # 2ples.  (index into sub_clique array, index into spot's hkl array)
    sub_clique = []
    for conn in most_connected_hkl.connections:
      print("SPOT %3d (% 3d, % 3d, % 3d) <--> SPOT %3d (% 3d, % 3d, % 3d) dObs: %6.6f, dCalc: %6.6f, diff: %6.6f"% \
            (conn.spot1.ID, conn.hkl1.ohkl[0],conn.hkl1.ohkl[1],conn.hkl1.ohkl[2],
             conn.spot2.ID, conn.hkl2.ohkl[0],conn.hkl2.ohkl[1],conn.hkl2.ohkl[2],
             conn.dobs, conn.dcalc, abs(conn.dobs - conn.dcalc)))

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

    print(mapping)
    for row in graph:
      print(row)

    graph_lines = []

    for j in range(len(mapping)):
      conn_count = 0
      spotA = sub_clique[mapping[j][0]]
      hklA = spotA.hkls[mapping[j][1]]
      line = "%d(%d,%d,%d) typeA "%(spotA.ID,hklA.ohkl[0],hklA.ohkl[1],hklA.ohkl[2])
      print("Spot %d %s is connected to "%(spotA.ID,hklA.ohkl.elems), end=' ')
      for i in range(j+1):
        if graph[j][i]:
          conn_count = conn_count + 1
          spotB = sub_clique[mapping[i][0]]
          hklB = spotB.hkls[mapping[i][1]]
          print("[%d, %s]"%(spotB.ID, hklB.ohkl.elems), end=' ')
          line += "%d(%d,%d,%d) "%(spotB.ID,hklB.ohkl[0],hklB.ohkl[1],hklB.ohkl[2])
      print("Conn count:", conn_count)
      graph_lines.append(line + "\n")

    print("converting to flex")
    graph_flex = flex.bool(flex.grid(len(graph),len(graph)))
    for j in range(len(graph)):
      for i in range(len(graph)):
        graph_flex[i,j] = bool(graph[i][j])

    #calcuate maximum size cliques using the Bron-Kerbosch algorithm:
    #http://en.wikipedia.org/wiki/Bron-Kerbosch_algorithm
    #code re-written from here (originally by Andy Hayden at #organizationName):
    #http://stackoverflow.com/questions/13904636/implementing-bronkerbosch-algorithm-in-python
    #choose the pivot to be the node with highest degree in the union of P and X, based on this paper:
    #http://www.sciencedirect.com/science/article/pii/S0304397508003903

    print("starting to find max clique of ", path)

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

    print("Total calls to bronk: ", total_calls)

    # map the cliques to the spots
    cliques = []
    for row in unmapped_cliques:
      new_row = []
      for column in row:
        new_row.append(mapping[column])
      cliques.append(new_row)
    print("cliques", end=' ')
    print(list(cliques))

    #find the biggest clique
    biggest = -1
    max_clique = []
    for clique in cliques:
      #print len(clique)
      # use >= here since the list should be sorted such that later entries have nodes with higher degree
      # spots, so in the case of a tie, picking a later entry will get a clique where the nodes are more
      # connected
      if len(clique) >= biggest:
        max_clique = clique
        biggest = len(clique)
    print("max clique:", max_clique)

    max_clique_spots = []
    max_clique_spots.append(most_connected_spot)

    for entry in max_clique:
      spot = sub_clique[entry[0]]
      if spot in max_clique_spots:
        print("Duplicate spot in the max_clique, can't continue.")
        print("IMAGE STATS %s: spots %5d, max clique: %5d, integrated %5d spots"%(path,len(reflections),max_clique_len,integrated_count))
        return

      assert spot.hkl is None
      spot.hkl = spot.hkls[entry[1]]
      max_clique_spots.append(spot)


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
      print("Resolving ambiguity in", key, ": ", len(matched[key]), "spots with the same hkl")
      best_spot = None
      best_dist = float("inf")
      for spot in matched[key]:
        avg_dist = 0
        num_conns = 0
        for conn in spot.hkl.connections:
          if conn.spot2.hkl is not None and conn.spot2.hkl.ohkl.elems not in ambig_keys:
            print(spot.ID, conn.hkl1.ohkl.elems,conn.hkl2.ohkl.elems, end=' ')
            dist = abs(conn.dobs - conn.dcalc)
            print(dist)
            avg_dist = avg_dist + dist
            num_conns = num_conns + 1
        avg_dist = avg_dist / num_conns
        print(spot.ID, "avgdist", avg_dist)
        assert avg_dist != best_dist  # i mean, i guess this could happen, but not really prepared to deal with it
        if avg_dist < best_dist:
          best_dist = avg_dist
          best_spot = spot
      assert best_spot is not None
      for spot in matched[key]:
        if spot is not best_spot:
          max_clique_spots.remove(spot)

    if len(max_clique_spots) > 4:
      print("############################")
    print("Final resolved clique spots:")
    for spot in max_clique_spots:
      print(spot.ID, spot.hkl.ohkl.elems)
    if len(max_clique_spots) > 4:
      print("############################")
      # Uncomment this to write a sif file, useful as input to graph display programs
      #gfile = open(os.path.splitext(os.path.basename(path))[0] + ".sif", 'w')
      #for line in graph_lines:
      #  gfile.write(line)
      #gfile.close()

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
      ori = get_crystal_orientation(working_set, sym, False, loop_count)
      if ori is None:
        print("Couldn't get basis vectors for max clique")
        break
      ok_to_integrate = True

      # here I should test the angles too
      if approx_equal(ori.unit_cell().reciprocal().parameters()[0], a, out=None, eps=1.e-2) and \
         approx_equal(ori.unit_cell().reciprocal().parameters()[1], b, out=None, eps=1.e-2) and \
         approx_equal(ori.unit_cell().reciprocal().parameters()[2], c, out=None, eps=1.e-2):

        print("cell parameters approx. equal")
      else:
        print("cell parameters NOT APPROX EQUAL")

      sym.unit_cell().show_parameters()
      ori.unit_cell().show_parameters()

      from dials.algorithms.refinement.prediction.managed_predictors import ExperimentsPredictorFactory
      from cctbx import miller
      import copy
      crystal = ori_to_crystal(ori, horiz_phil.small_cell.spacegroup)
      experiments = ExperimentListFactory.from_imageset_and_crystal(imageset, crystal)
      reflections['id'] = flex.int(len(reflections), -1)
      index_reflections(reflections, experiments)#, tolerance=0.15)
      reflections['miller_index_asymmetric'] = copy.deepcopy(reflections['miller_index'])
      miller.map_to_asu(crystal.get_space_group().type(), True, reflections['miller_index_asymmetric'])
      ref_predictor = ExperimentsPredictorFactory.from_experiments(experiments, force_stills=experiments.all_stills())
      reflections = ref_predictor(reflections)

      indexed = []
      for i, ref in enumerate(reflections):
        if ref['id'] < 0: continue
        spot = small_cell_spot(ref, i)
        spot.pred = ref['xyzcal.px'][0:2]
        spot.pred_panel_id = ref['panel']
        spot.hkl = small_cell_hkl(col(ref['miller_index_asymmetric']), col(ref['miller_index']))
        indexed.append(spot)

      indexed_rmsd = spots_rmsd(indexed)
      working_rmsd = spots_rmsd(working_set)

      print("Working set: %d spots, RMSD: %f"%(len(working_set),working_rmsd))
      print("Indexed set: %d spots, RMSD: %f"%(len(indexed),indexed_rmsd))
      print("Working set: ", end=' ')
      for s in working_set: print(s.ID, end=' ')
      print()
      print("Indexed set: ", end=' ')
      for s in indexed: print(s.ID, end=' ')
      print()

      #print "**** SHOWING DISTS ****"
      #for spot in indexed:
      #  if spot.pred is None:
      #    print "NO PRED"; continue
      #  s_x, s_y, _ = spot.spot_dict['xyzobs.px.value']
      #  _dist = measure_distance((s_x, s_y),(spot.pred[0],spot.pred[1]))
      #  print _dist
      #print "**** SHOWED DISTS ****"

      if len(working_set) < len(indexed):# and working_rmsd * 1.1 < indexed_rmsd: # allow a small increase in RMSD
        working_set = indexed
        print("Doing another round of unit cell refinement")
      else:
        print("Done refining unit cell.  No new spots to add.")
        break
      # end finding preds and crystal orientation matrix refinement loop

    if ori is not None and horiz_phil.small_cell.write_gnuplot_input:
      write_cell(ori,beam,indexed,horiz_phil)

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
      mapped_panels = flex.size_t()
      max_signal = flex.double()
      xyzobs = flex.vec3_double()
      xyzvar = flex.vec3_double()
      shoeboxes = flex.shoebox()
      s1 = flex.vec3_double()
      bbox = flex.int6()

      rmsd = 0
      rmsd_n = 0
      for spot in indexed:
        if spot.pred is None: continue
        peakpix = []
        peakvals = []
        tmp = []
        is_bad = False
        panel = detector[spot.pred_panel_id]
        panel_raw_data = raw_data[spot.pred_panel_id]
        for p in spot.peak_pixels:
          #if is_bad_pixel(panel_raw_data,p):
          #  is_bad = True
          #  break
          p = (p[0]+.5,p[1]+.5)
          peakpix.append(p)
          tmp.append(p)
          peakvals.append(panel_raw_data[int(p[1]),int(p[0])])
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
          try:
            i = panel_raw_data[int(p[1]),int(p[0])]
          except IndexError:
            continue
          if i is not None and i > 0:
            background.append(p)
            bg_vals.append(i)
            raw_bg_sum += i

        ret = reject_background_outliers(background, bg_vals)
        if ret is None:
          print("Not enough background pixels to integrate spot %d"%spot.ID)
          continue
        background, bg_vals = ret
        backgrounds[-1] = background

        bp_a,bp_b,bp_c = get_background_plane_parameters(bg_vals, background)

        intensity = 0
        bg_peak = 0
        for v,p in zip(peakvals,peakpix):
          intensity += v - (bp_a*p[0] + bp_b*p[1] + bp_c)
          bg_peak += bp_a*p[0] + bp_b*p[1] + bp_c

        gain = panel.get_gain()
        sigma = math.sqrt(gain * (intensity + bg_peak + ((len(peakvals)/len(bg_vals))**2) * raw_bg_sum))

        print("ID: %3d, ohkl: %s, ahkl: %s, I: %9.1f, sigI: %9.1f, RDiff: %9.6f"%( \
          spot.ID, spot.hkl.get_ohkl_str(), spot.hkl.get_ahkl_str(), intensity, sigma,
          (sqr(ori.reciprocal_matrix())*spot.hkl.ohkl - spot.xyz).length()))

        max_sig = panel_raw_data[int(spot.spot_dict['xyzobs.px.value'][1]),int(spot.spot_dict['xyzobs.px.value'][0])]

        s = "Orig HKL: % 4d % 4d % 4d "%(spot.hkl.ohkl.elems)
        s = s + "Asu HKL: % 4d % 4d % 4d "%(spot.hkl.ahkl.elems)
        s = s + "I: % 10.1f sigI: % 8.1f I/sigI: % 8.1f "%(intensity, sigma, intensity/sigma)
        s = s + "Size (pix): %3d Max pix val: %6d\n"%(len(spot.peak_pixels),max_sig)
        results.append(s)

        if spot.pred is None:
          mapped_predictions.append((spot.spot_dict['xyzobs.px.value'][0], spot.spot_dict['xyzobs.px.value'][1]))
          mapped_panels.append(spot.spot_dict['panel'])
        else:
          mapped_predictions.append((spot.pred[0],spot.pred[1]))
          mapped_panels.append(spot.pred_panel_id)
        xyzobs.append(spot.spot_dict['xyzobs.px.value'])
        xyzvar.append(spot.spot_dict['xyzobs.px.variance'])
        shoeboxes.append(spot.spot_dict['shoebox'])

        indexed_hkls.append(spot.hkl.ohkl.elems)
        indexed_intensities.append(intensity)
        indexed_sigmas.append(sigma)
        max_signal.append(max_sig)
        s1.append(s0+spot.xyz)
        bbox.append(spot.spot_dict['bbox'])

        if spot.pred is not None:
          rmsd_n += 1
          rmsd += measure_distance(col((spot.spot_dict['xyzobs.px.value'][0],spot.spot_dict['xyzobs.px.value'][1])),col(spot.pred))**2

      if len(results) >= horiz_phil.small_cell.min_spots_to_integrate:
        # Uncomment to get a text version of the integration results
        #f = open(os.path.splitext(os.path.basename(path))[0] + ".int","w")
        #for line in results:
        #  f.write(line)
        #f.close()

        if write_output:
          info = dict(
            xbeam = refined_bcx,
            ybeam = refined_bcy,
            distance = distance,
            wavelength = wavelength,
            pointgroup = horiz_phil.small_cell.spacegroup,
            observations = [cctbx.miller.set(sym,indexed_hkls).array(indexed_intensities,indexed_sigmas)],
            mapped_predictions = [mapped_predictions],
            mapped_panels = [mapped_panels],
            model_partialities = [None],
            sa_parameters = [None],
            max_signal = [max_signal],
            current_orientation = [ori],
            current_cb_op_to_primitive = [sgtbx.change_of_basis_op()], #identity.  only support primitive lattices.
            pixel_size = pixel_size,
          )
          G = open("int-" + os.path.splitext(os.path.basename(path))[0] +".pickle","wb")
          import pickle
          pickle.dump(info,G,pickle.HIGHEST_PROTOCOL)

        crystal = ori_to_crystal(ori, horiz_phil.small_cell.spacegroup)
        experiments = ExperimentListFactory.from_imageset_and_crystal(imageset, crystal)
        if write_output:
          dump = ExperimentListDumper(experiments)
          dump.as_json(os.path.splitext(os.path.basename(path).strip())[0]+"_integrated_experiments.json")

        refls = flex.reflection_table()
        refls['id'] = flex.int(len(indexed_hkls), 0)
        refls['panel'] = mapped_panels
        refls['integration.sum.value'] = indexed_intensities
        refls['integration.sum.variance'] = indexed_sigmas**2
        refls['xyzobs.px.value'] = xyzobs
        refls['xyzobs.px.variance'] = xyzvar
        refls['miller_index'] = indexed_hkls
        refls['xyzcal.px'] = flex.vec3_double(mapped_predictions.parts()[0], mapped_predictions.parts()[1], flex.double(len(mapped_predictions), 0))
        refls['shoebox'] = shoeboxes
        refls['entering'] = flex.bool(len(refls), False)
        refls['s1'] = s1
        refls['bbox'] = bbox

        refls.centroid_px_to_mm(detector)

        refls.set_flags(flex.bool(len(refls), True), refls.flags.indexed)
        if write_output:
          refls.as_pickle(os.path.splitext(os.path.basename(path).strip())[0]+"_integrated.pickle")

        print("cctbx.small_cell: integrated %d spots."%len(results), end=' ')
        integrated_count = len(results)
      else:
        raise RuntimeError("cctbx.small_cell: not enough spots to integrate (%d)."%len(results))

      if rmsd_n > 0:
        print(" RMSD: %f"%math.sqrt((1/rmsd_n)*rmsd))
      else:
        print(" Cannot calculate RMSD.  Not enough integrated spots or not enough clique spots near predictions.")

    print("IMAGE STATS %s: spots %5d, max clique: %5d, integrated %5d spots"%(path,len(all_spots),max_clique_len,integrated_count))
    return max_clique_len, experiments, refls

def spots_rmsd(spots):
  """ Calculate the rmsd for a series of small_cell_spot objects
  @param list of small_cell_spot objects
  @param RMSD (pixels) of each spot
  """
  rmsd = 0
  count = 0
  print('Spots with no preds', [spot.pred is None for spot in spots].count(True), 'of', len(spots))
  for spot in spots:
    if spot.pred is None:
      continue
    rmsd += measure_distance(col((spot.spot_dict['xyzobs.px.value'][0],spot.spot_dict['xyzobs.px.value'][1])),col(spot.pred))**2
    count += 1
  if count == 0: return 0
  return math.sqrt(rmsd/count)

def hkl_to_xyz(hkl,abasis,bbasis,cbasis):
  """ Compute reciprocal space coordinates of an hkl given a set of basis vectors
  @param hkl miller index (tuple)
  @param abasis vector a
  @param bbasis vector b
  @param cbasis vector c
  @return reciprocal space coordinates of the hkl
  """
  return (hkl[0]*abasis) + (hkl[1]*bbasis) + (hkl[2]*cbasis)

def get_crystal_orientation(spots, sym, use_minimizer=True, loop_count = 0):
  """ given a set of refelctions and input geometry, determine a crystal orientation using a set
  of linear equations, then refine it.
  @param spots list of small cell spot objects
  @param sym cctbx symmetry object
  @param use_minimizer if true, refine final orientation using a miminizer
  @param loop_count output during printout (debug use only)
  @return the orientation
  """

  # determine initial orientation matrix from the set of reflections
  miller_indices = flex.vec3_double([i.hkl.ohkl for i in spots])
  u_vectors = flex.vec3_double([i.xyz for i in spots])
  from xfel.small_cell.solve_orientation import small_cell_orientation
  solver = small_cell_orientation(miller_indices, u_vectors, sym)
  ori = solver.unrestrained_setting()

  det = sqr(ori.crystal_rotation_matrix()).determinant()
  print("Got crystal rotation matrix, deteriminant", det)
  if det <= 0:
    ori = ori.make_positive()
    for spot in spots: # fix the signs of the hkls in the clique using this new basis
      spot.hkl.flipped = True

  from cctbx import crystal_orientation
  F = sqr(sym.unit_cell().fractionalization_matrix()).transpose()

  try:
    Amat_start = sqr(ori.crystal_rotation_matrix()) * F
    ori_start = crystal_orientation.crystal_orientation(Amat_start, crystal_orientation.basis_type.reciprocal)

    print("powder  cell and residuals round %3d from %3d spots "%(loop_count,len(spots)),"[%.7f, %.7f, %.7f]"%(0,0,0), end=' ')         ;sym.unit_cell().show_parameters()
    print("derived cell and residuals round %3d from %3d spots "%(loop_count,len(spots)),"[%.7f, %.7f, %.7f]"%tuple(solver.residuals), end=' ');ori.unit_cell().show_parameters()
  except Exception as e: # can fail here w/ a corrupt metrical matrix or math domain error
    print(e.message)
    return None

  if not use_minimizer:
    return ori_start

  from scitbx.math import euler_angles_as_matrix
  ori_rot = ori_start

  """
  FIXME: remove scipy dependency
  """

  # minimize using scipy: http://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html

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

    return f

  x0 = np.array([0,0,0])
  res = minimize(simplex_uc_ori_only, x0, args=(sym,), method='nelder-mead',
                 options={'xtol': 1e-8, 'disp': True})

  rotation_final = euler_angles_as_matrix((res.x[0],res.x[1],res.x[2]),deg=True)

  sym_final = sym
  F_final = sqr(sym_final.unit_cell().fractionalization_matrix()).transpose()
  Amat_final = rotation_final * sqr(ori_rot.crystal_rotation_matrix()) * F_final

  ori_final = crystal_orientation.crystal_orientation(Amat_final, crystal_orientation.basis_type.reciprocal)

  m = ori_final.crystal_rotation_matrix()
  print("Final crystal rotation matrix: ", list(m))

  # not sure this is the correct conversion to euler angles...
  m = m.as_list_of_lists()
  print("Final euler angles: ", math.atan2(m[2][0],m[2][1])*180/math.pi,math.acos(m[2][2])*180/math.pi,math.atan2(m[0][2],m[1][2])*180/math.pi)

  return ori_final

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
  for j in range(2*square_radius+1):
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
  for i in range(len(fit_vals)):
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
  try:
    for p in pixels:
      i = raw_data[pix[1]+p[1],pix[0]+p[0]]
      if i == None or i <= 0:
        return True
  except IndexError:
    return True

  return False
