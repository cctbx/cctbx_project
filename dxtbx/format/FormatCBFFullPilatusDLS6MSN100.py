#!/usr/bin/env python
# FormatCBFFullPilatus.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Pilatus implementation of fullCBF format, for use with Dectris detectors.

from __future__ import absolute_import, division
from __future__ import print_function

#import pycbf

from builtins import zip
from builtins import range
from dxtbx.format.FormatCBFFullPilatus import FormatCBFFullPilatus
#from dxtbx.format.FormatPilatusHelpers import determine_pilatus_mask

class FormatCBFFullPilatusDLS6MSN100(FormatCBFFullPilatus):
  '''An image reading class for full CBF format images from Pilatus
  detectors.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an CBF format image, i.e. we can
    make sense of it.'''

    # this depends on DIALS for the goniometer shadow model; if missing
    # simply return False

    try:
      from dials.util.masking import GoniometerShadowMaskGenerator
    except ImportError as e:
      return False

    header = FormatCBFFullPilatus.get_cbf_header(image_file)

    for record in header.split('\n'):
      if '# Detector' in record and \
             'PILATUS' in record and 'S/N 60-0100 Diamond' in header:
        return True

    return False

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file.'''

    import libtbx
    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    self._dynamic_shadowing = kwargs.get('dynamic_shadowing', False)
    if self._dynamic_shadowing in (libtbx.Auto, 'Auto'):
      self._dynamic_shadowing = True
    FormatCBFFullPilatus.__init__(self, image_file, **kwargs)

    return

  def get_mask(self, goniometer=None):
    mask = super(FormatCBFFullPilatusDLS6MSN100, self).get_mask()
    if self._dynamic_shadowing:
      gonio_masker = self.get_goniometer_shadow_masker(goniometer=goniometer)
      scan = self.get_scan()
      detector = self.get_detector()
      shadow_mask = gonio_masker.get_mask(detector, scan.get_oscillation()[0])
      assert len(mask) == len(shadow_mask)
      for m, sm in zip(mask, shadow_mask):
        if sm is not None:
          m &= sm
    return mask

  def get_goniometer_shadow_masker(self, goniometer=None):
    if goniometer is None:
      goniometer = self.get_goniometer()

    assert goniometer is not None

    #avoid a module-level import from the DIALS namespace that kills LABELIT
    from dials.util.masking import GoniometerShadowMaskGenerator

    if goniometer.get_names()[1] == 'GON_CHI':
      # SmarGon

      class SmarGonShadowMaskGenerator(GoniometerShadowMaskGenerator):
        def __init__(SMG, goniometer):
          from scitbx.array_family import flex
          import math
          SMG.goniometer = goniometer

          coords = flex.vec3_double()
          axis = flex.size_t()

          # FACE A: Sample holder
          #   Defined as semi-circle of radius r(A) = 10 mm (centred on PHI axis)
          #   with rectangle of size a(A) = 12.8 mm (x 20 mm)

          offsetA = 33.0
          # semi-circle for phi=-90 ... +90
          radiusA = 10.0
          phi = flex.double_range(-90, 100, step=10) * math.pi/180
          x = flex.double(phi.size(), -offsetA)
          y = radiusA * flex.cos(phi)
          z = radiusA * flex.sin(phi)

          # corners of square
          sqdA = 12.8 # square depth
          nsteps = 10
          for i in range(nsteps+1):
            for sign in (+1, -1):
              x.append(-offsetA)
              y.append(i * -sqdA/nsteps)
              z.append(sign * radiusA)
          x.append(-offsetA)
          y.append(-sqdA)
          z.append(0)

          SMG.faceA = flex.vec3_double(x, y, z)

          # FACE B: Lower arm
          sx = -28.50
          sy = -4.90
          sz = 8.50
          mx = -13.80
          my = -26.00
          nx = -27.50
          ny = -29.50
          px = -65.50
          py = -29.50
          SMG.faceB = flex.vec3_double(((sx,sy,sz),(mx,my,0),(nx,ny,0),(px,py,0)))

          # FACE E: Rim of sample holder
          #   Defined as circle of radius r(E) = 6 mm (centred on PHI axis) at an
          #   offset o(E) = 19 mm

          offsetE = 19.0
          radiusE = 6.0
          phi = flex.double_range(0, 360, step=15) * math.pi/180
          x = flex.double(phi.size(), -offsetE)
          y = radiusE * flex.cos(phi)
          z = radiusE * flex.sin(phi)

          SMG.faceE = flex.vec3_double(x, y, z)

        def extrema_at_scan_angle(SMG, scan_angle):
          from scitbx.array_family import flex

          # Align end station coordinate system with ImgCIF coordinate system
          from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
          from scitbx import matrix
          R = align_reference_frame(matrix.col((-1,0,0)), matrix.col((1,0,0)),
                                    matrix.col((0,-1,0)), matrix.col((0,1,0)))
          faceA = R.elems * SMG.faceA
          faceE = R.elems * SMG.faceE

          axes = SMG.goniometer.get_axes()
          angles = SMG.goniometer.get_angles()
          scan_axis = SMG.goniometer.get_scan_axis()
          angles[scan_axis] = scan_angle

          extrema = flex.vec3_double()

          for coords in (faceA, faceE):
            coords = coords.deep_copy()
            for i, axis in enumerate(axes):
              if i == 0:
                continue # shadow doesn't change with phi setting
              sel = flex.bool(len(coords), True)
              rotation = matrix.col(
                axis).axis_and_angle_as_r3_rotation_matrix(angles[i], deg=True)
              coords.set_selected(sel, rotation.elems * coords.select(sel))
            extrema.extend(coords)

          s = matrix.col(SMG.faceB[0])
          mx, my, _ = SMG.faceB[1]
          nx, ny, _ = SMG.faceB[2]
          px, py, _ = SMG.faceB[3]

          Rchi = (R.inverse() * matrix.col(axes[1])).axis_and_angle_as_r3_rotation_matrix(angles[1], deg=True)
          sk = Rchi * s
          sxk, syk, szk = sk.elems
          coords = flex.vec3_double((
            (sxk, syk, 0),
            (sxk, syk, szk),
            (sxk+mx/2, syk+my/2, szk),
            (sxk+mx, syk+my, szk),
            (sxk+(mx+nx)/2, syk+(my+ny)/2, szk),
            (sxk+nx, syk+ny, szk),
            (sxk+(nx+px)/2, syk+(ny+py)/2, szk),
            (sxk+px, syk+py, szk),
            (sxk+px, syk+py, 0),
            (sxk+px, syk+py, -szk),
            (sxk+(nx+px)/2, syk+(ny+py)/2, -szk),
            (sxk+nx, syk+ny, -szk),
            (sxk+(mx+nx)/2, syk+(my+ny)/2, -szk),
            (sxk+mx, syk+my, -szk),
            (sxk+mx/2, syk+my/2, -szk),
            (sxk, syk, -szk),
          ))

          coords = R.elems * coords
          Romega = matrix.col(axes[2]).axis_and_angle_as_r3_rotation_matrix(angles[2], deg=True)
          coords = Romega.elems * coords
          extrema.extend(coords)

          return extrema

      #------------------ finished defining SmarGonShadowMaskGenerator
      return SmarGonShadowMaskGenerator(goniometer)

    elif goniometer.get_names()[1] == 'GON_KAPPA':
      # mini Kappa

      from dials.util.masking import GoniometerShadowMaskGenerator
      from scitbx.array_family import flex
      import math

      # Simple model of cone around goniometer phi axis
      # Exact values don't matter, only the ratio of height/radius
      height = 50 # mm
      radius = 20 # mm

      steps_per_degree = 1
      theta = flex.double([list(range(360*steps_per_degree))]) * math.pi/180 * 1/steps_per_degree
      y = radius * flex.cos(theta) # x
      z = radius * flex.sin(theta) # y
      x = flex.double(theta.size(), height) # z

      coords = flex.vec3_double(list(zip(x, y, z)))
      coords.insert(0, (0,0,0))

      if goniometer is None:
        goniometer = self.get_goniometer()
      return GoniometerShadowMaskGenerator(
        goniometer, coords, flex.size_t(len(coords), 0))

    else:
      raise RuntimeError(
        "Don't understand this goniometer: %s" %list(goniometer.get_names()))

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatCBFFullPilatus.understand(arg))
