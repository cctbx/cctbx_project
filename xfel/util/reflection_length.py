from __future__ import division

from libtbx import easy_pickle
from dials.array_family import flex
from scitbx import matrix
import math
from dials.algorithms.shoebox import MaskCode

class ReflectionsRadialLengths(object):
  """Compute the length of each spotfinder spot in the radial direction."""
  valid_code = MaskCode.Valid | MaskCode.Foreground
  def __init__(self, strong_spots, experiment=None, datablock=None):
    self.strong = strong_spots
    assert 'bbox' in self.strong and 'shoebox' in self.strong, \
      "Spotfinder shoeboxes are required for spot length calculations."
    assert experiment or datablock, \
      "Supply one experiment or datablock object."
    if datablock:
      imgset = datablock._imagesets[0]
      self.det = imgset.get_detector()
      self.beam = imgset.get_beam()
    else:
      self.det = experiment.detector
      self.beam = experiment.beam
    self.s0 = self.beam.get_unit_s0()
    self.panel_s0_intersections = flex.vec2_double(
      [self.det[i].get_ray_intersection_px(self.s0) for i in xrange(len(self.det))])
  def get_one_spot_length_width_angle(self, id):
    # the radial direction is along the vector from the beam center to
    # the spot centroid for each spot
    shoebox = self.strong['shoebox'][id]
    s0_position = self.panel_s0_intersections[shoebox.panel]
    centroid = shoebox.centroid_strong().px_xy
    dx, dy = [centroid[i] - s0_position[i] for i in (0, 1)]
    radial_abs = matrix.col((dx, dy))
    radial = radial_abs.normalize()
    transverse = self.s0.normalize().cross(radial)
    mask = flex.bool([(m & self.valid_code) != 0 for m in shoebox.mask])
    bbox = self.strong['bbox'][id]
    x_start, y_start = bbox[0], bbox[2]
    x_range, y_range = bbox[1] - bbox[0], bbox[3] - bbox[2]
    radial_distances = flex.double()
    transverse_distances = flex.double()
    for i, valid_foreground in enumerate(mask):
      if valid_foreground:
        position = matrix.col((x_start + i%x_range, y_start + i//y_range))
        projection_radial = position.dot(radial)
        projection_transverse = position.dot(transverse)
        radial_distances.append(projection_radial)
        transverse_distances.append(projection_transverse)
    length = flex.max(radial_distances) - flex.min(radial_distances)
    width = flex.max(transverse_distances) - flex.min(transverse_distances)
    # The angle subtended is centered at the spot centroid, spanning the
    # spot width. Half this angle makes a right triangle with legs of lengths
    # [distance to beam center] and [half the spot width]. Use the tangent.
    angle = 2*math.atan(width/(2*radial_abs.length()))
    return (length, width, angle)
  def get_spot_lengths_px(self):
    self.lengths, self.widths, self.angles = \
      [self.get_one_spot_length_width_angle(id) for id in xrange(len(self.strong))]
    return flex.double(self.lengths)
  def get_spot_width(self):
    if not hasattr(self, "widths"):
      self.get_spot_lengths_px()
    return flex.double(self.widths)
  def get_spot_subtended_angles_deg(self):
    if not hasattr(self, "angles"):
      self.get_spot_lengths_px()
    return flex.double(self.angles)*180/math.pi
  def get_intensities(self):
    return self.strong['intensity.sum.value']

class ReflectionsRadialLengthsFromFiles(ReflectionsRadialLengths):
  def __init__(self, files):
    from dials.util.options import Importer, flatten_reflections, flatten_experiments, flatten_datablocks
    importer = Importer(files, read_experiments=True, read_datablocks=True,
      read_reflections=True, check_format=False)
    if importer.unhandled:
      print "Unable to handle one or more files:", importer.unhandled
      return
    reflections = flatten_reflections(importer.reflections)
    assert len(reflections) == 1, "Implemented only for one reflection table at a time presently"
    datablock = None
    experiment = None
    if importer.experiments:
      experiments = flatten_experiments(importer.experiments)
      assert len(experiments) == 1, "Implemented only for one experiment at a time presently"
      experiment = experiments[0]
    if importer.datablocks:
      datablocks  = flatten_datablocks(importer.datablocks)
      assert len(datablocks) == 1, "Implemented only for one datablock at a time presently"
      datablock = datablocks[0]
    super(ReflectionsRadialLengthsFromFiles, self).__init__(
      reflections[0], datablock=datablock, experiment=experiment)

if __name__ == "__main__":
  import sys
  assert len(sys.argv) == 3
  strong_spot_lengths = ReflectionsRadialLengthsFromFiles(sys.argv[1:]).get_spot_lengths_px()
  easy_pickle.dump("spot_lengths_px.pickle", strong_spot_lengths)
  print list(strong_spot_lengths)
