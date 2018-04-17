from __future__ import division

from libtbx import easy_pickle
from dials.array_family import flex
from scitbx import matrix
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
    s0 = self.beam.get_unit_s0()
    self.panel_s0_intersections = flex.vec2_double(
      [self.det[i].get_ray_intersection_px(s0) for i in xrange(len(self.det))])
  def get_one_spot_length(self, id):
    # the radial direction is along the vector from the beam center to
    # the spot centroid for each spot
    shoebox = self.strong['shoebox'][id]
    s0_position = self.panel_s0_intersections[shoebox.panel]
    centroid = shoebox.centroid_strong().px_xy
    dx, dy = [centroid[i] - s0_position[i] for i in (0, 1)]
    radial = matrix.col((dx, dy)).normalize()
    mask = flex.bool([(m & self.valid_code) != 0 for m in shoebox.mask])
    bbox = self.strong['bbox'][id]
    x_start, y_start = bbox[0], bbox[2]
    x_range, y_range = bbox[1] - bbox[0], bbox[3] - bbox[2]
    radial_distances = flex.double()
    for i, valid_foreground in enumerate(mask):
      if valid_foreground:
        position = matrix.col((x_start + i%x_range, y_start + i//y_range))
        projection = position.dot(radial)
        radial_distances.append(projection)
    length = flex.max(radial_distances) - flex.min(radial_distances)
    return length
  def get_spot_lengths_px(self):
    return flex.double([self.get_one_spot_length(id) for id in xrange(len(self.strong))])

class ReflectionsRadialLengthsFromFiles(ReflectionsRadialLengths):
  def __init__(self, *files):
    from dials.util.options import Importer, flatten_reflections, flatten_experiments, flatten_datablocks
    importer = Importer(list(files), read_experiments=True, read_datablocks=True,
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
