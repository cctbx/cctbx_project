from __future__ import absolute_import, division, print_function
from spotfinder.array_family import flex
from spotfinder.core_toolbox import Distl


class heuristics_base(object):

  def __init__(self,phil_params):
    self.phil_params = phil_params
    self.pd = dict(resolution_inspection='100.0')
    self.NspotMin = self.phil_params.distl_minimum_number_spots_for_indexing
    self.NspotMax = self.phil_params.distl_maximum_number_spots_for_indexing
    self.BinMin = 25
    self.errormessage = None
    self.images = {}
    self.reporters = {}
    self.protocol = 'tnear2'
    self.overlapping = False #flag indicates whether the special procedure was used
    self.force_detail = False #flag indicates whether percent_overlap > force_detail cutoff

  def register_frames(self):
    nimages = len(self.phil_params.distl.image)

    return self.oneImage(0)

  def oneImage(self,framenumber):
    self.reporters[framenumber] = []

    import dxtbx.format.Registry
    filename = self.phil_params.distl.image[framenumber]
    reader = dxtbx.format.Registry.get_format_class_for_file(filename)
    img = reader(filename)

    detector = img.get_detector()
    beam = img.get_beam()
    S0 = beam.get_s0()
    data = img.get_raw_data()
    scan = img.get_scan()
    print(scan)
    if scan is None:
      print("No scan")
      RR = (0,1)
    else:
      print(scan.get_oscillation())
      RR = scan.get_oscillation_range()

    from spotfinder.dxtbx_toolbox import Distl

    sfall = Distl(params = self.phil_params, detector = detector, beam = beam, data = data)

    resolutions = flex.double()
    spotlist = []
    from dials.model.data import ReflectionList,Reflection
    reflections = ReflectionList()


    for ip,panel in enumerate(detector):
      for spot in sfall.finderlist[ip].spots:
        resolutions.append( panel.get_resolution_at_pixel(S0, (spot.ctr_mass_x(), spot.ctr_mass_y())) )
        spotlist.append(spot)
        refl = Reflection()
        refl.panel_number = ip
        refl.centroid_position = (spot.ctr_mass_x(), spot.ctr_mass_y(),0.0)
        refl.centroid_variance = (0.5,0.5,0.0)
        reflections.append(refl)


    selection = (resolutions>0.0)
    if self.phil_params.distl.res.outer is not None:
      selection = (selection and (resolutions>self.phil_params.distl.res.outer))
    if self.phil_params.distl.res.inner is not None:
      selection = (selection and (resolutions<self.phil_params.distl.res.inner))

    reflections = reflections.select(selection)

    return dict(detector=detector, beam=beam, reflections=reflections, scan = scan,
                gonio = img.get_goniometer())
