from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("rstbx_diffraction_fastbragg_ext")
from rstbx_diffraction_fastbragg_ext import *

@bp.inject_into(ext.detector)
class _():

  def show(self):
    print("Detector with %d x %d square pixels each with %5.3f mm edge"%(
       self.raw.focus()[0],self.raw.focus()[1],self.pixel_sz*1000.))

class standard_camera(camera):

  # Derived class; assumes camera is perpendicular to beam with
  # beam exactly hitting the center.
  def __init__(self,detector,mean_xtd_distance_m,mean_xray_wavelength_m):
    camera.__init__(self)
    self.distance = mean_xtd_distance_m
    self.lambda0 = mean_xray_wavelength_m
    self.Ybeam = detector.raw.focus()[0] * detector.pixel_sz/2.
    self.Zbeam = detector.raw.focus()[1] * detector.pixel_sz/2.

  def corner_resolution(self):
    import math
    corner_distance = self.Ybeam*math.sqrt(2.) # again, assumes a square
    theta = math.atan2(corner_distance,self.distance)/2.0;
    return self.lambda0/(2.0 * math.sin(theta));

  def show(self):
    import math
    print("Camera at sample distance %5.3f mm, wavelength %5.3f Angstroms"%(
    self.distance*1000.,self.lambda0*1.E10))

    assert self.Ybeam==self.Zbeam # everything assumes the detector is square.

    edge_distance = self.Ybeam
    theta = math.atan2(edge_distance,self.distance)/2.0;
    resolution = 1.E10*self.lambda0/(2.0 * math.sin(theta));
    print("The detector edge is at %5.3f Angstroms,"%resolution)

    print("   and the corner is at %5.3f Angstroms."%(self.corner_resolution()*1.E10))

class crystal_structure(crystal):

  def __init__(self,standard_camera,pdb_code=None,pdb_file=None):
    self.limiting_resolution = 1.E10*standard_camera.corner_resolution() # in Angstroms
    self.pdb_code = pdb_code
    from iotbx import pdb
    if pdb_file != None:
      pdb_inp = pdb.input(file_name = pdb_file)
    elif pdb_code != None:
      import iotbx.pdb.fetch
      pdb_url = iotbx.pdb.fetch.fetch(id=pdb_code)
      pdb_inp = pdb.input(source_info=None,lines=pdb_url.readlines())
    self.xray_structure = pdb_inp.xray_structure_simple()
    primitive_xray_structure = self.xray_structure.primitive_setting()
    self.P1_primitive_xray_structure = primitive_xray_structure.expand_to_p1()

  def p1_amplitudes(self,resolution,algorithm="fft"):
    fcalc = self.P1_primitive_xray_structure.structure_factors(
      d_min=resolution, anomalous_flag=True, algorithm=algorithm).f_calc()
    return fcalc.amplitudes()

  def p1_intensities(self):
    return self.p1_amplitudes(resolution=self.limiting_resolution
           ).as_intensity_array()

  def cell(self):
    return self.xray_structure.unit_cell()

  def symmetry(self):
    return self.xray_structure.crystal_symmetry()

  def p1_cell(self):
    return self.P1_primitive_xray_structure.unit_cell()

  def show(self):
    print("Structure with primitive unit cell %s diffracting out to %4.2f Angstrom"%(
      self.p1_cell(),self.limiting_resolution))
