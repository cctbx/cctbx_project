#!/usr/bin/env python
# FormatCBFMiniPilatus.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the CBF image reader for Pilatus images. Inherits from
# FormatCBFMini.

from __future__ import division
import os

from dxtbx.format.FormatCBFMini import FormatCBFMini
from dxtbx.format.FormatCBFMiniPilatusHelpers import \
     get_pilatus_timestamp
from dxtbx.format.FormatPilatusHelpers import determine_pilatus_mask
from dxtbx.model import ParallaxCorrectedPxMmStrategy

if 'DXTBX_OVERLOAD_SCALE' in os.environ:
  dxtbx_overload_scale = float(os.environ['DXTBX_OVERLOAD_SCALE'])
else:
  dxtbx_overload_scale = 1

class FormatCBFMiniPilatus(FormatCBFMini):
  '''A class for reading mini CBF format Pilatus images, and correctly
  constructing a model for the experiment from this.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an Pilatus mini CBF format image,
    i.e. we can make sense of it.'''

    if 'ENABLE_PHOTON_FACTORY_TWO_EIGER' in os.environ:
      return False

    header = FormatCBFMini.get_cbf_header(image_file)

    for record in header.split('\n'):
      if '# Detector' in record and \
             'EIGER' in record.upper():
        return False

    for record in header.split('\n'):
      if '_array_data.header_convention' in record and \
             'PILATUS' in record:
        return True
      if '_array_data.header_convention' in record and \
             'SLS' in record:
        return True
      if '_array_data.header_convention' in record and \
             '?' in record:
        return True
      if '# Detector' in record and \
             'PILATUS' in record:  #CBFlib v0.8.0 allowed
        return True

    return False

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    assert(self.understand(image_file))

    FormatCBFMini.__init__(self, image_file, **kwargs)

    self._raw_data = None

    return

  def _start(self):
    FormatCBFMini._start(self)

  def _goniometer(self):
    '''Return a model for a simple single-axis goniometer. This should
    probably be checked against the image header, though for miniCBF
    there are limited options for this.'''

    if 'Phi' in self._cif_header_dictionary:
      phi_value = float(self._cif_header_dictionary['Phi'].split()[0])

    return self._goniometer_factory.single_axis()

  def _detector(self):
    '''Return a model for a simple detector, presuming no one has
    one of these on a two-theta stage. Assert that the beam centre is
    provided in the Mosflm coordinate frame.'''

    distance = float(
        self._cif_header_dictionary['Detector_distance'].split()[0])

    beam_xy = self._cif_header_dictionary['Beam_xy'].replace(
        '(', '').replace(')', '').replace(',', '').split()[:2]

    wavelength = float(
        self._cif_header_dictionary['Wavelength'].split()[0])

    beam_x, beam_y = map(float, beam_xy)

    pixel_xy = self._cif_header_dictionary['Pixel_size'].replace(
        'm', '').replace('x', '').split()

    pixel_x, pixel_y = map(float, pixel_xy)

    thickness = float(
      self._cif_header_dictionary['Silicon'].split()[2]) * 1000.0

    nx = int(
        self._cif_header_dictionary['X-Binary-Size-Fastest-Dimension'])
    ny = int(
        self._cif_header_dictionary['X-Binary-Size-Second-Dimension'])

    overload = dxtbx_overload_scale * int(
        self._cif_header_dictionary['Count_cutoff'].split()[0])
    underload = -1

    # take into consideration here the thickness of the sensor also the
    # wavelength of the radiation (which we have in the same file...)
    from cctbx.eltbx import attenuation_coefficient
    table = attenuation_coefficient.get_table("Si")
    mu = table.mu_at_angstrom(wavelength) / 10.0
    t0 = thickness

    detector = self._detector_factory.simple(
        'PAD', distance * 1000.0, (beam_x * pixel_x * 1000.0,
                                   beam_y * pixel_y * 1000.0), '+x', '-y',
        (1000 * pixel_x, 1000 * pixel_y),
        (nx, ny), (underload, overload), [],
        ParallaxCorrectedPxMmStrategy(mu, t0))

    for f0, s0, f1, s1 in determine_pilatus_mask(detector):
      detector[0].add_mask(f0, s0, f1, s1)

    detector[0].set_thickness(thickness)
    detector[0].set_material('Si')
    detector[0].set_mu(mu)

    return detector

  def _beam(self):
    '''Return a simple model for the beam.'''

    wavelength = float(
        self._cif_header_dictionary['Wavelength'].split()[0])

    return self._beam_factory.simple(wavelength)

  def _scan(self):
    '''Return the scan information for this image.'''

    format = self._scan_factory.format('CBF')

    exposure_time = float(
        self._cif_header_dictionary['Exposure_period'].split()[0])

    osc_start = float(
        self._cif_header_dictionary['Start_angle'].split()[0])
    osc_range = float(
        self._cif_header_dictionary['Angle_increment'].split()[0])

    timestamp = get_pilatus_timestamp(
        self._cif_header_dictionary['timestamp'])

    return self._scan_factory.single(
        self._image_file, format, exposure_time,
        osc_start, osc_range, timestamp)

  def read_cbf_image(self, cbf_image):
    from cbflib_adaptbx import uncompress
    import binascii
    from scitbx.array_family import flex

    start_tag = binascii.unhexlify('0c1a04d5')

    data = self.open_file(cbf_image, 'rb').read()
    data_offset = data.find(start_tag) + 4
    cbf_header = data[:data_offset - 4]

    fast = 0
    slow = 0
    length = 0

    for record in cbf_header.split('\n'):
      if 'X-Binary-Size-Fastest-Dimension' in record:
        fast = int(record.split()[-1])
      elif 'X-Binary-Size-Second-Dimension' in record:
        slow = int(record.split()[-1])
      elif 'X-Binary-Number-of-Elements' in record:
        length = int(record.split()[-1])
      elif 'X-Binary-Size:' in record:
        size = int(record.split()[-1])

    assert(length == fast * slow)

    pixel_values = uncompress(packed = data[data_offset:data_offset + size],
                              fast = fast, slow = slow)

    return pixel_values

  def get_raw_data(self):
    if self._raw_data is None:
      data = self.read_cbf_image(self._image_file)
      self._raw_data = data

    return self._raw_data

  def get_mask(self, goniometer=None):
    from scitbx.array_family import flex
    detector = self.get_detector()
    mask = [flex.bool(flex.grid(reversed(p.get_image_size())), True)
             for p in detector]
    for i, p in enumerate(detector):
      untrusted_regions = p.get_mask()
      for j, (f0, s0, f1, s1) in enumerate(untrusted_regions):
        sub_array = flex.bool(flex.grid(s1-s0+1, f1-f0+1), False)
        mask[i].matrix_paste_block_in_place(sub_array, s0-1, f0-1)

    if len(detector) == 1:
      raw_data = [self.get_raw_data()]
    else:
      raw_data = self.get_raw_data()
      assert(len(raw_data) ==  len(detector))
    trusted_mask = [
      p.get_trusted_range_mask(im) for im, p in zip(raw_data, detector)]

    # returns merged untrusted pixels and active areas using bitwise AND (pixels are accepted
    # if they are inside of the active areas AND inside of the trusted range)
    return tuple([m & tm for m, tm in zip(mask, trusted_mask)])

  def detectorbase_start(self):

    from iotbx.detectors.pilatus_minicbf import PilatusImage
    self.detectorbase = PilatusImage(self._image_file)
    self.detectorbase.readHeader() # necessary for LABELIT

  @staticmethod
  def as_file(detector,beam,gonio,scan,data,path):
    import pycbf
    from dxtbx.format.FormatCBFMultiTile import cbf_wrapper

    cbf=cbf_wrapper()
    cbf_root=os.path.splitext(os.path.basename(path))[0]+".cbf"
    cbf.new_datablock(os.path.splitext(os.path.basename(path))[0])

    """ Data items in the ARRAY_DATA category are the containers for the array data
    items described in the category ARRAY_STRUCTURE. """
    cbf.add_category("array_data",["header_convention","header_contents","data"])

    "also need to add invalid pixels on sensor boundaries"

    panel = detector[0]
    pixel_xy = panel.get_pixel_size()
    pixel_x, pixel_y = map(float, pixel_xy)
    pixel_x_microns = pixel_x*1000
    pixel_y_microns = pixel_y*1000

    thickness = panel.get_thickness()

    exposure_period = scan.get_exposure_times()[0]
    phi_start = scan.get_oscillation()[0]
    osc = scan.get_oscillation()[1]

    # exposure_time = exposure_period+0.00203
    exposure_time = exposure_period  # simulation is a perfect detector

    tau = 0 # simulation is a perfect detector
    count_cutoff = 2**20  # not actually sure what this is

    wavelength = beam.get_wavelength()  # get the wavelength in the conventional way
    energy = 12398.4245/wavelength
    threshold = energy/2  # presume normal data collection convention

    bad_pixels = 0   # maybe get this from negative pixel values?

    assert len(detector)==1,"only single-panel detectors supported"
    distance_meters = detector[0].get_distance()/1000

    beam_center = detector[0].get_beam_centre_px(beam.get_s0())
    ORGX, ORGY = map(float, beam_center)

    cbf.add_row(["PILATUS_1.2", """
# Detector: PILATUS3 6M, S/N 60-0000
# 1972-01-01T00:00:00.000
# Pixel_size %(pixel_x_microns).0fe-6 m x %(pixel_x_microns).0fe-6 m
# Silicon sensor, thickness %(thickness)f m
# Exposure_time %(exposure_time)f s
# Exposure_period %(exposure_period)f s
# Tau = %(tau)f s
# Count_cutoff %(count_cutoff)d counts
# Threshold_setting: %(threshold)d eV
# Gain_setting: autog (vrf = 1.000)
# N_excluded_pixels = %(bad_pixels)d
# Excluded_pixels: badpix_mask.tif
# Flat_field: (nil)
# Trim_file: p6m0000_E%(energy)d_T%(threshold)d_vrf_m0p100.bin
# Image_path: /ramdisk/
# Wavelength %(wavelength).5f A
# Detector_distance %(distance_meters).5f m
# Beam_xy (%(ORGX).2f, %(ORGY).2f) pixels
# Start_angle %(phi_start).4f deg.
# Angle_increment %(osc).4f deg.
# Detector_2theta 0.0000 deg.
"""%locals()
])

    binary_id = 1
    focus = data.focus()
    data2 = data.copy_to_byte_str()
    elements = len(data)
    byteorder = "little_endian"
    dimfast = focus[1]
    dimmid = focus[0]
    dimslow = 1
    padding = 0

    elsize = 4
    elsigned = 1

    cbf.set_integerarray_wdims_fs(
        pycbf.CBF_PACKED,
        binary_id,
        data2,
        elsize,
        elsigned,
        elements,
        byteorder,
        dimfast,
        dimmid,
        dimslow,
        padding)

    cbf.write_widefile(cbf_root,pycbf.CBF,\
      pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K,0)

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatCBFMiniPilatus.understand(arg)
