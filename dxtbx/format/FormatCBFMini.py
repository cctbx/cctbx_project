#!/usr/bin/env python
# FormatCBFMini.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Base implementation of miniCBF format - as used with Dectris detectors -
# this will read the header and populate a dictionary of the keyword / value
# pairs.

from __future__ import absolute_import, division
from __future__ import print_function
from builtins import map
import os

from dxtbx.format.FormatCBF import FormatCBF
from dxtbx.model import ParallaxCorrectedPxMmStrategy

if 'DXTBX_OVERLOAD_SCALE' in os.environ:
  dxtbx_overload_scale = float(os.environ['DXTBX_OVERLOAD_SCALE'])
else:
  dxtbx_overload_scale = 1

class FormatCBFMini(FormatCBF):
  '''An image reading class for mini CBF format images i.e. those from
  Dectris, which will read the header into a dictionary.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an CBF format image, i.e. we can
    make sense of it.'''

    header = FormatCBF.get_cbf_header(image_file)

    if '_diffrn.id' in header and '_diffrn_source' in header:
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
      if '_array_data.header_convention' in record and \
             'XDS special' in record:
        return True
      if '_array_data.header_convention' in record and \
             'GENERIC_MINI' in record: # intended for simulated PAD data, non-Pilatus array size
        return True
      if '# Detector' in record and \
             'PILATUS' in record:  #CBFlib v0.8.0 allowed
        return True
      if '# Detector' in record and \
             'ADSC' in record and 'HF-4M' in header:
        return True

    return False

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file.'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    FormatCBF.__init__(self, image_file, **kwargs)

    self._raw_data = None

    return

  def _start(self):
    '''Open the image file, read the image header, copy it into a
    dictionary for future reference.'''

    FormatCBF._start(self)

    cif_header = FormatCBF.get_cbf_header(self._image_file)

    self._cif_header_dictionary = { }

    for record in cif_header.split('\n'):
      if not '#' in record[:1]:
        continue

      if len(record[1:].split()) <= 2 and record.count(':') == 2:
        self._cif_header_dictionary['timestamp'] = record[1:].strip()
        continue

      tokens = record.replace('=', '').replace(':', '').split()[1:]

      self._cif_header_dictionary[tokens[0]] = ' '.join(tokens[1:])

    for record in self._mime_header.split('\n'):
      if not record.strip():
        continue
      token, value = record.split(':')
      self._cif_header_dictionary[token.strip()] = value.strip()

    return

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

    beam_x, beam_y = list(map(float, beam_xy))

    pixel_xy = self._cif_header_dictionary['Pixel_size'].replace(
        'm', '').replace('x', '').split()

    pixel_x, pixel_y = list(map(float, pixel_xy))

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

    detector[0].set_thickness(thickness)
    detector[0].set_material('Si')
    detector[0].set_mu(mu)

    return detector

  def _goniometer(self):
    '''Return a model for a simple single-axis goniometer. This should
    probably be checked against the image header, though for miniCBF
    there are limited options for this.'''

    if 'Phi' in self._cif_header_dictionary:
      phi_value = float(self._cif_header_dictionary['Phi'].split()[0])

    return self._goniometer_factory.single_axis()

  def _beam(self):
    '''Return a simple model for the beam.'''

    wavelength = float(
        self._cif_header_dictionary['Wavelength'].split()[0])

    beam = self._beam_factory.simple(wavelength)

    try:
      flux = float(self._cif_header_dictionary['Flux'].split()[0])
      beam.set_flux(flux)
    except KeyError:
      pass

    try:
      transmission = float(self._cif_header_dictionary['Transmission'].split()[0])
      beam.set_transmission(transmission)
    except KeyError:
      pass

    return beam

  def read_cbf_image(self, cbf_image):
    from cbflib_adaptbx import uncompress
    import binascii

    start_tag = binascii.unhexlify('0c1a04d5')

    data = self.open_file(cbf_image, 'rb').read()
    data_offset = data.find(start_tag) + 4
    cbf_header = data[:data_offset - 4]

    fast = 0
    slow = 0
    length = 0
    byte_offset = False
    no_compression = False

    for record in cbf_header.split('\n'):
      if 'X-Binary-Size-Fastest-Dimension' in record:
        fast = int(record.split()[-1])
      elif 'X-Binary-Size-Second-Dimension' in record:
        slow = int(record.split()[-1])
      elif 'X-Binary-Number-of-Elements' in record:
        length = int(record.split()[-1])
      elif 'X-Binary-Size:' in record:
        size = int(record.split()[-1])
      elif 'conversions' in record:
        if 'x-CBF_BYTE_OFFSET' in record:
          byte_offset = True
        elif 'x-CBF_NONE' in record:
          no_compression = True

    assert(length == fast * slow)

    if byte_offset:
      pixel_values = uncompress(packed = data[data_offset:data_offset + size],
                                fast = fast, slow = slow)
    elif no_compression:
      from boost.python import streambuf
      from dxtbx import read_int32
      from scitbx.array_family import flex
      assert(len(self.get_detector()) == 1)
      f = self.open_file(self._image_file)
      f.read(data_offset)
      pixel_values = read_int32(streambuf(f), int(slow * fast))
      pixel_values.reshape(flex.grid(slow, fast))

    else:
      raise ValueError("Uncompression of type other than byte_offset or none is not supported (contact authors)")

    return pixel_values

  def get_raw_data(self):
    if self._raw_data is None:
      data = self.read_cbf_image(self._image_file)
      self._raw_data = data

    return self._raw_data

  def detectorbase_start(self):

    from iotbx.detectors.pilatus_minicbf import PilatusImage
    self.detectorbase = PilatusImage(self._image_file)
    self.detectorbase.readHeader() # necessary for LABELIT

  @staticmethod
  def as_file(detector,beam,gonio,scan,data,path,header_convention="GENERIC_MINI",det_type="GENERIC"):
    """Note to developers: first attempt to write a miniCBF given a dxtbx-style experiment,
       But fields are not filled rigorously as in cbflib/src/cbf_minicbf_header.c
       Present code does not account for:
         Pilatus model number and serial number
         Data collection date and time
         Sensor material
         Dead time (exposure time - exposure period)
         Highest trusted value
         Energy threshold for photon counting
         Gain
         Bad pixels
         Auxiliary files (bad pixel mask, flat field, trim, image path)
         Detector not normal to beam
    """
    import pycbf
    from dxtbx.format.FormatCBFMultiTile import cbf_wrapper

    cbf=cbf_wrapper()
    cbf_root=os.path.splitext(os.path.basename(path))[0]+".cbf"
    cbf.new_datablock(os.path.splitext(os.path.basename(path))[0])

    """ Data items in the ARRAY_DATA category are the containers for the array data
    items described in the category ARRAY_STRUCTURE. """
    cbf.add_category("array_data",["header_convention","header_contents","data"])

    panel = detector[0]
    pixel_xy = panel.get_pixel_size()
    pixel_x_microns = pixel_xy[0]*1000
    pixel_y_microns = pixel_xy[1]*1000

    thickness = max(0.000001,panel.get_thickness())

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

    # why the heck is this backwards?
    #print detector[0].get_fast_axis()
    #print detector[0].get_slow_axis()
    beam_center = detector[0].get_beam_centre_px(beam.get_s0())
    ORGY, ORGX = beam_center

    cbf.add_row([header_convention, """
# Detector: %(det_type)s, S/N 60-0000
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
# Trim_directory:
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
        pycbf.CBF_BYTE_OFFSET,
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
    print(FormatCBFMini.understand(arg))
