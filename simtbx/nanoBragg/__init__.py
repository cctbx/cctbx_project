from __future__ import absolute_import, division, print_function
import boost.python
import cctbx.uctbx # possibly implicit
ext = boost.python.import_ext("simtbx_nanoBragg_ext")
from simtbx_nanoBragg_ext import *

@boost.python.inject_into(ext.nanoBragg)
class _():

  def __getattr__(self,name):
    """assemble miller array of structure factors used to compute spot intensities from the internal C cube array
       how do we specify docstrings for individual overriden members? """
    if name == "Fhkl":
      from cctbx.crystal import symmetry
      cs = symmetry(unit_cell = self.unit_cell_Adeg,space_group="P 1")
      from cctbx.miller import set, array
      indices,data = self.Fhkl_tuple
      mset = set(crystal_symmetry=cs, anomalous_flag=True, indices=indices)
      return array(mset, data=data).set_observation_type_xray_amplitude()

  def __setattr__(self,name,value):
    """use a P1 anomalous=True miller array to initialize the internal C cube array with structure factors for the spot intensities
       how do we specify docstrings for individual overriden members? """
    if name in ["Fhkl"]:
      value=value.expand_to_p1()
      value=value.generate_bijvoet_mates()
      assert value.space_group_info().type().lookup_symbol() == "P 1"
      # handle exception by expanding to P1
      assert value.anomalous_flag() == True
      # handle exception by copying all F(hkl) to F(-h-k-l)
      #assert values are amplitudes # not sure how to guarantee this
      self.unit_cell_Adeg = value.unit_cell()
      #self.mock_up_group = value.space_group()
      #self.mock_up_anomalous_flag = value.anomalous_flag()
      self.Fhkl_tuple = (value.indices(),value.data())
    else:
      super(ext.nanoBragg,self).__setattr__(name,value)

  def to_smv_format_py(self,fileout,intfile_scale=0.0,debug_x=-1,debug_y=-1,
    rotmat=False,extra=None,verbose=False,gz=False):

    byte_order = "little_endian";

    #recast the image file write to Python to afford extra options: rotmat, extra, gz
    if gz:
      from libtbx.smart_open import for_writing
      outfile = for_writing(file_name=fileout+".gz", gzip_mode="wb")
    else:
      outfile = open(fileout,"w");

    outfile.write(("{\nHEADER_BYTES=1024;\nDIM=2;\nBYTE_ORDER=%s;\nTYPE=unsigned_short;\n"%byte_order).encode());
    outfile.write(b"SIZE1=%d;\nSIZE2=%d;\nPIXEL_SIZE=%g;\nDISTANCE=%g;\n"%(
      self.detpixels_fastslow[0],self.detpixels_fastslow[1],self.pixel_size_mm,self.distance_mm));
    outfile.write(b"WAVELENGTH=%g;\n"%self.wavelength_A);
    outfile.write(b"BEAM_CENTER_X=%g;\nBEAM_CENTER_Y=%g;\n"%self.beam_center_mm);
    outfile.write(b"ADXV_CENTER_X=%g;\nADXV_CENTER_Y=%g;\n"%self.adxv_beam_center_mm);
    outfile.write(b"MOSFLM_CENTER_X=%g;\nMOSFLM_CENTER_Y=%g;\n"%self.mosflm_beam_center_mm);
    outfile.write(b"DENZO_X_BEAM=%g;\nDENZO_Y_BEAM=%g;\n"%self.denzo_beam_center_mm);
    outfile.write(b"DIALS_ORIGIN=%g,%g,%g\n"%self.dials_origin_mm);
    outfile.write(b"XDS_ORGX=%g;\nXDS_ORGY=%g;\n"%self.XDS_ORGXY);
    outfile.write(b"CLOSE_DISTANCE=%g;\n"%self.close_distance_mm);
    outfile.write(b"PHI=%g;\nOSC_START=%g;\nOSC_RANGE=%g;\n"%(self.phi_deg,self.phi_deg,self.osc_deg));
    outfile.write(b"TIME=%g;\n"%self.exposure_s);
    outfile.write(b"TWOTHETA=%g;\n"%self.detector_twotheta_deg);
    outfile.write(b"DETECTOR_SN=000;\n");
    outfile.write(b"ADC_OFFSET=%g;\n"%self.adc_offset_adu);
    outfile.write(b"BEAMLINE=fake;\n");
    if rotmat:
      from scitbx.matrix import sqr
      RSABC = sqr(self.Amatrix).inverse().transpose()
      outfile.write( ("DIRECT_SPACE_ABC=%s;\n"%(",".join([repr(a) for a in RSABC.elems]))).encode() )
    if extra is not None:
      outfile.write(extra.encode())
    outfile.write(b"}\f");
    assert outfile.tell() < 1024, "SMV header too long, please edit this code and ask for more bytes."
    while ( outfile.tell() < 1024 ): outfile.write(b" ")
    from six import PY3
    if PY3:
      # Python3-compatible method for populating the output buffer.
      # Py2 implementation is more elegant in that the streambuf may be passed to C++,
      #   and the data are gzipped in chunks (default 1024). Py3 will not accept this method
      #   as it is PyString-based, with no converter mechanisms to bring data into PyBytes.
      # The Py3 method brings the full data in one chunk into PyBytes and then populates
      #   the output buffer in Python rather than C++.
      image_bytes = self.raw_pixels_unsigned_short_as_python_bytes(intfile_scale,debug_x,debug_y)
      outfile.write(image_bytes)
      outfile.close();
      return
    from boost.python import streambuf
    self.to_smv_format_streambuf(streambuf(outfile),intfile_scale,debug_x,debug_y)

    outfile.close();
