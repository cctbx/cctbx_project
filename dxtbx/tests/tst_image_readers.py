from __future__ import division

def get_image_examples():
  from os.path import join
  import libtbx.load_env
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError, e:
    print 'SKIP: dials_regression not configured'
    exit(0)


  images = {
    "./image_examples/ALS_501/als501_q4_1_001.img"                        : "smv",
    "./image_examples/SPring8_BL26B1_SaturnA200/A200_000001.img"          : "smv",
    "./image_examples/SPring8_BL26B1_SaturnA200/A200_000002.img"          : "smv",
    "./image_examples/APS_19ID/q315_unbinned_a.0001.img"                  : "smv",
    "./image_examples/ALS_821/q210_lyso_1_101.img"                        : "smv",
    "./image_examples/MLFSOM_simulation/fake_00001.img"                   : "smv",
    "./image_examples/ALS_831/q315r_lyso_001.img"                         : "smv",
    "./image_examples/DESY_ID141/q210_2_001.img"                          : "smv",
    "./image_examples/SSRL_bl91/q315_1_001.img"                           : "smv",
    "./image_examples/APS_24IDC/q315_1_001.img"                           : "smv",
    "./image_examples/ALS_1231/q315r_lyso_1_001.img"                      : "smv",
    "./image_examples/SRS_142/q4_1_001.img"                               : "smv",
    "./image_examples/ALS_422/lyso_041013a_1_001.img"                     : "smv",
    "./image_examples/APS_17ID/q210_1_001.img"                            : "smv",
    "./image_examples/saturn/lyso_00001.img"                              : "smv",

    # "./image_examples/SPring8_BL26B1_Raxis5/raxis5_000091.img"            : "raxis",
    # "./image_examples/SPring8_BL26B1_Raxis5/raxis5_000001.img"            : "raxis",

    "./image_examples/SPring8_BL32XU_MX225HS/ds_000045.img"               : "tiff",
    "./image_examples/SPring8_BL32XU_MX225HS/ds_000001.img"               : "tiff",
    "./image_examples/SPring8_BL44XU_MX300HE/bl44xu_lys_000002.img"       : "tiff",
    "./image_examples/SPring8_BL44XU_MX300HE/bl44xu_lys_000001.img"       : "tiff",
    "./image_examples/SPring8_BL32XU/rayonix225hs_0001.img"               : "tiff",
    "./image_examples/SPring8_BL32XU/rayonix225_0001.img"                 : "tiff",
    "./image_examples/SLS_X06SA/mar225_2_001.img"                         : "tiff",
    "./image_examples/CLS1_08ID1/mar225_2_E0_0001.img"                    : "tiff",
    "./image_examples/SPring8_BL38B1_MX225HE/bl38b1_001.img"              : "tiff",
    "./image_examples/SPring8_BL38B1_MX225HE/bl38b1_090.img"              : "tiff",
    "./image_examples/SRS_101/mar225_001.img"                             : "tiff",
    "./image_examples/SPring8_BL12B2_MX225HE/lys001_000091.img"           : "tiff",
    "./image_examples/SPring8_BL12B2_MX225HE/lys001_000001.img"           : "tiff",
    "./image_examples/SPring8_BL26B2_MX225/2sec_Al200um_000001.img"       : "tiff",
    "./image_examples/SPring8_BL26B2_MX225/2sec_Al200um_000090.img"       : "tiff",

    "./stills_test_data/hit-20111202210224984.cbf"                        : "cbf_multitile", # Multitile
    "./stills_test_data/hit-s00-20140306002935980.cbf"                    : "cbf_multitile", # Multitile
    "./stills_test_data/hit-s00-20140306002857363.cbf"                    : "cbf_multitile", # Multitile

     "./image_examples/ESRF_ID29/trypsin_1_0001.cbf"                       : "cbf",
     "./image_examples/xia2/merge2cbf_averaged_0001.cbf"                   : "cbf",
     "./image_examples/dials-190/whatev1_01_00002.cbf"                     : "cbf",
     "./image_examples/dials-190/whatev1_02_00001.cbf"                     : "cbf",
     "./image_examples/dials-190/whatev1_02_00002.cbf"                     : "cbf",
     "./image_examples/dials-190/whatev1_03_00001.cbf"                     : "cbf",
     "./image_examples/dials-190/whatev1_01_00001.cbf"                     : "cbf",
     "./image_examples/dials-190/whatev1_03_00002.cbf"                     : "cbf",
     "./image_examples/APS_24IDE_test/thaum-12_1_0005.cbf"                 : "cbf",
     "./image_examples/APS_24IDE_test/thaum-12_1_0004.cbf"                 : "cbf",
     "./image_examples/APS_24IDE_test/thaum-12_1_0006.cbf"                 : "cbf",
     "./image_examples/APS_24IDE_test/thaum-12_1_0008.cbf"                 : "cbf",
     "./image_examples/APS_24IDE_test/thaum-12_1_0009.cbf"                 : "cbf",
     "./image_examples/APS_24IDE_test/thaum-12_1_0010.cbf"                 : "cbf",
     "./image_examples/APS_24IDE_test/thaum-12_1_0007.cbf"                 : "cbf",
     "./image_examples/APS_24IDE_test/thaum-12_1_0002.cbf"                 : "cbf",
     "./image_examples/APS_24IDE_test/thaum-12_1_0001.cbf"                 : "cbf",
     "./image_examples/APS_24IDE_test/thaum-12_1_0003.cbf"                 : "cbf",
     "./image_examples/APS_24IDC/pilatus_1_0001.cbf"                       : "cbf",
     "./image_examples/SLS_Eiger_16M_as_CBF/insu_with_bs_labelit_0901.cbf" : "cbf",
     "./image_examples/SLS_Eiger_16M_as_CBF/insu_with_bs_labelit_0001.cbf" : "cbf",
     "./image_examples/SPring8_BL41XU_PILATUS3_6M/data1_000001.cbf"        : "cbf",
     "./image_examples/SPring8_BL41XU_PILATUS3_6M/data1_000901.cbf"        : "cbf",
     "./image_examples/DLS_I02/X4_wide_M1S4_1_0001.cbf"                    : "cbf",
     "./image_examples/DLS_I23/I23_P12M_alpha_0001.cbf"                    : "cbf",
     "./image_examples/DLS_I23/germ_13KeV_0001.cbf"                        : "cbf",
     "./image_examples/SLS_X06SA/pilatus6m_1_00001.cbf"                    : "cbf",
     "./image_examples/SPring8_ADSC_SN916/Xtal17-2phi_3_015.cbf"           : "cbf",
     "./image_examples/DLS_I19/I19_P300k_00001.cbf"                        : "cbf",
     "./image_examples/ED_From_TIFF/170112330001.cbf"                      : "cbf",

    "./image_examples/putative_imgCIF_HDF5_mapping/minicbf.h5"             : "hdf5",
  }

  return dict((join(dials_regression, k), v) for k, v in images.iteritems())

IMAGE_EXAMPLES = get_image_examples()

def get_smv_header(image_file):
  header_size = int(open(image_file, 'rb').read(45).split(
      '\n')[1].split('=')[1].replace(';', '').strip())
  header_text = open(image_file, 'rb').read(header_size)
  header_dictionary = { }

  # Check that we have the whole header, contained within { }.  Stop
  # extracting data once a record solely composed of a closing curly
  # brace is seen.  If there is no such character in header_text
  # either HEADER_BYTES caused a short read of the header or the
  # header is malformed.
  for record in header_text.split('\n'):
    if record == '}':
      break
    if not '=' in record:
      continue

    key, value = record.replace(';', '').split('=')

    header_dictionary[key.strip()] = value.strip()

  return header_size, header_dictionary

def read_smv_image(image_file):
  from boost.python import streambuf
  from dxtbx import read_uint16, read_uint16_bs, is_big_endian
  from scitbx.array_family import flex

  header_size, header_dictionary = get_smv_header(image_file)

  f = open(image_file, 'rb')
  f.read(header_size)

  if header_dictionary['BYTE_ORDER'] == 'big_endian':
    big_endian = True
  else:
    big_endian = False

  image_size = (int(header_dictionary['SIZE1']),
                int(header_dictionary['SIZE2']))

  if big_endian == is_big_endian():
    raw_data = read_uint16(streambuf(f), int(image_size[0] * image_size[1]))
  else:
    raw_data = read_uint16_bs(streambuf(f), int(image_size[0] * image_size[1]))

  raw_data.reshape(flex.grid(image_size[1], image_size[0]))

  return raw_data

def get_tiff_header(image_file):

  from dxtbx.format.FormatTIFFHelpers import read_basic_tiff_header
  width, height, depth, header, order = read_basic_tiff_header(
      image_file)

  header_bytes = open(image_file, 'rb').read(header)

  return width, height, depth // 8, order, header_bytes

def read_tiff_image(image_file):
  # currently have no non-little-endian machines...

  from boost.python import streambuf
  from dxtbx import read_uint16
  from scitbx.array_family import flex

  width, height, depth, order, header_bytes = get_tiff_header(image_file)
  image_size = (width, height)
  header_size = 4096

  f = open(image_file)
  f.read(header_size)
  raw_data = read_uint16(streambuf(f), int(image_size[0] * image_size[1]))
  raw_data.reshape(flex.grid(image_size[1], image_size[0]))

  return raw_data

def read_cbf_image(cbf_image):
  from cbflib_adaptbx import uncompress
  import binascii

  start_tag = binascii.unhexlify('0c1a04d5')

  data = open(cbf_image, 'rb').read()
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



def read_multitile_cbf_image(cbf_image):
  import numpy
  from scitbx.array_family import flex
  import pycbf

  raw_data = []
  cbf = pycbf.cbf_handle_struct()
  cbf.read_widefile(cbf_image, pycbf.MSG_DIGEST)
  cbf.find_category('array_structure')
  cbf.find_column('encoding_type')
  cbf.select_row(0)
  types = []
  for i in xrange(cbf.count_rows()):
    types.append(cbf.get_value())
    cbf.next_row()
  assert len(types) == cbf.count_rows()

  # read the data
  data = {}
  cbf.find_category("array_data")
  for i in xrange(cbf.count_rows()):
    cbf.find_column("array_id")
    name = cbf.get_value()

    cbf.find_column("data")
    assert cbf.get_typeofvalue().find('bnry') > -1

    if types[i] == 'signed 32-bit integer':
      array_string = cbf.get_integerarray_as_string()
      array = flex.int(numpy.fromstring(array_string, numpy.int32))
      parameters = cbf.get_integerarrayparameters_wdims_fs()
      array_size = (parameters[11], parameters[10], parameters[9])
    elif types[i] == 'signed 64-bit real IEEE':
      array_string = cbf.get_realarray_as_string()
      array = flex.double(numpy.fromstring(array_string, numpy.float))
      parameters = cbf.get_realarrayparameters_wdims_fs()
      array_size = (parameters[7], parameters[6], parameters[5])
    else:
      return None # type not supported

    array.reshape(flex.grid(*array_size))
    data[name] = array
    cbf.next_row()

  # extract the data for each panel
  try:
    cbf.find_category("array_structure_list_section")
    has_sections = True
  except Exception:
    has_sections = False
  if has_sections:
    section_shapes = {}
    for i in xrange(cbf.count_rows()):
      cbf.find_column("id")
      section_name = cbf.get_value()
      if not section_name in section_shapes:
        section_shapes[section_name] = {}
      cbf.find_column("array_id")
      if not "array_id" in section_shapes[section_name]:
        section_shapes[section_name]["array_id"] = cbf.get_value()
      else:
        assert section_shapes[section_name]["array_id"] == cbf.get_value()
      cbf.find_column("index"); axis_index = int(cbf.get_value())-1
      cbf.find_column("start"); axis_start = int(cbf.get_value())-1
      cbf.find_column("end");   axis_end   = int(cbf.get_value())

      section_shapes[section_name][axis_index] = slice(axis_start, axis_end)
      cbf.next_row()

    for section_name in sorted(section_shapes):
      section_shape  = section_shapes[section_name]
      section = data[section_shape["array_id"]][ \
        section_shape[2], section_shape[1], section_shape[0]]
      section.reshape(flex.grid(section.focus()[-2], section.focus()[-1]))
      raw_data.append(section)
  else:
    for key in sorted(data):
      data[key].reshape(flex.grid(data[key].focus()[-2],data[key].focus()[-1]))
      raw_data.append(data[key])


  return tuple(raw_data)




def tst_smv(filename):
  from dxtbx.format.image import SMVReader
  from dxtbx.datablock import DataBlockFactory
  from scitbx.array_family import flex

  image = SMVReader(filename).image()
  assert image.n_tiles() == 1
  data1 = image.tile(0).as_int()

  data2 = read_smv_image(filename)

  diff = flex.abs(data1 - data2)
  assert flex.max(diff) < 1e-7

  print 'OK'



def tst_tiff(filename):
  from dxtbx.datablock import DataBlockFactory
  from scitbx.array_family import flex
  from dxtbx.format.image import TIFFReader


  image = TIFFReader(filename).image()
  assert image.n_tiles() == 1
  data1 = image.tile(0).as_int()

  data2 = read_tiff_image(filename)

  diff = flex.abs(data1 - data2)
  assert flex.max(diff) < 1e-7

  print 'OK'


def tst_cbf_fast(filename):
  from dxtbx.datablock import DataBlockFactory
  from scitbx.array_family import flex
  from dxtbx.format.image import CBFFastReader

  image = CBFFastReader(filename).image()
  assert image.n_tiles() == 1
  data1 = image.tile(0).as_int()

  data2 = read_cbf_image(filename)

  diff = flex.abs(data1 - data2)
  assert flex.max(diff) < 1e-7

  print 'OK'


def tst_cbf(filename):
  from dxtbx.datablock import DataBlockFactory
  from scitbx.array_family import flex
  from dxtbx.format.image import CBFReader

  image = CBFReader(filename).image()
  data1 = tuple(image.tile(i).as_int() for i in range(image.n_tiles()))

  try:
    data2 = read_multitile_cbf_image(filename)
  except Exception, e:
    data2 = (read_cbf_image(filename),)

  assert len(data1) == len(data2)

  for d1, d2 in zip(data1, data2):
    diff = flex.abs(d1 - d2)
    assert flex.max(diff) < 1e-7

  print 'OK'

def tst_hdf5(filename):
  from dxtbx.format.image import HDF5Reader
  from dxtbx.format.nexus import dataset_as_flex_int
  from scitbx.array_family import flex
  import h5py

  handle = h5py.File(filename, "r")

  reader = HDF5Reader(
    handle.id.id,
    flex.std_string(["/entry/data/data"]))

  image = reader.image(0)

  data1 = image.tile(0).as_int()

  dataset = handle['/entry/data/data']
  N, height, width = dataset.shape
  data2 = dataset_as_flex_int(
    dataset.id.id,
      (slice(0, 1, 1),
       slice(0, height, 1),
       slice(0, width, 1)))
  data2.reshape(flex.grid(data2.all()[1:]))

  assert N == len(reader)
  assert data1.all()[0] == data2.all()[0]
  assert data1.all()[1] == data2.all()[1]
  diff = flex.abs(data1 - data2)
  assert flex.max(diff) < 1e-7

  print 'OK'



def tst_all():
  for k, v in IMAGE_EXAMPLES.iteritems():
    #if v == 'smv':
    #  tst_smv(k)
    #elif v == 'tiff':
    #  tst_tiff(k)
    #elif v == 'cbf':
    #  tst_cbf_fast(k)
    #  tst_cbf(k)
    #elif v == 'cbf_multitile':
    #  tst_cbf(k)
    if v == "hdf5":
      tst_hdf5(k)

if __name__ == '__main__':

  tst_all()
