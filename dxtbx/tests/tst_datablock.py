from __future__ import absolute_import, division
import libtbx.load_env
from dxtbx.datablock import DataBlockFactory
from dxtbx.imageset import ImageSweep
from os.path import join

class Test(object):

  def __init__(self):
    import os
    dials_regression = libtbx.env.dist_path('dials_regression')

    self.dials_regression = dials_regression
    self.centroid_test_data = os.path.join(dials_regression, 'centroid_test_data')
    self.image_examples = os.path.join(dials_regression, 'image_examples')

  def single_sweep_filenames(self):
    path = self.centroid_test_data
    filenames = []
    image_indices = range(1, 10)
    for i in image_indices:
      filenames.append(join(path, 'centroid_000{0}.cbf'.format(i)))
    return filenames

  def multiple_sweep_filenames(self):
    path = self.centroid_test_data
    filenames = []
    image_indices = list(range(1, 4)) + list(range(7, 10))
    for i in image_indices:
      filenames.append(join(path, 'centroid_000{0}.cbf'.format(i)))
    return filenames

  def all_image_examples(self):
    path = self.image_examples
    filenames = [
        ('ALS_1231', 'q315r_lyso_1_001.img'),
        ('ALS_501', 'als501_q4_1_001.img'),
        ('ALS_821', 'q210_lyso_1_101.img'),
        ('ALS_831', 'q315r_lyso_001.img'),
        ('APS_14BMC', 'q315_1_001.img'),
        ('APS_17ID', 'q210_1_001.img'),
        ('APS_19ID', 'q315_unbinned_a.0001.img'),
        ('APS_22ID', 'mar300.0001'),
        ('APS_23IDD', 'mar300_1_E1.0001'),
        ('APS_24IDC', 'pilatus_1_0001.cbf'),
        ('APS_24IDC', 'q315_1_001.img'),
        ('CLS1_08ID1', 'mar225_2_E0_0001.img'),
        ('DESY_ID141', 'q210_2_001.img'),
        ('ESRF_BM14', 'mar165_001.mccd'),
        ('ESRF_BM14', 'mar225_1_001.mccd'),
        ('ESRF_ID231', 'q315r_7_001.img'),
        ('RAXIS-HTC', 'test1_lysozyme_0111060001.osc'),
        ('SLS_X06SA', 'mar225_2_001.img'),
        ('SLS_X06SA', 'pilatus6m_1_00001.cbf'),
        ('SRS_101', 'mar225_001.img'),
        ('SRS_142', 'q4_1_001.img'),
        ('SSRL_bl111', 'mar325_1_001.mccd'),
        ('xia2', 'merge2cbf_averaged_0001.cbf'),
#        ('XDS', 'XPARM.XDS'),
#        ('XDS', 'INTEGRATE.HKL'),
#        ('XDS', 'XDS_ASCII.HKL')
        ]
    return [join(path, *f) for f in filenames]

  def multiple_block_filenames(self):
    return self.single_sweep_filenames() + self.all_image_examples()

  def pickle_then_unpickle(self, obj):
    '''Pickle to a temp file then un-pickle.'''
    import cPickle as pickle
    import cStringIO

    # Create a temporary "file"
    temp = cStringIO.StringIO()

    # Pickle the object
    pickle.dump(obj, temp)

    # Read the object
    temp.seek(0)
    return pickle.load(temp)

  def encode_json_then_decode(self, obj, check_format=True):
    import json
    string = json.dumps([db.to_dict() for db in obj], ensure_ascii=True)
    return DataBlockFactory.from_json(string, check_format=check_format)

  def run(self):
    self.tst_create_single_sweep()
    self.tst_create_multiple_sweeps()
    self.tst_create_multiple_blocks()
    self.tst_pickling()
    self.tst_json()
    self.tst_from_null_sweep()
    self.tst_with_external_lookup()
    self.tst_with_bad_external_lookup()

  def tst_create_single_sweep(self):
    filenames = self.single_sweep_filenames()
    blocks = DataBlockFactory.from_filenames(filenames)
    assert(len(blocks) == 1)
    assert(blocks[0].num_images() == 9)
    imageset = blocks[0].extract_imagesets()
    assert(len(imageset) == 1)
    assert(len(imageset[0]) == 9)
    sweeps = blocks[0].extract_sweeps()
    assert(len(sweeps) == 1)
    assert(len(sweeps[0]) == 9)
    print 'OK'

  def tst_create_multiple_sweeps(self):
    filenames = self.multiple_sweep_filenames()
    blocks = DataBlockFactory.from_filenames(filenames)
    assert(len(blocks) == 1)
    assert(blocks[0].num_images() == 6)
    imageset = blocks[0].extract_imagesets()
    assert(len(imageset) == 2)
    sweeps = blocks[0].extract_sweeps()
    assert(len(sweeps) == 2)
    assert(len(sweeps[0]) == 3)
    assert(len(sweeps[1]) == 3)
    print 'OK'

  def tst_create_multiple_blocks(self):
    filenames = self.multiple_block_filenames()
    blocks = DataBlockFactory.from_filenames(filenames, verbose=False)
    assert(len(blocks) == 22)

    # Block 1
    assert(blocks[0].num_images() == 9)
    imageset = blocks[0].extract_imagesets()
    assert(len(imageset) == 1)
    assert(len(imageset[0]) == 9)
    sweeps = blocks[0].extract_sweeps()
    assert(len(sweeps) == 1)
    assert(len(sweeps[0]) == 9)

    print 'OK'

  def tst_pickling(self):
    filenames = self.multiple_block_filenames()
    blocks1 = DataBlockFactory.from_filenames(filenames)
    blocks2 = self.pickle_then_unpickle(blocks1)
    assert(len(blocks2) == len(blocks1))
    for b1, b2 in zip(blocks1, blocks2):
      assert(b1.format_class() == b2.format_class())
      assert(b1 == b2)
    assert(blocks1 == blocks2)

    print 'OK'

  def tst_json(self):

    filenames = self.multiple_block_filenames()
    blocks1 = DataBlockFactory.from_filenames(filenames)
    blocks2 = self.encode_json_then_decode(blocks1)
    assert(len(blocks2) == len(blocks1))
    for b1, b2 in zip(blocks1, blocks2):
      assert(b1.format_class() == b2.format_class())
      assert(b1 == b2)
    assert(blocks1 == blocks2)

    filenames = self.multiple_block_filenames()
    blocks1 = DataBlockFactory.from_filenames(filenames)
    blocks2 = self.encode_json_then_decode(blocks1, check_format=False)
    assert(len(blocks2) == len(blocks1))
    for b1, b2 in zip(blocks1, blocks2):
      for im1, im2 in zip(b1.extract_imagesets(), b2.extract_imagesets()):
        assert(len(im1) == len(im2))
        if isinstance(im1, ImageSweep):
          assert(isinstance(im2, ImageSweep))
          assert(im1.get_beam() == im2.get_beam())
          assert(im1.get_detector() == im2.get_detector())
          assert(im1.get_goniometer() == im2.get_goniometer())
          assert(im1.get_scan() == im2.get_scan())
        else:
          assert(not isinstance(im2, ImageSweep))
          for i in xrange(len(im1)):
            assert(im1.get_beam(i) == im2.get_beam(i))
            assert(im1.get_detector(i) == im2.get_detector(i))

    print 'OK'

  def tst_from_null_sweep(self):
    from dxtbx.format.Format import Format
    from dxtbx.imageset import ImageSweep
    from dxtbx.model import Beam, Detector, Goniometer, Scan

    filenames = ["template_%2d.cbf" % (i+1) for i in range(0, 10)]
    sweep = Format.get_imageset(
      filenames,
      beam = Beam((0, 0, 1)),
      detector = Detector(),
      goniometer = Goniometer((1, 0, 0)),
      scan = Scan((1, 10), (0, 0.1)))

    # Create the datablock
    datablock = DataBlockFactory.from_imageset(sweep)
    assert(len(datablock) == 1)
    datablock = datablock[0]

    sweeps = datablock.extract_sweeps()
    assert(len(sweeps) == 1)
    assert(sweeps[0].get_beam() == sweep.get_beam())
    assert(sweeps[0].get_detector() == sweep.get_detector())
    assert(sweeps[0].get_goniometer() == sweep.get_goniometer())
    assert(sweeps[0].get_scan() == sweep.get_scan())

    print 'OK'

  def tst_with_bad_external_lookup(self):
    filename = join(self.dials_regression, "centroid_test_data",
                    "datablock_with_bad_lookup.json")
    blocks = DataBlockFactory.from_json_file(filename, check_format=False)
    assert(len(blocks) == 1)
    imageset = blocks[0].extract_imagesets()[0]
    assert imageset.external_lookup.mask.filename is not None
    assert imageset.external_lookup.gain.filename is not None
    assert imageset.external_lookup.pedestal.filename is not None
    assert imageset.external_lookup.mask.data is None
    assert imageset.external_lookup.gain.data is None
    assert imageset.external_lookup.pedestal.data is None

    blocks = self.encode_json_then_decode(blocks, check_format=False)
    assert(len(blocks) == 1)
    imageset = blocks[0].extract_imagesets()[0]
    assert imageset.external_lookup.mask.filename is not None
    assert imageset.external_lookup.gain.filename is not None
    assert imageset.external_lookup.pedestal.filename is not None
    assert imageset.external_lookup.mask.data is None
    assert imageset.external_lookup.gain.data is None
    assert imageset.external_lookup.pedestal.data is None

    print 'OK'

  def tst_with_external_lookup(self):
    filename = join(self.dials_regression, "centroid_test_data",
                    "datablock_with_lookup.json")
    blocks = DataBlockFactory.from_json_file(filename)
    assert(len(blocks) == 1)
    imageset = blocks[0].extract_imagesets()[0]
    assert imageset.external_lookup.mask.data is not None
    assert imageset.external_lookup.gain.data is not None
    assert imageset.external_lookup.pedestal.data is not None
    assert imageset.external_lookup.mask.filename is not None
    assert imageset.external_lookup.gain.filename is not None
    assert imageset.external_lookup.pedestal.filename is not None
    assert imageset.external_lookup.mask.data.all_eq(True)
    assert imageset.external_lookup.gain.data.all_eq(1)
    assert imageset.external_lookup.pedestal.data.all_eq(0)

    blocks = self.encode_json_then_decode(blocks)
    assert(len(blocks) == 1)
    imageset = blocks[0].extract_imagesets()[0]
    assert imageset.external_lookup.mask.data is not None
    assert imageset.external_lookup.gain.data is not None
    assert imageset.external_lookup.pedestal.data is not None
    assert imageset.external_lookup.mask.filename is not None
    assert imageset.external_lookup.gain.filename is not None
    assert imageset.external_lookup.pedestal.filename is not None
    assert imageset.external_lookup.mask.data.all_eq(True)
    assert imageset.external_lookup.gain.data.all_eq(1)
    assert imageset.external_lookup.pedestal.data.all_eq(0)

    print 'OK'

if __name__ == '__main__':
  if not libtbx.env.has_module("dials"):
    print "Skipping test: dials not present"
  elif not libtbx.env.has_module("dials_regression"):
    print "Skipping test: dials_regression not present"
  else:
    test = Test()
    test.run()
