
class Test(object):

  def __init__(self):
    import libtbx.load_env
    import os

    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      return

    self.centroid_test_data = os.path.join(dials_regression, 'centroid_test_data')
    self.image_examples = os.path.join(dials_regression, 'image_examples')

  def single_sweep_filenames(self):
    from os.path import join
    path = self.centroid_test_data
    filenames = []
    image_indices = range(1, 10)
    for i in image_indices:
      filenames.append(join(path, 'centroid_000{0}.cbf'.format(i)))
    return filenames

  def multiple_sweep_filenames(self):
    from os.path import join
    path = self.centroid_test_data
    filenames = []
    image_indices = list(range(1, 4)) + list(range(7, 10))
    for i in image_indices:
      filenames.append(join(path, 'centroid_000{0}.cbf'.format(i)))
    return filenames

  def all_image_examples(self):
    from os.path import join
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
    import pickle
    import tempfile

    # Create a tmp file
    temp = tempfile.TemporaryFile()

    # Pickle the object
    pickle.dump(obj, temp)

    # Read the object
    temp.flush()
    temp.seek(0)
    return pickle.load(temp)

  def run(self):
    self.tst_create_single_sweep()
    self.tst_create_multiple_sweeps()
    self.tst_create_multiple_blocks()
    self.tst_pickling()

  def tst_create_single_sweep(self):

    from dxtbx.datablock import DataBlockFactory

    filenames = self.single_sweep_filenames()
    blocks = DataBlockFactory.from_filenames(filenames)
    assert(len(blocks) == 1)
    assert(len(blocks[0]) == 9)
    assert(len(blocks[0].filenames()) == 9)
    assert(len(blocks[0].metadata()) == 9)
    imageset = blocks[0].extract_all()
    assert(len(imageset) == 9)
    imageset = blocks[0].extract_stills()
    assert(imageset == None)
    sweeps = blocks[0].extract_sweeps()
    assert(len(sweeps) == 1)
    assert(len(sweeps[0]) == 9)
    print 'OK'

  def tst_create_multiple_sweeps(self):

    from dxtbx.datablock import DataBlockFactory

    filenames = self.multiple_sweep_filenames()
    blocks = DataBlockFactory.from_filenames(filenames)
    assert(len(blocks) == 1)
    assert(len(blocks[0]) == 6)
    assert(len(blocks[0].filenames()) == 6)
    assert(len(blocks[0].metadata()) == 6)
    imageset = blocks[0].extract_all()
    assert(len(imageset) == 6)
    imageset = blocks[0].extract_stills()
    assert(imageset == None)
    sweeps = blocks[0].extract_sweeps()
    assert(len(sweeps) == 2)
    assert(len(sweeps[0]) == 3)
    assert(len(sweeps[0]) == 3)
    print 'OK'

  def tst_create_multiple_blocks(self):

    from dxtbx.datablock import DataBlockFactory

    filenames = self.multiple_block_filenames()
    blocks = DataBlockFactory.from_filenames(filenames)

    assert(len(blocks) == 16)

    # Block 1
    assert(len(blocks[0]) == 9)
    assert(len(blocks[0].filenames()) == 9)
    assert(len(blocks[0].metadata()) == 9)
    imageset = blocks[0].extract_all()
    assert(len(imageset) == 9)
    imageset = blocks[0].extract_stills()
    assert(imageset == None)
    sweeps = blocks[0].extract_sweeps()
    assert(len(sweeps) == 1)
    assert(len(sweeps[0]) == 9)

    # Block 2
    assert(len(blocks[1]) == 7)
    assert(len(blocks[1].filenames()) == 7)
    assert(len(blocks[1].metadata()) == 7)
    imageset = blocks[1].extract_all()
    assert(len(imageset) == 7)
    imageset = blocks[1].extract_stills()
    assert(imageset == None)
    sweeps = blocks[1].extract_sweeps()
    assert(len(sweeps) == 7)
    assert(all(len(s) == 1 for s in sweeps))

    # Block 3
    assert(len(blocks[2]) == 2)
    assert(len(blocks[2].filenames()) == 2)
    assert(len(blocks[2].metadata()) == 2)
    imageset = blocks[2].extract_all()
    assert(len(imageset) == 2)
    imageset = blocks[2].extract_stills()
    assert(imageset == None)
    sweeps = blocks[2].extract_sweeps()
    assert(len(sweeps) == 2)
    assert(all(len(s) == 1 for s in sweeps))

    # Block 4
    assert(len(blocks[3]) == 1)
    assert(len(blocks[3].filenames()) == 1)
    assert(len(blocks[3].metadata()) == 1)
    imageset = blocks[3].extract_all()
    assert(len(imageset) == 1)
    imageset = blocks[3].extract_stills()
    assert(imageset == None)
    sweeps = blocks[3].extract_sweeps()
    assert(len(sweeps) == 1)
    assert(all(len(s) == 1 for s in sweeps))

    print 'OK'

  def tst_pickling(self):

    from dxtbx.datablock import DataBlockFactory

    filenames = self.multiple_block_filenames()
    blocks1 = DataBlockFactory.from_filenames(filenames)
    blocks2 = self.pickle_then_unpickle(blocks1)
    assert(len(blocks2) == len(blocks1))
    for b1, b2 in zip(blocks1, blocks2):
      assert(all(f1 == f2 for f1, f2 in zip(b1.filenames(), b2.filenames())))
      assert(all(m1 == m2 for m1, m2 in zip(b1.metadata(), b2.metadata())))
      assert(b1.format_class() == b2.format_class())
      assert(b1 == b2)
    assert(blocks1 == blocks2)

    print 'OK'


if __name__ == '__main__':
  test = Test()
  test.run()
