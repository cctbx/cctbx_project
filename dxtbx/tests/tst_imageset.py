from __future__ import division
import libtbx.load_env

dials_regression = libtbx.env.dist_path('dials_regression')

class TestFormat(object):

  def get_images(self):
    images = [
      "./image_examples/ALS_501/als501_q4_1_001.img",
      "./image_examples/SPring8_BL26B1_SaturnA200/A200_000001.img",
      "./image_examples/SPring8_BL26B1_SaturnA200/A200_000002.img",
      "./image_examples/APS_19ID/q315_unbinned_a.0001.img",
      "./image_examples/ALS_821/q210_lyso_1_101.img",
      "./image_examples/MLFSOM_simulation/fake_00001.img",
      "./image_examples/ALS_831/q315r_lyso_001.img",
      "./image_examples/DESY_ID141/q210_2_001.img",
      "./image_examples/SSRL_bl91/q315_1_001.img",
      "./image_examples/APS_24IDC/q315_1_001.img",
      "./image_examples/ALS_1231/q315r_lyso_1_001.img",
      "./image_examples/SRS_142/q4_1_001.img",
      "./image_examples/ALS_422/lyso_041013a_1_001.img",
      "./image_examples/APS_17ID/q210_1_001.img",
      "./image_examples/saturn/lyso_00001.img",

      # "./image_examples/SPring8_BL26B1_Raxis5/raxis5_000091.img",
      # "./image_examples/SPring8_BL26B1_Raxis5/raxis5_000001.img",

      "./image_examples/SPring8_BL32XU_MX225HS/ds_000045.img",
      "./image_examples/SPring8_BL32XU_MX225HS/ds_000001.img",
      "./image_examples/SPring8_BL44XU_MX300HE/bl44xu_lys_000002.img",
      "./image_examples/SPring8_BL44XU_MX300HE/bl44xu_lys_000001.img",
      "./image_examples/SPring8_BL32XU/rayonix225hs_0001.img",
      "./image_examples/SPring8_BL32XU/rayonix225_0001.img",
      "./image_examples/SLS_X06SA/mar225_2_001.img",
      "./image_examples/CLS1_08ID1/mar225_2_E0_0001.img",
      "./image_examples/SPring8_BL38B1_MX225HE/bl38b1_001.img",
      "./image_examples/SPring8_BL38B1_MX225HE/bl38b1_090.img",
      "./image_examples/SRS_101/mar225_001.img",
      "./image_examples/SPring8_BL12B2_MX225HE/lys001_000091.img",
      "./image_examples/SPring8_BL12B2_MX225HE/lys001_000001.img",
      "./image_examples/SPring8_BL26B2_MX225/2sec_Al200um_000001.img",
      "./image_examples/SPring8_BL26B2_MX225/2sec_Al200um_000090.img",

      "./stills_test_data/hit-20111202210224984.cbf", # Multitile
      "./stills_test_data/hit-s00-20140306002935980.cbf", # Multitile
      "./stills_test_data/hit-s00-20140306002857363.cbf", # Multitile

       "./image_examples/ESRF_ID29/trypsin_1_0001.cbf",
       "./image_examples/xia2/merge2cbf_averaged_0001.cbf",
       "./image_examples/dials-190/whatev1_01_00002.cbf",
       "./image_examples/dials-190/whatev1_02_00001.cbf",
       "./image_examples/dials-190/whatev1_02_00002.cbf",
       "./image_examples/dials-190/whatev1_03_00001.cbf",
       "./image_examples/dials-190/whatev1_01_00001.cbf",
       "./image_examples/dials-190/whatev1_03_00002.cbf",
       "./image_examples/APS_24IDE_test/thaum-12_1_0005.cbf",
       "./image_examples/APS_24IDE_test/thaum-12_1_0004.cbf",
       "./image_examples/APS_24IDE_test/thaum-12_1_0006.cbf",
       "./image_examples/APS_24IDE_test/thaum-12_1_0008.cbf",
       "./image_examples/APS_24IDE_test/thaum-12_1_0009.cbf",
       "./image_examples/APS_24IDE_test/thaum-12_1_0010.cbf",
       "./image_examples/APS_24IDE_test/thaum-12_1_0007.cbf",
       "./image_examples/APS_24IDE_test/thaum-12_1_0002.cbf",
       "./image_examples/APS_24IDE_test/thaum-12_1_0001.cbf",
       "./image_examples/APS_24IDE_test/thaum-12_1_0003.cbf",
       "./image_examples/APS_24IDC/pilatus_1_0001.cbf",
       "./image_examples/SLS_Eiger_16M_as_CBF/insu_with_bs_labelit_0901.cbf",
       "./image_examples/SLS_Eiger_16M_as_CBF/insu_with_bs_labelit_0001.cbf",
       "./image_examples/SPring8_BL41XU_PILATUS3_6M/data1_000001.cbf",
       "./image_examples/SPring8_BL41XU_PILATUS3_6M/data1_000901.cbf",
       "./image_examples/DLS_I02/X4_wide_M1S4_1_0001.cbf",
       "./image_examples/DLS_I23/I23_P12M_alpha_0001.cbf",
       "./image_examples/DLS_I23/germ_13KeV_0001.cbf",
       "./image_examples/SLS_X06SA/pilatus6m_1_00001.cbf",
       "./image_examples/SPring8_ADSC_SN916/Xtal17-2phi_3_015.cbf",
       "./image_examples/DLS_I19/I19_P300k_00001.cbf",
       "./image_examples/ED_From_TIFF/170112330001.cbf",

      "./image_examples/LCLS_cspad_nexus/cxi78513_bslz4_r0014_subset4_master.h5"
    ]

    from os.path import join

    images = map(lambda f: join(dials_regression, f), images)

    return images

  def run(self):
    from dxtbx.format.Registry import Registry

    db_fail_count = 0
    for image in self.get_images():
      print image
      format_class = Registry.find(image)
      reader = format_class.get_reader()([image])
      masker = format_class.get_masker()([image])

      N = len(reader)

      for i in range(N):
        data = reader.read(i)
        mask = masker.get(i)

      iset = format_class.get_imageset([image])

    print 'OK'




class TestImageSet(object):
  def get_file_list(self):
    import os.path

    dials_regression = libtbx.env.dist_path('dials_regression')
    path = os.path.join(dials_regression, 'centroid_test_data')

    # Non-sequential Filenames and image indices
    filenames = []
    image_indices = range(1, 10)
    for i in image_indices:
      filenames.append(os.path.join(path, 'centroid_000{0}.cbf'.format(i)))

    return filenames

  def run(self):
    from dxtbx.imageset import ImageSet
    from dxtbx.format.Registry import Registry

    # Get the filenames
    filenames = self.get_file_list()

    # Create the format class
    format_class = Registry.find(filenames[0])

    # Create the reader
    imageset = format_class.get_imageset(filenames, as_imageset=True)

    # Run a load of tests
    self.tst_get_item(imageset)
    self.tst_len(imageset, len(filenames))
    self.tst_iter(imageset)
    self.tst_paths(imageset, filenames)
    self.tst_get_detectorbase(imageset, range(len(filenames)), 9)
    self.tst_get_models(imageset, range(len(filenames)), 9)

  def tst_get_item(self, imageset):
    image = imageset[0]
    try:
      image = imageset[9]
      assert(False)
    except IndexError:
      pass

    imageset2 = imageset[3:7]
    image = imageset2[0]
    try:
      image = imageset2[5]
      assert(False)
    except IndexError:
      pass

    self.tst_len(imageset2, 4)
    self.tst_get_detectorbase(imageset2, range(0, 4), 5)
    self.tst_get_models(imageset2, range(0, 4), 5)
    self.tst_paths(imageset2, imageset.paths()[3:7])
    self.tst_iter(imageset2)

    imageset2 = imageset[3:7:2]
    image = imageset2[0]
    try:
      image = imageset2[2]
      assert(False)
    except IndexError:
      pass

    self.tst_len(imageset2, 2)
    self.tst_get_detectorbase(imageset2, range(0, 2), 2)
    self.tst_get_models(imageset2, range(0, 2), 2)
    self.tst_paths(imageset2, imageset.paths()[3:7:2])
    self.tst_iter(imageset2)

    print 'OK'

  def tst_len(self, imageset, length):
    assert(len(imageset) == length)
    print 'OK'

  def tst_iter(self, imageset):
    for image in imageset:
      pass
    print 'OK'

  def tst_paths(self, imageset, filenames1):
    filenames2 = imageset.paths()
    for f1, f2 in zip(filenames1, filenames2):
      assert(f1 == f2)
    print 'OK'

  def tst_get_detectorbase(self, imageset, indices, outside_index):
    for i in indices:
      imageset.get_detectorbase(i)

    try:
      imageset.get_detectorbase(outside_index)
      assert(False)
    except IndexError:
      pass
    print 'OK'

  def tst_get_models(self, imageset, indices, outside_index):
    for i in indices:
      self.tst_get_models_index(imageset, i)

    try:
      self.tst_get_models_index(imageset, outside_index)
      assert(False)
    except IndexError:
      pass
    print 'OK'

  def tst_get_models_index(self, imageset, index=None):
    imageset.get_detector(index)
    imageset.get_beam(index)

  def tst_set_models(self, imageset):
    from dxtbx.model import Beam, Detector, Panel

    # Create some other models
    beam = Beam((1, 0, 0), 0.5)
    detector = Detector(Panel("UNKNOWN", "Panel",
                              (1, 0, 0), (0, 1, 0), (0, 0, 1),
                              (0.1, 0.1), (1000, 1000), (0, 1)))

    # Override sweep models
    imageset.set_beam(beam)
    imageset.set_detector(detector)

    # Ensure this doens't interfere with reading
    for i in imageset:
      pass

    # Get the models back and check they're ok
    beam2 = imageset.get_beam()
    detector2 = imageset.get_detector()
    assert(beam2 == beam)
    assert(detector2 == detector)

    # Get the models from an index back and check they're not the same
    beam2 = imageset.get_beam(0)
    detector2 = imageset.get_detector(0)
    assert(beam2 != beam)
    assert(detector2 != detector)


class TestImageSweep(object):
  def get_file_list(self):
    import os.path

    dials_regression = libtbx.env.dist_path('dials_regression')
    path = os.path.join(dials_regression, 'centroid_test_data')

    # Non-sequential Filenames and image indices
    template = os.path.join(path, 'centroid_%04d.cbf')
    array_range = (0, 9)

    filenames = [template % (i+1) for i in range(*array_range)]

    return filenames

  def run(self):
    from dxtbx.imageset import ImageSweep
    from dxtbx.format.Registry import Registry

    # Get the filenames
    filenames = self.get_file_list()

    # Create the format class
    format_class = Registry.find(filenames[0])

    # Create the sweep
    sweep = format_class.get_imageset(filenames)

    # Run a load of tests
    self.tst_get_item(sweep)
    self.tst_len(sweep, len(filenames))
    self.tst_iter(sweep)
    self.tst_paths(sweep, filenames)
    self.tst_get_detectorbase(sweep, range(len(filenames)), 9)
    self.tst_get_models(sweep, range(len(filenames)), 9)
    self.tst_get_array_range(sweep, (0, 9))
    self.tst_set_models(sweep)

  def tst_get_item(self, sweep):
    image = sweep[0]
    try:
      image = sweep[9]
      assert(False)
    except IndexError:
      pass

    sweep2 = sweep[3:7]
    image = sweep2[0]
    try:
      image = sweep2[5]
      assert(False)
    except IndexError:
      pass

    self.tst_len(sweep2, 4)
    self.tst_get_detectorbase(sweep2, range(0, 4), 5)
    self.tst_get_models(sweep2, range(0, 4), 5)
    self.tst_paths(sweep2, sweep.paths()[3:7])
    self.tst_iter(sweep2)
    self.tst_get_array_range(sweep2, (3, 7))

    try:
      sweep2 = sweep[3:7:2]
      assert(False)
    except IndexError:
      pass

    print 'OK'

  def tst_len(self, sweep, length):
    assert(len(sweep) == length)
    print 'OK'

  def tst_iter(self, sweep):
    for image in sweep:
      pass
    print 'OK'

  def tst_paths(self, sweep, filenames1):
    filenames2 = sweep.paths()
    for f1, f2 in zip(filenames1, filenames2):
      assert(f1 == f2)
    print 'OK'


  def tst_get_detectorbase(self, sweep, indices, outside_index):
    for i in indices:
      sweep.get_detectorbase(i)

    try:
      sweep.get_detectorbase(outside_index)
      assert(False)
    except IndexError:
      pass
    print 'OK'

  def tst_get_models(self, sweep, indices, outside_index):
    self.tst_get_models_index(sweep)
    for i in indices:
      self.tst_get_models_index(sweep, i)

    print 'OK'

  def tst_get_models_index(self, sweep, index=None):
    sweep.get_detector(index)
    sweep.get_beam(index)
    sweep.get_goniometer(index)
    sweep.get_scan(index)

    # Ensure state at zero
    sweep[0]
    scan1 = sweep.get_scan()
    # Put sweep to end
    sweep[len(sweep)-1]
    scan2 = sweep.get_scan()
    assert(scan1 == scan2)
    print 'OK'

  def tst_get_array_range(self, sweep, array_range):
    assert(sweep.get_array_range() == array_range)
    print 'OK'

  def tst_set_models(self, sweep):
    from dxtbx.model import Beam, Detector, Panel

    # Get some models
    beam = sweep.get_beam()
    gonio = sweep.get_goniometer()
    detector = sweep.get_detector()

    # Modify the geometry
    assert(len(detector) == 1)
    beam.set_direction((1, 0, 0))
    gonio.set_rotation_axis((0, 1, 0))
    detector[0].set_local_frame((1, 0, 0), (0, 1, 0), (0, 0, 1))

    # Override sweep models
    sweep.set_beam(beam)
    sweep.set_goniometer(gonio)
    sweep.set_detector(detector)

    # Ensure this doens't interfere with reading
    for i in sweep:
      pass

    # Get the models back and check they're ok
    beam2 = sweep.get_beam()
    gonio2 = sweep.get_goniometer()
    detector2 = sweep.get_detector()
    assert(beam2 == beam)
    assert(gonio2 == gonio)
    assert(detector2 == detector)

    # Get the models from an index back and check they're the same
    beam2 = sweep.get_beam(0)
    gonio2 = sweep.get_goniometer(0)
    detector2 = sweep.get_detector(0)
    assert(beam2 == beam)
    assert(gonio2 == gonio)
    assert(detector2 == detector)

    # Get a sub sweep
    sub_sweep = sweep[3:7]

    # Get the models back and check they're ok
    beam2 = sub_sweep.get_beam()
    gonio2 = sub_sweep.get_goniometer()
    detector2 = sub_sweep.get_detector()
    assert(beam2 == beam)
    assert(gonio2 == gonio)
    assert(detector2 == detector)

    # Get the models from an index back and check they're not the same
    beam2 = sub_sweep.get_beam(0)
    gonio2 = sub_sweep.get_goniometer(0)
    detector2 = sub_sweep.get_detector(0)
    assert(beam2 == beam)
    assert(gonio2 == gonio)
    assert(detector2 == detector)


class TestNexusFile(object):

  def __init__(self):
    from os.path import join
    self.filename = join(dials_regression,
                    "./image_examples/LCLS_cspad_nexus/cxi78513_bslz4_r0014_subset4_master.h5")

  def run(self):

    from dxtbx.format.Registry import Registry
    format_class = Registry.find(self.filename)

    iset = format_class.get_imageset([self.filename])

    assert len(iset) == 2
    for i in range(len(iset)):
      data = iset.get_raw_data(i)
      mask = iset.get_mask(i)
      b = iset.get_beam(i)
      d = iset.get_detector(i)
      g = iset.get_goniometer(i)
      s = iset.get_scan(i)

    print 'OK'

    iset = format_class.get_imageset([self.filename], single_file_indices=[1])
    assert len(iset) == 1

    for i in range(len(iset)):
      data = iset.get_raw_data(i)
      mask = iset.get_mask(i)
      b = iset.get_beam(i)
      d = iset.get_detector(i)
      g = iset.get_goniometer(i)
      s = iset.get_scan(i)

    print 'OK'



class TestImageSetFactory(object):
  def get_file_list(self):
    import os.path

    dials_regression = libtbx.env.dist_path('dials_regression')
    path = os.path.join(dials_regression, 'centroid_test_data')

    # Non-sequential Filenames and image indices
    filenames = []
    image_indices = range(1, 10)
    for i in image_indices:
      filenames.append(os.path.join(path, 'centroid_000{0}.cbf'.format(i)))

    return filenames

  def run(self):
    from dxtbx.imageset import ImageSetFactory, ImageSweep
    from os.path import join

    filenames = self.get_file_list()

    sweep = ImageSetFactory.new(filenames)

    assert(isinstance(sweep[0], ImageSweep) == True)

    print 'OK'

    template = join(dials_regression, "centroid_test_data", "centroid_####.cbf")
    image_range = (3, 6)

    sweep = ImageSetFactory.from_template(template, image_range)

    assert(isinstance(sweep[0], ImageSweep) == True)
    assert len(sweep[0]) == 4
    assert sweep[0].paths()[0].endswith("3.cbf")
    assert sweep[0].paths()[-1].endswith("6.cbf")

    print 'OK'

    imageset = ImageSetFactory.make_imageset(filenames)
    assert len(imageset) == 9

    print 'OK'

    imageset = ImageSetFactory.make_imageset(
      filenames,
      check_format=False)
    assert len(imageset) == 9

    print 'OK'

    sweep = ImageSetFactory.make_sweep(
      template,
      list(range(1, 9+1)))
    assert len(sweep) == 9

    print 'OK'

    sweep = ImageSetFactory.make_sweep(
      template,
      list(range(3, 6+1)))
    assert len(sweep) == 4

    print 'OK'


class TestPickleImageSet(object):

  def __init__(self):
    import os.path
    dials_regression = libtbx.env.dist_path('dials_regression')
    path = os.path.join(dials_regression, 'centroid_test_data')

    # Non-sequential Filenames and image indices
    filenames = []
    image_indices = range(1, 10)
    for i in image_indices:
      filenames.append(os.path.join(path, 'centroid_000{0}.cbf'.format(i)))

    from dxtbx.imageset import ImageSetFactory, ImageSweep
    self.sweep = ImageSetFactory.new(filenames)[0]

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

  def run(self):

    # Read the 5th image
    image = self.sweep[4]

    sweep2 = self.pickle_then_unpickle(self.sweep)

    assert(self.sweep.get_template() == sweep2.get_template())
    assert(self.sweep.get_array_range() == sweep2.get_array_range())
    assert(self.sweep.get_beam() == sweep2.get_beam())
    assert(self.sweep.get_goniometer() == sweep2.get_goniometer())
    assert(self.sweep.get_scan() == sweep2.get_scan())
    assert(self.sweep.paths() == sweep2.paths())
    assert(self.sweep == sweep2)

    print 'OK'

def run():

  TestFormat().run()
  TestImageSet().run()
  TestImageSweep().run()
  TestNexusFile().run()
  TestImageSetFactory().run()
  TestPickleImageSet().run()


if __name__ == '__main__':
  if not libtbx.env.has_module("dials"):
    print "Skipping test: dials not present"
  elif not libtbx.env.has_module("dials_regression"):
    print "Skipping test: dials_regression not present"
  else:
    run()
