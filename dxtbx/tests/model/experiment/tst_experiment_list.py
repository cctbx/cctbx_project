

from __future__ import absolute_import, division

from dxtbx.model import Experiment, ExperimentList
from dxtbx.model.experiment_list import ExperimentListFactory, \
  ExperimentListDumper

class TestExperiment(object):

  def __init__(self, path):
    self.path = path

  def run(self):
    self.tst_contains()
    self.tst_equality()
    self.tst_consistent()

  def tst_contains(self):
    from dxtbx.model import Beam, Detector, Goniometer, Scan
    from dxtbx.model import Crystal

    # Create a load of models
    b1 = Beam()
    d1 = Detector()
    g1 = Goniometer()
    s1 = Scan()
    c1 = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), space_group_symbol="P1")

    # Create an experiment
    e = Experiment(
      beam=b1, detector=d1, goniometer=g1,
      scan=s1, crystal=c1, imageset=None)

    # Check experiment contains model
    assert(b1 in e)
    assert(d1 in e)
    assert(g1 in e)
    assert(s1 in e)
    assert(c1 in e)

    # Create a load of models that look the same but aren't
    b2 = Beam()
    d2 = Detector()
    g2 = Goniometer()
    s2 = Scan()
    c2 = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), space_group_symbol="P1")

    # Check experiment doesn't contain model
    assert(b2 not in e)
    assert(d2 not in e)
    assert(g2 not in e)
    assert(s2 not in e)
    assert(c2 not in e)

    # Test passed
    print 'OK'

  def tst_equality(self):

    from dxtbx.model import Beam, Detector, Goniometer, Scan
    from dxtbx.model import Crystal

    # Create a load of models
    b1 = Beam()
    d1 = Detector()
    g1 = Goniometer()
    s1 = Scan()
    c1 = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), space_group_symbol="P1")

    # Create a load of models that look the same but aren't
    b2 = Beam()
    d2 = Detector()
    g2 = Goniometer()
    s2 = Scan()
    c2 = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), space_group_symbol="P1")

    # Create an experiment
    e1 = Experiment(
      beam=b1, detector=d1, goniometer=g1,
      scan=s1, crystal=c1, imageset=None)

    # Create an experiment
    e2 = Experiment(
      beam=b1, detector=d1, goniometer=g1,
      scan=s1, crystal=c1, imageset=None)

    # Create an experiment
    e3 = Experiment(
      beam=b2, detector=d2, goniometer=g2,
      scan=s2, crystal=c2, imageset=None)

    # Check e1 equals e2 but not e3
    assert(e1 == e2)
    assert(e1 != e3)
    assert(e2 != e3)

    # Test passed
    print 'OK'

  def tst_consistent(self):

    from dxtbx.imageset import ImageSetFactory
    from glob import glob
    from os.path import join
    from dxtbx.model import Scan

    # Create a sweep
    sweep_filenames = join(self.path, 'centroid_test_data', 'centroid*.cbf')
    sweep = ImageSetFactory.new(sorted(glob(sweep_filenames)))[0]

    # Create experiment with sweep and good scan
    e = Experiment(imageset=sweep, scan=sweep.get_scan())
    assert(e.is_consistent())

    # Create experiment with sweep and defective scan
    scan = sweep.get_scan()
    scan.set_image_range((1, 1))
    e = Experiment(imageset=sweep, scan=scan)
    #assert(not e.is_consistent()) # FIXME

    ## Create experiment with imageset and good scan
    #assert(e.is_consistent())

    ## Create experiment with imageset and non-still scan
    #assert(not e.is_consistent())

    ## Create experiment with imageset and scan with more than 1 image
    #assert(not e.is_consistent())

    ## Create experiment with imageset and defective scan
    #assert(not e.is_consistent())

    # Test passed
    print 'OK'

class TestExperimentList(object):

  def __init__(self, path):
    self.path = path
    self.el = self.generate()

  def run(self):
    self.tst_contains()
    # self.tst_index()
    self.tst_replace()
    self.tst_indices()
    self.tst_models()
    self.tst_to_dict()
    self.tst_where()

  def tst_contains(self):

    from dxtbx.model import Beam, Detector, Goniometer, Scan

    # Check all the models are found
    for e in self.el:
      assert(e.beam in self.el)
      assert(e.detector in self.el)
      assert(e.goniometer in self.el)
      assert(e.scan in self.el)

    # Create some more models
    b = Beam()
    d = Detector()
    g = Goniometer()
    s = Scan()

    # Check that models not in are not found
    assert(b not in self.el)
    assert(d not in self.el)
    assert(g not in self.el)
    assert(s not in self.el)

    # Test passed
    print 'OK'

  # def tst_index(self):

  #   # Check the indices of exisiting experiments
  #   assert(self.el.index(self.el[0]) is 0)
  #   assert(self.el.index(self.el[1]) is 1)
  #   assert(self.el.index(self.el[2]) is 2)
  #   assert(self.el.index(self.el[3]) is 1)
  #   assert(self.el.index(self.el[4]) is 0)

  #   # Check index of non exisiting experiment
  #   try:
  #     self.el.index(Experiment())
  #     assert(False)
  #   except ValueError:
  #     pass

  #   # Test passed
  #   print 'OK'

  def tst_replace(self):

    # Get the models
    b = [e.beam for e in self.el]
    d = [e.detector for e in self.el]
    g = [e.goniometer for e in self.el]
    s = [e.scan for e in self.el]

    # Replace some models
    self.el.replace(b[0], b[1])
    assert(self.el[0].beam is b[1])
    assert(self.el[4].beam is b[1])

    # Replace again
    self.el[0].beam = b[0]
    self.el[4].beam = b[4]

    # Test passed
    print 'OK'

  def tst_indices(self):
    from dxtbx.model import Beam, Detector, Goniometer, Scan

    # Get the models
    b = [e.beam for e in self.el]
    d = [e.detector for e in self.el]
    g = [e.goniometer for e in self.el]
    s = [e.scan for e in self.el]

    # Check indices of beams
    assert(list(self.el.indices(b[0])) == [0, 4])
    assert(list(self.el.indices(b[1])) == [1, 3])
    assert(list(self.el.indices(b[2])) == [2])
    assert(list(self.el.indices(b[3])) == [1, 3])
    assert(list(self.el.indices(b[4])) == [0, 4])

    # Check indices of detectors
    assert(list(self.el.indices(d[0])) == [0, 4])
    assert(list(self.el.indices(d[1])) == [1, 3])
    assert(list(self.el.indices(d[2])) == [2])
    assert(list(self.el.indices(d[3])) == [1, 3])
    assert(list(self.el.indices(d[4])) == [0, 4])

    # Check indices of goniometer
    assert(list(self.el.indices(g[0])) == [0, 4])
    assert(list(self.el.indices(g[1])) == [1, 3])
    assert(list(self.el.indices(g[2])) == [2])
    assert(list(self.el.indices(g[3])) == [1, 3])
    assert(list(self.el.indices(g[4])) == [0, 4])

    # Check indices of scans
    assert(list(self.el.indices(s[0])) == [0, 4])
    assert(list(self.el.indices(s[1])) == [1, 3])
    assert(list(self.el.indices(s[2])) == [2])
    assert(list(self.el.indices(s[3])) == [1, 3])
    assert(list(self.el.indices(s[4])) == [0, 4])

    # Check some models not in the list
    assert(len(self.el.indices(Beam())) == 0)
    assert(len(self.el.indices(Detector())) == 0)
    assert(len(self.el.indices(Goniometer())) == 0)
    assert(len(self.el.indices(Scan())) == 0)

    # Test passed
    print 'OK'

  def tst_models(self):

    # Get all the unique models
    b = self.el.beams()
    d = self.el.detectors()
    g = self.el.goniometers()
    s = self.el.scans()

    # Check we have the expected number
    assert(len(b) == 3)
    assert(len(d) == 3)
    assert(len(g) == 3)
    assert(len(s) == 3)

    # Check we have the expected order
    assert(b[0] == self.el[0].beam)
    assert(b[1] == self.el[1].beam)
    assert(b[2] == self.el[2].beam)

    assert(d[0] == self.el[0].detector)
    assert(d[1] == self.el[1].detector)
    assert(d[2] == self.el[2].detector)

    assert(g[0] == self.el[0].goniometer)
    assert(g[0] == self.el[0].goniometer)
    assert(g[1] == self.el[1].goniometer)

    assert(s[2] == self.el[2].scan)
    assert(s[1] == self.el[1].scan)
    assert(s[2] == self.el[2].scan)

    # Test passed
    print 'OK'

  def tst_to_dict(self):

    # Convert the list to a dictionary
    obj = self.el.to_dict()

    # Check this is the right object
    assert(obj['__id__'] == 'ExperimentList')

    # Check length of items
    assert(len(obj['experiment']) == 5)
    assert(len(obj['beam']) == 3)
    assert(len(obj['detector']) == 3)
    assert(len(obj['goniometer']) == 3)
    assert(len(obj['scan']) == 3)

    # The expected models
    b = [0, 1, 2, 1, 0]
    d = [0, 1, 2, 1, 0]
    g = [0, 1, 2, 1, 0]
    s = [0, 1, 2, 1, 0]

    # Check all the experiments
    for i, eobj in enumerate(obj['experiment']):
      assert(eobj['__id__'] == 'Experiment')
      assert(eobj['beam'] == b[i])
      assert(eobj['detector'] == d[i])
      assert(eobj['goniometer'] == g[i])
      assert(eobj['scan'] == s[i])

    # Test passed
    print 'OK'

  def tst_where(self):
    for beam in self.el.beams():
      assert beam is not None
      for i in self.el.where(beam=beam):
        assert self.el[i].beam is beam
    for goniometer in self.el.goniometers():
      assert goniometer is not None
      for i in self.el.where(goniometer=goniometer):
        assert self.el[i].goniometer is goniometer
    for scan in self.el.scans():
      assert scan is not None
      for i in self.el.where(scan=scan):
        assert self.el[i].scan is scan
    for detector in self.el.detectors():
      assert detector is not None
      for i in self.el.where(detector=detector):
        assert self.el[i].detector is detector

    print 'OK'


  def generate(self):
    from dxtbx.model import Beam, Detector, Goniometer, Scan

    # Initialise a list of experiments
    experiments = ExperimentList()

    # Create a few beams
    b1 = Beam()
    b2 = Beam()
    b3 = Beam()

    # Create a few detectors
    d1 = Detector()
    d2 = Detector()
    d3 = Detector()

    # Create a few goniometers
    g1 = Goniometer()
    g2 = Goniometer()
    g3 = Goniometer()

    # Create a few scans
    s1 = Scan()
    s2 = Scan()
    s3 = Scan()

    # Create a list of models
    b = [b1, b2, b3, b2, b1]
    d = [d1, d2, d3, d2, d1]
    g = [g1, g2, g3, g2, g1]
    s = [s1, s2, s3, s2, s1]

    # Populate with various experiments
    for i in range(5):
      experiments.append(Experiment(
        beam=b[i],
        detector=d[i],
        goniometer=g[i],
        scan=s[i]))

    # Return the list of experiments
    return experiments


class TestExperimentListFactory(object):

  def __init__(self, path):
    self.path = path

  def run(self):
    self.tst_from_json()
    self.tst_from_pickle()
    self.tst_from_args()
    self.tst_from_imageset()
    self.tst_from_sweep()
    self.tst_from_datablock()

  def tst_from_json(self):
    from os.path import join
    import os

    os.environ['DIALS_REGRESSION'] = self.path

    # Get all the filenames
    filename1 = join(self.path, 'experiment_test_data', 'experiment_1.json')
    filename2 = join(self.path, 'experiment_test_data', 'experiment_2.json')
    filename3 = join(self.path, 'experiment_test_data', 'experiment_3.json')
    filename4 = join(self.path, 'experiment_test_data', 'experiment_4.json')

    # Read all the experiment lists in
    el1 = ExperimentListFactory.from_json_file(filename1)
    #el2 = ExperimentListFactory.from_json_file(filename2)
    el3 = ExperimentListFactory.from_json_file(filename3)
    el4 = ExperimentListFactory.from_json_file(filename4)

    # All the experiment lists should be the same length
    assert(len(el1) == 1)
    #assert(len(el1) == len(el2))
    assert(len(el1) == len(el3))
    assert(len(el1) == len(el4))

    # Check all the models are the same
    for e in zip(el1, el3, el4):
      e1 = e[0]
      assert(e1.imageset is not None)
      assert(e1.beam is not None)
      assert(e1.detector is not None)
      assert(e1.goniometer is not None)
      assert(e1.scan is not None)
      assert(e1.crystal is not None)
      for ee in e[1:]:
        assert(e1.imageset == ee.imageset)
        assert(e1.beam == ee.beam)
        assert(e1.detector == ee.detector)
        assert(e1.goniometer == ee.goniometer)
        assert(e1.scan == ee.scan)
        assert(e1.crystal == ee.crystal)

    # test passed
    print 'OK'

  def tst_from_pickle(self):
    from os.path import join
    import os

    os.environ['DIALS_REGRESSION'] = self.path

    # Get all the filenames
    filename1 = join(self.path, 'experiment_test_data', 'experiment_1.json')

    # Read all the experiment lists in
    el1 = ExperimentListFactory.from_json_file(filename1)

    # Pickle then load again
    el2 = self.pickle_then_unpickle(el1)

    # All the experiment lists should be the same length
    assert(len(el1) == 1)
    assert(len(el1) == len(el2))

    # Check all the models are the same
    for e1, e2 in zip(el1, el2):
      assert(e1.imageset is not None)
      assert(e1.beam is not None)
      assert(e1.detector is not None)
      assert(e1.goniometer is not None)
      assert(e1.scan is not None)
      assert(e1.crystal is not None)
      assert(e1.imageset == e2.imageset)
      assert(e1.beam == e2.beam)
      assert(e1.detector == e2.detector)
      assert(e1.goniometer == e2.goniometer)
      assert(e1.scan == e2.scan)
      assert(e1.crystal == e2.crystal)

    # test passed
    print 'OK'

  def tst_from_args(self):
    from os.path import join
    from glob import glob

    # Get all the filenames
    filenames = [
      join(self.path, 'experiment_test_data', 'experiment_1.json'),
      #join(self.path, 'experiment_test_data', 'experiment_2.json'),
      join(self.path, 'experiment_test_data', 'experiment_3.json'),
      join(self.path, 'experiment_test_data', 'experiment_4.json')]

    # Get the experiments from a list of filenames
    experiments = ExperimentListFactory.from_args(filenames)

    # Have 4 experiment
    assert(len(experiments) == 3)
    for i in range(3):
      assert(experiments[i].imageset is not None)
      assert(experiments[i].beam is not None)
      assert(experiments[i].detector is not None)
      assert(experiments[i].goniometer is not None)
      assert(experiments[i].scan is not None)

    # Test passed
    print 'OK'

  def tst_from_imageset(self):
    from dxtbx.imageset import ImageSet
    from dxtbx.model import Beam, Detector, Goniometer, Scan
    from dxtbx.model import Crystal
    from dxtbx.format.Format import Format

    imageset = Format.get_imageset(["filename.cbf"], as_imageset=True)
    imageset.set_beam(Beam(), 0)
    imageset.set_detector(Detector(), 0)

    crystal = Crystal(
      (1, 0, 0), (0, 1, 0), (0, 0, 1), space_group_symbol="P1")

    experiments = ExperimentListFactory.from_imageset_and_crystal(
      imageset, crystal)


    assert(len(experiments) == 1)
    assert(experiments[0].imageset is not None)
    assert(experiments[0].beam is not None)
    assert(experiments[0].detector is not None)
    assert(experiments[0].crystal is not None)

    print 'OK'

  def tst_from_sweep(self):
    from dxtbx.imageset import ImageSweep
    from dxtbx.model import Beam, Detector, Goniometer, Scan
    from dxtbx.model import Crystal
    from dxtbx.format.Format import Format

    filenames = ["filename_%01d.cbf" % (i+1) for i in range(0, 2)]

    imageset = Format.get_imageset(
      filenames,
      beam = Beam(),
      detector = Detector(),
      goniometer = Goniometer(),
      scan = Scan((1,2), (0,1)),
      as_sweep=True)

    crystal = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), space_group_symbol="P1")

    experiments = ExperimentListFactory.from_imageset_and_crystal(
      imageset, crystal)

    assert(len(experiments) == 1)
    assert(experiments[0].imageset is not None)
    assert(experiments[0].beam is not None)
    assert(experiments[0].detector is not None)
    assert(experiments[0].goniometer is not None)
    assert(experiments[0].scan is not None)
    assert(experiments[0].crystal is not None)

    print 'OK'

  def tst_from_datablock(self):
    from dxtbx.imageset import ImageSweep
    from dxtbx.model import Beam, Detector, Goniometer, Scan
    from dxtbx.datablock import DataBlockFactory
    from dxtbx.model import Crystal
    from dxtbx.format.Format import Format

    filenames = ["filename_%01d.cbf" % (i+1) for i in range(0, 2)]

    imageset = Format.get_imageset(
      filenames,
      beam = Beam(),
      detector = Detector(),
      goniometer = Goniometer(),
      scan = Scan((1,2), (0,1)),
      as_sweep=True)


    crystal = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), space_group_symbol="P1")

    datablock = DataBlockFactory.from_imageset(imageset)

    experiments = ExperimentListFactory.from_datablock_and_crystal(
      datablock, crystal)

    assert(len(experiments) == 1)
    assert(experiments[0].imageset is not None)
    assert(experiments[0].beam is not None)
    assert(experiments[0].detector is not None)
    assert(experiments[0].goniometer is not None)
    assert(experiments[0].scan is not None)
    assert(experiments[0].crystal is not None)

    print 'OK'

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


class TestExperimentListDumper(object):

  def __init__(self, path):
    self.path = path

  def run(self):
    self.tst_dump_formats()
    self.tst_dump_empty_sweep()
    self.tst_dump_with_lookup()
    self.tst_dump_with_bad_lookup()

  def tst_dump_formats(self):
    from uuid import uuid4
    from os.path import join
    import os

    os.environ['DIALS_REGRESSION'] = self.path

    # Get all the filenames
    filename1 = join(self.path, 'experiment_test_data', 'experiment_1.json')

    # Read all the experiment lists in
    elist1 = ExperimentListFactory.from_json_file(filename1)

    # Create the experiment list dumper
    dump = ExperimentListDumper(elist1)

    # Dump as JSON file and reload
    filename = 'temp%s.json' % uuid4().hex
    dump.as_json(filename)
    elist2 = ExperimentListFactory.from_json_file(filename)
    self.check(elist1, elist2)

    # Dump as split JSON file and reload
    filename = 'temp%s.json' % uuid4().hex
    dump.as_json(filename, split=True)
    elist2 = ExperimentListFactory.from_json_file(filename)
    self.check(elist1, elist2)

    # Dump as pickle and reload
    filename = 'temp%s.pickle' % uuid4().hex
    dump.as_pickle(filename)
    elist2 = ExperimentListFactory.from_pickle_file(filename)
    self.check(elist1, elist2)

  def tst_dump_empty_sweep(self):
    from dxtbx.imageset import ImageSweep
    from dxtbx.model import Beam, Detector, Goniometer, Scan
    from dxtbx.model import Crystal
    from uuid import uuid4
    from dxtbx.format.Format import Format

    filenames = ["filename_%01d.cbf" % (i+1) for i in range(0, 2)]

    imageset = Format.get_imageset(
      filenames,
      beam = Beam((1, 0, 0)),
      detector = Detector(),
      goniometer = Goniometer(),
      scan = Scan((1,2), (0.0, 1.0)),
      as_sweep=True)

    crystal = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), space_group_symbol="P1")

    experiments = ExperimentListFactory.from_imageset_and_crystal(
      imageset, crystal)

    dump = ExperimentListDumper(experiments)
    filename = 'temp%s.json' % uuid4().hex
    dump.as_json(filename)
    experiments2 = ExperimentListFactory.from_json_file(filename,
                                                        check_format=False)
    self.check(experiments, experiments2)

    print 'OK'

  def tst_dump_with_lookup(self):
    from dxtbx.imageset import ImageSweep
    from dxtbx.model import Beam, Detector, Goniometer, Scan
    from dxtbx.model import Crystal
    from uuid import uuid4
    import libtbx.load_env
    import os
    from os.path import join

    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    filename = join(dials_regression, "centroid_test_data",
                    "experiments_with_lookup.json")

    experiments = ExperimentListFactory.from_json_file(
      filename,
      check_format=True)

    imageset = experiments[0].imageset
    assert imageset.external_lookup.mask.data is not None
    assert imageset.external_lookup.gain.data is not None
    assert imageset.external_lookup.pedestal.data is not None
    assert imageset.external_lookup.mask.filename is not None
    assert imageset.external_lookup.gain.filename is not None
    assert imageset.external_lookup.pedestal.filename is not None
    assert imageset.external_lookup.mask.data.all_eq(True)
    assert imageset.external_lookup.gain.data.all_eq(1)
    assert imageset.external_lookup.pedestal.data.all_eq(0)

    dump = ExperimentListDumper(experiments)
    filename = 'temp%s.json' % uuid4().hex
    dump.as_json(filename)

    experiments = ExperimentListFactory.from_json_file(
      filename,
      check_format=True)

    imageset = experiments[0].imageset
    assert imageset.external_lookup.mask.data is not None
    assert imageset.external_lookup.gain.data is not None
    assert imageset.external_lookup.pedestal.data is not None
    assert imageset.external_lookup.mask.filename is not None
    assert imageset.external_lookup.gain.filename is not None
    assert imageset.external_lookup.pedestal.filename is not None
    assert imageset.external_lookup.mask.data.all_eq(True)
    assert imageset.external_lookup.gain.data.all_eq(1)
    assert imageset.external_lookup.pedestal.data.all_eq(0)

  def tst_dump_with_bad_lookup(self):
    from dxtbx.imageset import ImageSweep
    from dxtbx.model import Beam, Detector, Goniometer, Scan
    from dxtbx.model import Crystal
    from uuid import uuid4
    import libtbx.load_env
    import os
    from os.path import join

    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    filename = join(dials_regression, "centroid_test_data",
                    "experiments_with_bad_lookup.json")

    experiments = ExperimentListFactory.from_json_file(
      filename,
      check_format=False)

    imageset = experiments[0].imageset
    assert imageset.external_lookup.mask.data is None
    assert imageset.external_lookup.gain.data is None
    assert imageset.external_lookup.pedestal.data is None
    assert imageset.external_lookup.mask.filename is not None
    assert imageset.external_lookup.gain.filename is not None
    assert imageset.external_lookup.pedestal.filename is not None

    dump = ExperimentListDumper(experiments)
    filename = 'temp%s.json' % uuid4().hex
    dump.as_json(filename)

    experiments = ExperimentListFactory.from_json_file(
      filename,
      check_format=False)

    imageset = experiments[0].imageset
    assert imageset.external_lookup.mask.data is None
    assert imageset.external_lookup.gain.data is None
    assert imageset.external_lookup.pedestal.data is None
    assert imageset.external_lookup.mask.filename is not None
    assert imageset.external_lookup.gain.filename is not None
    assert imageset.external_lookup.pedestal.filename is not None

    print "OK"

  def check(self, el1, el2):

    # All the experiment lists should be the same length
    assert(len(el1) == 1)
    assert(len(el1) == len(el2))

    # Check all the models are the same
    for e1, e2 in zip(el1, el2):
      assert(e1.imageset is not None)
      assert(e1.beam is not None)
      assert(e1.detector is not None)
      assert(e1.goniometer is not None)
      assert(e1.scan is not None)
      assert(e1.crystal is not None)
      assert(e1.imageset == e2.imageset)
      assert(e1.beam == e2.beam)
      assert(e1.detector == e2.detector)
      assert(e1.goniometer == e2.goniometer)
      assert(e1.scan == e2.scan)
      assert(e1.crystal == e2.crystal)
    print 'OK'

class Test(object):
  def __init__(self):
    import libtbx

    if not libtbx.env.has_module("dials"):
      print 'Skipping: dials_regresson not configured'
      exit(0)
    if not libtbx.env.has_module("dials_regression"):
      print 'Skipping: dials_regresson not configured'
      exit(0)

    dials_regression = libtbx.env.dist_path('dials_regression')
    self.tst_experiment = TestExperiment(dials_regression)
    self.tst_list = TestExperimentList(dials_regression)
    self.tst_factory = TestExperimentListFactory(dials_regression)
    self.tst_dumper = TestExperimentListDumper(dials_regression)

  def run(self):
    self.tst_experiment.run()
    self.tst_list.run()
    self.tst_factory.run()
    self.tst_dumper.run()

if __name__ == '__main__':
  test = Test()
  test.run()
