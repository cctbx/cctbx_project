from __future__ import absolute_import, division, print_function

import cPickle as pickle
from glob import glob
import os

import pytest

from dxtbx.model import Experiment, ExperimentList
from dxtbx.model.experiment_list import ExperimentListFactory, \
  ExperimentListDumper, ExperimentListDict

def test_experiment_contains():
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
  assert b1 in e
  assert d1 in e
  assert g1 in e
  assert s1 in e
  assert c1 in e

  # Create a load of models that look the same but aren't
  b2 = Beam()
  d2 = Detector()
  g2 = Goniometer()
  s2 = Scan()
  c2 = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), space_group_symbol="P1")

  # Check experiment doesn't contain model
  assert b2 not in e
  assert d2 not in e
  assert g2 not in e
  assert s2 not in e
  assert c2 not in e

def test_experiment_equality():
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
  assert e1 == e2
  assert e1 != e3
  assert e2 != e3

def test_experiment_consistent(dials_regression):
  from dxtbx.imageset import ImageSetFactory
  from dxtbx.model import Scan

  # Create a sweep
  sweep_filenames = os.path.join(dials_regression, 'centroid_test_data', 'centroid*.cbf')
  sweep = ImageSetFactory.new(sorted(glob(sweep_filenames)))[0]

  # Create experiment with sweep and good scan
  e = Experiment(imageset=sweep, scan=sweep.get_scan())
  assert e.is_consistent()

  # Create experiment with sweep and defective scan
  scan = sweep.get_scan()
  scan.set_image_range((1, 1))
  e = Experiment(imageset=sweep, scan=scan)
  #assert not e.is_consistent()) # FIXME

  ## Create experiment with imageset and good scan
  #assert e.is_consistent()

  ## Create experiment with imageset and non-still scan
  #assert not e.is_consistent()

  ## Create experiment with imageset and scan with more than 1 image
  #assert not e.is_consistent()

  ## Create experiment with imageset and defective scan
  #assert not e.is_consistent()

def test_experimentlist_contains(experiment_list):
  from dxtbx.model import Beam, Detector, Goniometer, Scan

  # Check all the models are found
  for e in experiment_list:
    assert e.beam in experiment_list
    assert e.detector in experiment_list
    assert e.goniometer in experiment_list
    assert e.scan in experiment_list

  # Create some more models
  b = Beam()
  d = Detector()
  g = Goniometer()
  s = Scan()

  # Check that models not in are not found
  assert b not in experiment_list
  assert d not in experiment_list
  assert g not in experiment_list
  assert s not in experiment_list

# def test_experimentlist_index(experiment_list):

#   # Check the indices of exisiting experiments
#   assert experiment_list.index(experiment_list[0]) is 0
#   assert experiment_list.index(experiment_list[1]) is 1
#   assert experiment_list.index(experiment_list[2]) is 2
#   assert experiment_list.index(experiment_list[3]) is 1
#   assert experiment_list.index(experiment_list[4]) is 0

#   # Check index of non exisiting experiment
#   try:
#     experiment_list.index(Experiment())
#     assert False
#   except ValueError:
#     pass

def test_experimentlist_replace(experiment_list):
  # Get the models
  b = [e.beam for e in experiment_list]
  d = [e.detector for e in experiment_list]
  g = [e.goniometer for e in experiment_list]
  s = [e.scan for e in experiment_list]

  # Replace some models
  experiment_list.replace(b[0], b[1])
  assert experiment_list[0].beam is b[1]
  assert experiment_list[4].beam is b[1]

  # Replace again
  experiment_list[0].beam = b[0]
  experiment_list[4].beam = b[4]

def test_experimentlist_indices(experiment_list):
  from dxtbx.model import Beam, Detector, Goniometer, Scan

  # Get the models
  b = [e.beam for e in experiment_list]
  d = [e.detector for e in experiment_list]
  g = [e.goniometer for e in experiment_list]
  s = [e.scan for e in experiment_list]

  # Check indices of beams
  assert list(experiment_list.indices(b[0])) == [0, 4]
  assert list(experiment_list.indices(b[1])) == [1, 3]
  assert list(experiment_list.indices(b[2])) == [2]
  assert list(experiment_list.indices(b[3])) == [1, 3]
  assert list(experiment_list.indices(b[4])) == [0, 4]

  # Check indices of detectors
  assert list(experiment_list.indices(d[0])) == [0, 4]
  assert list(experiment_list.indices(d[1])) == [1, 3]
  assert list(experiment_list.indices(d[2])) == [2]
  assert list(experiment_list.indices(d[3])) == [1, 3]
  assert list(experiment_list.indices(d[4])) == [0, 4]

  # Check indices of goniometer
  assert list(experiment_list.indices(g[0])) == [0, 4]
  assert list(experiment_list.indices(g[1])) == [1, 3]
  assert list(experiment_list.indices(g[2])) == [2]
  assert list(experiment_list.indices(g[3])) == [1, 3]
  assert list(experiment_list.indices(g[4])) == [0, 4]

  # Check indices of scans
  assert list(experiment_list.indices(s[0])) == [0, 4]
  assert list(experiment_list.indices(s[1])) == [1, 3]
  assert list(experiment_list.indices(s[2])) == [2]
  assert list(experiment_list.indices(s[3])) == [1, 3]
  assert list(experiment_list.indices(s[4])) == [0, 4]

  # Check some models not in the list
  assert len(experiment_list.indices(Beam())) == 0
  assert len(experiment_list.indices(Detector())) == 0
  assert len(experiment_list.indices(Goniometer())) == 0
  assert len(experiment_list.indices(Scan())) == 0

def test_experimentlist_models(experiment_list):
  # Get all the unique models
  b = experiment_list.beams()
  d = experiment_list.detectors()
  g = experiment_list.goniometers()
  s = experiment_list.scans()

  # Check we have the expected number
  assert len(b) == 3
  assert len(d) == 3
  assert len(g) == 3
  assert len(s) == 3

  # Check we have the expected order
  assert b[0] == experiment_list[0].beam
  assert b[1] == experiment_list[1].beam
  assert b[2] == experiment_list[2].beam

  assert d[0] == experiment_list[0].detector
  assert d[1] == experiment_list[1].detector
  assert d[2] == experiment_list[2].detector

  assert g[0] == experiment_list[0].goniometer
  assert g[0] == experiment_list[0].goniometer
  assert g[1] == experiment_list[1].goniometer

  assert s[2] == experiment_list[2].scan
  assert s[1] == experiment_list[1].scan
  assert s[2] == experiment_list[2].scan

def test_experimentlist_to_dict(experiment_list):
  # Convert the list to a dictionary
  obj = experiment_list.to_dict()

  # Check this is the right object
  assert obj['__id__'] == 'ExperimentList'

  # Check length of items
  assert len(obj['experiment']) == 5
  assert len(obj['beam']) == 3
  assert len(obj['detector']) == 3
  assert len(obj['goniometer']) == 3
  assert len(obj['scan']) == 3

  # The expected models
  b = [0, 1, 2, 1, 0]
  d = [0, 1, 2, 1, 0]
  g = [0, 1, 2, 1, 0]
  s = [0, 1, 2, 1, 0]

  # Check all the experiments
  for i, eobj in enumerate(obj['experiment']):
    assert eobj['__id__'] == 'Experiment'
    assert eobj['beam'] == b[i]
    assert eobj['detector'] == d[i]
    assert eobj['goniometer'] == g[i]
    assert eobj['scan'] == s[i]

def test_experimentlist_where(experiment_list):
  for beam in experiment_list.beams():
    assert beam is not None
    for i in experiment_list.where(beam=beam):
      assert experiment_list[i].beam is beam
  for goniometer in experiment_list.goniometers():
    assert goniometer is not None
    for i in experiment_list.where(goniometer=goniometer):
      assert experiment_list[i].goniometer is goniometer
  for scan in experiment_list.scans():
    assert scan is not None
    for i in experiment_list.where(scan=scan):
      assert experiment_list[i].scan is scan
  for detector in experiment_list.detectors():
    assert detector is not None
    for i in experiment_list.where(detector=detector):
      assert experiment_list[i].detector is detector

@pytest.fixture
def experiment_list():
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
  ident = ["sausage", "eggs", "bacon", "toast", "beans"]

  # Populate with various experiments
  for i in range(5):
    experiments.append(Experiment(
      beam=b[i],
      detector=d[i],
      goniometer=g[i],
      scan=s[i],
      identifier=ident[i]))

  # Return the list of experiments
  return experiments

def test_experimentlist_factory_from_json(dials_regression):
  os.environ['DIALS_REGRESSION'] = dials_regression

  # Get all the filenames
  filename1 = os.path.join(dials_regression, 'experiment_test_data', 'experiment_1.json')
  filename2 = os.path.join(dials_regression, 'experiment_test_data', 'experiment_2.json')
  filename3 = os.path.join(dials_regression, 'experiment_test_data', 'experiment_3.json')
  filename4 = os.path.join(dials_regression, 'experiment_test_data', 'experiment_4.json')

  # Read all the experiment lists in
  el1 = ExperimentListFactory.from_json_file(filename1)
  #el2 = ExperimentListFactory.from_json_file(filename2)
  el3 = ExperimentListFactory.from_json_file(filename3)
  el4 = ExperimentListFactory.from_json_file(filename4)

  # All the experiment lists should be the same length
  assert len(el1) == 1
  #assert len(el1) == len(el2)
  assert len(el1) == len(el3)
  assert len(el1) == len(el4)

  # Check all the models are the same
  for e in zip(el1, el3, el4):
    e1 = e[0]
    assert e1.imageset is not None
    assert e1.beam is not None
    assert e1.detector is not None
    assert e1.goniometer is not None
    assert e1.scan is not None
    assert e1.crystal is not None
    for ee in e[1:]:
      assert e1.imageset == ee.imageset
      assert e1.beam == ee.beam
      assert e1.detector == ee.detector
      assert e1.goniometer == ee.goniometer
      assert e1.scan == ee.scan
      assert e1.crystal == ee.crystal

def test_experimentlist_factory_from_pickle(dials_regression):
  os.environ['DIALS_REGRESSION'] = dials_regression

  # Get all the filenames
  filename1 = os.path.join(dials_regression, 'experiment_test_data', 'experiment_1.json')

  # Read all the experiment lists in
  el1 = ExperimentListFactory.from_json_file(filename1)

  # Pickle then load again
  el2 = pickle.loads(pickle.dumps(el1))

  # All the experiment lists should be the same length
  assert len(el1) == 1
  assert len(el1) == len(el2)

  # Check all the models are the same
  for e1, e2 in zip(el1, el2):
    assert e1.imageset and e1.imageset == e2.imageset
    assert e1.beam and e1.beam == e2.beam
    assert e1.detector and e1.detector == e2.detector
    assert e1.goniometer and e1.goniometer == e2.goniometer
    assert e1.scan and e1.scan == e2.scan
    assert e1.crystal and e1.crystal == e2.crystal

def test_experimentlist_factory_from_args(dials_regression):
  pytest.importorskip('dials')
  os.environ['DIALS_REGRESSION'] = dials_regression

  # Get all the filenames
  filenames = [
    os.path.join(dials_regression, 'experiment_test_data', 'experiment_1.json'),
    #os.path.join(dials_regression, 'experiment_test_data', 'experiment_2.json'),
    os.path.join(dials_regression, 'experiment_test_data', 'experiment_3.json'),
    os.path.join(dials_regression, 'experiment_test_data', 'experiment_4.json')]

  # Get the experiments from a list of filenames
  experiments = ExperimentListFactory.from_args(filenames, verbose=True)

  # Have 4 experiment
  assert len(experiments) == 3
  for i in range(3):
    assert experiments[i].imageset is not None
    assert experiments[i].beam is not None
    assert experiments[i].detector is not None
    assert experiments[i].goniometer is not None
    assert experiments[i].scan is not None

def test_experimentlist_factory_from_imageset():
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


  assert len(experiments) == 1
  assert experiments[0].imageset is not None
  assert experiments[0].beam is not None
  assert experiments[0].detector is not None
  assert experiments[0].crystal is not None

def test_experimentlist_factory_from_sweep():
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

  assert len(experiments) == 1
  assert experiments[0].imageset is not None
  assert experiments[0].beam is not None
  assert experiments[0].detector is not None
  assert experiments[0].goniometer is not None
  assert experiments[0].scan is not None
  assert experiments[0].crystal is not None

def test_experimentlist_factory_from_datablock():
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

  assert len(experiments) == 1
  assert experiments[0].imageset is not None
  assert experiments[0].beam is not None
  assert experiments[0].detector is not None
  assert experiments[0].goniometer is not None
  assert experiments[0].scan is not None
  assert experiments[0].crystal is not None

def test_experimentlist_dumper_dump_formats(dials_regression, tmpdir):
  tmpdir.chdir()
  os.environ['DIALS_REGRESSION'] = dials_regression

  # Get all the filenames
  filename1 = os.path.join(dials_regression, 'experiment_test_data', 'experiment_1.json')

  # Read all the experiment lists in
  elist1 = ExperimentListFactory.from_json_file(filename1)

  # Create the experiment list dumper
  dump = ExperimentListDumper(elist1)

  # Dump as JSON file and reload
  filename = 'temp1.json'
  dump.as_json(filename)
  elist2 = ExperimentListFactory.from_json_file(filename)
  check(elist1, elist2)

  # Dump as split JSON file and reload
  filename = 'temp2.json'
  dump.as_json(filename, split=True)
  elist2 = ExperimentListFactory.from_json_file(filename)
  check(elist1, elist2)

  # Dump as pickle and reload
  filename = 'temp.pickle'
  dump.as_pickle(filename)
  elist2 = ExperimentListFactory.from_pickle_file(filename)
  check(elist1, elist2)

def test_experimentlist_dumper_dump_scan_varying(dials_regression, tmpdir):
  tmpdir.chdir()
  os.environ['DIALS_REGRESSION'] = dials_regression

  # Get all the filenames
  filename1 = os.path.join(dials_regression, 'experiment_test_data', 'experiment_1.json')

  # Read the experiment list in
  elist1 = ExperimentListFactory.from_json_file(filename1)

  # Make trivial scan-varying models
  crystal = elist1[0].crystal
  beam = elist1[0].beam
  crystal.set_A_at_scan_points([crystal.get_A()] * 5)
  beam.set_s0_at_scan_points([beam.get_s0()] * 5)

  # Create the experiment list dumper
  dump = ExperimentListDumper(elist1)

  # Dump as JSON file and reload
  filename = 'temp.json'
  dump.as_json(filename)
  elist2 = ExperimentListFactory.from_json_file(filename)
  check(elist1, elist2)

def test_experimentlist_dumper_dump_empty_sweep(tmpdir):
  tmpdir.chdir()
  from dxtbx.model import Beam, Detector, Goniometer, Scan
  from dxtbx.model import Crystal
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
  filename = 'temp.json'
  dump.as_json(filename)
  experiments2 = ExperimentListFactory.from_json_file(filename,
                                                      check_format=False)
  check(experiments, experiments2)

def test_experimentlist_dumper_dump_with_lookup(dials_regression, tmpdir):
  tmpdir.chdir()
  from dxtbx.model import Beam, Detector, Goniometer, Scan
  from dxtbx.model import Crystal

  filename = os.path.join(dials_regression, "centroid_test_data",
                  "experiments_with_lookup.json")

  experiments = ExperimentListFactory.from_json_file(
    filename, check_format=True)

  imageset = experiments[0].imageset
  assert not imageset.external_lookup.mask.data.empty()
  assert not imageset.external_lookup.gain.data.empty()
  assert not imageset.external_lookup.pedestal.data.empty()
  assert imageset.external_lookup.mask.filename is not None
  assert imageset.external_lookup.gain.filename is not None
  assert imageset.external_lookup.pedestal.filename is not None
  assert imageset.external_lookup.mask.data.tile(0).data().all_eq(True)
  assert imageset.external_lookup.gain.data.tile(0).data().all_eq(1)
  assert imageset.external_lookup.pedestal.data.tile(0).data().all_eq(0)

  dump = ExperimentListDumper(experiments)
  filename = 'temp.json'
  dump.as_json(filename)

  experiments = ExperimentListFactory.from_json_file(
    filename,
    check_format=True)

  imageset = experiments[0].imageset
  assert not imageset.external_lookup.mask.data.empty()
  assert not imageset.external_lookup.gain.data.empty()
  assert not imageset.external_lookup.pedestal.data.empty()
  assert imageset.external_lookup.mask.filename is not None
  assert imageset.external_lookup.gain.filename is not None
  assert imageset.external_lookup.pedestal.filename is not None
  assert imageset.external_lookup.mask.data.tile(0).data().all_eq(True)
  assert imageset.external_lookup.gain.data.tile(0).data().all_eq(1)
  assert imageset.external_lookup.pedestal.data.tile(0).data().all_eq(0)

def test_experimentlist_dumper_dump_with_bad_lookup(dials_regression, tmpdir):
  tmpdir.chdir()
  from dxtbx.model import Beam, Detector, Goniometer, Scan
  from dxtbx.model import Crystal

  filename = os.path.join(dials_regression, "centroid_test_data",
                  "experiments_with_bad_lookup.json")

  experiments = ExperimentListFactory.from_json_file(
    filename, check_format=False)

  imageset = experiments[0].imageset
  assert imageset.external_lookup.mask.data.empty()
  assert imageset.external_lookup.gain.data.empty()
  assert imageset.external_lookup.pedestal.data.empty()
  assert imageset.external_lookup.mask.filename is not None
  assert imageset.external_lookup.gain.filename is not None
  assert imageset.external_lookup.pedestal.filename is not None

  dump = ExperimentListDumper(experiments)
  filename = 'temp.json'
  dump.as_json(filename)

  experiments = ExperimentListFactory.from_json_file(
    filename, check_format=False)

  imageset = experiments[0].imageset
  assert imageset.external_lookup.mask.data.empty()
  assert imageset.external_lookup.gain.data.empty()
  assert imageset.external_lookup.pedestal.data.empty()
  assert imageset.external_lookup.mask.filename is not None
  assert imageset.external_lookup.gain.filename is not None
  assert imageset.external_lookup.pedestal.filename is not None

def test_experimentlist_with_identifiers():
  from dxtbx.model import Beam, Detector, Goniometer, Scan

  # Initialise a list of experiments
  experiments = ExperimentList()

  experiments.append(Experiment(
    beam=Beam(s0=(0,0,-1)),
    detector=Detector(),
    identifier="bacon"))

  experiments.append(Experiment(
    beam=Beam(s0=(0,0,-1)),
    detector=Detector(),
    identifier="sausage"))

  with pytest.raises(Exception):
    experiments.append(Experiment(
      beam=Beam(),
      detector=Detector(),
      identifier="bacon"))

  d = experiments.to_dict()
  e2 = ExperimentListDict(d).decode()

  assert experiments[0].identifier == e2[0].identifier
  assert experiments[1].identifier == e2[1].identifier

def check(el1, el2):
  # All the experiment lists should be the same length
  assert len(el1) == 1
  assert len(el1) == len(el2)

  # Check all the models are the same
  for e1, e2 in zip(el1, el2):
    assert e1.imageset and e1.imageset == e2.imageset
    assert e1.beam and e1.beam == e2.beam
    assert e1.detector is not None and e1.detector == e2.detector
    assert e1.goniometer and e1.goniometer == e2.goniometer
    assert e1.scan and e1.scan == e2.scan
    assert e1.crystal and e1.crystal == e2.crystal
    assert e1.identifier == e2.identifier
