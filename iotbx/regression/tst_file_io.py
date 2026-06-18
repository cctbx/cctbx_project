import os

def test_maps():
  from iotbx.file_io import (
    any_file_type, data_manager_type, extension_to_datatype, CIF_SENTINEL)
  # datatype -> any_file type
  assert any_file_type['model'] == 'pdb'
  assert any_file_type['miller_array'] == 'hkl'
  assert any_file_type['real_map'] == 'ccp4_map'
  assert any_file_type['restraint'] == 'cif'
  # reverse, with hkl/mtz aliases pointing at the miller_array parent
  assert data_manager_type['pdb'] == 'model'
  assert data_manager_type['hkl'] == 'miller_array'
  assert data_manager_type['mtz'] == 'miller_array'
  assert data_manager_type['cif'] == 'restraint'
  # extension map
  assert extension_to_datatype['pdb'] == 'model'
  assert extension_to_datatype['mtz'] == 'miller_array'
  assert extension_to_datatype['ccp4'] == 'real_map'
  assert extension_to_datatype['json'] == 'json'
  assert extension_to_datatype['yaml'] == 'yaml'
  assert extension_to_datatype['ncs'] == 'ncs_spec'
  assert extension_to_datatype['cif'] == CIF_SENTINEL
  print('test_maps OK')

def test_backward_compat_reexport():
  # existing importers must keep working (enabled in a later task once
  # iotbx.data_manager re-exports the maps)
  from iotbx.data_manager import data_manager_type as dmt
  from iotbx.file_io import data_manager_type as fmt
  assert dmt is fmt
  print('test_backward_compat_reexport OK')

def _make_fixtures(tmpdir):
  '''Create one real file of each datatype; return {key: path}. Keys that are
  not datatype names (model_cif) are noted in the tests.'''
  import iotbx.pdb
  from cctbx import crystal, miller
  from cctbx.array_family import flex
  from iotbx.map_model_manager import map_model_manager
  paths = {}

  pdb_str = '''CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1      2
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
END
'''
  p = os.path.join(tmpdir, 'm.pdb')
  with open(p, 'w') as f: f.write(pdb_str)
  paths['model'] = p

  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  cif_str = pdb_inp.construct_hierarchy().as_mmcif_string(
    crystal_symmetry=pdb_inp.crystal_symmetry())
  p = os.path.join(tmpdir, 'm.cif')
  with open(p, 'w') as f: f.write(cif_str)
  paths['model_cif'] = p

  cs = crystal.symmetry((10, 10, 10, 90, 90, 90), 'P1')
  ms = miller.set(cs, flex.miller_index([(1, 0, 0), (0, 1, 0), (0, 0, 1)]),
                  anomalous_flag=False)
  ma = miller.array(ms, data=flex.double([1., 2., 3.]))
  p = os.path.join(tmpdir, 'r.mtz')
  ma.as_mtz_dataset(column_root_label='F').mtz_object().write(p)
  paths['miller_array'] = p

  mmm = map_model_manager()
  mmm.generate_map()
  p = os.path.join(tmpdir, 'm.map')
  mmm.map_manager().write_map(p)
  paths['real_map'] = p

  restraint_cif = '''data_comp_list
loop_
_chem_comp.id
_chem_comp.name
TST Test
data_comp_TST
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
 TST C1 C
 TST C2 C
'''
  p = os.path.join(tmpdir, 'r.cif')
  with open(p, 'w') as f: f.write(restraint_cif)
  paths['restraint'] = p

  p = os.path.join(tmpdir, 's.seq')
  with open(p, 'w') as f: f.write('>1SAR\nDVSGTVCLSALPPEATDTLNLIASDGPFPYSQDG\n')
  paths['sequence'] = p

  p = os.path.join(tmpdir, 'p.phil')
  with open(p, 'w') as f: f.write('foo {\n  bar = 1\n    .type = int\n}\n')
  paths['phil'] = p

  import json as _json
  p = os.path.join(tmpdir, 'd.json')
  with open(p, 'w') as f: f.write(_json.dumps({'a': [1, 2, 3]}))
  paths['json'] = p

  p = os.path.join(tmpdir, 'd.yaml')
  with open(p, 'w') as f: f.write('a:\n  - 1\n  - 2\n')
  paths['yaml'] = p

  # real ncs_spec layout: free-text header lines BEFORE new_ncs_group, so the
  # whole-prefix scan is exercised (not a synthetic file starting at the token).
  ncs_spec_str = (
    '\nSummary of NCS information\nThu May 21 14:48:36 2020\n\n\n'
    'new_ncs_group\nnew_operator\n\n'
    'rota_matrix    1.0000    0.0000    0.0000\n'
    'rota_matrix    0.0000    1.0000    0.0000\n'
    'rota_matrix    0.0000    0.0000    1.0000\n'
    'tran_orth     0.0000    0.0000    0.0000\n\n'
    'center_orth   14.4035    7.4690    0.2556\n'
    'CHAIN A\nRMSD 0\nMATCHING 3\n  RESSEQ 1:3\n\n')
  p = os.path.join(tmpdir, 'n.ncs_spec')
  with open(p, 'w') as f: f.write(ncs_spec_str)
  paths['ncs_spec'] = p

  return paths

def test_detection_by_extension():
  from iotbx.file_io import get_file_type
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    assert get_file_type(paths['model']) == 'model'
    assert get_file_type(paths['model_cif']) == 'model'
    assert get_file_type(paths['miller_array']) == 'miller_array'
    assert get_file_type(paths['real_map']) == 'real_map'
    assert get_file_type(paths['restraint']) == 'restraint'
    assert get_file_type(paths['sequence']) == 'sequence'
    assert get_file_type(paths['phil']) == 'phil'
    assert get_file_type(paths['json']) == 'json'
    assert get_file_type(paths['yaml']) == 'yaml'
  finally:
    shutil.rmtree(tmp)
  print('test_detection_by_extension OK')

def test_detection_valid_types_and_missing():
  from iotbx.file_io import get_file_type
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    assert get_file_type(paths['json'], valid_types=['model']) is None
    assert get_file_type(paths['json'], valid_types=['json']) == 'json'
    assert get_file_type(os.path.join(tmp, 'nope.pdb')) is None
  finally:
    shutil.rmtree(tmp)
  print('test_detection_valid_types_and_missing OK')

def test_detection_misnamed_binary():
  from iotbx.file_io import get_file_type
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    import shutil as _sh
    fake = os.path.join(tmp, 'fake.pdb')
    _sh.copyfile(paths['miller_array'], fake)
    assert get_file_type(fake) == 'miller_array'
    fake_cif = os.path.join(tmp, 'fake.cif')
    _sh.copyfile(paths['miller_array'], fake_cif)
    assert get_file_type(fake_cif) != 'model'
  finally:
    shutil.rmtree(tmp)
  print('test_detection_misnamed_binary OK')

def test_detection_never_raises():
  from iotbx.file_io import get_file_type
  assert get_file_type('/definitely/does/not/exist.pdb') is None
  print('test_detection_never_raises OK')

def test_read_file_attribute_contracts():
  from iotbx.file_io import read_file
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    # model: file_object exposes .input for mmtbx.model.manager
    r = read_file(paths['model'], file_type='model')
    assert r.data_type == 'model'
    assert hasattr(r.file_object, 'input')
    assert r.file_content is r.file_object  # any_file-compatible alias
    # model from mmCIF
    r = read_file(paths['model_cif'], file_type='model')
    assert hasattr(r.file_object, 'input')
    # restraint: file_object.model() returns the cif model
    r = read_file(paths['restraint'], file_type='restraint')
    assert r.data_type == 'restraint'
    assert r.file_object.model() is not None
    # miller_array / real_map: parsed object present
    assert read_file(paths['miller_array'], file_type='miller_array').file_object is not None
    assert read_file(paths['real_map'], file_type='real_map').file_object is not None
    # json / yaml return parsed data
    assert read_file(paths['json'], file_type='json').file_object == {'a': [1, 2, 3]}
    assert read_file(paths['yaml'], file_type='yaml').file_object == {'a': [1, 2]}
  finally:
    shutil.rmtree(tmp)
  print('test_read_file_attribute_contracts OK')

def test_read_file_restraint_engines():
  from iotbx.file_io import read_file
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    for engine in ['xcif', 'ucif']:
      r = read_file(paths['restraint'], file_type='restraint', cif_engine=engine)
      assert r.file_object.model() is not None
  finally:
    shutil.rmtree(tmp)
  print('test_read_file_restraint_engines OK')

def test_read_file_detects_when_no_type():
  from iotbx.file_io import read_file
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    r = read_file(paths['model'])  # no file_type -> detect
    assert r.data_type == 'model'
  finally:
    shutil.rmtree(tmp)
  print('test_read_file_detects_when_no_type OK')

def test_parity_with_any_file():
  from iotbx.file_io import get_file_type, data_manager_type
  from iotbx.file_reader import any_file
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    for key in ['model', 'model_cif', 'miller_array', 'real_map',
                'restraint', 'sequence', 'phil', 'ncs_spec']:
      f = paths[key]
      expected = data_manager_type.get(any_file(f).file_type)
      got = get_file_type(f)
      assert got == expected, (key, got, expected)
  finally:
    shutil.rmtree(tmp)
  print('test_parity_with_any_file OK')

def test_read_file_missing():
  '''A missing file reports "Couldn't find the file", consistent across
  json/yaml/model (not a misleading "not a valid ..." message).'''
  from iotbx.file_io import read_file
  from libtbx.utils import Sorry
  for ft in ['json', 'yaml', 'model']:
    try:
      read_file('/no/such/dir/missing_zzz.' + ft, file_type=ft)
      raise AssertionError('expected Sorry for missing %s' % ft)
    except Sorry as e:
      assert "Couldn't find the file" in str(e), (ft, str(e))
  print('test_read_file_missing OK')

def test_dm_process_file_and_get_file_type():
  from iotbx.data_manager import DataManager
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    dm = DataManager(['model', 'miller_array', 'real_map', 'restraint',
                      'sequence', 'phil', 'map_coefficients', 'ncs_spec',
                      'json', 'yaml'])
    assert dm.process_file(paths['model']) == ['model']
    assert paths['model'] in dm.get_model_names()
    assert dm.process_file(paths['miller_array']) == ['miller_array']
    assert dm.process_file(paths['restraint']) == ['restraint']
    assert dm.process_file(paths['json']) == ['json']
    assert dm.process_file(paths['yaml']) == ['yaml']
    # idempotent
    assert dm.process_file(paths['model']) == ['model']
    # get_file_type scoped to supported datatypes
    assert dm.get_file_type(paths['model']) == 'model'
  finally:
    shutil.rmtree(tmp)
  print('test_dm_process_file_and_get_file_type OK')

def test_dm_process_file_cif_engine():
  '''process_file forwards cif_engine to the restraint reader (both engines).'''
  from iotbx.data_manager import DataManager
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    for engine in ['xcif', 'ucif']:
      dm = DataManager(['restraint'])
      assert dm.process_file(paths['restraint'], cif_engine=engine) == ['restraint']
      assert dm.get_restraint(paths['restraint']) is not None
  finally:
    shutil.rmtree(tmp)
  print('test_dm_process_file_cif_engine OK')

def test_dm_json_opt_in():
  from iotbx.data_manager import DataManager
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    dm = DataManager()  # default_datatypes (no json/yaml)
    assert dm.get_file_type(paths['json']) is None
    assert dm.process_file(paths['json']) == []
  finally:
    shutil.rmtree(tmp)
  print('test_dm_json_opt_in OK')

def test_process_files_loop():
  # exercise the same detect+dispatch loop cli_parser.process_files now uses,
  # covering BOTH the recognized (-> found) and unsupported (-> unused) branches
  from iotbx.data_manager import DataManager
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    # this DataManager intentionally does NOT support 'restraint'
    dm = DataManager(['model', 'miller_array', 'json', 'yaml',
                      'map_coefficients'])
    found = {}
    unused = []
    for f in [paths['model'], paths['miller_array'], paths['json'],
              paths['yaml'], paths['restraint']]:
      dts = dm.process_file(f)
      if dts:
        for dt in dts: found[dt] = f
      else:
        unused.append(f)
    assert 'model' in found and 'miller_array' in found
    assert 'json' in found and 'yaml' in found
    # the restraint file is not a supported datatype here -> unused
    assert unused == [paths['restraint']], unused
  finally:
    shutil.rmtree(tmp)
  print('test_process_files_loop OK')

def test_detection_short_binary():
  '''A short (<1000-byte) binary file with a recognized text extension must be
  caught by the text gate, not trusted by extension. detect_binary_file only
  renders a verdict after its monitor window, which _is_binary caps at the
  prefix length.'''
  from iotbx.file_io import get_file_type
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    p = os.path.join(tmp, 'garbage.pdb')
    with open(p, 'wb') as f:
      f.write(b'\x00\x01\x02\x03not a real file\xff\xfe')
    assert get_file_type(p) != 'model', get_file_type(p)
  finally:
    shutil.rmtree(tmp)
  print('test_detection_short_binary OK')

def test_compression_roundtrip():
  '''smart_open round-trips for the new compressors (.xz/.lzma, and .zst when
  zstandard is installed), and a compressed file is still detected via the
  compression-aware splitext + decompressed prefix peek.'''
  from libtbx import smart_open
  from iotbx.file_io import get_file_type
  import tempfile, shutil, gzip
  tmp = tempfile.mkdtemp()
  try:
    data = b'line 1\nline 2\nthe end\n'
    for ext in ['.xz', '.lzma']:
      fn = os.path.join(tmp, 'rt' + ext)
      with smart_open.for_writing(fn, mode='wb') as fh:
        fh.write(data)
      with smart_open.for_reading(fn, mode='rb') as fh:
        assert fh.read() == data, ext
    try:
      import zstandard  # noqa: F401
      fn = os.path.join(tmp, 'rt.zst')
      with smart_open.for_writing(fn, mode='wb') as fh:
        fh.write(data)
      with smart_open.for_reading(fn, mode='rb') as fh:
        assert fh.read() == data
    except ImportError:
      print('  (zstd round-trip skipped: zstandard not installed)')
    # a gzipped model is still detected as 'model' (compression-aware splitext +
    # decompressed prefix binary check)
    paths = _make_fixtures(tmp)
    gz = os.path.join(tmp, 'm.pdb.gz')
    with open(paths['model'], 'rb') as src, gzip.open(gz, 'wb') as dst:
      dst.write(src.read())
    assert get_file_type(gz) == 'model', get_file_type(gz)
  finally:
    shutil.rmtree(tmp)
  print('test_compression_roundtrip OK')

def test_dm_gzipped_model_population():
  '''A gzipped model loads end-to-end through DataManager.process_file
  (detection is compression-aware and iotbx.pdb reads .gz natively).'''
  from iotbx.data_manager import DataManager
  import tempfile, shutil, gzip
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    gz = os.path.join(tmp, 'm.pdb.gz')
    with open(paths['model'], 'rb') as src, gzip.open(gz, 'wb') as dst:
      dst.write(src.read())
    dm = DataManager(['model'])
    assert dm.process_file(gz) == ['model']
    assert gz in dm.get_model_names()
    assert dm.get_model(gz) is not None
  finally:
    shutil.rmtree(tmp)
  print('test_dm_gzipped_model_population OK')

def test_compressed_population_all_formats():
  '''Compressed files load end-to-end for every datatype: text readers handle
  .gz/.Z natively and .xz/.lzma/.zst via the smart_open gzip_mode fix; .bz2 and
  the binary MTZ/MRC readers go through read_file's temp-decompress fallback
  (which also fixes the pre-existing .mtz.gz gap).'''
  from iotbx.data_manager import DataManager
  import tempfile, shutil, gzip, gc
  try:
    import bz2, lzma
  except ImportError:
    print('Skipping test_compressed_population_all_formats (bz2/lzma not available)')
    return
  tmp = tempfile.mkdtemp()
  def _compress(src, dst, opener):
    with open(src, 'rb') as s, opener(dst, 'wb') as d:
      d.write(s.read())
  try:
    paths = _make_fixtures(tmp)
    # model: .xz (native via gzip_mode fix) and .bz2 (temp-decompress fallback)
    for ext, opener in [('.xz', lzma.open), ('.bz2', bz2.open)]:
      mp = paths['model'] + ext
      _compress(paths['model'], mp, opener)
      dm = DataManager(['model'])
      assert dm.process_file(mp) == ['model'], (ext, dm.process_file(mp))
      assert dm.get_model(mp) is not None
    # miller_array: .gz (MTZ reader handles NO compression natively -> temp) + .xz
    for ext, opener in [('.gz', gzip.open), ('.xz', lzma.open)]:
      rp = paths['miller_array'] + ext
      _compress(paths['miller_array'], rp, opener)
      dm = DataManager(['miller_array'])
      assert dm.process_file(rp) == ['miller_array'], (ext, dm.process_file(rp))
      assert dm.get_miller_array(rp) is not None
    # real_map: .xz (mrcfile handles .gz but not .xz -> temp)
    mp = paths['real_map'] + '.xz'
    _compress(paths['real_map'], mp, lzma.open)
    dm = DataManager(['real_map'])
    assert dm.process_file(mp) == ['real_map'], dm.process_file(mp)
    assert dm.get_real_map(mp) is not None
    # json: .gz (json uses plain open -> temp-decompress fallback)
    jp = paths['json'] + '.gz'
    _compress(paths['json'], jp, gzip.open)
    dm = DataManager(['json'])
    assert dm.process_file(jp) == ['json'], dm.process_file(jp)
    assert dm.get_json(jp) == {'a': [1, 2, 3]}
  finally:
    # The compressed-map read leaves an mrcfile handle open until the next
    # collection (mrcfile cannot read .xz directly, so the failed direct read's
    # MrcFile object is held by a reference cycle on its traceback). Force the
    # collection so the handle is closed before shutil.rmtree deletes m.map.xz;
    # otherwise Windows raises "the process cannot access the file". POSIX is
    # unaffected (an open handle does not block unlink there).
    gc.collect()
    shutil.rmtree(tmp)
  print('test_compressed_population_all_formats OK')

def test_process_files_loop_skips_unparseable():
  '''A stray text file with a recognized extension (e.g. notes.pdb whose content
  is not a model) must not abort the detect+dispatch loop. Content verification
  recognizes it is not a model and returns None, so dm.process_file reports it as
  unused ([]); the loop (as cli_parser.process_files does) skips it and keeps
  going. (A file that detects as a type but then fails to parse still raises
  Sorry, which the loop also catches -- see the try/except below.)'''
  from iotbx.data_manager import DataManager
  from libtbx.utils import Sorry
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    stray = os.path.join(tmp, 'notes.pdb')
    with open(stray, 'w') as f:
      f.write('just some notes, not a model\n')
    dm = DataManager(['model', 'miller_array'])
    # content verification recognizes the stray file is not a model, so
    # process_file reports it as unused ([]) rather than raising
    assert dm.process_file(stray) == []
    # the cli_parser loop catches that and continues
    found, unused = {}, []
    for f in [stray, paths['model'], paths['miller_array']]:
      try:
        dts = dm.process_file(f)
      except Sorry:
        dts = []
      if dts:
        for dt in dts: found[dt] = f
      else:
        unused.append(f)
    assert 'model' in found and 'miller_array' in found  # valid files still loaded
    assert stray in unused
  finally:
    shutil.rmtree(tmp)
  print('test_process_files_loop_skips_unparseable OK')

def test_cif_datatypes():
  '''_cif_datatypes parses a CIF (xcif) and reports every datatype present;
  a model's bare _chem_comp / _reflns stats must not be read as restraint /
  miller_array.'''
  from iotbx.file_io.detection import _cif_datatypes
  from iotbx.file_io import get_file_type
  import iotbx.pdb, iotbx.cif, tempfile, shutil
  from cctbx import miller
  from cctbx.array_family import flex
  tmp = tempfile.mkdtemp()
  try:
    pdb = ('CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1      2\n'
           'ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N\n'
           'ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C\n'
           'END\n')
    pin = iotbx.pdb.input(source_info=None, lines=pdb)
    cs = pin.crystal_symmetry()
    mcif = pin.construct_hierarchy().as_mmcif_string(crystal_symmetry=cs)
    ma = miller.array(
      miller.set(cs, flex.miller_index([(1,0,0),(0,1,0),(0,0,1),(1,1,0)]),
                 anomalous_flag=False),
      data=flex.double([1., 2., 3., 4.])).set_observation_type_xray_amplitude()
    refl_block = ma.as_cif_block(array_type='meas')
    def write(name, text):
      p = os.path.join(tmp, name)
      with open(p, 'w') as f: f.write(text)
      return p
    # model-only (carries a bare _chem_comp list)
    pm = write('model.cif', mcif)
    assert _cif_datatypes(pm) == {'model'}, _cif_datatypes(pm)
    assert get_file_type(pm) == 'model'
    # model + injected _reflns statistics -> still just model (not miller_array)
    pm2 = write('model_reflns.cif',
                mcif + '\nloop_\n_reflns.d_resolution_high\n_reflns.number_obs\n 1.80 1234\n')
    assert _cif_datatypes(pm2) == {'model'}, _cif_datatypes(pm2)
    # reflections-only
    cm = iotbx.cif.model.cif(); cm['refl'] = refl_block
    pr = write('refl.cif', str(cm))
    assert _cif_datatypes(pr) == {'miller_array'}, _cif_datatypes(pr)
    # restraint (monomer library)
    pres = write('restr.cif',
                 'data_comp_list\nloop_\n_chem_comp.id\nTST\ndata_comp_TST\n'
                 'loop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n'
                 '_chem_comp_atom.type_symbol\n TST C1 C\n')
    assert _cif_datatypes(pres) == {'restraint'}, _cif_datatypes(pres)
    assert get_file_type(pres) == 'restraint'
    # combined model + reflections
    combo = iotbx.cif.reader(input_string=mcif).model()
    combo[list(combo.keys())[0]].update(refl_block)
    pc = write('combined.cif', str(combo))
    assert _cif_datatypes(pc) == {'model', 'miller_array'}, _cif_datatypes(pc)
    # get_file_type returns the precedence primary (miller_array)
    assert get_file_type(pc) == 'miller_array', get_file_type(pc)
  finally:
    shutil.rmtree(tmp)
  print('test_cif_datatypes OK')

def test_dm_combined_cif():
  '''A CIF with model + reflections loads as every supported type, each object.'''
  from iotbx.data_manager import DataManager
  import iotbx.pdb, iotbx.cif, tempfile, shutil
  from cctbx import miller
  from cctbx.array_family import flex
  tmp = tempfile.mkdtemp()
  try:
    pdb = ('CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1      2\n'
           'ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N\n'
           'ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C\n'
           'END\n')
    pin = iotbx.pdb.input(source_info=None, lines=pdb)
    cs = pin.crystal_symmetry()
    mcif = pin.construct_hierarchy().as_mmcif_string(crystal_symmetry=cs)
    ma = miller.array(
      miller.set(cs, flex.miller_index([(1,0,0),(0,1,0),(0,0,1),(1,1,0)]),
                 anomalous_flag=False),
      data=flex.double([1.,2.,3.,4.])).set_observation_type_xray_amplitude()
    combo = iotbx.cif.reader(input_string=mcif).model()
    combo[list(combo.keys())[0]].update(ma.as_cif_block(array_type='meas'))
    p = os.path.join(tmp, 'combined.cif')
    with open(p, 'w') as f: f.write(str(combo))
    # supports both -> loads as both, each with its object
    dm = DataManager(['model', 'miller_array'])
    assert sorted(dm.process_file(p)) == ['miller_array', 'model']
    assert dm.get_model(p) is not None
    assert dm.get_miller_array(p) is not None
    # model-only DM -> only the model loads (reflections ignored)
    dm2 = DataManager(['model'])
    assert dm2.process_file(p) == ['model']
    assert dm2.get_model(p) is not None
  finally:
    shutil.rmtree(tmp)
  print('test_dm_combined_cif OK')

def test_text_reclassification():
  '''A text file whose content is PHIL but whose extension says ncs_spec is
  detected as phil (content wins). This is the tst_ncs_1.py .ncs-PHIL bug.'''
  from iotbx.file_io import get_file_type
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    phil_body = ('refinement.pdb_interpretation.ncs_group {\n'
                 '  reference = chain A and resseq 20:28\n'
                 '  selection = chain B and resseq 20:28\n}\n')
    ncs = os.path.join(tmp, 'override.ncs')
    with open(ncs, 'w') as f: f.write(phil_body)
    assert get_file_type(ncs) == 'phil', get_file_type(ncs)
    # a PHIL body misnamed .pdb also reclassifies to phil
    pdb = os.path.join(tmp, 'override.pdb')
    with open(pdb, 'w') as f: f.write(phil_body)
    assert get_file_type(pdb) == 'phil', get_file_type(pdb)
  finally:
    shutil.rmtree(tmp)
  print('test_text_reclassification OK')

def test_text_valid_types():
  '''Reclassification respects valid_types: a real ncs_spec is dropped on a
  refine-like datatype set (no ncs_spec), while a .ncs of PHIL is kept as phil.'''
  from iotbx.file_io import get_file_type
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  refine = ['model', 'phil', 'miller_array', 'restraint']
  try:
    paths = _make_fixtures(tmp)
    assert get_file_type(paths['ncs_spec'], valid_types=refine) is None
    ncs = os.path.join(tmp, 'override.ncs')
    with open(ncs, 'w') as f:
      f.write('refinement.pdb_interpretation.ncs_group {\n'
              '  reference = chain A\n  selection = chain B\n}\n')
    assert get_file_type(ncs, valid_types=refine) == 'phil'
  finally:
    shutil.rmtree(tmp)
  print('test_text_valid_types OK')

def test_dat_sequence_behavior_change():
  '''Documented change: a numeric .dat is no longer typed sequence-by-extension.
  The sniffer misses (not residue-only) and detection defers to any_file, which
  does not recognize it -> None. A residue-only .dat is still sequence.'''
  from iotbx.file_io import get_file_type
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    dat = os.path.join(tmp, 'cols.dat')
    with open(dat, 'w') as f: f.write('1.0 2.0 3.0\n4.0 5.0 6.0\n')
    assert get_file_type(dat) is None, get_file_type(dat)
    seqdat = os.path.join(tmp, 'seq.dat')
    with open(seqdat, 'w') as f: f.write('DVSGTVCLSALPPEATDT\n')
    assert get_file_type(seqdat) == 'sequence'
  finally:
    shutil.rmtree(tmp)
  print('test_dat_sequence_behavior_change OK')

def test_dm_ncs_phil_end_to_end():
  '''DataManager.process_file on a .ncs file of PHIL loads it as phil (not a
  dropped/Sorry ncs_spec) under a refine-like datatype set.'''
  from iotbx.data_manager import DataManager
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    ncs = os.path.join(tmp, 'override.ncs')
    with open(ncs, 'w') as f:
      f.write('refinement.pdb_interpretation.ncs_group {\n'
              '  reference = chain A\n  selection = chain B\n}\n')
    dm = DataManager(['model', 'phil', 'miller_array', 'restraint'])
    result = dm.process_file(ncs)
    assert result == ['phil'], result
    assert ncs in dm.get_phil_names()
  finally:
    shutil.rmtree(tmp)
  print('test_dm_ncs_phil_end_to_end OK')

def test_process_file_returns_empty_on_read_failure():
  '''A detection false-positive -- a file that sniffs to a datatype but whose
  type-specific reader then rejects it -- is reported as unused ([]), not raised,
  matching the combined-CIF branch and the process_file contract.'''
  from iotbx.data_manager import DataManager
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    # an ATOM record (so it sniffs as model) but invalid coordinate content
    # that the model reader rejects (any_file does not see a model)
    weird = os.path.join(tmp, 'weird.pdb')
    with open(weird, 'w') as f:
      f.write('ATOM  this is not valid pdb coordinate data at all\n')
    dm = DataManager(['model', 'phil'])
    assert dm.process_file(weird) == []
  finally:
    shutil.rmtree(tmp)
  print('test_process_file_returns_empty_on_read_failure OK')

def test_cif_compression_detection():
  '''A compressed CIF (.cif.gz/.xz/...) is detected via iotbx.cif routing through
  smart_open, matching the uncompressed CIF -- not silently dropped. Regression
  for the new compressors missing from iotbx.cif._XCIF_COMPRESSED_SUFFIXES.'''
  from iotbx.file_io import get_file_type
  from iotbx.data_manager import DataManager
  import tempfile, shutil, gzip, lzma
  tmp = tempfile.mkdtemp()
  try:
    cif = os.path.join(tmp, 'lig.cif')
    with open(cif, 'w') as f:
      f.write('data_comp_list\nloop_\n_chem_comp.id\n_chem_comp.name\nTST Test\n'
              'data_comp_TST\nloop_\n_chem_comp_atom.comp_id\n'
              '_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n'
              ' TST C1 C\n TST C2 C\n')
    with open(cif, 'rb') as f:
      raw = f.read()
    # .gz already worked; .xz exercises the fix (same code path as .lzma/.zst)
    for ext, opener in [('.gz', gzip.open), ('.xz', lzma.open)]:
      p = cif + ext
      with opener(p, 'wb') as f: f.write(raw)
      assert get_file_type(p) == 'restraint', (ext, get_file_type(p))
      assert DataManager(['restraint']).process_file(p) == ['restraint'], ext
  finally:
    shutil.rmtree(tmp)
  print('test_cif_compression_detection OK')

def test_corrupt_compressed_no_warning():
  '''A corrupt compressed file is detected as None WITHOUT the unexpected-error
  warning, for both a bad header (LZMAError) and a valid-header-but-truncated
  stream (EOFError) -- decompression errors are caught, not mislabeled as bugs.'''
  from iotbx.file_io import get_file_type
  import tempfile, shutil, warnings, gzip
  tmp = tempfile.mkdtemp()
  try:
    cases = []
    # bad header
    bad = os.path.join(tmp, 'bad.pdb.xz')
    with open(bad, 'wb') as f: f.write(b'this is not valid xz compressed data')
    cases.append(bad)
    # valid gzip header but a truncated body -> EOFError on the prefix read
    trunc = os.path.join(tmp, 'trunc.pdb.gz')
    with gzip.open(trunc, 'wb') as f: f.write(b'ATOM line\n' * 5000)
    with open(trunc, 'r+b') as f: f.truncate(30)
    cases.append(trunc)
    for p in cases:
      with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        result = get_file_type(p)
      assert result is None, (p, result)
      assert len(w) == 0, (p, [str(x.message) for x in w])
  finally:
    shutil.rmtree(tmp)
  print('test_corrupt_compressed_no_warning OK')

def test_reader_decompress_error_messages():
  '''read_file gives clear errors when decompression cannot proceed: a .zip with
  an explicit type is reported as unsupported (not "could not be decompressed"),
  and a corrupt .Z surfaces a decompression failure naming the file (#7, #8).'''
  from iotbx.file_io import read_file
  from libtbx.utils import Sorry
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    z = os.path.join(tmp, 'lig.cif.zip')
    with open(z, 'wb') as f: f.write(b'PK\x03\x04 not really a zip archive')
    try:
      read_file(z, file_type='restraint')
      raise AssertionError('expected Sorry for .zip')
    except Sorry as e:
      assert 'ZIP' in str(e), str(e)
      assert 'could not be decompressed' not in str(e), str(e)
    cz = os.path.join(tmp, 'bad.cif.Z')
    with open(cz, 'wb') as f: f.write(b'not valid LZW .Z data at all')
    try:
      read_file(cz, file_type='restraint')
      raise AssertionError('expected Sorry for corrupt .Z')
    except Sorry as e:
      assert 'bad.cif.Z' in str(e), str(e)
  finally:
    shutil.rmtree(tmp)
  print('test_reader_decompress_error_messages OK')

def test_read_file_map_coefficients_and_ncs_spec():
  '''read_file populates the less-common datatypes: map_coefficients (forced hkl)
  and ncs_spec, and DataManager loads an ncs_spec end-to-end.'''
  from iotbx.file_io import read_file
  from iotbx.data_manager import DataManager
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    mc = read_file(paths['miller_array'], file_type='map_coefficients')
    assert mc.data_type == 'map_coefficients' and mc.file_object is not None
    ncs = read_file(paths['ncs_spec'], file_type='ncs_spec')
    assert ncs.data_type == 'ncs_spec' and ncs.file_object is not None
    assert DataManager(['ncs_spec']).process_file(paths['ncs_spec']) == ['ncs_spec']
  finally:
    shutil.rmtree(tmp)
  print('test_read_file_map_coefficients_and_ncs_spec OK')

def test_detection_real_map_magic_mismatch():
  '''A .map file without the CCP4 "MAP " magic is not detected as real_map: the
  magic-byte sniff rejects it and detection falls back.'''
  from iotbx.file_io import get_file_type
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    p = os.path.join(tmp, 'notreally.map')
    with open(p, 'w') as f: f.write('this is plainly not a CCP4 map file\n')
    assert get_file_type(p) != 'real_map', get_file_type(p)
  finally:
    shutil.rmtree(tmp)
  print('test_detection_real_map_magic_mismatch OK')

def test_dot_z_population():
  '''A valid .Z file reads end-to-end through the gunzip-subprocess branch.
  Skipped when the 'compress' utility is unavailable to create a .Z.'''
  from iotbx.file_io import read_file
  import tempfile, shutil, subprocess
  tmp = tempfile.mkdtemp()
  try:
    src = _make_fixtures(tmp)['model']
    z = src + '.Z'
    try:
      with open(z, 'wb') as out:
        subprocess.run(['compress', '-c', src], stdout=out,
                       stderr=subprocess.DEVNULL, check=True)
    except (OSError, subprocess.CalledProcessError):
      print('test_dot_z_population (skipped: compress unavailable)')
      return
    assert hasattr(read_file(z, file_type='model').file_object, 'input')
  finally:
    shutil.rmtree(tmp)
  print('test_dot_z_population OK')

def test_direct_read_phil():
  '''read_file parses phil with the direct reader (iotbx.phil), not any_file:
  FileIOResult.reader is None (the any_file fallback would set it).'''
  from iotbx.file_io import read_file
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    r = read_file(paths['phil'], file_type='phil')
    assert r.reader is None, 'expected the direct reader, not the any_file fallback'
    assert r.data_type == 'phil'
    assert len(r.file_object.objects) > 0      # a parsed phil scope
  finally:
    shutil.rmtree(tmp)
  print('test_direct_read_phil OK')

def test_direct_read_phil_malformed_falls_back():
  '''A malformed phil must not leak the parser's RuntimeError: the direct reader's
  parse error is caught and routed to the any_file fallback, which raises a clean
  Sorry.'''
  from iotbx.file_io import read_file
  from libtbx.utils import Sorry
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    bad = os.path.join(tmp, 'bad.phil')
    with open(bad, 'w') as f: f.write('bar 1\n')   # missing '=' -> phil syntax error
    try:
      read_file(bad, file_type='phil')
      raise AssertionError('expected Sorry for malformed phil')
    except Sorry:
      pass            # clean Sorry, not a leaked RuntimeError
  finally:
    shutil.rmtree(tmp)
  print('test_direct_read_phil_malformed_falls_back OK')

def test_direct_read_model():
  '''read_file parses model with iotbx.pdb.hierarchy directly; the .input the
  mmtbx.model.manager needs is preserved, and reader is None (direct path).'''
  from iotbx.file_io import read_file
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    r = read_file(paths['model'], file_type='model')
    assert r.reader is None, 'expected the direct reader, not the any_file fallback'
    assert hasattr(r.file_object, 'input')     # mmtbx.model.manager(model_input=...)
    rc = read_file(paths['model_cif'], file_type='model')
    assert rc.reader is None and hasattr(rc.file_object, 'input')
  finally:
    shutil.rmtree(tmp)
  print('test_direct_read_model OK')

def test_direct_read_model_wrong_type_falls_back():
  '''Reading a non-model file as model raises a clean Sorry (the direct reader's
  failure is caught and routed to the any_file fallback), not a leaked exception.'''
  from iotbx.file_io import read_file
  from libtbx.utils import Sorry
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    try:
      read_file(paths['phil'], file_type='model')   # a phil file is not a model
      raise AssertionError('expected Sorry')
    except Sorry:
      pass
  finally:
    shutil.rmtree(tmp)
  print('test_direct_read_model_wrong_type_falls_back OK')

def test_direct_read_real_map():
  '''read_file reads a CCP4/MRC map directly via iotbx.map_manager.'''
  from iotbx.file_io import read_file
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    r = read_file(paths['real_map'], file_type='real_map')
    assert r.reader is None, 'expected the direct reader, not the any_file fallback'
    assert r.file_object.map_data() is not None    # a map_manager
  finally:
    shutil.rmtree(tmp)
  print('test_direct_read_real_map OK')

def test_direct_read_real_map_wrong_type_falls_back():
  '''Reading a non-map file as real_map raises a clean Sorry (failure routed to
  the any_file fallback), not a leaked exception.'''
  from iotbx.file_io import read_file
  from libtbx.utils import Sorry
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    try:
      read_file(paths['phil'], file_type='real_map')
      raise AssertionError('expected Sorry')
    except Sorry:
      pass
  finally:
    shutil.rmtree(tmp)
  print('test_direct_read_real_map_wrong_type_falls_back OK')

def test_direct_read_reflections():
  '''read_file reads reflections directly via any_reflection_file, for both the
  miller_array and map_coefficients datatypes (both -> the reflection reader).'''
  from iotbx.file_io import read_file
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    for dt in ['miller_array', 'map_coefficients']:
      r = read_file(paths['miller_array'], file_type=dt)
      assert r.reader is None, (dt, 'expected the direct reader')
      assert r.data_type == dt
      assert r.file_object.file_type() is not None   # any_reflection_file object
  finally:
    shutil.rmtree(tmp)
  print('test_direct_read_reflections OK')

def test_direct_read_reflections_wrong_type_falls_back():
  '''Reading a non-reflection file as miller_array raises a clean Sorry (failure
  routed to the any_file fallback), not a leaked exception.'''
  from iotbx.file_io import read_file
  from libtbx.utils import Sorry
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    try:
      read_file(paths['phil'], file_type='miller_array')
      raise AssertionError('expected Sorry')
    except Sorry:
      pass
  finally:
    shutil.rmtree(tmp)
  print('test_direct_read_reflections_wrong_type_falls_back OK')

def test_direct_read_ncs_spec():
  '''read_file reads an ncs_spec directly via mmtbx.ncs.'''
  from iotbx.file_io import read_file
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    r = read_file(paths['ncs_spec'], file_type='ncs_spec')
    assert r.reader is None, 'expected the direct reader, not the any_file fallback'
    assert r.file_object.max_operators() > 0     # an mmtbx.ncs.ncs object
  finally:
    shutil.rmtree(tmp)
  print('test_direct_read_ncs_spec OK')

def test_direct_read_ncs_spec_wrong_type_falls_back():
  '''Reading a non-ncs file as ncs_spec raises a clean Sorry (failure routed to
  the any_file fallback), not a leaked exception.'''
  from iotbx.file_io import read_file
  from libtbx.utils import Sorry
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    try:
      read_file(paths['phil'], file_type='ncs_spec')
      raise AssertionError('expected Sorry')
    except Sorry:
      pass
  finally:
    shutil.rmtree(tmp)
  print('test_direct_read_ncs_spec_wrong_type_falls_back OK')

def test_direct_read_sequence():
  '''read_file reads sequences directly via any_sequence_format.'''
  from iotbx.file_io import read_file
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    paths = _make_fixtures(tmp)
    r = read_file(paths['sequence'], file_type='sequence')
    assert r.reader is None, 'expected the direct reader, not the any_file fallback'
    assert len(r.file_object) > 0                # a list of sequence objects
    assert hasattr(r.file_object[0], 'sequence')
  finally:
    shutil.rmtree(tmp)
  print('test_direct_read_sequence OK')

def test_direct_read_sequence_gaps_fall_back():
  '''A gapped (alignment) sequence is rejected by the direct reader and routes to
  the any_file fallback as a clean Sorry (any_file reads gaps as "aln", not seq).'''
  from iotbx.file_io import read_file
  from libtbx.utils import Sorry
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  try:
    gapped = os.path.join(tmp, 'aln.fasta')
    with open(gapped, 'w') as f: f.write('>aln\nDVSG--TVCLS--ALPP\n')
    try:
      read_file(gapped, file_type='sequence')
      raise AssertionError('expected Sorry')
    except Sorry:
      pass
  finally:
    shutil.rmtree(tmp)
  print('test_direct_read_sequence_gaps_fall_back OK')

def test_any_file_fallback_engages():
  '''When a direct reader fails, read_file falls back to any_file (the
  deprecation-window safety net): the object is still correct and
  FileIOResult.reader is the any_file_input (non-None).'''
  from iotbx.file_io import read_file
  import iotbx.file_io.reader as _reader
  from libtbx.utils import Sorry
  import tempfile, shutil
  tmp = tempfile.mkdtemp()
  orig = _reader._read_direct
  def _boom(filename, datatype, force=False):
    raise Sorry('forced direct-reader failure')
  try:
    paths = _make_fixtures(tmp)
    _reader._read_direct = _boom
    r = read_file(paths['model'], file_type='model')
    assert r.reader is not None, 'expected the any_file fallback to be used'
    assert hasattr(r.file_object, 'input')      # object still correct via fallback
  finally:
    _reader._read_direct = orig
    shutil.rmtree(tmp)
  print('test_any_file_fallback_engages OK')

if __name__ == '__main__':
  test_maps()
  test_backward_compat_reexport()
  test_cif_datatypes()
  test_detection_by_extension()
  test_detection_valid_types_and_missing()
  test_detection_misnamed_binary()
  test_detection_short_binary()
  test_detection_never_raises()
  test_corrupt_compressed_no_warning()
  test_text_reclassification()
  test_text_valid_types()
  test_dat_sequence_behavior_change()
  test_dm_ncs_phil_end_to_end()
  test_read_file_attribute_contracts()
  test_direct_read_phil()
  test_direct_read_phil_malformed_falls_back()
  test_direct_read_model()
  test_direct_read_model_wrong_type_falls_back()
  test_direct_read_real_map()
  test_direct_read_real_map_wrong_type_falls_back()
  test_direct_read_reflections()
  test_direct_read_reflections_wrong_type_falls_back()
  test_direct_read_ncs_spec()
  test_direct_read_ncs_spec_wrong_type_falls_back()
  test_direct_read_sequence()
  test_direct_read_sequence_gaps_fall_back()
  test_any_file_fallback_engages()
  test_read_file_restraint_engines()
  test_read_file_detects_when_no_type()
  test_parity_with_any_file()
  test_read_file_missing()
  test_dm_process_file_and_get_file_type()
  test_process_file_returns_empty_on_read_failure()
  test_dm_process_file_cif_engine()
  test_dm_combined_cif()
  test_cif_compression_detection()
  test_dm_json_opt_in()
  test_process_files_loop()
  test_process_files_loop_skips_unparseable()
  test_compression_roundtrip()
  test_dm_gzipped_model_population()
  test_compressed_population_all_formats()
  test_reader_decompress_error_messages()
  test_read_file_map_coefficients_and_ncs_spec()
  test_detection_real_map_magic_mismatch()
  test_dot_z_population()
  print('OK')
