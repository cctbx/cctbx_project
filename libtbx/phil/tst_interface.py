
from __future__ import absolute_import, division, print_function
from libtbx.test_utils import show_diff, Exception_expected
from libtbx.phil import interface
import libtbx.load_env
import libtbx.phil
from six.moves import cStringIO as StringIO
import sys

def exercise():
  master_phil = libtbx.phil.parse("""
refinement {
  input {
    pdb {
      file_name = None
        .type = path
        .multiple = True
    }
    sequence = None
      .type = path
      .style = seq_file
  }
  refine {
    strategy = *individual_sites *individual_adp *occupancies tls rigid_body
      .type = choice(multi=True)
    adp {
      tls = None
        .type = str
        .multiple = True
        .help = Selection for TLS group
    }
  }
  main {
    ncs = False
      .type = bool
      .help = This turns on NCS restraints
    ordered_solvent = False
      .type = bool
    number_of_macro_cycles = 3
      .type = int
    ias = False
      .type = bool
  }
  developer
    .expert_level = 3
  {
    place_elemental_ions = False
      .type = bool
  }
  gui {
    include scope libtbx.phil.interface.tracking_params
    output_dir = None
      .type = path
      .style = output_dir
  }
}
""", process_includes=True)
  refine_phil1 = libtbx.phil.parse("""
refinement {
  input {
    pdb {
      file_name = protein.pdb
      file_name = ligand.pdb
    }
  }
  refine {
    adp {
      tls = "chain A"
      tls = "chain B"
    }
  }
  main {
    ncs = True
    ordered_solvent = True
  }
}
""")
  refine_phil2_str = """
refinement {
  input {
    pdb {
      file_name = model1.pdb
    }
  }
  main {
    ncs = True
    ordered_solvent = False
  }
}"""
  refine_phil3 = libtbx.phil.parse("refinement.main.number_of_macro_cycles=5")
  refine_phil4_str = """
refinement.refine.adp.tls = None
"""
  i = libtbx.phil.interface.index(master_phil=master_phil,
    working_phil=refine_phil1,
    fetch_new=True)
  params = i.get_python_object()
  i.update(refine_phil2_str)
  # object retrieval
  pdb_phil = i.get_scope_by_name("refinement.input.pdb.file_name")
  assert len(pdb_phil) == 1
  os_phil = i.get_scope_by_name("refinement.main.ordered_solvent")
  assert os_phil.full_path() == "refinement.main.ordered_solvent"
  os_param = os_phil.extract()
  assert os_param == False
  params = i.get_python_object()
  assert len(params.refinement.refine.adp.tls) == 2
  # more updating, object extraction
  i.merge_phil(phil_object=refine_phil3)
  params = i.get_python_object()
  assert params.refinement.main.ncs == True
  assert params.refinement.main.ordered_solvent == False
  assert params.refinement.main.number_of_macro_cycles == 5
  assert params.refinement.input.pdb.file_name == ["model1.pdb"]
  i.merge_phil(phil_string=refine_phil4_str)
  params = i.get_python_object()
  assert len(params.refinement.refine.adp.tls) == 0
  phil1 = libtbx.phil.parse("""refinement.refine.strategy = *tls""")
  phil2 = libtbx.phil.parse("""refinement.input.pdb.file_name = ligand2.pdb""")
  i.save_param_file(
    file_name="tst_params.eff",
    sources=[phil1, phil2],
    extra_phil="refinement.main.ias = True",
    diff_only=True)
  params = i.get_python_from_file("tst_params.eff")
  assert params.refinement.refine.strategy == ["tls"]
  assert params.refinement.input.pdb.file_name == ["model1.pdb","ligand2.pdb"]
  assert params.refinement.main.ias == True
  i2 = i.copy(preserve_changes=False)
  params2 = i2.get_python_object()
  assert not params2.refinement.main.ncs
  i3 = i.copy(preserve_changes=True)
  params3 = i3.get_python_object()
  assert params3.refinement.main.ncs == True
  seq_file_def = i.get_seq_file_def_name()
  assert (seq_file_def == "refinement.input.sequence")

  # text searching (we can assume this will break quickly, but easily checked
  # by uncommenting the print statements)
  names = i.search_phil_text("macro_cycles", phil_name_only=True)
  assert len(names) == 1
  names = i.search_phil_text("elemental")
  assert (len(names) == 1)
  i.update("refinement.gui.output_dir=/var/tmp")
  i.update("refinement.gui.job_title=\"Hello, world!\"")
  assert (i.get_output_dir() == "/var/tmp")
  assert (i.get_job_title() == "Hello, world!")

  assert (libtbx.phil.interface.get_adjoining_phil_path(
    "refinement.input.xray_data.file_name", "labels") ==
    "refinement.input.xray_data.labels")

  master_phil = libtbx.phil.parse("""
pdb_in = None
  .type = path
  .short_caption = Input model
  .style = file_type:pdb input_file
pdb_out = None
  .type = path
  .style = file_type:pdb new_file
""")
  working_phil = master_phil.fetch(source=libtbx.phil.parse("""
pdb_in = foo.pdb
pdb_out = foo.modified.pdb
"""))
  i = libtbx.phil.interface.index(master_phil=master_phil,
    working_phil=working_phil,
    fetch_new=False)
  pdb_map = i.get_file_type_map("pdb")
  assert (pdb_map.get_param_names() == ['pdb_in'])
  assert (i.get_input_files() == [('foo.pdb', 'Input model', 'pdb_in')])
  # test captions
  master_phil = libtbx.phil.parse("""
my_options {
opt1 = *foo bar
  .type = choice
  .caption = Foo Bar
opt2 = *two_fofc fofc
  .type = choice(multi=True)
  .caption = 2mFo-DFc
}
""")
  try :
    libtbx.phil.interface.validate_choice_captions(master_phil)
  except AssertionError as e :
    assert (str(e) == "my_options.opt2")
  else :
    raise Exception_expected

# XXX sorry about the cross-import here, but I really need to test this on
# something large and complex
def exercise_2(verbose=False):
  if (not libtbx.env.has_module(name="phenix")):
    print("phenix module not available: skipping advanced tests")
    return
  from phenix.refinement import runtime
  import iotbx.phil
  from time import time
  phil_str = """
refinement.pdb_interpretation.secondary_structure.protein {
  helix {
    selection = "chain A and resseq 10:20"
  }
  helix {
    selection = "chain A and resseq 30:40"
  }
  helix {
    selection = "chain A and resseq 50:60"
  }
}
"""

  phil_str_2 = """
refinement.pdb_interpretation.secondary_structure.protein {
  helix {
    selection = "chain B and resseq 10:20"
  }
  helix {
    selection = "chain B and resseq 30:40"
  }
  helix {
    selection = "chain B and resseq 50:60"
  }
}
"""

  master_phil = runtime.master_phil()
  for phil_object in master_phil.objects:
    if phil_object.name == 'data_manager':
      master_phil.objects.remove(phil_object)
  i = interface.index(master_phil=master_phil,
    parse=iotbx.phil.parse)
  t1 = time()
  i.merge_phil(phil_string=phil_str)
  t2 = time()
  params = i.get_python_object()
  assert (params.refinement.pdb_interpretation.secondary_structure.\
      protein.helix[0].selection == "chain A and resseq 10:20")
  t3 = time()
  i.merge_phil(phil_string=phil_str_2,
    only_scope="refinement.pdb_interpretation.secondary_structure")
  t4 = time()
  params = i.get_python_object()
  assert (params.refinement.pdb_interpretation.secondary_structure.\
      protein.helix[0].selection == "chain B and resseq 10:20")
  scope = i.get_scope_by_name("refinement.pdb_interpretation.secondary_structure")
  params2 = scope.extract()
  assert (params2.protein.helix[0].selection == "chain B and resseq 10:20")
  if verbose :
    print("Merge with global fetch: %6.1fms" % ((t2-t1) * 1000))
    print("Merge with local fetch:  %6.1fms" % ((t4-t3) * 1000))
  i.merge_phil(phil_string="""
refinement.gui.migration.refinement.input.pdb.file_name = protein.pdb
refinement.gui.migration.refinement.input.pdb.file_name = ligand.pdb
refinement.gui.migration.refinement.input.monomers.file_name = ligand.cif
refinement.output.job_title = Test refinement run
""")
  names = i.search_phil_text("CIF")
  assert (set(names) == {
    'refinement.output.write_model_cif_file',
    'refinement.output.write_reflection_cif_file',
    'refinement.gui.migration.refinement.input.monomers.file_name',}
)

  expected_result = [
    ('protein.pdb', 'Input model (X-ray)', 'refinement.gui.migration.refinement.input.pdb.file_name'),
    ('ligand.pdb', 'Input model (X-ray)', 'refinement.gui.migration.refinement.input.pdb.file_name'),
    ('ligand.cif', 'Restraints (CIF)', 'refinement.gui.migration.refinement.input.monomers.file_name')]
  for f in i.get_input_files():
    assert f in expected_result
  assert len(i.get_input_files()) == len(expected_result)
  assert (i.get_job_title() == "Test refinement run")
  #
  # .style processing
  style = i.get_scope_style("refinement.refine.strategy")
  #assert (style.auto_launch_dialog == [
  #  'refinement.refine.sites.individual', 'refinement.refine.sites.individual',
  #  'refinement.refine.sites.rigid_body', 'refinement.refine.adp.individual',
  #  'refinement.refine.adp.group', 'refinement.refine.adp.tls',
  #  'refinement.refine.occupancies', 'refinement.refine.anomalous_scatterers'])
  assert (style.file_type is None)
  style = i.get_scope_style("refinement.gui.migration.refinement.input.xray_data.file_name")
  assert (style.get_list("file_type") == ["hkl"])
  assert (style.get_child_params() == {'fobs': 'labels',
    'd_max': 'low_resolution', 'd_min': 'high_resolution',
    'rfree_file': 'r_free_flags.file_name'})
  assert i.is_list_type("refinement.gui.migration.refinement.input.xray_data.labels")
  style = i.get_scope_style("refinement.gui.migration.refinement.input.xray_data.labels")
  assert (style.get_parent_params() == {"file_name" : "file_name"})
  file_map = i.get_file_type_map("pdb")
  expected_result = ['refinement.gui.migration.refinement.input.pdb.file_name',
                     'refinement.gui.migration.refinement.input.pdb.electron_file_name',
                     'refinement.gui.migration.refinement.input.pdb.neutron_file_name',
                     'refinement.reference_model.file']
  for p in file_map.get_multiple_params():
    assert p in expected_result, (p, expected_result)
  assert len(file_map.get_multiple_params()) == len(expected_result)
  assert (file_map.get_default_param() == "refinement.gui.migration.refinement.input.pdb.file_name")
  file_map = i.get_file_type_map("hkl")
  assert (file_map.get_overall_max_count() == 7)
  assert (len(file_map.get_multiple_params()) == 0)
  assert (file_map.get_max_count("refinement.gui.migration.refinement.input.xray_data.file_name") == 1)
  menu = i.get_menu_db()
  assert (len(menu.get_items()) > 15) # XXX ballpark (currently 17)
  submenu = menu.get_submenu("Atom_selections")
  assert (str(submenu.get_items()[0]) == "refinement.refine.sites")

def exercise_3():
  if (not libtbx.env.has_module(name="phaser")):
    print("phaser module not available: skipping advanced tests")
    return
  import iotbx.phil
  import phaser.phenix_interface
  master_phil = phaser.phenix_interface.master_phil()
  i = interface.index(master_phil=master_phil,
    parse=iotbx.phil.parse)
  i.merge_phil(phil_string="""\
phaser {
  hklin = "/Users/nat/Documents/beta-blip/beta_blip_P3221.mtz"
  labin = Fobs,Sigma
  composition {
    chain {
      sequence_file = "/Users/nat/Documents/beta-blip/beta.seq"
    }
    chain {
      sequence_file = "/Users/nat/Documents/beta-blip/blip.seq"
    }
  }
  ensemble {
    model_id = "beta"
    coordinates {
      pdb = "/Users/nat/Documents/beta-blip/beta.pdb"
    }
  }
  ensemble {
    model_id = "blip"
    coordinates {
      pdb = "/Users/nat/Documents/beta-blip/blip.pdb"
    }
  }
}
""")
  files_in = [
    ('/Users/nat/Documents/beta-blip/beta_blip_P3221.mtz', 'Data file',
      'phaser.hklin'),
    ('/Users/nat/Documents/beta-blip/blip.seq', 'Sequence file',
      'phaser.composition.chain.sequence_file'),
    ('/Users/nat/Documents/beta-blip/blip.pdb', 'Ensemble model',
      'phaser.ensemble.coordinates.pdb')
  ]
  assert (i.get_input_files() == files_in)
  i.save_param_file(
    file_name="phaser.eff",
    extra_phil="""
phaser.search {
  ensembles = beta
  copies = 1
}
phaser.search {
  ensembles = blip
  copies = 1
}""",
    replace_path="/Users/nat/Documents/beta-blip")
  i = interface.index(master_phil=master_phil,
    parse=iotbx.phil.parse)
  i.merge_phil(phil_file="phaser.eff")
  p = i.get_python_object().phaser
  assert (i.get_input_files() == files_in)
  assert (p.hklin == "/Users/nat/Documents/beta-blip/beta_blip_P3221.mtz")
  assert (len(p.search) == 2)
  assert (p.search[0].ensembles == ["beta"])
  # update file in-place (with variable substitution)
  interface.update_phil_file_paths(
    master_phil=master_phil,
    file_name="phaser.eff",
    old_path="/Users/nat/Documents/beta-blip",
    new_path="/Users/nat/Documents/projects/beta-blip",
    use_iotbx_parser=True)
  i = interface.index(master_phil=master_phil,
    parse=iotbx.phil.parse)
  i.merge_phil(phil_file="phaser.eff")
  p = i.get_python_object().phaser
  assert (p.hklin ==
    "/Users/nat/Documents/projects/beta-blip/beta_blip_P3221.mtz")
  # update file in-place, by modifying phil objects directly
  i.save_param_file(file_name="phaser2.eff")
  interface.update_phil_file_paths(
    master_phil=master_phil,
    file_name="phaser2.eff",
    old_path="/Users/nat/Documents/projects/beta-blip",
    new_path="/home/nat/projects/beta-blip",
    use_iotbx_parser=True)
  i = interface.index(master_phil=master_phil,
    parse=iotbx.phil.parse)
  i.merge_phil(phil_file="phaser2.eff")
  p = i.get_python_object().phaser
  assert (p.hklin == "/home/nat/projects/beta-blip/beta_blip_P3221.mtz")
  i.set_prefix("phaser")
  assert (i.get_full_path(".hklin") == "phaser.hklin")
  assert (i.get_scope_by_name(".keywords") is not None)
  # and now with Windows-style paths
  interface.update_phil_file_paths(
    master_phil=master_phil,
    file_name="phaser2.eff",
    old_path="/home/nat/projects/beta-blip",
    new_path="C:\\projects\\xtal\\beta-blip", # \x and \b are key here
    use_iotbx_parser=True)
  i = interface.index(master_phil=master_phil,
    parse=iotbx.phil.parse)
  i.merge_phil(phil_file="phaser2.eff")
  p = i.get_python_object().phaser
  # XXX obviously these are not entirely transferrable between Unix and
  # Windows - we need to caution users against this
  assert (p.hklin == "C:\\projects\\xtal\\beta-blip/beta_blip_P3221.mtz")

def exercise_adopt_phil():
  master_phil = libtbx.phil.parse("""\
scope1 {
  a = 1
    .type = int
  b = 2
    .type = int
}
""")
  working_phil = libtbx.phil.parse("""\
scope1.a = 3
scope1.b = 4
""")
  i = libtbx.phil.interface.index(master_phil=master_phil,
                                  working_phil=working_phil,
                                  fetch_new=True)
  params = i.get_python_object()
  other_master_phil = libtbx.phil.parse("""\
scope1 {
  c = 3
    .type = int
}
scope2 {
  subscope2 {
    d = 4
      .type = int
  }
  e = 5
    .type = int
}
""")
  i.adopt_phil(phil_object=other_master_phil)
  scope1 = i.get_scope_by_name("scope1")
  s = StringIO()
  scope1.show(out=s)
  assert not show_diff(s.getvalue(), """\
scope1 {
  a = 3
  b = 4
  c = 3
}
""")
  s = StringIO()
  i.working_phil.show(out=s)
  assert not show_diff(s.getvalue(), """\
scope1 {
  a = 3
  b = 4
  c = 3
}
scope2 {
  subscope2 {
    d = 4
  }
  e = 5
}
""")


if __name__ == "__main__" :
  exercise()
  exercise_2(verbose=("-v" in sys.argv[1:] or "--verbose" in sys.argv[1:]))
  exercise_3()
  exercise_adopt_phil()
  print("OK")

#---end
