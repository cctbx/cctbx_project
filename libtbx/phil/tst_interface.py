
import os, sys
import libtbx.phil
from libtbx.phil import interface
from libtbx.utils import Sorry
from libtbx import easy_pickle

def exercise () :
  master_phil = libtbx.phil.parse("""
refinement {
  input {
    pdb {
      file_name = None
        .type = path
        .multiple = True
    }
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
  ncs {
    restraint_group
      .multiple = True
      .optional = True
      .short_caption = Restraint group
    {
      reference = None
        .type = str
        .help = Reference selection for restraint group
      selection = None
        .type = str
        .multiple = True
        .optional = False
        .help = Restrained selection
    }
  }
}
""")
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
  ncs {
    restraint_group {
      reference = "chain A"
      selection = "chain B"
      selection = "chain C"
      selection = "chain D"
    }
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
  ncs {
    restraint_group {
      reference = "chain A"
      selection = "chain B"
    }
    restraint_group {
      reference = "chain C"
      selection = "chain D"
    }
  }
}"""
  refine_phil3 = libtbx.phil.parse("refinement.main.number_of_macro_cycles=5")
  refine_phil4_str = """
refinement.refine.adp.tls = None
refinement.ncs.restraint_group {
  reference = "chain C and resseq 1:100"
  selection = "chain D and resseq 1:100"
}"""
  i = libtbx.phil.interface.index(master_phil=master_phil,
                                  working_phil=refine_phil1)
  params = i.get_python_object()
  assert len(params.refinement.ncs.restraint_group) == 1
  i.update(refine_phil2_str)
  # object retrieval
  ncs_phil = i.get_scope_by_name("refinement.ncs.restraint_group")
  assert len(ncs_phil) == 2
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
  assert len(params.refinement.ncs.restraint_group) == 2
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

  # text searching (we can assume this will break quickly, but easily checked
  # by uncommenting the print statements)
  names = i.search_phil_text("Restraint group", match_all=True,
    labels_only=True)
  assert len(names) == 0
  names = i.search_phil_text("Restraint group", match_all=True,
    labels_only=False)
  assert len(names) == 3
  names = i.search_phil_text("selection group", match_all=True,
    labels_only=False)
  assert len(names) == 3

  assert (libtbx.phil.interface.get_adjoining_phil_path(
    "refinement.input.xray_data.file_name", "labels") ==
    "refinement.input.xray_data.labels")

# XXX sorry about the cross-import here, but I
def exercise_2 (verbose=False) :
  try :
    from phenix.refinement import runtime
    import iotbx.phil
  except ImportError :
    print "PHENIX sources not found, skipping advanced tests"
    return False
  from time import time
  phil_str = """
refinement.secondary_structure {
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
refinement.secondary_structure {
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

  phil_str_3 = """
refinement.ncs.restraint_group {
  reference = "chain A"
  selection = "chain B"
}
"""
  phil_str_4 = """
refinement.ncs.restraint_group {
  reference = "chain C"
  selection = "chain D"
}"""

  master_phil = runtime.master_phil
  i = interface.index(master_phil=master_phil,
    parse=iotbx.phil.parse)
  t1 = time()
  i.merge_phil(phil_string=phil_str)
  t2 = time()
  params = i.get_python_object()
  assert (params.refinement.secondary_structure.helix[0].selection ==
    "chain A and resseq 10:20")
  t3 = time()
  i.merge_phil(phil_string=phil_str_2,
    only_scope="refinement.secondary_structure")
  t4 = time()
  params = i.get_python_object()
  assert (params.refinement.secondary_structure.helix[0].selection ==
    "chain B and resseq 10:20")
  scope = i.get_scope_by_name("refinement.secondary_structure")
  params2 = scope.extract()
  assert (params2.helix[0].selection == "chain B and resseq 10:20")
  if verbose :
    print "Merge with global fetch: %6.1fms" % ((t2-t1) * 1000)
    print "Merge with local fetch:  %6.1fms" % ((t4-t3) * 1000)
  i.merge_phil(phil_string=phil_str_3,
    only_scope="refinement.ncs.restraint_group")
  params = i.get_python_object()
  assert (params.refinement.ncs.restraint_group[0].reference == "chain A")
  i.merge_phil(phil_string=phil_str_4,
    only_scope="refinement.ncs.restraint_group")
  params = i.get_python_object()
  assert (params.refinement.ncs.restraint_group[0].reference == "chain C")

if __name__ == "__main__" :
  exercise()
  exercise_2(verbose=("-v" in sys.argv[1:] or "--verbose" in sys.argv[1:]))
  print "OK"

#---end
