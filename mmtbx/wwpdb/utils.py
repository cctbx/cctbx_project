
"""
Utilities for re-refining (or otherwise working with) structures downloaded
directly from the PDB.
"""

# XXX MARKED_FOR_DELETION_OLEG
# Reason: outdated unused code.


from __future__ import absolute_import, division, print_function
from libtbx import easy_run
from libtbx.utils import null_out, Sorry
import os
import sys

def get_program(pdb_file):
  assert 0
  from iotbx.pdb import remark_3_interpretation
  with open(pdb_file) as f:
    lines = f.readlines()
  program = program_full = None
  for line in lines :
    if (line.startswith("REMARK   3")) and ("PROGRAM" in line):
      program = remark_3_interpretation.get_program(line)
      if (program is not None):
        program_full = line.split(":")[1].strip()
        break
  if (program == "PHENIX"):
    program = "PHENIX.REFINE"
  return program, program_full

def fetch_pdb_data(
    pdb_id,
    pdb_dir=None,
    sf_dir=None,
    log=None,
    verbose=False):
  """
  Copy data from local repository if defined and available, or download it
  from the PDB, and run cif_as_mtz.
  """
  assert 0
  from mmtbx.command_line import fetch_pdb
  from mmtbx.command_line import cif_as_mtz
  if (log is None):
    if (verbose) : log = sys.stdout
    else : log = null_out()
  pdb_file = "%s.pdb" % pdb_id
  cif_file = "%s-sf.cif" % pdb_id
  mtz_file = "%s.mtz" % pdb_id
  fetch_pdb.run2(args=[pdb_id], log=log)
  assert (os.path.isfile(pdb_file))
  fetch_pdb.run2(args=["-x", pdb_id], log=log)
  if (not os.path.isfile("%s-sf.cif" % pdb_id)):
    raise Sorry("Structure factors are not available for %s." % pdb_id)
  cif_as_mtz.run(args=[
      cif_file,
      "--symmetry=%s" % pdb_file,
      "--merge",
      "--output_file_name=%s" % mtz_file])
  if (not os.path.isfile(mtz_file)):
    raise RuntimeError("Missing %s!\ncif_as_mtz stderr:\n%s" %
      (mtz_file, "\n".join(import_out.stderr_lines)))
  return os.path.abspath(pdb_file), os.path.abspath(mtz_file)

def find_data_arrays(mtz_file, log=None, merge_anomalous=False,
    preferred_labels=None, crystal_symmetry=None):
  """
  Guess an appropriate data array to use for refinement, plus optional
  Hendrickson-Lattman coefficients and R-free flags if present.
  """
  assert 0
  from iotbx import reflection_file_utils
  from iotbx.file_reader import any_file
  if (log is None) : log = sys.stdout
  phases = data = flags = flag_value = None
  hkl_in = any_file(mtz_file, force_type="hkl")
  hkl_server = hkl_in.file_server
  # always use anomalous data if available!  also, prefer amplitudes over
  # intensities if possible, as they may already be on an absolute scale
  data_arrays = hkl_server.get_xray_data(
    file_name               = None,
    labels                  = None,
    ignore_all_zeros        = False,
    parameter_scope         = "",
    return_all_valid_arrays = True,
    minimum_score           = 4,
    prefer_amplitudes       = True,
    prefer_anomalous        = True)
  if (len(data_arrays) == 0):
    raise Sorry("No data arrays found in %s." % mtz_file)
  data = data_arrays[0]
  if (preferred_labels is not None):
    for array in data_arrays :
      array_labels = array.info().label_string()
      if (array_label== preferred_labels):
        data = array
        break
    else :
      raise Sorry("Can't find label string '%s'!" % preferred_labels)
  else :
    print("Defaulting to using %s" % data.info().label_string(), file=log)
  hl_arrays = hkl_server.get_experimental_phases(
    file_name               = None,
    labels                  = None,
    ignore_all_zeros        = True,
    parameter_scope         = "",
    return_all_valid_arrays = True,
    minimum_score           = 1)
  if (len(hl_arrays) > 0):
    phases = hl_arrays[0]
  flags_and_values = hkl_server.get_r_free_flags(
    file_name=None,
    label=None,
    test_flag_value=None,
    disable_suitability_test=False,
    parameter_scope="",
    return_all_valid_arrays=True,
    minimum_score=1)
  if (len(flags_and_values) > 0):
    flags, flag_value = flags_and_values[0]
  if (crystal_symmetry is not None):
    data = data.customized_copy(crystal_symmetry=crystal_symmetry)
    if (flags is not None):
      flags = flags.customized_copy(crystal_symmetry=crystal_symmetry)
    if (phases is not None):
      phases = phases.customized_copy(crystal_symmetry=crystal_symmetry)
  return reflection_file_utils.process_raw_data(
    obs=data,
    r_free_flags=flags,
    test_flag_value=flag_value,
    phases=phases,
    log=log,
    merge_anomalous=merge_anomalous)

def combine_split_structure(
    pdb_file,
    pdb_id,
    base_dir=None,
    log=None):
  """
  Assembles complete structures from split PDB files (e.g. ribosomes),
  preserving the original file name.  Return value is a list of IDs which
  were added to the current model (or None).
  """
  assert 0
  from mmtbx.command_line import fetch_pdb
  from iotbx import pdb
  if (log is None) : log = sys.stdout
  pdb_in = pdb.input(file_name=pdb_file)
  title = pdb_in.title_section()
  other_ids = None
  for line in title :
    if (line.startswith("SPLIT")):
      fields = line.strip().lower().split()
      other_ids = fields[1:]
      assert (len(other_ids) > 0)
  if (other_ids is not None):
    pdb_files = [pdb_file]
    combined_ids = []
    for other_id in other_ids :
      if (other_id.lower() == pdb_id.lower()):
        continue
      dest_dir_2 = os.path.join(base_dir, other_id)
      if (not os.path.isdir(dest_dir_2)):
        dest_dir_2 = os.getcwd()
      pdb_file_2 = os.path.join(dest_dir_2, "%s.pdb" % other_id)
      if (not os.path.isfile(pdb_file_2)):
        fetch_pdb.run2(args=[other_id])
      if (not os.path.isfile(pdb_file_2)):
        break
      pdb_files.append(pdb_file_2)
      combined_ids.append(other_id)
    if (len(pdb_files) > 1):
      pdb_all = os.path.join(base_dir, "%s_new.pdb" % pdb_id)
      print("Joining multi-part structure: %s %s" % (pdb_id,
        " ".join(other_ids)), file=log)
      easy_run.call("iotbx.pdb.join_fragment_files %s > %s" %
        (" ".join(pdb_files), pdb_all))
      os.remove(pdb_file)
      os.rename(pdb_all, pdb_file)
    return combined_ids
  return None

class filter_pdb_file(object):
  """
  Processing of PDB files to remove common pathologies and enable automatic
  refinement behavior.  In particular, delete unknown atoms and ligands,
  reduce the occupancy of Se atoms from 1 to trigger occupancy refinement,
  and optionally remove zero-occupancy atoms.
  """
  def __init__(self,
                pdb_file,
                output_file=None,
                log=None,
                quiet=False,
                set_se_occ=True,
                remove_atoms_with_zero_occupancy=False):
    assert 0
    import iotbx.pdb
    if (log is None):
      log = null_out()
    pdb_in = iotbx.pdb.input(pdb_file)
    hierarchy = pdb_in.construct_hierarchy()
    if (len(hierarchy.models()) > 1):
      raise Sorry("Multi-MODEL PDB files are not supported.")
    n_unknown = 0
    all_atoms = hierarchy.atoms()
    cache = hierarchy.atom_selection_cache()
    # resname UNK is now okay (with some restrictions)
    known_sel = cache.selection("not (element X or resname UNX or resname UNL)")
    semet_sel = cache.selection("element SE and resname MSE")
    zero_occ_sel = all_atoms.extract_occ() == 0
    self.n_unknown = known_sel.count(False)
    self.n_semet = semet_sel.count(True)
    self.n_zero_occ = zero_occ_sel.count(True)
    keep_sel = known_sel
    modified = False
    if ((self.n_unknown > 0) or
        ((self.n_semet > 0) and (set_se_occ)) or
        (self.n_zero_occ > 0) and (remove_atoms_with_zero_occupancy)):
      modified = True
      if (output_file is None):
        output_file = pdb_file
    if (self.n_unknown > 0) and (not quiet):
      print("Warning: %d unknown atoms or ligands removed:" % \
        self.n_unknown, file=log)
      for i_seq in (~known_sel).iselection():
        print("  %s" % all_atoms[i_seq].id_str(), file=log)
    if (self.n_zero_occ > 0):
      msg = "Warning: %d atoms with zero occupancy present in structure:"
      if (remove_atoms_with_zero_occupancy):
        msg = "Warning: %d atoms with zero occupancy removed:"
        keep_sel &= ~zero_occ_sel
      if (not quiet):
        print(msg % self.n_zero_occ, file=log)
        for i_seq in zero_occ_sel.iselection():
          print("  %s" % all_atoms[i_seq].id_str(), file=log)
    hierarchy_filtered = hierarchy.select(keep_sel)
    if (self.n_semet > 0) and (set_se_occ):
      for atom in hierarchy_filtered.atoms():
        if (atom.element == "SE") and (atom.fetch_labels().resname == "MSE"):
          if (atom.occ == 1.0):
            if (not quiet):
              print("Set occupancy of %s to 0.99" % atom.id_str(), file=log)
            atom.occ = 0.99 # just enough to trigger occupancy refinement
    if (modified):
      f = open(output_file, "w")
      # if the input file is actually from the PDB, we need to preserve the
      # header information for downstream code.
      print("\n".join(pdb_in.title_section()), file=f)
      print("\n".join(pdb_in.remark_section()), file=f)
      print(iotbx.pdb.format_cryst1_record(
        crystal_symmetry=pdb_in.crystal_symmetry()), file=f)
      print(hierarchy_filtered.as_pdb_string(), file=f)
      f.close()

# XXX END_MARKED_FOR_DELETION_OLEG
