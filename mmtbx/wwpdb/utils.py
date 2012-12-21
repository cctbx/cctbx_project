
"""
Utilities for re-refining (or otherwise working with) structures downloaded
directly from the PDB.
"""

from __future__ import division
from libtbx import easy_run
from libtbx.utils import null_out, Sorry
from libtbx import group_args
import os
import sys

def get_program (pdb_file) :
  from iotbx.pdb import remark_3_interpretation
  lines = open(pdb_file).readlines()
  program = program_full = None
  for line in lines :
    if (line.startswith("REMARK   3")) and ("PROGRAM" in line) :
      program = remark_3_interpretation.get_program(line)
      if (program is not None) :
        program_full = line.split(":")[1].strip()
        break
  if (program == "PHENIX") :
    program = "PHENIX.REFINE"
  return program, program_full

def fetch_pdb_data (
    pdb_id,
    pdb_dir=None,
    sf_dir=None,
    log=None,
    verbose=False) :
  """
  Copy data from local repository if defined and available, or download it
  from the PDB, and run cif_as_mtz.
  """
  from mmtbx.command_line import fetch_pdb
  from mmtbx.command_line import cif_as_mtz
  if (log is None) :
    if (verbose) : log = sys.stdout
    else : log = null_out()
  pdb_file = "%s.pdb" % pdb_id
  cif_file = "%s-sf.cif" % pdb_id
  mtz_file = "%s.mtz" % pdb_id
  fetch_pdb.run2(args=[pdb_id], log=log)
  assert (os.path.isfile(pdb_file))
  fetch_pdb.run2(args=["-x", pdb_id], log=log)
  if (not os.path.isfile("%s-sf.cif" % pdb_id)) :
    raise Sorry("Structure factors are not available for %s." % pdb_id)
  cif_as_mtz.run(args=[
      cif_file,
      "--symmetry=%s" % pdb_file,
      "--merge",
      "--output_file_name=%s" % mtz_file])
  if (not os.path.isfile(mtz_file)) :
    raise RuntimeError("Missing %s!\ncif_as_mtz stderr:\n%s" %
      (mtz_file, "\n".join(import_out.stderr_lines)))
  return os.path.abspath(pdb_file), os.path.abspath(mtz_file)

def find_data_array (mtz_file, log=None) :
  """
  Guess an appropriate data array to use for refinement, plus optional
  Hendrickson-Lattman coefficients if present.
  """
  if (log is None) : log = sys.stdout
  from iotbx.file_reader import any_file
  hl = data = None
  hkl_in = any_file(mtz_file, force_type="hkl")
  hkl_server = hkl_in.file_server
  data_arrays = hkl_server.get_xray_data(
    file_name               = None,
    labels                  = None,
    ignore_all_zeros        = False,
    parameter_scope         = "",
    return_all_valid_arrays = True,
    minimum_score           = 4)
  # XXX always use anomalous data if available!
  if (len(data_arrays) > 0) :
    data_labels = [ array.info().label_string() for array in data_arrays ]
    for array_label in ["FOBS(+),SIGFOBS(+),FOBS(-),SIGFOBS(-)",
                        "F(+),SIGF(+),F(-),SIGF(-)",
                        "I(+),SIGI(+),I(-),SIGI(-)",
                        "FOBS,SIGFOBS",
                        "IOBS,SIGIOBS"] :
      if (array_label in data_labels) :
        data = array_label
        break
    else :
      data = data_labels[0]
  hl_arrays = hkl_server.get_experimental_phases(
    file_name               = None,
    labels                  = None,
    ignore_all_zeros        = True,
    parameter_scope         = "",
    return_all_valid_arrays = True,
    minimum_score           = 1)
  if (len(hl_arrays) > 0) :
    hl_labels = [ array.info().label_string() for array in data_arrays ]
    if ("HLA,HLB,HLC,HLD" in hl_labels) :
      hl = "HLA,HLB,HLC,HLD"
  # very important - filter experimental data to include only those
  # reflections with R-free flags (otherwise phenix.refine will crash)
  filtered_data = filter_raw_data(
    file_name=mtz_file,
    data_label=data,
    hl_label=hl,
    log=log)
  os.remove(mtz_file)
  os.rename(filtered_data.mtz_file, mtz_file)
  if (filtered_data.new_flags) :
    new_flags = os.path.join(os.path.dirname(mtz_file), "NEW_FLAGS")
    open(new_flags, "w").write(mtz_file)
  return filtered_data

def filter_raw_data (
    file_name,
    data_label,
    hl_label=None,
    log=None) :
  """
  Delete any data which do not have R-free flags (unless no flags are present,
  in which case they will be created).
  """
  from iotbx import reflection_file_utils
  from iotbx import file_reader
  if (log is None) : log = sys.stdout
  hkl_in = file_reader.any_file(file_name, force_type="hkl")
  miller_arrays = hkl_in.file_server.miller_arrays
  data = r_free = hl_coeffs = None
  for array in miller_arrays :
    labels = array.info().label_string()
    if (labels == data_label) :
      data = array
    elif (labels == "R-free-flags") :
      r_free = array
    elif (labels == hl_label) :
      hl_coeffs = array
  assert (data is not None)
  processed = reflection_file_utils.process_raw_data(
    obs=data,
    r_free_flags=r_free,
    test_flag_value=0, # XXX this is what cif_as_mtz does
    phases=hl_coeffs,
    log=log,
    merge_reconstructed_amplitudes=False)
  new_mtz_file = os.path.splitext(file_name)[0] + "_new.mtz"
  processed.write_mtz_file(new_mtz_file, single_dataset=True)
  hkl_new = file_reader.any_file(new_mtz_file, force_type="hkl")
  new_arrays = hkl_new.file_server.miller_arrays
  assert (len(new_arrays) >= 2)
  data_label = new_arrays[0].info().label_string()
  r_free_label = new_arrays[1].info().label_string()
  if (len(new_arrays) == 3) :
    hl_label = new_arrays[2].info().label_string()
  else :
    hl_label = None
  return group_args(
    mtz_file=new_mtz_file,
    data_label=data_label,
    r_free_label=r_free_label,
    hl_label=hl_label,
    new_flags=processed.flags_are_new(),
    anomalous_flag=data.anomalous_flag())

def combine_split_structure (
    pdb_file,
    pdb_id,
    base_dir=None,
    log=None) :
  """
  Assembles complete structures from split PDB files (e.g. ribosomes),
  preserving the original file name.  Return value is a list of IDs which
  were added to the current model (or None).
  """
  from mmtbx.command_line import fetch_pdb
  from iotbx import pdb
  if (log is None) : log = sys.stdout
  pdb_in = pdb.input(file_name=pdb_file)
  title = pdb_in.title_section()
  other_ids = None
  for line in title :
    if (line.startswith("SPLIT")) :
      fields = line.strip().lower().split()
      other_ids = fields[1:]
      assert (len(other_ids) > 0)
  if (other_ids is not None) :
    pdb_files = [pdb_file]
    combined_ids = []
    for other_id in other_ids :
      if (other_id.lower() == pdb_id.lower()) :
        continue
      dest_dir_2 = os.path.join(base_dir, other_id)
      if (not os.path.isdir(dest_dir_2)) :
        dest_dir_2 = os.getcwd()
      pdb_file_2 = os.path.join(dest_dir_2, "%s.pdb" % other_id)
      if (not os.path.isfile(pdb_file_2)) :
        fetch_pdb.run2(args=[other_id])
      if (not os.path.isfile(pdb_file_2)) :
        break
      pdb_files.append(pdb_file_2)
      combined_ids.append(other_id)
    if (len(pdb_files) > 1) :
      pdb_all = os.path.join(base_dir, "%s_new.pdb" % pdb_id)
      print >> log, "Joining multi-part structure: %s %s" % (pdb_id,
        " ".join(other_ids))
      easy_run.call("iotbx.pdb.join_fragment_files %s > %s" %
        (" ".join(pdb_files), pdb_all))
      os.remove(pdb_file)
      os.rename(pdb_all, pdb_file)
    return combined_ids
  return None

class filter_pdb_file (object) :
  """
  Get rid of those pesky unknown atoms, and make other modifications that we
  deem prudent, such as reducing the occupancy of Se atoms from 1.
  """
  def __init__ (self,
                pdb_file,
                output_file=None,
                log=None,
                set_se_occ=True) :
    from iotbx.file_reader import any_file
    import iotbx.pdb
    if (log is None) :
      log = null_out()
    pdb_in = any_file(pdb_file, force_type="pdb")
    pdb_in.assert_file_type("pdb")
    hierarchy = pdb_in.file_object.construct_hierarchy()
    if (len(hierarchy.models()) > 1) :
      raise Sorry("Multi-MODEL PDB files are not supported.")
    n_unknown = 0
    cache = hierarchy.atom_selection_cache()
    known_sel = cache.selection("not (element X or resname UNX or resname UNL)")
    semet_sel = cache.selection("element SE and resname MSE")
    self.n_unknown = known_sel.count(False)
    self.n_semet = semet_sel.count(True)
    if (self.n_unknown > 0) or (self.n_semet > 0) :
      hierarchy_filtered = hierarchy.select(known_sel)
      if (output_file is None) :
        output_file = pdb_file
      if (self.n_semet > 0) and (set_se_occ) :
        for atom in hierarchy_filtered.atoms() :
          if (atom.element == "SE") and (atom.fetch_labels().resname == "MSE") :
            if (atom.occ == 1.0) :
              atom.occ = 0.99 # just enough to trigger occupancy refinement
      f = open(output_file, "w")
      # if the input file is actually from the PDB, we need to preserve the
      # header information for downstream code.
      print >> f, "\n".join(pdb_in.file_object.title_section())
      print >> f, "\n".join(pdb_in.file_object.remark_section())
      print >> f, iotbx.pdb.format_cryst1_record(
        crystal_symmetry=pdb_in.file_object.crystal_symmetry())
      print >> f, hierarchy_filtered.as_pdb_string()
      f.close()
      pdb_base = os.path.basename(pdb_file)
      print >> log, "WARNING: removed %d unknown atoms from %s" % (n_unknown,
        pdb_base)
