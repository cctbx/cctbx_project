
import iotbx.pdb.secondary_structure
import libtbx.phil
from libtbx import smart_open
from libtbx.utils import Sorry
from math import sqrt
import sys, os

ss_restraint_params_str = """
  restrain_base_pairs = False
    .type = bool
  remove_outliers = True
    .type = bool
  restrain_initial_values = False
    .type = bool
  slack = 0.1
    .type = float
  sigma = 0.05
    .type = float
  n_o_distance_ideal = 3.0
    .type = float
  n_o_outlier_cutoff = 3.5
    .type = float
  h_o_distance_ideal = 2.0
    .type = float
  h_o_outlier_cutoff = 2.5
    .type = float
"""

sec_str_master_phil_str = """
%s
%s
""" % (iotbx.pdb.secondary_structure.ss_input_params_str,
       ss_restraint_params_str)

sec_str_master_phil = libtbx.phil.parse(sec_str_master_phil_str)

def analyze_h_bonds (pdb_hierarchy,
                     bonded_selections,
                     params,
                     bond_i_seqs=None,
                     selection_cache=None,
                     log=sys.stdout) :
  if bond_i_seqs is not None :
    assert len(bond_i_seqs) == len(bonded_selections)
  sites = pdb_hierarchy.atoms().extract_xyz()
  if selection_cache is None and bond_i_seqs is None :
    selection_cache = pdb_hierarchy.atom_selection_cache()
  isel = selection_cache.iselection
  distances = []
  approved_bonds = []
  approved_selections = []
  if params.use_hydrogens :
    cutoff = params.h_o_outlier_cutoff
  else :
    cutoff = params.n_o_outlier_cutoff
  for k, (donor_sele, acceptor_sele) in enumerate(bonded_selections) :
    if bond_i_seqs is not None :
      (donor_i_seqs, acceptor_i_seqs) = bond_i_seqs[k]
    else :
      donor_i_seqs = isel(donor_sele)
      acceptor_i_seqs = isel(acceptor_sele)
    n_donor_sel = donor_i_seqs.size()
    n_acceptor_sel = acceptor_i_seqs.size()
    if n_donor_sel == 0 or n_acceptor_sel == 0 :
      print >> log, "analyze_h_bonds: one or more atoms missing"
      print >> log, "  %s (%d atoms)" % (donor_sele, n_donor_sel)
      print >> log, "  %s (%d atoms)" % (acceptor_sele, n_acceptor_sel)
      continue
    elif n_donor_sel > 1 or n_acceptor_sel > 1 :
      print >> log, "analyze_h_bonds: multiple atoms matching a selection"
      print >> log, "  %s (%d atoms)" % (donor_sele, n_donor_sel)
      print >> log, "  %s (%d atoms)" % (acceptor_sele, n_acceptor_sel)
      continue
    (x1, y1, z1) = sites[donor_i_seqs[0]]
    (x2, y2, z2) = sites[acceptor_i_seqs[0]]
    dist = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
    if dist > cutoff and params.remove_outliers :
      print >> log, "analyze_h_bonds: distance exceeds cutoff (%.2fA, %.2fA)"%\
        (dist, cutoff)
      print >> log, "  %s" % donor_sele
      print >> log, "  %s" % acceptor_sele
      continue
    distances.append(dist)
    approved_selections.append((donor_sele, acceptor_sele))
    approved_bonds.append((donor_i_seqs[0], acceptor_i_seqs[0]))
  return (approved_bonds, approved_selections, distances)

# XXX: in practice, this will be much too slow for normal use, but it may
# be useful for testing purposes.
def generate_custom_bond_restraints (params, sec_str, pdb_hierarchy,
    log=sys.stdout) :
  bonded_selections = sec_str.as_bond_selections(params)
  (approved_bonds, approved_selections, distances) = analyze_h_bonds(
    pdb_hierarchy=pdb_hierarchy,
    bonded_selections=bonded_selections,
    params=params,
    log=log)
  bond_params = []
  for (selections, dist) in zip(approved_selections, distances) :
    if params.restrain_initial_values :
      distance = dist
    elif params.use_hydrogens :
      distance = params.h_o_distance_ideal
    else :
      distance = params.n_o_distance_ideal
    bond_params.append(
"""bond {
  atom_selection_1 = "%s"
  atom_selection_2 = "%s"
  distance_ideal = %.3f
  sigma = %.4f
  slack = %.4f
}""" % (selections[0], selections[1], distance, params.sigma, params.slack))
  return "\n".join(bond_params)

def process_files (params, records=None) :
  from iotbx import file_reader
  assert len(params.file_name) > 0 or records is not None
  if records is None :
    records = []
    for file_name in params.file_name :
      lines = smart_open.for_reading(file_name).readlines()
      for line in lines :
        if line.startswith("HELIX") or line.startswith("SHEET") :
          records.append(line)
  secondary_structure = iotbx.pdb.secondary_structure.annotation(
    records=records)
  return secondary_structure

def get_pdb_hierarchy (file_names) :
  import iotbx.pdb
  from scitbx.array_family import flex
  pdb_combined = iotbx.pdb.combine_unique_pdb_files(file_names=file_names)
  pdb_structure = iotbx.pdb.input(source_info=None,
    lines=flex.std_string(pdb_combined.raw_records))
  return pdb_structure.construct_hierarchy()

def run (args, out=sys.stdout, log=sys.stderr) :
  import iotbx.pdb
  master_phil = libtbx.phil.parse("""
  show_restraints = True
    .type = bool
  echo_pdb_records = False
    .type = bool
  echo_pymol_cmds = False
    .type = bool
  %s
""" % sec_str_master_phil_str)
  user_phil = []
  for arg in args :
    if os.path.isfile(arg) :
      file_name = os.path.abspath(arg)
      assert iotbx.pdb.is_pdb_file(file_name)
      user_phil.append(libtbx.phil.parse("file_name=\"%s\"" % file_name))
    else :
      if arg.startswith("--") :
        arg = arg[2:] + "=True"
      try :
        cmdline_phil = libtbx.phil.parse(arg)
      except RuntimeError, e :
        print >> log, str(e)
      else :
        user_phil.append(cmdline_phil)
  working_phil = master_phil.fetch(sources=user_phil)
  params = working_phil.extract()
  if len(params.file_name) == 0 :
    raise Usage("Please supply at least one PDB file.")
  secondary_structure = process_files(params)
  if params.echo_pdb_records :
    print >> out, secondary_structure.as_pdb_str()
  elif params.echo_pymol_cmds :
    print >> out, secondary_structure.as_pymol_dashes()
  elif params.show_restraints :
    print >> out, generate_custom_bond_restraints(params, secondary_structure,
      get_pdb_hierarchy(params.file_name), log=sys.stderr)

if __name__ == "__main__" :
  run(sys.argv[1:])
