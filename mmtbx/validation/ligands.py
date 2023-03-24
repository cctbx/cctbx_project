
from __future__ import absolute_import, division, print_function
from libtbx.math_utils import percentile_based_spread
from libtbx import slots_getstate_setstate
from libtbx.utils import Sorry
from libtbx import str_utils
import iotbx.pdb
import json
import os.path
import sys
from six.moves import zip

def extract_ligand_residue(hierarchy, ligand_code):
  copies = []
  for chain in hierarchy.only_model().chains():
    main_conf = chain.conformers()[0]
    for residue in main_conf.residues():
      if (residue.resname == ligand_code):
        copies.append(residue)
  return copies

def extract_ligand_atom_group(hierarchy, ligand_code, only_segid=None):
  copies = []
  for chain in hierarchy.only_model().chains():
    for residue_group in chain.residue_groups():
      for atom_group in residue_group.atom_groups():
        if (atom_group.resname == ligand_code):
          if (only_segid is not None):
            have_segid = True
            for atom in atom_group.atoms():
              if (atom.segid.strip() != only_segid.strip()):
                have_segid = False
                break
            if (not have_segid) : continue
          copies.append(atom_group)
  return copies

def compare_ligands(ligand_code,
    hierarchy_1=None,
    hierarchy_2=None,
    pdb_file_1=None,
    pdb_file_2=None,
    max_distance_between_centers_of_mass=8.0,
    exclude_hydrogens=True,
    verbose=False,
    implicit_matching=False,
    out=sys.stdout):
  assert (ligand_code.isalnum())
  assert (([hierarchy_1, hierarchy_2] == [None, None]) or
          ([pdb_file_1, pdb_file_2] == [None, None]))
  if (hierarchy_1 is None):
    assert (hierarchy_2 is None)
    hierarchy_1 = iotbx.pdb.input(pdb_file_1).construct_hierarchy()
    hierarchy_2 = iotbx.pdb.input(pdb_file_2).construct_hierarchy()
  # XXX should this use residues or atom_groups?
  ligands_1 = extract_ligand_residue(hierarchy_1, ligand_code)
  ligands_2 = extract_ligand_residue(hierarchy_2, ligand_code)
  if (len(ligands_1) == 0) or (len(ligands_2) == 0):
    raise Sorry("One or both models missing residue '%s'!" % ligand_code)
  print("%d copies in 1st model, %d copies in 2nd model" % \
    (len(ligands_1), len(ligands_2)), file=out)
  rmsds = []
  pbss = []
  for ligand_1 in ligands_1 :
    rmsds_curr, pbss_curr = compare_ligands_impl(
      ligand=ligand_1,
      reference_ligands=ligands_2,
      verbose=verbose,
      exclude_hydrogens=exclude_hydrogens,
      implicit_matching=False,
      quiet=False,
      out=out)
    rmsds.append(rmsds_curr)
    pbss.append(pbss_curr)
  return rmsds, pbss

def compare_ligands_impl(ligand,
    reference_ligands,
    max_distance_between_centers_of_mass=8.0,
    exclude_hydrogens=True,
    implicit_matching=False,
    verbose=False,
    quiet=False,
    raise_sorry_if_no_matching_atoms=True,
    out=sys.stdout):
  """
  Given a target ligand and a list of reference ligands, return the RMSD(s)
  for any ligand determined to be approximately equivalent.  (Usually there
  will be just one of these, but this allows for alternate conformations.)
  """
  from scitbx.array_family import flex
  from scitbx.matrix import col
  matching = []
  atoms_1 = ligand.atoms()
  sites_1 = atoms_1.extract_xyz()
  xyz_mean_1 = sites_1.mean()
  for ligand_2 in reference_ligands :
    sites_2 = ligand_2.atoms().extract_xyz()
    xyz_mean_2 = sites_2.mean()
    dxyz = abs(col(xyz_mean_1) - col(xyz_mean_2))
    if (dxyz < max_distance_between_centers_of_mass):
      matching.append(ligand_2)
  rmsds = []
  pbss = []
  for ligand_2 in matching :
    atoms_2 = ligand_2.atoms()
    isel_1 = flex.size_t()
    isel_2 = flex.size_t()
    for i_seq, atom_1 in enumerate(ligand.atoms()):
      if (atom_1.element.strip() in ["H","D"]) and (exclude_hydrogens):
        continue
      for j_seq, atom_2 in enumerate(ligand_2.atoms()):
        if (atom_1.name == atom_2.name):
          isel_1.append(i_seq)
          isel_2.append(j_seq)
          break
    if (len(isel_1) == 0):
      if (implicit_matching):
        print("  warning: no atom name matches found - will guess equivalence from sites", file=out)
        # XXX this is embarrassing... needs to be much smarter
        for i_seq, atom_1 in enumerate(ligand.atoms()):
          if (atom_1.element.strip() in ["H","D"]) and (exclude_hydrogens):
            continue
          j_seq_best = None
          name_best = None
          dxyz_best = sys.maxsize
          for j_seq, atom_2 in enumerate(ligand_2.atoms()):
            if (atom_1.element == atom_2.element):
              dxyz = abs(col(atom_1.xyz) - col(atom_2.xyz))
              if (dxyz < dxyz_best):
                j_seq_best = j_seq
                name_best = atom_2.name
                dxyz_best = dxyz
          if (j_seq_best is not None):
            print("    '%s' : '%s' (distance = %.2f)" % (atom_1.name,
              name_best, dxyz_best), file=out)
            isel_1.append(i_seq)
            isel_2.append(j_seq_best)
      if (len(isel_1) == 0):
        if (raise_sorry_if_no_matching_atoms):
          raise Sorry("No matching atoms found!")
        else :
          print("  WARNING: no matching atoms found!", file=out)
          return None
    sites_1 = sites_1.select(isel_1)
    sites_2 = ligand_2.atoms().extract_xyz().select(isel_2)
    rmsd = sites_1.rms_difference(sites_2)
    pbs = percentile_based_spread((sites_2 - sites_1).norms())
    if (not quiet):
      print("  '%s' matches '%s': atoms=%d rmsd=%.3f" % (
        ligand.id_str(), ligand_2.id_str(), sites_1.size(), rmsd), file=out)
    rmsds.append(rmsd)
    pbss.append(pbs)
    if (verbose) and (not quiet):
      atoms = ligand.atoms()
      dxyz = (sites_2 - sites_1).norms()
      for i_seq, j_seq in zip(isel_1, isel_2):
        print("    %s: dxyz=%.2f" % (atoms_1[i_seq].id_str(),
          dxyz[i_seq]), file=out)
  return rmsds, pbss

class ligand_validation(slots_getstate_setstate):
  __slots__ = [
    "cc", "two_fofc_min", "two_fofc_max", "two_fofc_mean", "fofc_min",
    "fofc_max", "fofc_mean", "n_below_two_fofc_cutoff", "n_below_fofc_cutoff",
    "b_iso_mean", "occupancy_mean", "rmsds", "pbss", "id_str",
    "atom_selection", "xyz_center",
  ]
  def __init__(self,
      ligand,
      pdb_hierarchy,
      xray_structure,
      two_fofc_map,
      fofc_map,
      fmodel_map,
      reference_ligands=None,
      two_fofc_map_cutoff=1.5,
      fofc_map_cutoff=-3.0):
    from mmtbx import real_space_correlation
    from cctbx import adptbx
    from scitbx.array_family import flex
    atom_selection = ligand.atoms().extract_i_seq()
    assert (len(atom_selection) == 1) or (not atom_selection.all_eq(0))
    manager = real_space_correlation.selection_map_statistics_manager(
      atom_selection=atom_selection,
      xray_structure=xray_structure,
      fft_m_real=two_fofc_map.all(),
      fft_n_real=two_fofc_map.focus(),
      exclude_hydrogens=True)
    stats_two_fofc = manager.analyze_map(
      map=two_fofc_map,
      model_map=fmodel_map,
      min=1.5)
    stats_fofc = manager.analyze_map(
      map=fofc_map,
      model_map=fmodel_map,
      min=-3.0)
    self.atom_selection = manager.atom_selection # XXX non-hydrogens only!
    sites_cart = xray_structure.sites_cart().select(self.atom_selection)
    self.xyz_center = sites_cart.mean()
    self.id_str = ligand.id_str()
    self.cc = stats_two_fofc.cc
    self.two_fofc_min = stats_two_fofc.min
    self.two_fofc_max = stats_two_fofc.max
    self.two_fofc_mean = stats_two_fofc.mean
    self.fofc_min = stats_fofc.min
    self.fofc_max = stats_fofc.max
    self.fofc_mean = stats_fofc.mean
    self.n_below_two_fofc_cutoff = stats_two_fofc.n_below_min
    self.n_below_fofc_cutoff = stats_fofc.n_below_min
    u_iso = xray_structure.extract_u_iso_or_u_equiv().select(
      self.atom_selection)
    u_iso_mean = flex.mean(u_iso)
    self.b_iso_mean = adptbx.u_as_b(u_iso_mean)
    occ = xray_structure.scatterers().extract_occupancies().select(
      self.atom_selection)
    self.occupancy_mean = flex.mean(occ)
    self.rmsds = self.pbss = None
    if (reference_ligands is not None) and (len(reference_ligands) > 0):
      self.rmsds, self.pbss = compare_ligands_impl(ligand=ligand,
        reference_ligands=reference_ligands,
        max_distance_between_centers_of_mass=8.0,
        raise_sorry_if_no_matching_atoms=False,
        verbose=False,
        quiet=True)

  def show(self, out=sys.stdout):
    box = str_utils.framed_output(out=out, width=80,
      title="Residue: %s" % self.id_str)
    print("""\
Number of non-H atoms = %d   Mean B_iso = %.2f   Mean occ. = %.2f
2mFo-DFc map:  min = %6.2f  max = %6.2f  mean = %6.2f   CC = %5.3f
 mFo-DFc map:  min = %6.2f  max = %6.2f  mean = %6.2f""" % (
      len(self.atom_selection), self.b_iso_mean, self.occupancy_mean,
      self.two_fofc_min, self.two_fofc_max, self.two_fofc_mean, self.cc,
      self.fofc_min, self.fofc_max, self.fofc_mean), file=box)
    if (self.rmsds is not None):
      rmsd_formatted = ", ".join([ "%.3f" % r for r in self.rmsds ])
      if (len(self.rmsds) == 0):
        rmsd_formatted = "[none found]"
      print("RMSD to reference ligand(s): %s" % rmsd_formatted, file=box)
      if (len(self.pbss) > 0):
        pbs_formatted = ", ".join([ "%.3f" % s for s in self.pbss ])
        print("    percentile-based spread: %s" % pbs_formatted, file=box)
    self._show_warnings(box)
    del box

  def _show_warnings(self, out, prefix="    "):
    if (self.cc < 0.8):
      print(prefix+"!!! warning: CC to 2mFo-DFc is poor (< 0.8)", file=out)
    elif (self.cc < 0.9):
      print(prefix+"!!! warning: CC to 2mFo-DFc is sub-optimal (< 0.9)", file=out)
    if (self.n_below_fofc_cutoff > 0):
      print(prefix + \
        "!!! warning: %d atoms have mFo-DFc density < -3sigma" \
        % self.n_below_fofc_cutoff, file=out)

  def show_simple(self, out=sys.stdout, warnings=True):
    print("%-12s %6.2f %6.2f  %5.3f %5.1f %5.1f %5.1f  %5.1f %5.1f %5.1f" % \
      (self.id_str, self.b_iso_mean, self.occupancy_mean,
       self.cc, self.two_fofc_min, self.two_fofc_max, self.two_fofc_mean,
       self.fofc_min, self.fofc_max, self.fofc_mean), file=out)
    if (warnings):
      self._show_warnings(out=out, prefix=" "*8)

# TODO refactor to avoid having 3 maps in memory
def validate_ligands(
    pdb_hierarchy,
    fmodel,
    ligand_code,
    only_segid=None,
    reference_structure=None,
    export_for_web=False,
    output_dir=None):
  assert (0 < len(ligand_code) <= 3)
  two_fofc_map = fmodel.map_coefficients(map_type="2mFo-DFc",
    fill_missing=False,
    merge_anomalous=True).fft_map(
      resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()
  fofc_map = fmodel.map_coefficients(map_type="mFo-DFc",
    fill_missing=False,
    merge_anomalous=True).fft_map(
      resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()
  fmodel_map = fmodel.map_coefficients(map_type="Fmodel",
    fill_missing=False,
    merge_anomalous=True).fft_map(
      resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()
  ligands = extract_ligand_atom_group(pdb_hierarchy, ligand_code,
    only_segid=only_segid)
  if (len(ligands) == 0):
    return None
    #raise Sorry("No residues named '%s' found." % ligand_code)
  if (output_dir is None) : # for web export only
    output_dir = os.getcwd()
  reference_ligands = None
  if (reference_structure is not None):
    import iotbx.pdb
    reference_hierarchy = iotbx.pdb.input(reference_structure).construct_hierarchy()
    # XXX currently using the residue objects for reference ligands, but the
    # atom_group object for the target ligands.  This is because we assume
    # that the target is newly placed and doesn't already have alternate
    # conformations, while the reference may or may not have these.
    reference_ligands = extract_ligand_residue(reference_hierarchy,
      ligand_code)
  unit_cell = fmodel.xray_structure.unit_cell()
  validations = []
  json_block = []
  for i_lig, ligand in enumerate(ligands):
    v = ligand_validation(
      ligand=ligand,
      pdb_hierarchy=pdb_hierarchy,
      xray_structure=fmodel.xray_structure,
      reference_ligands=reference_ligands,
      two_fofc_map=two_fofc_map,
      fofc_map=fofc_map,
      fmodel_map=fmodel_map)
    validations.append(v)
    if (export_for_web):
      import iotbx.map_tools
      two_fofc_file_name = os.path.join(output_dir, "final_%d_2Fo-Fc.dsn6" %
        (i_lig+1))
      fofc_file_name = os.path.join(output_dir, "final_%d_Fo-Fc.dsn6" %
        (i_lig+1))
      sites_cart = ligand.atoms().extract_xyz()
      iotbx.map_tools.write_dsn6_map(
        sites_cart=sites_cart,
        unit_cell=unit_cell,
        map_data=two_fofc_map,
        n_real=two_fofc_map,
        file_name=two_fofc_file_name,
        buffer=3)
      iotbx.map_tools.write_dsn6_map(
        sites_cart=sites_cart,
        unit_cell=unit_cell,
        map_data=fofc_map,
        n_real=fofc_map,
        file_name=fofc_file_name,
        buffer=3)
      ligand_dict = {
        "id_str" : v.id_str,
        "cc" : v.cc,
        "xyz" : v.xyz_center,
        "two_fofc" : os.path.basename(two_fofc_file_name),
        "fofc" : os.path.basename(fofc_file_name),
      }
      json_block.append(ligand_dict)
  if (export_for_web):
    json_file = os.path.join(output_dir, "final_ligands.json")
    open(json_file, "w").write(json.dumps(json_block))
  return validations

def show_validation_results(validations, out, verbose=True):
  if (not verbose):
    box = str_utils.framed_output(out=out, title="Ligand summary", width=80)
    print("%30s |-----2mFo-DFc-----|    |---mFo-Dfc---|" % "", file=box)
    print("%-12s %6s %6s  %5s %5s %5s %5s  %5s %5s %5s" % (
      "ID", "B_iso", "Occ", "CC", "min", "max", "mean", "min", "max", "mean"), file=box)
    box.add_separator()
    out = box
  for v in validations :
    if (verbose):
      v.show(out=out)
      print("", file=out)
    else :
      v.show_simple(out=out)
