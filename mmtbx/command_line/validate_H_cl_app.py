"""validate hydrogens"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.validate_H
import sys, time
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
import iotbx.phil
import mmtbx.model
from mmtbx.monomer_library import pdb_interpretation
from mmtbx.hydrogens.validate_H import validate_H

master_params_str = """\
cif_file = None
  .type = path
  .multiple = True
use_neutron_distances = True
  .type = bool
output_prefix = "output"
  .type = str
input_model_fname = None
  .type = path
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)


class cl_validate_H():
  def __init__(self, cl_args):
    self.help_message = """\
This tool will validate hydrogen atoms in a model file.
Usage:
  %s model.pdb """ % 'validate_H_cl_app.py'
    self.log = sys.stdout
    self.cl_args = cl_args
    self.master_params_str = master_params_str
    self.pdbf_def = 'input_model_fname'
    self.ciff_def = 'cif_file'

  def parse_cl(self):
    self.master_par = master_params()
    self.params = self.master_par.extract()
    if len(self.cl_args) < 1 or "--help" in self.cl_args or "-h" in self.cl_args:
      self.show_help()
      return 1
    self.input_objects = iotbx.phil.process_command_line_with_files(
      args         = self.cl_args[0:],
      master_phil  = self.master_par,
      pdb_file_def = self.pdbf_def,
      cif_file_def = self.ciff_def
      )
    self.work_params = self.input_objects.work.extract()
    # If everything looks fine:
    return 0

  def show_help(self):
    print(self.help_message, file=self.log)
    self.master_par.show(self.log)

  def read_model_file(self):
    self.pdb_inp = iotbx.pdb.input(
      file_name=getattr(self.work_params, self.pdbf_def))
    return 0

  def print_overall_results(self, overall_counts_hd):
    oc = overall_counts_hd
    print('*'*79, file=self.log)
    print('H/D ATOMS IN THE INPUT MODEL:\n', file=self.log)
    print('Total number of hydrogen atoms: ',  oc.count_h, file=self.log)
    print('Total number of deuterium atoms: ', oc.count_d, file=self.log)
    print('Number of H atoms (protein): ', oc.count_h_protein, file=self.log)
    print('Number of D atoms (protein): ', oc.count_d_protein, file=self.log)
    print('Number of H atoms (water): ', oc.count_h_water, file=self.log)
    print('Number of D atoms (water): ', oc.count_d_water, file=self.log)
    print('*'*79, file=self.log)
    print('WATER MOLECULES:\n', file=self.log)
    print('Number of water: ', oc.count_water, file=self.log)
    print('Number of water with 0 H (or D): ', oc.count_water_0h, file=self.log)
    print('Number of water with 1 H (or D): ', oc.count_water_1h, file=self.log)
    print('Number of water with 2 H (or D): ', oc.count_water_2h, file=self.log)
    print('Number of water in alternative conformation: ', \
      oc.count_water_altconf, file=self.log)
    print('Number of water without oxygen atom: ', \
      oc.count_water_no_oxygen, file=self.log)

  def print_renamed(self, renamed):
    print('*'*79, file=self.log)
    print('The following atoms were renamed:', file=self.log)
    for entry in renamed:
      id_str = entry[0]
      oldname = entry[2]
      newname = entry[1]
      print('%s atom %s --> %s' % (id_str, oldname, newname), file=self.log)

  def print_atoms_occ_lt_1(self, hd_atoms_with_occ_0, single_hd_atoms_occ_lt_1):
    if hd_atoms_with_occ_0:
      print('*'*79, file=self.log)
      print('H (or D) atoms with zero occupancy:', file=self.log)
      for item in hd_atoms_with_occ_0:
        print(item[0])
    if single_hd_atoms_occ_lt_1:
      print('H (or D) atoms with occupancy < 1:', file=self.log)
      for item in single_hd_atoms_occ_lt_1:
        print('%s with occupancy %s' %(item[0], item[1]), file=self.log)

  def print_results_hd_sites(
        self, count_exchanged_sites, hd_sites_analysis, overall_counts_hd):
    sites_different_xyz = hd_sites_analysis.sites_different_xyz
    sites_different_b   = hd_sites_analysis.sites_different_b
    sites_sum_occ_not_1 = hd_sites_analysis.sites_sum_occ_not_1
    sites_occ_sum_no_scattering = hd_sites_analysis.sites_occ_sum_no_scattering
    print('*'*79, file=self.log)
    print('H/D EXCHANGED SITES:\n', file=self.log)
    print('Number of H/D exchanged sites: ', count_exchanged_sites, file=self.log)
    print('Number of atoms modelled only as H: ', \
     overall_counts_hd.count_h_protein - count_exchanged_sites, file=self.log)
    print('Number of atoms modelled only as D: ', \
     overall_counts_hd.count_d_protein - count_exchanged_sites, file=self.log)
    if sites_different_xyz:
      print('\nH/D pairs not at identical positions:', file=self.log)
      for item in sites_different_xyz:
        print('%s and  %s at distance %s' % \
          (item[0][5:-1], item[1][5:-1], item[2]), file=self.log)
    if sites_different_b:
      print('\nH/D pairs without identical ADPs:', file=self.log)
      for item in sites_different_b:
        print('%s and %s ' % (item[0][5:-1], item[1][5:-1]), file=self.log)
    if sites_sum_occ_not_1:
      print('\nH/D pairs with occupancy sum != 1:', file=self.log)
      for item in sites_sum_occ_not_1:
        print(' %s  and %s with occupancy sum %s' % \
          (item[0][5:-1], item[1][5:-1], item[2]), file=self.log)
    if sites_occ_sum_no_scattering:
      print('\nRotatable H/D pairs with zero scattering occupancy sum:', file=self.log)
      for item in sites_occ_sum_no_scattering:
        print(' %s with occ %s and  %s with occ %s' %\
          (item[0][5:-1], item[2], item[1][5:-1], item[3]), file=self.log)

  def print_missing_HD_atoms(self, missing_HD_atoms):
    print('*'*79, file=self.log)
    print('MISSING H or D atoms:\n', file=self.log)
    for item in missing_HD_atoms:
      print('%s conformer %s : %s ' % (item[0][8:-1], item[3], ", ".join(item[1])), file=self.log)

  def print_outliers_bonds_angles(self, outliers_bonds, outliers_angles):
    print('*'*79, file=self.log)
    if outliers_bonds:
      print('BOND OUTLIERS:\n', file=self.log)
      for item in outliers_bonds:
        print('Bond %s, observed: %s, delta from target: %s' % \
          (item[1], item[2], item[3]), file=self.log)
    if outliers_angles:
      print('ANGLE OUTLIERS:\n', file=self.log)
      for item in outliers_angles:
        print('Angle %s, observed: %s, delta from target: %s' % \
          (item[1], item[2], item[3]), file=self.log)

  def print_xray_distance_warning(self):
    print('*'*79, file=self.log)
    print('WARNING: Model has a majority of X-H bonds with X-ray bond lengths.\n \
          Input was to use neutron distances. Please check your model carefully.',
          file=self.log)

  def print_results(self, results):
    overall_counts_hd  = results.overall_counts_hd
    count_exchanged_sites = results.count_exchanged_sites
    renamed            = results.renamed
    hd_sites_analysis  = results.hd_sites_analysis
    missing_HD_atoms   = results.missing_HD_atoms
    hd_atoms_with_occ_0 = overall_counts_hd.hd_atoms_with_occ_0
    single_hd_atoms_occ_lt_1 = overall_counts_hd.single_hd_atoms_occ_lt_1
    outliers_bonds     = results.outliers_bonds
    outliers_angles    = results.outliers_angles
    bond_results       = results.bond_results
    if overall_counts_hd:
      self.print_overall_results(overall_counts_hd)
    if renamed:
      self.print_renamed(renamed)
    if hd_atoms_with_occ_0 or single_hd_atoms_occ_lt_1:
      self.print_atoms_occ_lt_1(hd_atoms_with_occ_0, single_hd_atoms_occ_lt_1)
    if count_exchanged_sites is not None:
      self.print_results_hd_sites(
        count_exchanged_sites, hd_sites_analysis, overall_counts_hd)
    if missing_HD_atoms:
      self.print_missing_HD_atoms(missing_HD_atoms)
    if outliers_bonds or outliers_angles:
      self.print_outliers_bonds_angles(outliers_bonds, outliers_angles)
    if bond_results.xray_distances_used:
      self.print_xray_distance_warning()

  def get_pdb_interpretation_params(self):
    pdb_interpretation_phil = iotbx.phil.parse(
      input_string = grand_master_phil_str, process_includes = True)
    pi_params = pdb_interpretation_phil.extract()
    pi_params.pdb_interpretation.use_neutron_distances = \
      self.work_params.use_neutron_distances
    #pi_params.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
    #pi_params.pdb_interpretation.restraints_library.cdl=False
    return pi_params

  def run(self):
    r = self.parse_cl()
    if r != 0:
      return
    r = self.read_model_file()
    if r != 0:
      return

    pi_params = self.get_pdb_interpretation_params()
    model = mmtbx.model.manager(
      model_input       = self.pdb_inp,
      stop_for_unknowns = False,
      restraint_objects = self.input_objects.cif_objects)
    model.process(pdb_interpretation_params=pi_params,
      make_restraints=True)
    print("Model object created from file %s:" % \
      getattr(self.work_params, self.pdbf_def), file=self.log)

    # If needed, this could be wrapped in try...except to catch errors.
    c = validate_H(model = model,
                   use_neutron_distances = pi_params.pdb_interpretation.use_neutron_distances)
    c.validate_inputs()
    c.run()
    results = c.get_results()
    self.print_results(results)

    #results.pdb_hierarchy_curated.write_pdb_file(file_name="%s.pdb" % 'bla2')
if __name__ == "__main__":

  t0 = time.time()
  validate_H_app = cl_validate_H(
    cl_args=sys.argv[1:])
  validate_H_app.run()
  print("Finished. Time: %8.3f"%(time.time()-t0))

