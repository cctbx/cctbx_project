from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.validate_H
import sys, time
#from libtbx.utils import Sorry
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
import iotbx.phil
import mmtbx.model
from mmtbx.monomer_library import pdb_interpretation
#from mmtbx import monomer_library
from mmtbx.hydrogens.validate_H import validate_H
#import mmtbx.monomer_library.server

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
    print >> self.log, self.help_message
    self.master_par.show(self.log)

  def read_model_file(self):
    self.pdb_inp = iotbx.pdb.input(
      file_name=getattr(self.work_params, self.pdbf_def))
    return 0

  def print_overall_results(self, overall_counts_hd):
    oc = overall_counts_hd
    print >> self.log, '*'*79
    print >> self.log, 'H/D ATOMS IN THE INPUT MODEL:\n'
    print >> self.log, 'Total number of hydrogen atoms: ',  oc.count_h
    print >> self.log, 'Total number of deuterium atoms: ', oc.count_d
    print >> self.log, 'Number of H atoms (protein): ', oc.count_h_protein
    print >> self.log, 'Number of D atoms (protein): ', oc.count_d_protein
    print >> self.log, 'Number of H atoms (water): ', oc.count_h_water
    print >> self.log, 'Number of D atoms (water): ', oc.count_d_water
    print >> self.log, '*'*79
    print >> self.log, 'WATER MOLECULES:\n'
    print >> self.log, 'Number of water: ', oc.count_water
    print >> self.log, 'Number of water with 0 H (or D): ', oc.count_water_0h
    print >> self.log, 'Number of water with 1 H (or D): ', oc.count_water_1h
    print >> self.log, 'Number of water with 2 H (or D): ', oc.count_water_2h
    print >> self.log, 'Number of water in alternative conformation: ', \
      oc.count_water_altconf
    print >> self.log, 'Number of water without oxygen atom: ', \
      oc.count_water_no_oxygen

  def print_renamed(self, renamed):
    print >> self.log, '*'*79
    print >> self.log, 'The following atoms were renamed:'
    for entry in renamed:
      atom = entry['atom']
      oldname = entry['oldname']
      newname = atom.name
      print >> self.log, 'Residue %s atom %s --> %s' % \
        (atom.id_str()[10:-1], oldname, newname)

  def print_atoms_occ_0(self, hd_atoms_with_occ_0):
    print >> self.log, '*'*79
    print >> self.log, 'The following H (or D) atoms have zero occupancy:'
    for item in hd_atoms_with_occ_0:
      print item[0][5:-1]

  def print_results_hd_sites(
        self, hd_exchanged_sites, hd_sites_analysis, overall_counts_hd):
    count_exchanged_sites = len(hd_exchanged_sites.keys())
    sites_different_xyz = hd_sites_analysis.sites_different_xyz
    sites_different_b   = hd_sites_analysis.sites_different_b
    sites_sum_occ_not_1 = hd_sites_analysis.sites_sum_occ_not_1
    sites_occ_sum_no_scattering = hd_sites_analysis.sites_occ_sum_no_scattering
    print >> self.log, '*'*79
    print >> self.log, 'H/D EXCHANGED SITES:\n'
    print >> self.log, 'Number of H/D exchanged sites: ', count_exchanged_sites
    print >> self.log, 'Number of atoms modelled only as H: ', \
     overall_counts_hd.count_h_protein - count_exchanged_sites
    print >> self.log, 'Number of atoms modelled only as D: ', \
     overall_counts_hd.count_d_protein - count_exchanged_sites
    if sites_different_xyz:
      print >> self.log, '\nH/D pairs not at identical positions:'
      for item in sites_different_xyz:
        print >> self.log, '%s and  %s at distance %s' % \
          (item[0][5:-1], item[1][5:-1], item[2])
    if sites_different_b:
      print >> self.log, '\nH/D pairs without identical ADPs:'
      for item in sites_different_b:
        print >> self.log, '%s and %s ' % (item[0][5:-1], item[1][5:-1])
    if sites_sum_occ_not_1:
      print >> self.log, '\nH/D pairs with occupancy sum != 1:'
      for item in sites_sum_occ_not_1:
        print >> self.log, ' %s  and %s with occupancy sum %s' % \
          (item[0][5:-1], item[1][5:-1], item[2])
    if sites_occ_sum_no_scattering:
      print >> self.log, \
        '\nRotatable H/D pairs with zero scattering occupancy sum:'
      for item in sites_occ_sum_no_scattering:
        print >> self.log, ' %s with occ %s and  %s with occ %s' %\
          (item[0][5:-1], item[2], item[1][5:-1], item[3])

  def print_missing_HD_atoms(self, missing_HD_atoms):
    print >> self.log, '*'*79
    print >> self.log, 'MISSING H or D atoms:\n'
    for item in missing_HD_atoms:
      print >> self.log, '%s : %s ' % (item[0][8:-1], ", ".join(item[1]))

  def print_results(self, results):
    overall_counts_hd  = results.overall_counts_hd
    hd_exchanged_sites = results.hd_exchanged_sites
    renamed            = results.renamed
    hd_sites_analysis  = results.hd_sites_analysis
    missing_HD_atoms   = results.missing_HD_atoms
    hd_atoms_with_occ_0 = overall_counts_hd.hd_atoms_with_occ_0
    if overall_counts_hd:
      self.print_overall_results(overall_counts_hd)
    if renamed:
      self.print_renamed(renamed)
    if hd_atoms_with_occ_0:
      self.print_atoms_occ_0(hd_atoms_with_occ_0)
    if len(hd_exchanged_sites.keys()) != 0:
      self.print_results_hd_sites(
        hd_exchanged_sites, hd_sites_analysis, overall_counts_hd)
    if missing_HD_atoms:
      self.print_missing_HD_atoms(missing_HD_atoms)

  def get_pdb_interpretation_params(self):
    pdb_interpretation_phil = iotbx.phil.parse(
      input_string = grand_master_phil_str, process_includes = True)
    pi_params = pdb_interpretation_phil.extract()
    pi_params.pdb_interpretation.use_neutron_distances = \
      self.params.use_neutron_distances
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
        model_input               = self.pdb_inp,
        process_input             = True,
        pdb_interpretation_params = pi_params,
        restraint_objects         = self.input_objects.cif_objects)

    print >> self.log, "Model object created from file %s:" % \
      getattr(self.work_params, self.pdbf_def)

    # If needed, this could be wrapped in try...except to catch errors.
    c = validate_H(
      model  = model,
      params = self.params,
      log    = self.log)
    c.validate_inputs()
    c.run()
    results = c.get_results()
    self.print_results(results)

    #results.pdb_hierarchy_curated.write_pdb_file(file_name="%s.pdb" % 'bla2')

t0 = time.time()
validate_H_app = cl_validate_H(
    cl_args=sys.argv[1:])
validate_H_app.run()
print "Finished. Time: %8.3f"%(time.time()-t0)
