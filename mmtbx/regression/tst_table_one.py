from __future__ import absolute_import, division, print_function

import os

from mmtbx.regression import model_1yjp
from mmtbx.command_line import table_one
import iotbx.pdb
from cctbx import sgtbx
from scitbx.array_family import flex
from libtbx.utils import null_out
import random
import libtbx.load_env

def exercise():
  flex.set_random_seed(123456)
  random.seed(123456)
  base = "tst_table_one"
  pdb_in = iotbx.pdb.input(source_info=None, lines=model_1yjp)
  xrs = pdb_in.xray_structure_simple()
  xrs.set_inelastic_form_factors(
    photon=1.54,
    table="sasaki")
  fc = abs(xrs.structure_factors(d_min=1.5).f_calc()).average_bijvoet_mates()
  fc.set_observation_type_xray_amplitude()
  flags = fc.generate_r_free_flags()
  mtz = fc.as_mtz_dataset(column_root_label="F",
    wavelength=1.54)
  mtz.add_miller_array(flags, column_root_label="FreeR_flag")
  mtz.mtz_object().write(base + ".mtz")
  xrs_p1 = xrs.expand_to_p1()
  xrs_p1.shake_sites_in_place(rms_difference=0.1)
  fc_p1 = xrs_p1.structure_factors(d_min=1.4).f_calc()
  fc_p1_extra = fc_p1.randomize_amplitude_and_phase(amplitude_error=1.0,
    phase_error_deg=0,
    random_seed=123456)
  fc_p1 = abs(fc_p1.concatenate(other=fc_p1_extra)).sort(
    by_value="packed_indices")
  fc_p1.set_observation_type_xray_amplitude()
  sg_p2 = sgtbx.space_group_info("P2")
  ic = fc_p1.f_as_f_sq().customized_copy(
    space_group_info=sg_p2,
    sigmas=flex.double(fc_p1.size(), 10.0))
  ic.export_as_scalepack_unmerged(file_name=base + ".sca")
  with open(base + ".pdb", "w") as f:
    f.write(model_1yjp)
  args = [
    base + ".mtz",
    base + ".pdb",
    "unmerged_data=%s.sca" % base,
    "prefix=tst_table_one_1",
  ]
# BROKEN
#  table_one.run(args=args, out=null_out(),
#    use_current_directory_if_not_specified=True)
#  # now with unmerged data in SHELX format
#  f = open(base + ".hkl", "w")
#  ic.export_as_shelx_hklf(file_object=f)
#  f.close()
#  args = [
#    base + ".mtz",
#    base + ".pdb",
#    "unmerged_data=%s.hkl=hklf4" % base,
#    "prefix=tst_table_one_2",
#  ]
#  table_one.run(args=args, out=null_out(),
#    use_current_directory_if_not_specified=True)
#  # now with phil file
#  f = open("tst_table_one_3.eff", "w")
#  f.write("""\
#table_one {
#  structure {
#    name = %(base)s
#    pdb_file = %(base)s.pdb
#    mtz_file = %(base)s.mtz
#    unmerged_data = %(base)s.hkl=hklf4
#  }
#  output {
#    directory = os.getcwd()
#    base_name = %(base)s_3
#  }
#}""" % {"base" : base })
#  args = [ "tst_table_one_3.eff" ]
#  table_one.run(args=args, out=null_out(),
#    use_current_directory_if_not_specified=True)

def exercise_counts():
  pdb_file = libtbx.env.under_dist('mmtbx', 'regression/pdbs/two_chains_ligand_water.pdb')
  base = "tst_table_one_counts"
  mtz_file = base + ".mtz"
  pdb_in = iotbx.pdb.input(pdb_file)
  xrs = pdb_in.xray_structure_simple()
  xrs.set_inelastic_form_factors(photon=1.54, table="sasaki")
  fc = abs(xrs.structure_factors(d_min=4.0).f_calc()).average_bijvoet_mates()
  fc.set_observation_type_xray_amplitude()
  fcs = fc.customized_copy(sigmas=flex.double(fc.size(), 10.0))
  flags = fcs.generate_r_free_flags()
  mtz = fcs.as_mtz_dataset(column_root_label="Fobs", wavelength=1.54)
  mtz.add_miller_array(flags, column_root_label="R-free-flags")
  mtz.mtz_object().write(mtz_file)
  params = """
table_one {
  structure {
    name = "abc"
    pdb_file = "%s"
    mtz_file = "./%s"
    data_labels = "Fobs,SIGFobs"
    r_free_flags_label = "R-free-flags"
    wavelength = None
    cif_file = None
    cif_directory = None
    data_type = *xray neutron
    unmerged_data = None
    unmerged_labels = None
    use_internal_variance = False
    count_anomalous_pairs_separately = False
    enable_twinning = False
    twin_law = Auto
  }
  processing {
    re_compute_r_factors = True
    n_bins = 10
    ligand_selection = None
  }
  multiprocessing {
    nproc = 1
    technology = *multiprocessing sge lsf pbs condor slurm
    qsub_command = None
  }
  output {
    directory = "%s"
    job_title = None
    show_missing_fields = True
    format = *txt *csv *rtf
    base_name = "%s"
    verbose = "True"
    text_field_separation = 2
  }
}
""" % (pdb_file, mtz_file, os.getcwd(), base)
  eff_file = base + ".eff"
  with open(eff_file, "w") as f:
    f.write(params)
  args = [ eff_file ]
  table_one.run(args=args, out=null_out(),
    use_current_directory_if_not_specified=True)

  count_stats = """\
  Number of non-hydrogen atoms                             65
                macromolecules                             26
                       ligands                             35
                       solvent                              4
              Protein residues                              3\
"""
  adp_stats = """\
              Average B-factor                          16.01
                macromolecules                           8.07
                       ligands                          23.16
                       solvent                           5.19\
"""
  with open(base + ".txt", "r") as f:
    lines = f.read()
  assert count_stats in lines
  assert adp_stats in lines

if (__name__ == "__main__"):
  exercise()
  exercise_counts()
  print("OK")
