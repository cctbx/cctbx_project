from __future__ import division, print_function
import sys, os, time

from libtbx.test_utils import approx_equal
import iotbx.phil

import libtbx.load_env
data_dir = libtbx.env.under_dist(
  module_name="mmtbx",
  path=os.path.join("regression","pdbs"),
  test=os.path.isdir)

pae_data_dir = libtbx.env.under_dist(
  module_name="mmtbx",
  path="regression",
  test=os.path.isdir)

from iotbx.data_manager import DataManager
from mmtbx.process_predicted_model import split_model_into_compact_units, \
   get_cutoff_b_value, \
   get_b_values_from_lddt, get_rmsd_from_lddt, \
   process_predicted_model, master_phil_str, get_lddt_from_b

from mmtbx.domains_from_pae import parse_pae_file
master_phil = iotbx.phil.parse(master_phil_str)
params = master_phil.extract()


model_file=os.path.join(data_dir,'fibronectin_af_ca_1358_1537.pdb')
pae_model_file=os.path.join(data_dir,'pae_model.pdb')
pae_file=os.path.join(pae_data_dir,'pae.json')

def tst_01(log = sys.stdout):


  # Check calculations of conversion between rmsd, lddt , and B values
  print("\nChecking conversions between rmsd, lddt and B-values", file = log)
  for maximum_rmsd, minimum_lddt, target_b in [
       (1.5, None, 59.2175263686),
       (None,0.7,59.2175263686),
       (1.5,0.7,59.2175263686),
       (1.0, None, 26.3189006083),
       (None,0.5,293.306328196),]:
    print()
    cutoff_b = get_cutoff_b_value(
      maximum_rmsd,
      minimum_lddt,
      log = log)
    print("maximum_rmsd: %s min lddt %s Cutoff B:  %.2f" %(
     maximum_rmsd, minimum_lddt,
     cutoff_b), file = log)
    assert approx_equal(cutoff_b, target_b)


  # Read in alphafold model and get LDDT from B-value field
  print("\nReading in alphafold model with lddt values in B-value field",
    file = log)

  dm = DataManager()
  dm.set_overwrite(True)
  m = dm.get_model(model_file)
  pae_m = dm.get_model(pae_model_file)
  pae_matrix = parse_pae_file(pae_file)

  lddt_values = m.get_hierarchy().atoms().extract_b().deep_copy()
  print("\nLDDT mean:",lddt_values.min_max_mean().mean)
  assert approx_equal(lddt_values.min_max_mean().mean, 82.5931111111)

  # Multiply lddt_values by 0.01 (fractional)
  fractional_lddt = lddt_values * 0.01

  #  Convert lddt to b
  b_values = get_b_values_from_lddt(lddt_values)
  print("B-value mean:",b_values.min_max_mean().mean)
  assert approx_equal(b_values.min_max_mean().mean, 24.7254093338)

  # Convert b to lddt
  lddt = get_lddt_from_b(b_values)
  assert approx_equal(lddt,fractional_lddt)
  lddt = get_lddt_from_b(b_values, input_lddt_is_fractional=False)
  assert approx_equal(lddt,lddt_values)

  # Convert  lddt to rmsd
  rmsd_values = get_rmsd_from_lddt(lddt_values)
  print("RMSD mean:",rmsd_values.min_max_mean().mean)
  assert approx_equal(rmsd_values.min_max_mean().mean, 0.93559254135)

  # use process_predicted_model to convert lddt or rmsd to B return with
  #  mark_atoms_to_ignore_with_occ_zero

  print("\nConverting lddt to B values and using mark_atoms_to_ignore_with_occ_zero", file = log)
  params.process_predicted_model.maximum_fraction_close = 0.5
  params.process_predicted_model.b_value_field_is = 'lddt'
  params.process_predicted_model.remove_low_confidence_residues = True
  params.process_predicted_model.maximum_rmsd = 1.5
  params.process_predicted_model.split_model_by_compact_regions = True
  params.process_predicted_model.maximum_domains = 3

  model_info = process_predicted_model(m, params, mark_atoms_to_keep_with_occ_one= True)
  models = model_info.model_list
  for mm,vrms,target_vrms,n1,n2 in zip(models,model_info.vrms_list,[1.1506528458663525,1.1506528458663525],[84,88],[88,84]):
    model_occ_values = mm.get_hierarchy().atoms().extract_occ()
    assert model_occ_values.count(1) == n1
    assert model_occ_values.count(0) == n2
    assert approx_equal(vrms, target_vrms, eps=0.01)

  # use process_predicted_model to convert lddt or rmsd to B

  print("\nConverting lddt to B values", file = log)
  params.process_predicted_model.maximum_fraction_close = 0.5
  params.process_predicted_model.b_value_field_is = 'lddt'
  params.process_predicted_model.remove_low_confidence_residues = False
  params.process_predicted_model.split_model_by_compact_regions = False
  params.process_predicted_model.input_lddt_is_fractional = None

  model_info = process_predicted_model(m, params)
  model = model_info.model
  model_b_values = model.get_hierarchy().atoms().extract_b()
  assert approx_equal(b_values, model_b_values, eps = 0.02) # come back rounded


  print("\nConverting fractional lddt to B values", file = log)
  ph = model.get_hierarchy().deep_copy()
  ph.atoms().set_b(fractional_lddt)
  test_model = model.as_map_model_manager().model_from_hierarchy(ph,
     return_as_model = True)
  params.process_predicted_model.maximum_fraction_close = 0.5
  params.process_predicted_model.b_value_field_is = 'lddt'
  params.process_predicted_model.remove_low_confidence_residues = False
  params.process_predicted_model.split_model_by_compact_regions = False
  params.process_predicted_model.input_lddt_is_fractional = None
  model_info = process_predicted_model(test_model, params)
  model = model_info.model
  model_b_values = model.get_hierarchy().atoms().extract_b()
  assert approx_equal(b_values, model_b_values, eps = 3) # come back very rounded

  ph = model.get_hierarchy().deep_copy()
  ph.atoms().set_b(rmsd_values)
  test_model = model.as_map_model_manager().model_from_hierarchy(ph,
     return_as_model = True)

  print("\nConverting rmsd to B values", file = log)
  params.process_predicted_model.maximum_fraction_close = 0.5
  params.process_predicted_model.b_value_field_is = 'rmsd'
  params.process_predicted_model.remove_low_confidence_residues = False
  params.process_predicted_model.split_model_by_compact_regions = False
  params.process_predicted_model.input_lddt_is_fractional = None
  model_info = process_predicted_model(test_model, params)
  model = model_info.model
  model_b_values = model.get_hierarchy().atoms().extract_b()
  assert approx_equal(b_values, model_b_values, eps = 0.5) # come back rounded

  print("B-values > 59: %s of %s" %(
     (model_b_values > 59).count(True), model_b_values.size()), file = log)

  print("\nConverting rmsd to B values and selecting rmsd < 1.5", file = log)
  params.process_predicted_model.maximum_fraction_close = 0.5
  params.process_predicted_model.b_value_field_is = 'rmsd'
  params.process_predicted_model.remove_low_confidence_residues = True
  params.process_predicted_model.maximum_rmsd = 1.5
  params.process_predicted_model.split_model_by_compact_regions = False
  params.process_predicted_model.input_lddt_is_fractional = None

  model_info = process_predicted_model(test_model, params)
  model = model_info.model
  print("Residues before: %s   After: %s " %(
    test_model.get_hierarchy().overall_counts().n_residues,
    model.get_hierarchy().overall_counts().n_residues,), file = log)

  # Check splitting model into domains
  print("\nSplitting model into domains", file = log)
  model_info = split_model_into_compact_units(model,
      maximum_fraction_close = 0.5, log = log)

  chainid_list = model_info.chainid_list
  print("Segments found: %s" %(" ".join(chainid_list)), file = log)
  assert len(chainid_list) == 2

  # Check processing and splitting model into domains
  print("\nProcessing and splitting model into domains", file = log)

  params.process_predicted_model.maximum_fraction_close = 0.5
  params.process_predicted_model.b_value_field_is = 'lddt'
  params.process_predicted_model.remove_low_confidence_residues = True
  params.process_predicted_model.maximum_rmsd = 1.5
  params.process_predicted_model.split_model_by_compact_regions = True
  params.process_predicted_model.maximum_domains = 3
  model_info = process_predicted_model(m,  params, log = log)

  chainid_list = model_info.chainid_list
  print("Segments found: %s" %(" ".join(chainid_list)), file = log)
  assert len(chainid_list) == 2


  mmm = model_info.model.as_map_model_manager()
  mmm.write_model('model_with_groupings.pdb')
  residue_count = []
  expected_residue_count = [84, 88]
  for chainid in chainid_list:
    selection_string = "chain %s" %(chainid)
    ph = model_info.model.get_hierarchy()
    asc1 = ph.atom_selection_cache()
    sel = asc1.selection(selection_string)
    m1 = model_info.model.select(sel)
    n = m1.get_hierarchy().overall_counts().n_residues
    print("Residues in %s: %s" %(
      selection_string, n),
       file = log)
    residue_count.append(n)
  assert expected_residue_count == residue_count

  # Now process and use pae model and pae model file
  print("\nProcessing and splitting model into domains with pae", file = log)


  params.process_predicted_model.maximum_fraction_close = 0.5
  params.process_predicted_model.b_value_field_is = 'lddt'
  params.process_predicted_model.remove_low_confidence_residues = True
  params.process_predicted_model.maximum_rmsd = 0.7
  params.process_predicted_model.split_model_by_compact_regions = True
  params.process_predicted_model.maximum_domains = 3
  params.process_predicted_model.pae_power= 2
  model_info = process_predicted_model(pae_m,  params, pae_matrix = pae_matrix,
     log = log)


if __name__ == "__main__":

  t0 = time.time()
  tst_01()
  print ("Time:", time.time()-t0)
  print ("OK")
