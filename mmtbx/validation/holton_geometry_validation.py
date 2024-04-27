from __future__ import absolute_import, division, print_function
from libtbx.utils import null_out, Sorry
from libtbx import group_args
from cctbx.geometry_restraints.linking_class import linking_class
from libtbx import adopt_init_args
from scitbx.array_family import flex
from scipy.special import erf, erfinv, gammainc
import math
import os, sys

''' Port of James Holton's untangle_score.py in Python. Header from original
file follows.


throw every validation we can at it                       -James Holton 1-18-24

  assemble an overall geometry score
  "energy" = (deviate/sigma)**2

  fraction of Gaussian-random events observed with at least one
     event above x=deviate/sigma
  in any of Nrep trials is:
  fracabove = 1-erf(deviate/sigma/sqrt(2))**Nreps

i.e. you expect to see > 2 sigma deviates ~5% of the time,
unless you look at 10 at a time, in which case 37% of the time you see at least one > 2 sigma deviate

so, probability that a x-sigma deviate is not noise, given N samples is:
Pnotnoise(deviate,N) = erf(abs(deviate)/sigma/sqrt(2))**N

final score will be sum of average square sigma deviations for each
validation metric, plus the sum of:
         Pnotnoise*(worstdeviate/sigma)^2
 over worst outlier of each validation metric
 to avoid overwhelming all other considerations, these outlier quantities
   are softened
 so that energy values above 10 are substitued with 10+log(energy)-log(10)
 also, the worst-clash energy is multiplied by the number of clashes

sigma level where probability is 0.5 is given by:
   sigma_P50 = sqrt(2)*inverf((0.5)**(1./exponent))
   which is reasonably approximated by:
   a2 = -0.0341977; a1 = 0.852225 ; a0 = 1.10939
   sigma_P50(log10e) = a2*log10e**2+a1*log10e+a0
   log10e(exponent) = log(exponent)/log(10)
   softer scoring function for large exponents is:
   Pnotnoise = 1-2**((deviate/sigma_P50)**5) , naturally clipped at [-1:1]

for rama, rota etc.
convert probability/frequency back to sigma with:
  deviate/sigma = inverf(1-probOK)*sqrt(2)

cbetadev - use 0.05 A as "sigma"

nonbonds:  use Leonard-Jones to convert to energy [-1:inf], but dont
let "worst" or avg be negative

omega twist: energy=((sin(omega)/0.07)^2+(1+cos(omega))^10)/(proxPRO*2+1)
where proxPRO means a neighboring residue is proline

  '''


def holton_geometry_validation(dm = None,
     filename = None,
     model = None,
     get_individual_residue_scores = None,
     round_numbers = True,
     worst_clash_is_one_plus_n_times_worst_clash = True,
     clash_energy_add_n = True,
     minimum_nonbond_score_to_be_worst = -0.1,
     minimum_nonbond_score_to_be_included_in_average = 0,
     keep_hydrogens = False,  # redo the hydrogens
     ignore_cis_peptides = False,
     ignore_h_except_in_nonbond = True,
     ignore_arg_h_nonbond = True,
     ignore_bond_lengths_with_h = False,
     ignore_water_h_bonds = False,
     rotalyze_max_energy = 99,
     overall_max_energy = 400,
     omega_angle_sigma = 4,  # Sigma for omega angle
     cbetadev_sigma = 0.05,  # Sigma for CB position
     clashscore_ideal_dist = 3,  # Ideal distance in LJ for clashscore result
     lj_dist_that_yields_zero = 6, # Distance for modified LJ to cross zero
     softPnna_params = group_args(group_args_type = 'softPnna_params',
       y0= 1,
       a2= -0.0192266,
       a1 = 0.751694,
       a0 = 1.12482,
       mx = 0.21805,
       my = 0.736621),
     verbose = False,
     log = sys.stdout,
    ):

  # Capture arguments and save in info
  O = group_args(group_args_type = 'info')
  adopt_init_args(O,locals())
  info = O

  # Read in model and get sequence
  get_model(info)

  # Get basic geometry restraints and filter them
  get_geometry_results(info)

  add_omega_results(info)

  filter_geometry_results(info)

  # Get CBdev analysis
  add_cbetadev_results(info)

  # Get Ramachandran analyses
  add_rama_results(info)

  # Get Rotamer analyses
  add_rotamer_results(info)

  # Get Molprobity results
  add_clashscore_results(info)

  # Sort them all
  sort_geometry_results(info)

  # Get worst and average values for each and rescale if desired
  analyze_geometry_values(info)

  # Get residue-based scores
  get_residue_scores(info)

  print_results(info)

  return info


def add_clashscore_results(info):
  clashscore_result = group_args(group_args_type = 'CLASHSCORE result',
    name = 'CLASH',
    value_list = [])
  info.geometry_results[clashscore_result.name] = clashscore_result

  from mmtbx.validation.clashscore import clashscore
  clashes = clashscore(
      info.model.get_hierarchy(),
      fast = False,
      keep_hydrogens=True,
      time_limit=120,
      save_modified_hierarchy=False,
      verbose=False,
      do_flips=False,
      out=null_out())

  for r in clashes.results:
    if info.round_numbers:
      delta = float("%.3f" %(r.overlap))
    else:
      delta = r.overlap

    # Convert mmtbx.validation.atom_info to atom_label objects
    labels = [atom_label(ai) for ai in r.atoms_info]

    dist = info.clashscore_ideal_dist + delta
    energy = lj(dist, info.clashscore_ideal_dist,
           dist_that_yields_zero = info.lj_dist_that_yields_zero,
           round_numbers = info.round_numbers)

    v = group_args(group_args_type = 'clashscore result as standard value ',
      clashscore_result = r, # contains all results
      as_string = "CLASH: Energy = %.4f dev = %.3f \n   Atom 1:   %s \n   Atom 2:   %s " %(
        tuple([energy, delta] + [l.as_string() for l in labels])),
      resseq = None,
      delta = delta,
      residual = energy,
      labels = labels,
      )
    clashscore_result.value_list.append(v)

def add_rotamer_results(info):
  rotamer_result = group_args(group_args_type = 'ROTA result',
    name = 'ROTA',
    value_list = [])
  info.geometry_results[rotamer_result.name] = rotamer_result
  args = [info.filename]
  from iotbx.cli_parser import run_program
  from mmtbx.programs import rotalyze
  result = run_program(program_class=rotalyze.Program,args=args,
     logger = null_out())

  for r in result.results:
    if info.round_numbers:
      prob = float("%.1f" %(r.score))/100
    else:
      prob = r.score/100

    prob = min(1.0, max(0.0, prob))
    if prob == 1:
      energy = info.rotalyze_max_energy
    else:
      prob = float("%.35g" %(prob)) -1.0e-16
      energy = energy_from_probability(prob)
      energy = min(energy, info.rotalyze_max_energy)
    v = group_args(group_args_type = 'rotamer result as standard value ',
      rotamer_result = r, # contains resseq, resseq_as_int, resname, chain_id
      as_string = "ROTA: Energy = %.4f \n   Residue:  %s %s %s %s" %(
        energy, r.altloc, r.resname, r.chain_id, r.resseq),
      resseq= r.resseq,
      delta = None,
      residual = energy,
      labels = [atom_label_from_info(
          r.chain_id, r.resseq,  r.altloc, r.resname, 'CA' )],
      )
    rotamer_result.value_list.append(v)

def add_rama_results(info):
  rama_result = group_args(group_args_type = 'RAMA result',
    name = 'RAMA',
    value_list = [])
  info.geometry_results[rama_result.name] = rama_result
  args = [info.filename]
  from iotbx.cli_parser import run_program
  from mmtbx.programs import ramalyze
  result = run_program(program_class=ramalyze.Program,args=args,
     logger = null_out())
  for r in result.results:
    if info.round_numbers:
      prob = float("%.2f" %(r.score))/100
    else:
      prob = r.score/100
    energy = energy_from_probability(prob)
    v = group_args(group_args_type = 'rama result as standard value ',
      rama_result = r, # contains resseq, resseq_as_int, resname, chain_id
      as_string = "RAMA: Energy = %.4f \n   Residue:  %s %s %s %s" %(
        energy, r.altloc, r.resname, r.chain_id, r.resseq),
      resseq= r.resseq,
      delta = None,
      residual = energy,
      labels = [atom_label_from_info(
          r.chain_id, r.resseq,  r.altloc, r.resname, 'CA' )],
      )
    rama_result.value_list.append(v)

def add_cbetadev_results(info):
  cbetadev_result = group_args(group_args_type = 'CBETADEV result',
    name = 'CBETADEV',
    value_list = [])
  info.geometry_results[cbetadev_result.name] = cbetadev_result
  args = [info.filename]
  from iotbx.cli_parser import run_program
  from mmtbx.programs import cbetadev
  result = run_program(program_class=cbetadev.Program,args=args,
     logger = null_out())
  for r in result.results:
    if info.round_numbers:
      delta = float("%.3f" %(r.deviation))
    else:
      delta = r.deviation
    energy = (delta/info.cbetadev_sigma)**2
    v = group_args(group_args_type = 'cbetadev result as standard value ',
      cbetadev_result = r, # contains resseq, resseq_as_int, resname, chain_id
      as_string =
        "CBETADEV: Energy = %.4f dev = %.3f A \n   Residue: %s %s %s %s" %(
        energy, r.deviation, r.altloc, r.resname, r.chain_id, r.resseq),
      resseq= r.resseq,
      delta = r.deviation,
      residual = energy,
      labels = [atom_label_from_info(
          r.chain_id, r.resseq,  r.altloc, r.resname, 'CA' )],
      )
    cbetadev_result.value_list.append(v)

def select_geometry_result(info, name = None, skip_name = None):
  result_list = []
  for key in info.geometry_results.keys():
    if name and info.geometry_results[key].name == name:
      return info.geometry_results[key]
    elif skip_name and info.geometry_results[key].name != skip_name:
      result_list.append(info.geometry_results[key])
  return result_list

def get_residue_scores(info):
  ''' run through all residues, select those results that involve each residue,
      recalculate overall score.
  '''
  info.residue_scores = []
  if not info.get_individual_residue_scores:
    return # nothing to do

  # For now, only one chain_id
  base_info = get_base_info(info)

  for chain_id in list(info.chain_dict.keys()):
    for resseq in info.chain_dict[chain_id].sequence_dict.keys():
      working_info = select_residue_info(info, base_info, chain_id, resseq)
      analyze_geometry_values(working_info)
      info.residue_scores.append(working_info.sum_energy)
      if info.verbose:
        print("Residue score for %s: %.2f" %(resseq, working_info.sum_energy),
          file = info.log)


def get_base_info(info):
  ''' just get everything except geometry'''

  log = info.log
  dm = info.dm
  geometry_results = info.geometry_results
  chain_dict = info.chain_dict
  model = info.model

  info.log = None
  info.dm= None
  info.geometry_results = {}
  info.chain_dict = {}
  info.model = None

  from copy import deepcopy
  base_info = deepcopy(info)

  info.log = log
  info.dm = dm
  info.geometry_results = geometry_results
  info.chain_dict = chain_dict
  info.model = model
  return base_info

def select_residue_info(info, base_info, chain_id, resseq):
  new_info = base_info # overwrite base_info each time

  keys = list(info.geometry_results.keys())
  for key in keys:
    orig = info.geometry_results[key]
    new_info.geometry_results[key] = group_args(
       group_args_type = orig.group_args_type,
       name = orig.name,
       value_list = [])

    for value in info.geometry_results[key].value_list:
      if labels_contain(value.labels, chain_id = chain_id, resseq = resseq):
        new_info.geometry_results[key].value_list.append(value)
  sort_geometry_results(new_info)
  return new_info

def labels_contain(labels, chain_id = None, resseq = None):
  for label in labels:
    if ((chain_id is None) or (label.chain_id().strip() == chain_id.strip())) and \
       ((resseq is None) or   (label.resseq().strip() == resseq.strip())):
      return True

  return False

def analyze_geometry_values(info):
  keys = list(info.geometry_results.keys())
  keys.sort()
  sum_energy = 0.0
  for key in keys:
    result = info.geometry_results[key]
    result.pnna = None
    result.energy = None
    result.chisq = None
    result.energy_using_mean = None
    result.worst_residual = None
    result.mean_residual = None

    if not result.value_list:
      continue

    worst_value = result.value_list[0]

    # Get worst residual value and save it:

    # Some special cases allowed in choosing worst value

    if info.worst_clash_is_one_plus_n_times_worst_clash and \
        result.name == 'CLASH': # Special case
      result.worst_residual = \
          1 + (len(result.value_list) * worst_value.residual)

    elif info.minimum_nonbond_score_to_be_worst is not None and \
        result.name == 'NONBOND' and (
          worst_value.residual < info.minimum_nonbond_score_to_be_worst):
      result.worst_residual = None
    else: # usual
      result.worst_residual = worst_value.residual

    # Get average value and save it too:
    value_list = result.value_list

    # Special case: non bond scores with low values may be excluded
    if info.minimum_nonbond_score_to_be_included_in_average is not None and \
       result.name == 'NONBOND':
      value_list = []
      for v in result.value_list:
        if v.residual > info.minimum_nonbond_score_to_be_included_in_average:
          value_list.append(v)

    values = flex.double()
    for v in value_list:
      values.append(v.residual)
    result.mean_residual = values.min_max_mean().mean
    result.n = values.size()
    if result.n < 1:

      continue # nothing to do

    energy = max(0, min(info.overall_max_energy, result.worst_residual))
    delta = energy**0.5
    result.pnna = softPnna(delta, result.n, info.softPnna_params)

    energy = filtered_energy(energy)

    # Special case for CLASH
    if result.name == "CLASH" and info.clash_energy_add_n:
      energy += result.n

    result.energy = result.pnna * energy

    sum_energy += result.energy

    # Repeat using mean instead of worst, multiplying mean energy * Chisq
    energy_mean = max(0, min(info.overall_max_energy, result.mean_residual))
    ssd = result.mean_residual * result.n
    result.chisq = chisq(ssd, result.n)

    energy_mean = filtered_energy(energy_mean)

    # Special case for CLASH
    if result.name == "CLASH" and info.clash_energy_add_n:
      energy_mean += result.n

    result.energy_using_mean = result.chisq * energy_mean

    sum_energy += result.energy_using_mean

  info.sum_energy = sum_energy

def sort_geometry_results(info):
  keys = list(info.geometry_results.keys())
  keys.sort()
  for key in keys:
    result = info.geometry_results[key]
    result.value_list = sorted(
      result.value_list, key = lambda v:
       v.residual, reverse = True)


def print_results(info):

  print("\nSUMMARY of Holton geometry validation scoring for %s" %(
    info.filename), file = info.log)
  print(file = info.log)

  info.result_table = {}
  info.result_header_row = [
       'Category','N','Mean','Worst','Chisq','Pnna','Energy','Using mean']
  fmt = len(info.result_header_row) * "%8s "
  print(file = info.log)
  print(fmt  %(tuple(info.result_header_row)), file = info.log)
  print(file = info.log)

  keys = list(info.geometry_results.keys())
  keys.sort()
  for key in keys:
    result = info.geometry_results[key]
    value_list = result.value_list
    info.result_table[result.name] = \
        [result.name,
         "%s" %str(len(value_list)),
         "%.2f" %(result.mean_residual) if result.mean_residual \
           is not None else "None",
         "%.2f" %(result.worst_residual) if result.worst_residual \
           is not None else "None",
         "%.4f" %(result.chisq) if result.chisq \
           is not None else "None",
         "%.4f" %(result.pnna) if result.pnna \
           is not None else "None",
         "%.4f" %(result.energy) if result.energy\
           is not None else "None",
         "%.4f" %(result.energy_using_mean) if result.energy_using_mean\
           is not None else "None" ]
    print(fmt  %(tuple(info.result_table[result.name])),
         file = info.log)

  info.worst_table = {}
  info.worst_header_row = ['   -----  Worst deviation in each category -----']
  fmt = "%s"
  print(file = info.log)
  print(fmt %(tuple(info.worst_header_row)), file = info.log)
  print(file = info.log)
  for key in keys:
    result = info.geometry_results[key]
    value_list = result.value_list
    info.worst_table[result.name] = [
       value_list[0].as_string if value_list else ""]
    print(fmt  %(tuple(info.worst_table[result.name])),
         file = info.log)

  if info.ignore_h_except_in_nonbond:
    print("\nTotal H interactions (not in non-bonded) removed: %s" %(
       info.ignore_h_except_in_nonbond_removed), file = info.log)
  if info.ignore_arg_h_nonbond:
    print("Total ARG H-nonbond removed: %s" %(
       info.ignore_arg_h_nonbond_removed), file = info.log)
  if info.ignore_water_h_bonds:
    print("Total water H-nonbond removed: %s" %(
       info.ignore_water_h_bonds_removed), file = info.log)
  if info.ignore_bond_lengths_with_h:
    print("Total bonds with H removed: %s" %(
       info.ignore_bond_lengths_with_h_removed), file = info.log)
  print("\nOverall geometry energy for %s: %.4f\n" %(info.filename,
     info.sum_energy), file = info.log)

def filter_geometry_results(info):

  info.ignore_arg_h_nonbond_removed = 0
  info.ignore_water_h_bonds_removed = 0
  info.ignore_bond_lengths_with_h_removed = 0
  info.ignore_h_except_in_nonbond_removed = 0

  if not True in [info.ignore_arg_h_nonbond, info.ignore_bond_lengths_with_h,
      info.ignore_water_h_bonds]:
    return # nothing to do

  if info.ignore_h_except_in_nonbond:
    # Remove anything with element H except in NONBOND
    for result in select_geometry_result(info,skip_name = 'NONBOND'):
      for v in result.value_list:
        elements = [l.atom.element for l in v.labels]
        if 'H' in elements:
          result.value_list.remove(v)
          info.ignore_h_except_in_nonbond_removed += 1
          continue

  if info.ignore_arg_h_nonbond or info.ignore_water_h_bonds:
    nonbond_result = select_geometry_result(info,'NONBOND')
    for v in nonbond_result.value_list:
      elements = [l.atom.element for l in v.labels]
      residues = [l.atom.resname for l in v.labels]

      if info.ignore_arg_h_nonbond:
        if elements == ['H','H'] and residues == ['ARG','ARG']:
          nonbond_result.value_list.remove(v)
          info.ignore_arg_h_nonbond_removed += 1
          continue

      if info.ignore_water_h_bonds:
        found = False
        for element, res in zip(elements, residues):
          if element == 'H' and res in ['WAT','HOH']:
            nonbond_result.value_list.remove(v)
            found = True
        if found:
          info.ignore_water_h_bonds_removed += 1
          continue

  if info.ignore_bond_lengths_with_h:
    bond_result = select_geometry_result(info,'BOND')
    for v in bond_result.value_list:
      elements = [l.atom.element for l in v.labels]
      if 'H' in elements:
        bond_result.value_list.remove(v)
        info.ignore_bond_lengths_with_h_removed += 1
        continue

def add_omega_results(info):
  torsion_result = select_geometry_result(info, 'TORSION')
  omega_result = group_args(group_args_type = 'OMEGA result',
    name = 'OMEGA',
    value_list = [])
  info.geometry_results[omega_result.name] = omega_result

  degtorad = math.atan2(1,1)/45.
  for v in torsion_result.value_list:
    atom_names = [l.atom.name for l in v.labels]
    if atom_names[0].strip() != 'CA' or atom_names[3].strip() != 'CA': continue
    resseq = v.labels[0].atom.resseq.strip()
    chain_id = v.labels[0].atom.chain_id
    altloc = v.labels[0].atom.altloc
    resnam = info.chain_dict[chain_id].sequence_dict.get(resseq)
    omega = v.model * degtorad
    n_pro = get_n_pro(info, chain_id, resseq)
    new_v = group_args(group_args_type = 'OMEGA derived result',
     )
    new_v.model = v.model
    new_v.n_pro = n_pro
    new_v.labels = [v.labels[0]]
    # Score function assuming sigma of 4 degrees: 0.07 is sin(4 deg)
    sigma = math.sin(info.omega_angle_sigma * degtorad)
    if info.round_numbers:
      sigma = float("%.3f" %(sigma))
    new_v.residual=((math.sin(omega)/sigma)**2+
         (1+math.cos(omega))**10)/(n_pro*2+1)
    at = new_v.labels[0].atom
    new_v.as_string = "OMEGA: Energy = "+\
     "%.6f (%s Proline) Angle: %.2f deg\n   Residue:  %s %s %s %s" %(
        new_v.residual,  new_v.n_pro, new_v.model,
         at.altloc, at.resname, at.chain_id, at.resseq)
    if info.ignore_cis_peptides and math.cos(omega) > 0:
      continue # skip it
    omega_result.value_list.append(new_v)


def get_n_pro(info, chain_id, resseq):
  n_pro = 0
  resseq_as_int = int(resseq)
  for i in [resseq_as_int - 1, resseq_as_int + 1]:
    if info.chain_dict[chain_id].sequence_dict.get(str(i),None) == 'PRO':
      n_pro += 1
  return n_pro

def get_geometry_results(info):
  geometry = info.model.get_restraints_manager().geometry
  sites_cart= info.model.get_sites_cart()
  site_labels = get_site_labels(info.model)
  origin_ids = linking_class()
  pair_proxies = geometry.pair_proxies(sites_cart=sites_cart)

  name_dict = {'nonbonded': 'NONBOND',
        'angle_proxies': 'ANGLE',
        'bond_proxies': 'BOND',
        'chirality_proxies': 'CHIR',
        'dihedral_proxies': 'TORSION',
        'planarity_proxies': 'PLANE',}
  geometry_results = {}
  info.geometry_results = geometry_results
  for key in name_dict.keys():
    geometry_results[key] = group_args(group_args_type = '%s result' %(key),
      name = key,
      value_list = [])

  # Non-bonded
  if pair_proxies.nonbonded_proxies is not None:
    proxy_name = 'nonbonded'
    result = pair_proxies.nonbonded_proxies.show_sorted(
          by_value="delta",
          sites_cart=sites_cart,
          site_labels=site_labels,
          f=null_out(),
          suppress_model_minus_vdw_greater_than = None,
          prefix="  ",
          max_items=None,
          return_result = True)
    # Calculate residual from Lennard-Jones potential
    for v in result.value_list:
      v.residual = lj(v.model,v.ideal,
           dist_that_yields_zero = info.lj_dist_that_yields_zero,
            round_numbers = info.round_numbers)
      v.delta = v.ideal - v.model
      v.group_args_type += " residual is LJ(model, ideal)"
      v.as_string = "NONBOND: Energy = %.6f dev = %.3f A  obs = %.3f target = %.3f\n   Atom 1:  %s\n   Atom 2:  %s" %(
        v.residual, v.delta, v.model, v.ideal,
         v.labels[0].as_string(), v.labels[1].as_string())
    geometry_results[proxy_name] = result

  for proxy_name in ['bond_proxies', 'angle_proxies', 'dihedral_proxies',
       'chirality_proxies', 'planarity_proxies']:
    proxies = getattr(pair_proxies, proxy_name,
                getattr(geometry, proxy_name,
                  getattr(geometry, proxy_name, None)))
    if proxy_name in ['bond_proxies', 'angle_proxies']:
      origin_id=origin_ids.get_origin_id('covalent geometry')
    else:
      origin_id=None
    if proxies:
      result = proxies.show_sorted(
          by_value="residual",
          sites_cart=sites_cart,
          site_labels=site_labels,
          f=null_out(),
          origin_id=origin_id,
          return_result = True)
      geometry_results[proxy_name] = result
      if result.group_args_type == 'Bond restraints':
        for v in result.value_list:
          v.as_string = "BOND: Energy = %.2f dev = %.3f A  obs = %.3f target = %.3f sigma = %.2f\n   Atom 1:  %s\n   Atom 2:  %s " %(
            v.residual, v.delta, v.model, v.ideal, v.sigma,
             v.labels[0].as_string(), v.labels[1].as_string())
      elif result.group_args_type == 'Bond angle restraints':
        for v in result.value_list:
          v.as_string = "ANGLE: Energy = "+\
           "%.2f dev = %.2f deg  obs = %.2f target = %.2f sigma = %.1f\n   Atom 1:  %s\n   Atom 2:  %s\n   Atom 3:  %s " %(
            v.residual, v.delta, v.model, v.ideal, v.sigma,
             v.labels[0].as_string(), v.labels[1].as_string(),
             v.labels[2].as_string())
      elif result.group_args_type == 'Dihedral angle restraints':
        for v in result.value_list:
          v.as_string = \
          "TORSION: Energy = %.4f dev = %.1f deg  obs = %.1f target = %.0f sigma = %.1f\n   Atom 1:  %s\n   Atom 2:  %s\n   Atom 3:  %s\n   Atom 4:  %s" %(
            v.residual, v.delta, v.model, v.ideal, v.sigma,
             v.labels[0].as_string(), v.labels[1].as_string(),
             v.labels[2].as_string(), v.labels[3].as_string())
      elif result.group_args_type == 'Chirality restraints':
        for v in result.value_list:
          v.as_string = \
          "CHIR: Energy = %.3f delta = %.2f A**3  obs = %.2f target = %.2f sigma = %.1f\n   Atom 1:  %s\n   Atom 2:  %s\n   Atom 3:  %s\n   Atom 4:  %s" %(
            v.residual, v.delta, v.model, v.ideal, v.sigma,
             v.labels[0].as_string(), v.labels[1].as_string(),
             v.labels[2].as_string(), v.labels[3].as_string())
      elif result.group_args_type == 'Planarity restraints':
        for v in result.value_list:
          v.residual = energy(v.delta, v.sigma,
             round_numbers = info.round_numbers)
          v.as_string = \
          "PLANE: Energy = %.4f dev = %.3f A  sigma = %.2f" %(
            v.residual, -v.delta, v.sigma)
          i_at = 0
          for x in v.labels:
            i_at += 1
            v.as_string += "\n   Atom %s:  %s" %(i_at, x.as_string())
  for key in name_dict.keys():
    geometry_results[key].name = name_dict[key]

def get_model(info):
  if not info.filename:
    raise Sorry("Please supply a filename even if a model object is supplied")
  if not info.model:
    if not os.path.isfile(info.filename):
      raise Sorry("The filename %s is missing" %(info.filename))
  if not info.dm:
    from iotbx.data_manager import DataManager
    info.dm = DataManager()
  if not info.model:
    info.model = info.dm.get_model(info.filename)
  if not os.path.isfile(info.filename):  # write to it
    info.dm.write_model_file(info.model, info.filename)
  info.model.set_log(null_out())
  if (not info.keep_hydrogens):
    info.model.add_crystal_symmetry_if_necessary()
    if info.model.has_hd():
      info.model.get_hierarchy().remove_hd(reset_i_seq=True)
  info.model.process(make_restraints=True)
  info.chain_dict = {}
  for chain_id in info.model.chain_ids(unique_only = True):
    entry = group_args(group_args_type = 'chain entry', )
    info.chain_dict[chain_id] = entry

    entry.chain = info.model.apply_selection_string(
        "chain %s and not water" %(chain_id))
    entry.sequence = entry.chain.as_sequence(as_string = True)
    entry.sequence_dict = get_sequence_dict(entry.chain)
    entry.residue_numbers = list(entry.sequence_dict.keys())
    entry.residue_numbers.sort()

def get_sequence_dict(m): # Assumes 1 model, one chain, any conformers
  sequence_dict = {}
  for model in m.get_hierarchy().models()[:1]:
    for chain in model.chains()[:1]:
      for conformer in chain.conformers()[:1]:
        for residue in conformer.residues():
          sequence_dict[residue.resseq.strip()] = residue.resname
  return sequence_dict

def get_site_labels(model):
  labels = []
  for at in model.get_hierarchy().atoms_with_labels():
    labels.append(atom_label(at))
  return labels

class pseudo_atom:
  def __init__(self, chain_id, resseq,  altloc, resname, name):
     adopt_init_args(self,locals())

def atom_label_from_info(chain_id, resseq,  altloc, resname, name):
  return atom_label(pseudo_atom(chain_id, resseq,  altloc, resname, name))

class atom_label:
  ''' Class to hold an atom_with_labels object and deliver
      the same string as model.get_site_labels() and allow right-ward addition
      of str representation and length of the str representation.
      Used to pass through pdb_interpretation and return the full atom so
      that all characteristics are preserved
      NOTE: Can be used for mmtbx.validation.atom_info as well
  '''
  def __init__(self, atom):
    self.atom = atom

  def as_string(self):
    ''' Print a longer string that has all fields separated '''
    atom = self.atom
    return "%3s %s %3s %2s %4s%s" %(atom.name,atom.altloc,atom.resname,
       atom.chain_id,atom.resseq,atom.icode)

  def __repr__(self):
    #  pdb=" CG1AVAL A   1 "  Mimics model.get_site_labels()
    atom = self.atom
    return 'pdb="%3s%1s%3s%2s%4s "' %(
      atom.name, atom.altloc, atom.resname, atom.chain_id, atom.resseq)
  def __len__(self):
    return len(self.__repr__())
  def __add__(self, a):
    return self.__repr__() + a
  def resseq(self):
    return self.atom.resseq
  def chain_id(self):
    return self.atom.chain_id
  def icode(self):
    return self.atom.icode

def energy(delta, sigma, round_numbers = None):
  if round_numbers:
    delta = float("%.3f" %(delta))
    sigma = float("%.3f" %(sigma))
  energy = (delta/max(1.e-10, sigma))**2
  return energy

def lj0(r, r0, round_numbers = None):
   ''' Lennard-Jones potential'''
   if round_numbers:
     r = float("%.3f" %(r))
     r0 = float("%.3f" %(r0))
   if( r == 0 ):
     return 1e40
   else:
     return 4*((r0*2**(-1./6)/r)**12-(r0*2**(-1./6)/r)**6)

def lj(r, r0, dist_that_yields_zero = 6, round_numbers = None):
   ''' Modified Lennard-Jones potential with value of dist_that_yields_zero at 6 (A)'''
   return lj0(r,r0, round_numbers = round_numbers)- \
       lj0(dist_that_yields_zero,r0, round_numbers = round_numbers)

def filtered_energy(energy):
  if energy <= 0:
    return 0
  elif energy <= 10:
    return energy
  else:
    return 10 + math.log(energy/10)


def softPnna(delta, n, params = None):
  return  1 - 2.0**clip(
     -abs(delta/asigma_Pnn50(safelog(n), params))**exponent(
         delta,safelog(n), params) ,1000)

def Pnn(delta,n):
  return erf(abs(delta)/math.sqrt(2))**n

def sigma_Pnn50(n):
  # sigma deviates that have Pnn=0.5
  return math.sqrt(2)*erfinv((0.5)**(1./n))

def asigma_Pnn50(x, params = None):
 return (params.a2*x**2+params.a1*x+params.a0)

def exponent(x, y, params = None):
  # exponent to use in approximation
  return params.my*y + params.mx*x + params.y0

def safelog(x):
 if x <= 0:
   return 0
 else:
   return math.log10(x)

def s_thresh(x):
  return math.sqrt(2)*erfinv((1.e-30)**(1./x))

def safePnn(delta,n):
  if delta < s_thresh(n):
    return 0
  else:
    return Pnn(delta,n)

def clip(x,y):
  if x > y:
    return y
  elif x <= -y:
    return -y
  else:
    return x

def chisq(x,k):
  if k<1:
    return 0
  elif k<3200:
    return gammainc(k/2.0,x/2.0)
  else:
    return chisqhi(x,k)

def chisqhi(x,k):
  return 0.5 + 0.5 * erf((x**0.5)-(k**.5))

def energy_from_probability(prob):
  prob = max(0, min(1, prob))
  return (erfinv(1-prob))**2

if __name__=="__main__":
  import sys
  filename = sys.argv[1] if len(sys.argv[1:]) > 0 else None
  holton_geometry_validation(filename = filename)
