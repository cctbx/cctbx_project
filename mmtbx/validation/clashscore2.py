
"""
All-atom contact analysis.  Calls mmtbx.reduce functions and mmtbx.probe2 program.
This is a rewrite of the original clashscore that used stand-alone
external reduce and probe programs.
"""

from __future__ import absolute_import, division, print_function
from mmtbx.validation.clashscore import clash
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.reduce import Optimizers
from mmtbx.programs import probe2
from mmtbx.validation import validation, atoms, atom_info
from libtbx.utils import Sorry, null_out
import iotbx.pdb
import iotbx.cli_parser
import mmtbx
import os
import re
import sys
import six
import json
import copy
import tempfile

def remove_models_except_index(model_manager, model_index):
    hierarchy = model_manager.get_hierarchy()
    models = hierarchy.models()
    if model_index < len(models):
        selected_model = models[model_index]
        for model in models:
            if model != selected_model:
                hierarchy.remove_model(model=model)
    return model_manager

class clashscore2(validation):
  __slots__ = validation.__slots__ + [
    "clashscore",
    "clashscore_b_cutoff",
    "clash_dict",
    "clash_dict_b_cutoff",
    "list_dict",
    "b_factor_cutoff",
    "fast",
    "condensed_probe",
    "probe_file",
    "probe_clashscore_manager"
  ]
  program_description = "Analyze clashscore for protein model"
  gui_list_headers = ["Atom 1", "Atom 2", "Overlap"]
  gui_formats = ["%s", "%s", ".3f"]
  wx_column_widths = [150, 150, 150] #actually set in GUI's Molprobity/Core.py

  def get_result_class(self) : return clash

  def __init__(self,
      probe_parameters,
      data_manager,
      fast = False, # do really fast clashscore, produce only the number
      condensed_probe = False, # Use -CON for probe. Reduces output 10x.
      keep_hydrogens=True,
      nuclear=False,
      force_unique_chain_ids=False,
      time_limit=120,
      b_factor_cutoff=None,
      verbose=False,
      do_flips=False,
      out=sys.stdout):
    validation.__init__(self)
    self.b_factor_cutoff = b_factor_cutoff
    self.fast = fast
    self.condensed_probe = condensed_probe
    self.clashscore = None
    self.clashscore_b_cutoff = None
    self.clash_dict = {}
    self.clash_dict_b_cutoff = {}
    self.list_dict = {}
    self.probe_file = None
    if verbose:
      if not nuclear:
        print("\nUsing electron cloud x-H distances and vdW radii")
      else:
        print("\nUsing nuclear cloud x-H distances and vdW radii")
    import iotbx.pdb
    from scitbx.array_family import flex
    from mmtbx.validation import utils

    data_manager_model = data_manager.get_model()
    # Fix up bogus unit cell when it occurs by checking crystal symmetry.
    data_manager_model.add_crystal_symmetry_if_necessary()

    # If we've been asked to, add hydrogens to all of the models in the PDB hierarchy
    # associated with our data_manager_model.
    data_manager_model,_ = check_and_add_hydrogen(
      probe_parameters=probe_parameters,
      data_manager_model=data_manager_model,
      nuclear=nuclear,
      verbose=verbose,
      keep_hydrogens=keep_hydrogens,
      do_flips = do_flips,
      log=out)

    # First we must rebuild the model from the new hierarchy so that the copy can succeed.
    # Make a copy of the original model to use for submodel processing, we'll trim atoms out
    # of it for each submodel.
    ro = data_manager_model.get_restraint_objects()
    data_manager_model = mmtbx.model.manager(
      model_input       = None,
      pdb_hierarchy     = data_manager_model.get_hierarchy(),
      stop_for_unknowns = False,
      crystal_symmetry  = data_manager_model.crystal_symmetry(),
      restraint_objects = ro,
      log               = None)
    original_model = data_manager_model.deep_copy()

    pdb_hierarchy = data_manager_model.get_hierarchy()
    n_models = len(pdb_hierarchy.models())
    use_segids = utils.use_segids_in_place_of_chainids(
                   hierarchy=pdb_hierarchy)
    for i_mod, model in enumerate(pdb_hierarchy.models()):

      # Select only the current submodel from the hierarchy
      submodel = original_model.deep_copy()
      remove_models_except_index(submodel, i_mod)

      # Construct a hierarchy for the current submodel
      r = iotbx.pdb.hierarchy.root()
      mdc = submodel.get_hierarchy().models()[0].detached_copy()
      r.append_model(mdc)

      occ_max = flex.max(r.atoms().extract_occ())

      # Make yet another model for the new hierarchy
      subset_model_manager = mmtbx.model.manager(
        model_input       = None,
        pdb_hierarchy     = r,
        stop_for_unknowns = False,
        crystal_symmetry  = submodel.crystal_symmetry(),
        restraint_objects = ro,
        log               = None)

      if verbose:
        print("\nFinding clashes with mmtbx.probe2...\n")
      self.probe_clashscore_manager = probe_clashscore_manager(
        nuclear=nuclear,
        fast=self.fast,
        condensed_probe=self.condensed_probe,
        largest_occupancy=occ_max,
        b_factor_cutoff=b_factor_cutoff,
        use_segids=use_segids,
        verbose=verbose,
        model_id=model.id)
      self.probe_clashscore_manager.run_probe_clashscore(data_manager, subset_model_manager)

      self.clash_dict[model.id] = self.probe_clashscore_manager.clashscore
      self.clash_dict_b_cutoff[model.id] = self.probe_clashscore_manager.\
                                           clashscore_b_cutoff
      self.list_dict[model.id] = self.probe_clashscore_manager.bad_clashes
      if (n_models == 1) or (self.clashscore is None):
        self.results = self.probe_clashscore_manager.bad_clashes
        self.n_outliers = len(self.results)
        self.clashscore = self.probe_clashscore_manager.clashscore
        self.clashscore_b_cutoff = self.probe_clashscore_manager.\
                                   clashscore_b_cutoff

  def get_clashscore(self):
    return self.clashscore

  def get_clashscore_b_cutoff(self):
    return self.clashscore_b_cutoff

  def show_old_output(self, out=sys.stdout, verbose=False):
    self.print_clashlist_old(out=out)
    self.show_summary(out=out)

  def show_summary(self, out=sys.stdout, prefix=""):
    if self.clashscore is None:
      raise Sorry("PROBE output is empty. Model is not compatible with PROBE.")
    elif (len(self.clash_dict) == 1):
      #FIXME indexing keys can break py2/3 compat if more than 1 key
      k = list(self.clash_dict.keys())[0]
      #catches case where file has 1 model, but also has model/endmdl cards
      print(prefix + "clashscore = %.2f" % self.clash_dict[k], file=out)
      if self.clash_dict_b_cutoff[k] is not None and self.b_factor_cutoff is not None:
        print("clashscore (B factor cutoff = %d) = %f" % \
          (self.b_factor_cutoff,
           self.clash_dict_b_cutoff[k]), file=out)
    else:
      for k in sorted(self.clash_dict.keys()):
        print(prefix + "MODEL %s clashscore = %.2f" % (k,
          self.clash_dict[k]), file=out)
        if self.clash_dict_b_cutoff[k] is not None and self.b_factor_cutoff is not None:
          print("MODEL%s clashscore (B factor cutoff = %d) = %f" % \
            (k, self.b_factor_cutoff, self.clash_dict_b_cutoff[k]), file=out)

  def print_clashlist_old(self, out=sys.stdout):
    if self.fast:
      print("Bad Clashes >= 0.4 Angstrom - not available in fast=True mode", file=out)
      return
    for k in self.list_dict.keys():
      if k == '':
        print("Bad Clashes >= 0.4 Angstrom:", file=out)
        for result in self.list_dict[k] :
          print(result.format_old(), file=out)
      else:
        print("Bad Clashes >= 0.4 Angstrom MODEL%s" % k, file=out)
        for result in self.list_dict[k] :
          print(result.format_old(), file=out)

  def show(self, out=sys.stdout, prefix="", outliers_only=None, verbose=None):
    if (len(self.clash_dict) == 1):
      for result in self.list_dict[''] :
        print(prefix + str(result), file=out)
    else :
      for k in self.list_dict.keys():
        for result in self.list_dict[k] :
          print(prefix + str(result), file=out)
    self.show_summary(out=out, prefix=prefix)

  def as_JSON(self, addon_json={}):
    if not addon_json:
      addon_json = {}
    addon_json["validation_type"] = "clashscore"
    data = addon_json
    flat_results = []
    hierarchical_results = {}
    residue_clash_list = []
    summary_results = {}
    for k in sorted(self.list_dict.keys()):
      for result in self.list_dict[k]:
        flat_results.append(json.loads(result.as_JSON()))
        hier_result = json.loads(result.as_hierarchical_JSON())
        hierarchical_results = self.merge_dict(hierarchical_results, hier_result)

    data['flat_results'] = flat_results
    data['hierarchical_results'] = hierarchical_results

    for k in sorted(self.clash_dict.keys()):
      summary_results[k] = {"clashscore": self.clash_dict[k],
                            "num_clashes": len(self.list_dict[k])}
    data['summary_results'] = summary_results
    return json.dumps(data, indent=2)

  def as_coot_data(self):
    data = []
    for result in self.results :
      if result.is_outlier():
        data.append((result.atoms_info[0].id_str(),
          result.atoms_info[1].id_str(), result.overlap, result.xyz))
    return data

class probe_line_info(object): # this is parent
  def __init__(self, line, model_id=""):
    self.overlap_value = None
    self.model_id = model_id

  def is_similar(self, other):
    assert type(self) is type(other)
    return (self.srcAtom == other.srcAtom and self.targAtom == other.targAtom)

  def as_clash_obj(self, use_segids):
    assert self.overlap_value is not None
    atom1 = self.srcAtom
    atom2 = self.targAtom
    if (self.srcAtom < self.targAtom):
      atoms = [ atom1, atom2 ]
    else:
      atoms = [ atom2, atom1 ]
    clash_obj = clash(
      atoms_info=atoms,
      overlap=self.overlap_value,
      probe_type=self.type,
      outlier=self.overlap_value <= -0.4,
      max_b_factor=max(self.sBval, self.tBval),
      xyz=(self.x,self.y,self.z))
    return clash_obj

def encode_atom_info(json_dict, model_id):
  return atom_info(
    model_id=model_id,
    chain_id=json_dict['chainID'],
    resseq=json_dict['resID'],
    icode=json_dict['iCode'],
    resname=json_dict['resName'],
    altloc=json_dict['alt'],
    name=json_dict['atomName']
    )

def atom_info_to_string(info):
  return "%2s%4s%1s%3s%4s%1s" % (
    info.chain_id,
    info.resseq,
    info.icode,
    info.resname,
    info.name,
    info.altloc)

class json_probe_line_info(probe_line_info):
  def __init__(self, line, model_id=""):
    super(json_probe_line_info, self).__init__(line, model_id)
    self.name = line['master']
    self.pat = line['group']
    self.type = line['type']
    self.stype = line['srcClass']
    self.ttype = line['targetClass']
    self.sBval = line['srcBFactor']
    self.tBval = line['targetBFactor']
    self.srcAtom = encode_atom_info(line['src'], model_id)
    if 'target' in line:
      self.targAtom = encode_atom_info(line['target'], model_id)
      self.min_gap = line['gap']
      self.x = line['loc'][0]
      self.y = line['loc'][1]
      self.z = line['loc'][2]
      if 'dotCount' in line:    # Present in condensed output
        self.dotCount = line['dotCount']
        # Replace the values with condensed output values
        self.gap = self.min_gap
      else:      # Not present in condensed output but present in full output
        self.gap = line['dotGap']
        self.spike = line['spike']
        self.spikeLen = line['spikeLen']
        self.score = line['scoreOverDensity']
      self.overlap_value = self.gap

class probe_clashscore_manager(object):
  def __init__(self,
               fast = False,
               condensed_probe=False,
               nuclear=False,
               largest_occupancy=10,
               b_factor_cutoff=None,
               use_segids=False,
               verbose=False,
               model_id=""):
    """
    Calculate probe (MolProbity) clashscore

    Args:
      nuclear (bool): When True use nuclear cloud x-H distances and vdW radii,
        otherwise use electron cloud x-H distances and vdW radii
      largest_occupancy (int)
      b_factor_cutoff (float)
      use_segids (bool)
      verbose (bool): verbosity of printout
      model_id (str): model ID number, used in json output
    """
    if fast and not condensed_probe:
      raise Sorry("Incompatible parameters: fast=True and condensed=False:\n"+\
        "There's no way to work fast without using condensed output.")

    self.b_factor_cutoff = b_factor_cutoff
    self.fast = fast
    self.condensed_probe = condensed_probe
    self.use_segids=use_segids
    self.model_id = model_id
    self.occupancy_frac = 0.1
    if largest_occupancy / 100 < self.occupancy_frac:
      self.occupancy_frac = largest_occupancy / 100
    self.nuclear = nuclear

  def put_group_into_dict(self, line_info, clash_hash, hbond_hash):
    srcString = atom_info_to_string(line_info.srcAtom)
    targString = atom_info_to_string(line_info.targAtom)
    key = targString+srcString
    if (srcString < targString):
      key = srcString+targString
    if line_info.type == "bo":
      if (line_info.overlap_value <= -0.4):
        if (key in clash_hash):
          if (line_info.overlap_value < clash_hash[key].overlap_value):
            clash_hash[key] = line_info
        else :
          clash_hash[key] = line_info
    elif (line_info.type == "hb"):
      if self.condensed_probe:
        hbond_hash[key] = line_info
      else: # not condensed
        if (key in hbond_hash):
          if (line_info.gap < hbond_hash[key].gap):
            hbond_hash[key] = line_info
        else :
          hbond_hash[key] = line_info

  def filter_dicts(self, new_clash_hash, new_hbond_hash):
    temp = []
    for k,v in six.iteritems(new_clash_hash):
      if k not in new_hbond_hash:
        temp.append(v.as_clash_obj(self.use_segids))
    return temp

  def process_json_probe_output(self, probe_json):
    new_clash_hash = {}
    new_hbond_hash = {}
    # Parse the json string into a Python object
    data = json.loads(probe_json)["flat_results"]
    if self.condensed_probe:
      for line in data:
        try:
          line_storage = json_probe_line_info(line, self.model_id)
        except KeyboardInterrupt: raise
        except ValueError:
          continue # something else (different from expected) got into output
        self.put_group_into_dict(line_storage, new_clash_hash, new_hbond_hash)
    else: # not condensed
      previous_line = None
      for line in data:
        processed=False
        try:
          line_storage = json_probe_line_info(line, self.model_id)
        except KeyboardInterrupt: raise
        except ValueError:
          continue # something else (different from expected) got into output

        if previous_line is not None:
          if line_storage.is_similar(previous_line):
            # modify previous line to store this one if needed
            previous_line.overlap_value = min(previous_line.overlap_value, line_storage.overlap_value)
          else:
            # seems like new group of lines, then dump previous and start new
            # one
            self.put_group_into_dict(previous_line, new_clash_hash, new_hbond_hash)
            previous_line = line_storage
        else:
          previous_line = line_storage
      if previous_line is not None:
        self.put_group_into_dict(previous_line, new_clash_hash, new_hbond_hash)
    return self.filter_dicts(new_clash_hash, new_hbond_hash)

  def get_condensed_clashes(self, probe_json):
    clashes = set() # [(src, targ), (src, targ)]
    hbonds = [] # (src, targ), (targ, src)
    lines = json.loads(probe_json)["flat_results"]
    for l in lines:
      info = json_probe_line_info(l, self.model_id)
      if info.type == 'bo':
        gap = info.gap
        if gap <= -0.4:
          # print "good gap, saving", gap
          srcString = atom_info_to_string(info.srcAtom)
          targString = atom_info_to_string(info.targAtom)
          if (srcString, targString) not in clashes and (targString, srcString) not in clashes:
            clashes.add((srcString, targString))
            # print (srcAtom, targAtom)
      elif info.type == 'hb':
        srcString = atom_info_to_string(info.srcAtom)
        targString = atom_info_to_string(info.targAtom)
        hbonds.append((srcString, targString))
        hbonds.append((targString, srcString))
        prev_line = l
    hbonds_set = set(hbonds)
    n_clashes = 0
    # print "clashes", len(clashes)
    # print "hbonds", len(hbonds)
    for clash in clashes:
      if clash not in hbonds_set:
        n_clashes += 1
      # else:
        # print "skipping", clash
    return n_clashes

  # We have to take both the original data manager, which is from the model
  # without hydrogens, and the hydrogenated modified model because we need one
  # to construct a Probe2 program and the other to replace its model to run on.
  def run_probe_clashscore(self, data_manager, hydrogenated_model):
    self.n_clashes = 0
    self.n_clashes_b_cutoff = 0
    self.clashscore_b_cutoff = None
    self.bad_clashes = []
    self.clashscore = None
    self.n_atoms = 0
    self.natoms_b_cutoff = 0

    # Construct override parameters and then run probe2 using them and delete the resulting
    # temporary file.
    tempName = tempfile.mktemp()
    parser = iotbx.cli_parser.CCTBXParser(program_class=probe2.Program, logger=null_out())
    args = [
      "source_selection='(occupancy > {}) and not water'".format(self.occupancy_frac),
      "target_selection='occupancy > {}'".format(self.occupancy_frac),
      "use_neutron_distances={}".format(self.nuclear),
      "approach=once",
      "output.filename='{}'".format(tempName),
      "output.format=json",
      "output.condensed={}".format(self.condensed_probe),
      "output.report_vdws=False",
      "ignore_lack_of_explicit_hydrogens=True",
    ]
    parser.parse_args(args)
    p2 = probe2.Program(data_manager, parser.working_phil.extract(),
                       master_phil=parser.master_phil, logger=null_out())
    p2.overrideModel(hydrogenated_model)
    dots, output = p2.run()
    probe_json = output
    os.unlink(tempName)

    # Debugging facility, do not remove!
    # import random
    # pdb_string = hydrogenated_model.get_hierarchy().as_pdb_string()
    # tempdir = "tmp_for_probe_debug_%d" % random.randint(1000,9999)
    # while os.path.isdir(tempdir):
    #   tempdir = "tmp_for_probe_debug_%d" % random.randint(1000,9999)
    # os.mkdir(tempdir)
    # print ("Dumping info to %s" % tempdir)
    # with open(tempdir + os.sep + 'model.pdb', 'w') as f:
    #   f.write(pdb_string)
    # with open(tempdir + os.sep + 'probe_out.txt', 'w') as f:
    #   f.write('\n'.join(probe_json))

    if not self.fast:
      temp = self.process_json_probe_output(probe_json)
      self.n_clashes = len(temp)

      # If we're not running fast, call probe2 again to get non-condensed output
      # including VDW contacts.  We make another temporary file for the output and
      # then delete it.  This option is for printing to file for coot usage.  The
      # no-VDWOUT option is used to speed up the parsing of the output.
      tempName = tempfile.mktemp()
      parser = iotbx.cli_parser.CCTBXParser(program_class=probe2.Program, logger=null_out())
      args = [
        "source_selection='(occupancy > {}) and not water'".format(self.occupancy_frac),
        "target_selection='occupancy > {}'".format(self.occupancy_frac),
        "use_neutron_distances={}".format(self.nuclear),
        "approach=once",
        "output.filename='{}'".format(tempName),
        "output.format=json",
        "ignore_lack_of_explicit_hydrogens=True",
      ]
      parser.parse_args(args)
      p2 = probe2.Program(data_manager, parser.working_phil.extract(),
                         master_phil=parser.master_phil, logger=null_out())
      p2.overrideModel(hydrogenated_model)
      dots, output = p2.run()
      self.probe_json = output
      os.unlink(tempName)
    else:
      self.n_clashes = self.get_condensed_clashes(probe_json)

    # Find the number of non-water atoms that pass the occupancy test
    # and then use it to compute the clashscore.
    self.n_atoms = 0
    for a in hydrogenated_model.get_hierarchy().atoms():
      isWater = iotbx.pdb.common_residue_names_get_class(name=a.parent().resname) == "common_water"
      if a.occ > self.occupancy_frac and not isWater:
        self.n_atoms += 1
    if self.n_atoms == 0:
      self.clashscore = 0.0
    else:
      self.clashscore = (self.n_clashes * 1000) / self.n_atoms

    if not self.fast:
      # The rest is not necessary, we already got clashscore
      if self.b_factor_cutoff is not None:
        clashes_b_cutoff = 0
        for clash_obj in temp:
          if clash_obj.max_b_factor < self.b_factor_cutoff:
            clashes_b_cutoff += 1
        self.n_clashes_b_cutoff = clashes_b_cutoff
      used = []

      for clash_obj in sorted(temp):
        test_key = clash_obj.id_str_no_atom_name()
        test_key = clash_obj.id_str()
        if test_key not in used:
          used.append(test_key)
          self.bad_clashes.append(clash_obj)

      if self.b_factor_cutoff is not None:
        # Find the number of non-water atoms who pass the B factor test and the occupancy test
        self.natoms_b_cutoff = 0
        for a in hydrogenated_model.get_hierarchy().atoms():
          isWater = iotbx.pdb.common_residue_names_get_class(name=a.parent().resname) == "common_water"
          if  (a.b < self.b_factor_cutoff) and (a.occ > self.occupancy_frac) and not isWater:
            self.natoms_b_cutoff += 1
      self.clashscore_b_cutoff = None
      if self.natoms_b_cutoff == 0:
        self.clashscore_b_cutoff = 0.0
      else :
        self.clashscore_b_cutoff = \
          (self.n_clashes_b_cutoff*1000) / self.natoms_b_cutoff

def check_and_add_hydrogen(
        probe_parameters=None,
        data_manager_model=None,
        nuclear=False,
        keep_hydrogens=True,
        verbose=False,
        n_hydrogen_cut_off=0,
        do_flips=False,
        log=None):
  """
  If no hydrogens present, force addition for clashscore calculation.
  Use REDUCE to add the hydrogen atoms.

  Args:
    data_manager_model : Model from the data_manager
    nuclear (bool): When True use nuclear cloud x-H distances and vdW radii,
      otherwise use electron cloud x-H distances and vdW radii
    keep_hydrogens (bool): when True, if there are hydrogen atoms, keep them
    verbose (bool): verbosity of printout
    n_hydrogen_cut_off (int): when number of hydrogen atoms < n_hydrogen_cut_off
      force keep_hydrogens tp True

  Returns:
    (model): Model with hydrogens added
    (bool): True when the model was modified/replaced
  """
  if not log: log = sys.stdout
  assert probe_parameters
  assert data_manager_model
  if keep_hydrogens:
    elements = data_manager_model.get_hierarchy().atoms().extract_element()
    # strangely the elements can have a space when coming from phenix.clashscore
    # but no space when coming from phenix.molprobity
    h_count = elements.count('H')
    if h_count <= n_hydrogen_cut_off: h_count += elements.count(' H')
    if h_count <= n_hydrogen_cut_off: h_count += elements.count('D')
    if h_count <= n_hydrogen_cut_off: h_count += elements.count(' D')
    if h_count > n_hydrogen_cut_off:
      has_hd = True
    else:
      has_hd = False
    if not has_hd:
      if verbose:
        print("\nNo H/D atoms detected - forcing hydrogen addition!\n", file=log)
      keep_hydrogens = False

  # add hydrogen if needed
  if not keep_hydrogens:
    # Remove hydrogens and add them back in
    if verbose:
      print("\nTrimming and adding hydrogens...\n")
    reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
      model = data_manager_model,
      use_neutron_distances=nuclear,
      n_terminal_charge="residue_one",
      exclude_water=True,
      stop_for_unknowns=False,
      keep_existing_H=False
    )
    reduce_add_h_obj.run()
    reduce_add_h_obj.show(log)
    missed_residues = set(reduce_add_h_obj.no_H_placed_mlq)
    if len(missed_residues) > 0:
      bad = ""
      for res in missed_residues:
        bad += " " + res
      raise Sorry("Restraints were not found for the following residues:"+bad)
    data_manager_model = reduce_add_h_obj.get_model()

    # Optimize H atoms with mmtbx.reduce
    if verbose:
      print("\nOptimizing H atoms with mmtbx.reduce2:Optimizer...\n")
    opt = Optimizers.Optimizer(probe_parameters, do_flips, data_manager_model, modelIndex=None,
      fillAtomDump = False)

    # Re-process the model because we have removed some atoms that were previously
    # bonded.  Don't make restraints during the reprocessing.
    # We had to do this to keep from crashing on a call to pair_proxies when generating
    # mmCIF files, so we always do it for safety.
    data_manager_model.get_hierarchy().sort_atoms_in_place()
    data_manager_model.get_hierarchy().atoms().reset_serial()
    p = reduce_hydrogen.get_reduce_pdb_interpretation_params(nuclear)
    # We need to turn this on because without it 1zz0.txt kept flipping the ring
    # in A TYR 214 every time we re-interpreted. The original interpretation done
    # by Hydrogen placement will have flipped them, so we don't need to do it again.
    p.pdb_interpretation.flip_symmetric_amino_acids=False
    p.pdb_interpretation.disable_uc_volume_vs_n_atoms_check=True
    p.pdb_interpretation.allow_polymer_cross_special_position=True
    p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
    p.pdb_interpretation.proceed_with_excessive_length_bonds=True
    #p.pdb_interpretation.sort_atoms=True
    data_manager_model.process(make_restraints=False, pdb_interpretation_params=p)

    return data_manager_model, True
  else:
    if verbose:
      print("\nUsing input model H/D atoms...\n")
    return data_manager_model, False

  def show(self, out=sys.stdout, prefix=""):
    if (self.n_outliers == 0):
      print(prefix+"No backwards Asn/Gln/His sidechains found.", file=out)
    else :
      for flip in self.results :
        print(prefix+flip.as_string(), file=out)
