
"""
All-atom contact analysis.  Requires Reduce and Probe (installed separately).
"""

from __future__ import absolute_import, division, print_function
from mmtbx.validation import validation, atoms, atom_info, residue
from mmtbx.utils import run_reduce_with_timeout
from libtbx.math_utils import cmp
from libtbx.utils import Sorry
from libtbx import easy_run
import libtbx.load_env
import iotbx.pdb
import os
import re
import sys
import six
import json

class clash(atoms):
  __clash_attr__ = [
    "overlap",
    "probe_type",
    "max_b_factor",
  ]
  __slots__ = atoms.__slots__ + __clash_attr__

  @staticmethod
  def header():
    return "%-20s %-20s  %7s" % ("Atom 1", "Atom 2", "Overlap")

  def id_str(self, spacer=" "):
    return "%s%s%s" % (self.atoms_info[0].id_str(), spacer,
      self.atoms_info[1].id_str())

  def id_str_no_atom_name(self):
    return "%s %s" % (self.atoms_info[0].id_str()[0:11],
      self.atoms_info[1].id_str()[0:11])

  def id_str_src_atom_no_atom_name(self):
    return self.atoms_info[0].id_str()[0:11]

  def format_old(self):
    return "%s :%.3f" % (self.id_str(), abs(self.overlap))

  def as_JSON(self):
    atom0_slots_list = [s for s in self.atoms_info[0].__slots__]
    atom0_slots_as_dict = ({s: getattr(self.atoms_info[0], s) for s in atom0_slots_list if s != 'xyz' })
    atom1_slots_list = [s for s in self.atoms_info[1].__slots__]
    atom1_slots_as_dict = ({"target_"+s: getattr(self.atoms_info[1], s) for s in atom1_slots_list if s != 'xyz' })
    atoms_dict = self.merge_two_dicts(atom0_slots_as_dict, atom1_slots_as_dict)
    #atom0_slots_as_dict["resid"] = atom0_slots_as_dict['resseq']+atom0_slots_as_dict['icode']
    #print("atom0_dict: " + str(atom0_slots_as_dict))
    serializable_slots = [s for s in self.__slots__ if s != 'atoms_info' and hasattr(self, s) ]
    slots_as_dict = ({s: getattr(self, s) for s in serializable_slots})
    #print("slots_dict: " + str(slots_as_dict))
    #print("combo:")
    #print({**slots_as_dict, **atom0_slots_as_dict})
    return json.dumps(self.merge_two_dicts(slots_as_dict, atoms_dict), indent=2)
    #return json.dumps({**slots_as_dict, **atom0_slots_as_dict, **atom1_slots_as_dict}, indent=2)

  def as_hierarchical_JSON(self):
    hierarchical_dict = {}
    hierarchy_nest_list = ['model_id', 'chain_id', 'resid', 'altloc']
    return json.dumps(self.nest_dict(hierarchy_nest_list, hierarchical_dict), indent=2)

  def as_string(self):
    return "%-20s %-20s  %7.3f" % (self.atoms_info[0].id_str(),
      self.atoms_info[1].id_str(), abs(self.overlap))

  def as_table_row_phenix(self):
    return [ self.atoms_info[0].id_str(), self.atoms_info[1].id_str(),
             abs(self.overlap) ]

  def __cmp__(self, other) : # sort in descending order
    return cmp(self.overlap, other.overlap)

  def __eq__(self, other):
    return self.overlap == other.overlap

  def __ne__(self, other):
    return self.overlap != other.overlap

  def __lt__(self, other):
    return self.overlap < other.overlap

  def __le__(self, other):
    return self.overlap <= other.overlap

  def __gt__ (self, other):
    return self.overlap > other.overlap

  def __ge__(self, other):
    return self.overlap >= other.overlap

class clashscore(validation):
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
      pdb_hierarchy,
      fast = False, # do really fast clashscore, produce only the number
      condensed_probe = False, # Use -CON for probe. Reduces output 10x.
      keep_hydrogens=True,
      nuclear=False,
      force_unique_chain_ids=False,
      time_limit=120,
      b_factor_cutoff=None,
      save_modified_hierarchy=False,
      verbose=False,
      do_flips=False,
      out=sys.stdout):
    if (not pdb_hierarchy.fits_in_pdb_format()):
      from iotbx.pdb.forward_compatible_pdb_cif_conversion \
        import forward_compatible_pdb_cif_conversion
      conversion_info = forward_compatible_pdb_cif_conversion(
        hierarchy = pdb_hierarchy)
      conversion_info.\
       convert_hierarchy_to_forward_compatible_pdb_representation(pdb_hierarchy)
      if verbose:
        print(
        "Converted model to forward_compatible PDB for clashscore",file = out)
    else:
      conversion_info = None
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
    if (not libtbx.env.has_module(name="probe")):
      raise RuntimeError(
        "Probe could not be detected on your system.  Please make sure "+
        "Probe is in your path.\nProbe is available at "+
        "http://kinemage.biochem.duke.edu/")
    if verbose:
      if not nuclear:
        print("\nUsing electron cloud x-H distances and vdW radii")
      else:
        print("\nUsing nuclear cloud x-H distances and vdW radii")
    import iotbx.pdb
    from scitbx.array_family import flex
    from mmtbx.validation import utils
    n_models = len(pdb_hierarchy.models())
    use_segids = utils.use_segids_in_place_of_chainids(
                   hierarchy=pdb_hierarchy)
    for i_mod, model in enumerate(pdb_hierarchy.models()):
      input_str,_ = check_and_add_hydrogen(
        pdb_hierarchy=pdb_hierarchy,
        model_number=i_mod,
        nuclear=nuclear,
        verbose=verbose,
        time_limit=time_limit,
        keep_hydrogens=keep_hydrogens,
        do_flips = do_flips,
        log=out)
      r = iotbx.pdb.hierarchy.root()
      mdc = model.detached_copy()
      r.append_model(mdc)
      occ_max = flex.max(r.atoms().extract_occ())
      self.probe_clashscore_manager = probe_clashscore_manager(
        h_pdb_string=input_str,
        nuclear=nuclear,
        fast=self.fast,
        condensed_probe=self.condensed_probe,
        largest_occupancy=occ_max,
        b_factor_cutoff=b_factor_cutoff,
        use_segids=use_segids,
        verbose=verbose,
        model_id=model.id)
      self.probe_clashscore_manager.run_probe_clashscore(input_str)
      if (save_modified_hierarchy):
        self.pdb_hierarchy = iotbx.pdb.input(
          pdb_string=self.probe_clashscore_manager.h_pdb_string).construct_hierarchy()
        if conversion_info:
          conversion_info.convert_hierarchy_to_full_representation(
             self.pdb_hierarchy)

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

    if conversion_info:
      if verbose:
        print("Converted model back to full representation", file = out)
      conversion_info.convert_hierarchy_to_full_representation(pdb_hierarchy)

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
    atom1 = decode_atom_string(self.srcAtom,  use_segids, self.model_id)
    atom2 = decode_atom_string(self.targAtom, use_segids, self.model_id)
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

class condensed_probe_line_info(probe_line_info):
  def __init__(self, line, model_id=""):
    super(condensed_probe_line_info, self).__init__(line, model_id)
    # What is in line:
    # name:pat:type:srcAtom:targAtom:dot-count:mingap:gap:spX:spY:spZ:spikeLen:score:stype:ttype:x:y:z:sBval:tBval:
    sp = line.split(":")
    self.type = sp[2]
    self.srcAtom = sp[3]
    self.targAtom = sp[4]
    self.min_gap = float(sp[6])
    self.gap = float(sp[7])
    self.x = float(sp[-5])
    self.y = float(sp[-4])
    self.z = float(sp[-3])
    self.sBval = float(sp[-2])
    self.tBval = float(sp[-1])
    self.overlap_value = self.gap

class raw_probe_line_info(probe_line_info):
  def __init__(self, line, model_id=""):
    super(raw_probe_line_info, self).__init__(line, model_id)
    self.name, self.pat, self.type, self.srcAtom, self.targAtom, self.min_gap, \
    self.gap, self.kissEdge2BullsEye, self.dot2BE, self.dot2SC, self.spike, \
    self.score, self.stype, self.ttype, self.x, self.y, self.z, self.sBval, \
    self.tBval = line.split(":")
    self.gap = float(self.gap)
    self.x = float(self.x)
    self.y = float(self.y)
    self.z = float(self.z)
    self.sBval = float(self.sBval)
    self.tBval = float(self.tBval)
    self.overlap_value = self.gap


class probe_clashscore_manager(object):
  def __init__(self,
               h_pdb_string,
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
      h_pdb_string (str): PDB string that contains hydrogen atoms
      nuclear (bool): When True use nuclear cloud x-H distances and vdW radii,
        otherwise use electron cloud x-H distances and vdW radii
      largest_occupancy (int)
      b_factor_cutoff (float)
      use_segids (bool)
      verbose (bool): verbosity of printout
      model_id (str): model ID number, used in json output
    """
    assert libtbx.env.has_module(name="probe")
    if fast and not condensed_probe:
      raise Sorry("Incompatible parameters: fast=True and condensed=False:\n"+\
        "There's no way to work fast without using condensed output.")

    self.b_factor_cutoff = b_factor_cutoff
    self.fast = fast
    self.condensed_probe = condensed_probe
    self.use_segids=use_segids
    self.model_id = model_id
    ogt = 10
    blt = self.b_factor_cutoff
    if largest_occupancy < ogt:
      ogt = largest_occupancy

    self.probe_atom_b_factor = None
    probe_command = libtbx.env.under_build(os.path.join('probe', 'exe', 'probe'))
    if os.getenv('CCP4'):
      ccp4_probe = os.path.join(os.environ['CCP4'],'bin','probe')
      if (os.path.isfile(ccp4_probe) and not os.path.isfile(probe_command)):
        probe_command = ccp4_probe
    probe_command = '"%s"' % probe_command   # in case of spaces in path
    self.probe_command = probe_command
    nuclear_flag = ""
    condensed_flag = ""
    if nuclear:
      nuclear_flag = "-nuclear"
    if self.condensed_probe:
      condensed_flag = "-CON"
    self.probe_txt = \
      '%s -u -q -mc -het -once -NOVDWOUT %s %s' % (probe_command, condensed_flag, nuclear_flag) +\
        ' "ogt%d not water" "ogt%d" -' % (ogt, ogt)
    #The -NOVDWOUT probe run above is faster for clashscore to parse,
    # the full_probe_txt version below is for printing to file for coot usage
    self.full_probe_txt = \
      '%s -u -q -mc -het -once %s' % (probe_command, nuclear_flag) +\
        ' "ogt%d not water" "ogt%d" -' % (ogt, ogt)
    self.probe_atom_txt = \
      '%s -q -mc -het -dumpatominfo %s' % (probe_command, nuclear_flag) +\
        ' "ogt%d not water" -' % ogt
    if blt is not None:
      self.probe_atom_b_factor = \
        '%s -q -mc -het -dumpatominfo %s' % (probe_command, nuclear_flag) +\
          ' "blt%d ogt%d not water" -' % (blt, ogt)

    self.h_pdb_string = h_pdb_string
    #self.run_probe_clashscore(self.h_pdb_string)

  def put_group_into_dict(self, line_info, clash_hash, hbond_hash):
    key = line_info.targAtom+line_info.srcAtom
    if (line_info.srcAtom < line_info.targAtom):
      key = line_info.srcAtom+line_info.targAtom
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

  def process_raw_probe_output(self, probe_unformatted):
    new_clash_hash = {}
    new_hbond_hash = {}
    if self.condensed_probe:
      for line in probe_unformatted:
        try:
          line_storage = condensed_probe_line_info(line, self.model_id)
        except KeyboardInterrupt: raise
        except ValueError:
          continue # something else (different from expected) got into output
        self.put_group_into_dict(line_storage, new_clash_hash, new_hbond_hash)
    else: # not condensed
      previous_line = None
      for line in probe_unformatted:
        processed=False
        try:
          line_storage = raw_probe_line_info(line, self.model_id)
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

  def get_condensed_clashes(self, lines):
    # Standalone faster parsing of output when only clashscore is needed.
    def parse_line(line):
      sp = line.split(':')
      return sp[3], sp[4], float(sp[7])
    def parse_h_line(line):
      sp = line.split(':')
      return sp[3], sp[4]

    clashes = set() # [(src, targ), (src, targ)]
    hbonds = [] # (src, targ), (targ, src)
    for l in lines:
      rtype = l[6:8]
      if rtype == 'bo':
        srcAtom, targAtom, gap = parse_line(l)
        if gap <= -0.4:
          # print l[:43], "good gap, saving", gap
          if (srcAtom, targAtom) not in clashes and (targAtom, srcAtom) not in clashes:
            clashes.add((srcAtom, targAtom))
            # print (srcAtom, targAtom)
      elif rtype == 'hb':
        srcAtom, targAtom = parse_h_line(l)
        hbonds.append((srcAtom, targAtom))
        hbonds.append((targAtom, srcAtom))
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

  def run_probe_clashscore(self, pdb_string):
    self.n_clashes = 0
    self.n_clashes_b_cutoff = 0
    self.clashscore_b_cutoff = None
    self.bad_clashes = []
    self.clashscore = None
    self.n_atoms = 0
    self.natoms_b_cutoff = 0

    probe_out = easy_run.fully_buffered(self.probe_txt,
      stdin_lines=pdb_string)
    if (probe_out.return_code != 0):
      raise RuntimeError("Probe crashed - dumping stderr:\n%s" %
        "\n".join(probe_out.stderr_lines))
    probe_unformatted = probe_out.stdout_lines

    # Debugging facility, do not remove!
    # import random
    # tempdir = "tmp_for_probe_debug_%d" % random.randint(1000,9999)
    # while os.path.isdir(tempdir):
    #   tempdir = "tmp_for_probe_debug_%d" % random.randint(1000,9999)
    # os.mkdir(tempdir)
    # print ("Dumping info to %s" % tempdir)
    # with open(tempdir + os.sep + 'model.pdb', 'w') as f:
    #   f.write(pdb_string)
    # with open(tempdir + os.sep + 'probe_out.txt', 'w') as f:
    #   f.write('\n'.join(probe_unformatted))

    if not self.fast:
      temp = self.process_raw_probe_output(probe_unformatted)
      self.n_clashes = len(temp)
      # XXX Warning: one more probe call here
      printable_probe_out = easy_run.fully_buffered(self.full_probe_txt,
                                                    stdin_lines=pdb_string)
      self.probe_unformatted = "\n".join(printable_probe_out.stdout_lines)
    else:
      self.n_clashes = self.get_condensed_clashes(probe_unformatted)

    # getting number of atoms from probe
    probe_info = easy_run.fully_buffered(self.probe_atom_txt,
      stdin_lines=pdb_string) #.raise_if_errors().stdout_lines
    err = probe_info.format_errors_if_any()
    if err is not None and err.find("No atom data in input.")>-1:
      return
    #if (len(probe_info) == 0):
    #  raise RuntimeError("Empty PROBE output.")
    n_atoms = 0
    for line in probe_info.stdout_lines:
      try:
        dump, n_atoms = line.split(":")
      except KeyboardInterrupt: raise
      except ValueError:
        pass # something else (different from expected) got into output
    self.n_atoms = int(n_atoms)
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

      if self.probe_atom_b_factor is not None:
        probe_info_b_factor = easy_run.fully_buffered(self.probe_atom_b_factor,
          stdin_lines=pdb_string).raise_if_errors().stdout_lines
        for line in probe_info_b_factor :
          dump_b, natoms_b_cutoff = line.split(":")
        self.natoms_b_cutoff = int(natoms_b_cutoff)
      self.clashscore_b_cutoff = None
      if self.natoms_b_cutoff == 0:
        self.clashscore_b_cutoff = 0.0
      else :
        self.clashscore_b_cutoff = \
          (self.n_clashes_b_cutoff*1000) / self.natoms_b_cutoff

def decode_atom_string(atom_str, use_segids=False, model_id=""):
  # Example:
  # ' A  49 LEU HD11B'
  if (not use_segids) or (len(atom_str) == 16):
    return atom_info(
      model_id=model_id,
      chain_id=atom_str[0:2],
      resseq=atom_str[2:6],
      icode=atom_str[6],
      resname=atom_str[7:10],
      altloc=atom_str[15],
      name=atom_str[11:15])
  else:
    return atom_info(
      model_id=model_id,
      chain_id=atom_str[0:4],
      resseq=atom_str[4:8],
      icode=atom_str[8],
      resname=atom_str[9:12],
      altloc=atom_str[17],
      name=atom_str[13:17])

def check_and_report_reduce_failure(fb_object, input_lines, output_fname):
  if (fb_object.return_code != 0):
    with open(output_fname, 'w') as f:
      f.write(input_lines)
    msg_str = "Reduce crashed with command '%s'.\nDumping stdin to file '%s'.\n" +\
        "Return code: %d\nDumping stderr:\n%s"
    raise Sorry(msg_str % (fb_object.command, output_fname, fb_object.return_code,
                           "\n".join(fb_object.stderr_lines)))

def check_and_add_hydrogen(
        pdb_hierarchy=None,
        file_name=None,
        nuclear=False,
        keep_hydrogens=True,
        verbose=False,
        model_number=0,
        n_hydrogen_cut_off=0,
        time_limit=120,
        allow_multiple_models=True,
        crystal_symmetry=None,
        do_flips=False,
        log=None):
  """
  If no hydrogens present, force addition for clashscore calculation.
  Use REDUCE to add the hydrogen atoms.

  Args:
    pdb_hierarchy : pdb hierarchy
    file_name (str): pdb file name
    nuclear (bool): When True use nuclear cloud x-H distances and vdW radii,
      otherwise use electron cloud x-H distances and vdW radii
    keep_hydrogens (bool): when True, if there are hydrogen atoms, keep them
    verbose (bool): verbosity of printout
    model_number (int): the number of model to use
    time_limit (int): limit the time it takes to add hydrogen atoms
    n_hydrogen_cut_off (int): when number of hydrogen atoms < n_hydrogen_cut_off
      force keep_hydrogens tp True
    allow_multiple_models (bool): Allow models that contain more than one model
    crystal_symmetry : must provide crystal symmetry when using pdb_hierarchy

  Returns:
    (str): PDB string
    (bool): True when PDB string was updated
  """
  if not log: log = sys.stdout
  if file_name:
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    cryst_sym = pdb_inp.crystal_symmetry()
  elif not allow_multiple_models:
    assert crystal_symmetry
    cryst_sym = crystal_symmetry
  else:
    cryst_sym = None
  assert pdb_hierarchy
  assert model_number < len(pdb_hierarchy.models())
  models = pdb_hierarchy.models()
  if (len(models) > 1) and (not allow_multiple_models):
    raise Sorry("When using CCTBX clashscore, provide only a single model.")
  model = models[model_number]
  r = iotbx.pdb.hierarchy.root()
  mdc = model.detached_copy()
  r.append_model(mdc)
  if keep_hydrogens:
    elements = r.atoms().extract_element()
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
  import libtbx.load_env
  has_reduce = libtbx.env.has_module(name="reduce")
  # add hydrogen if needed
  if has_reduce and (not keep_hydrogens):
    # set reduce running parameters
    if verbose:
      print("\nAdding H/D atoms with reduce...\n")
    build = "-oh -his -flip -keep -allalt -limit{}"
    if not do_flips : build += " -pen9999"
    if nuclear:
      build += " -nuc -"
    else:
      build += " -"
    build = build.format(time_limit)
    trim = " -quiet -trim -"
    stdin_lines = r.as_pdb_string(cryst_sym)
    clean_out = run_reduce_with_timeout(
        parameters=trim,
        stdin_lines=stdin_lines)
    stdin_fname = "reduce_fail.pdb"
    check_and_report_reduce_failure(
        fb_object=clean_out,
        input_lines=stdin_lines,
        output_fname="reduce_fail.pdb")
    build_out = run_reduce_with_timeout(
        parameters=build,
        stdin_lines=clean_out.stdout_lines)
    check_and_report_reduce_failure(
        fb_object=build_out,
        input_lines=stdin_lines,
        output_fname="reduce_fail.pdb")
    reduce_str = '\n'.join(build_out.stdout_lines)
    return reduce_str,True
  else:
    if not has_reduce:
      msg = 'molprobity.reduce could not be detected on your system.\n'
      msg += 'Cannot add hydrogen to PDB file'
      print(msg, file=log)
    if verbose:
      print("\nUsing input model H/D atoms...\n")
    return r.as_pdb_string(cryst_sym),False

#-----------------------------------------------------------------------
# this isn't really enough code to justify a separate module...
#
class nqh_flip(residue):
  """
  Backwards Asn/Gln/His sidechain, identified by Reduce's hydrogen-bond
  network optimization.
  """
  def as_string(self):
    return self.id_str()

  def as_table_row_phenix(self):
    if self.chain_id:
      return [ self.chain_id, "%1s%s %s" % (self.altloc,self.resname,self.resid) ]
    elif self.segid:
      return [ self.segid, "%1s%s %s" % (self.altloc,self.resname,self.resid) ]
    else:
      raise Sorry("no chain_id or segid found for nqh flip table row")

  #alternate residue class methods for segid compatibility
  #more method overrides may be necessary
  #a more robust propagation of segid would preferable, eventually
  def atom_selection_string(self):
    if self.chain_id:
      return "(chain '%s' and resid '%s' and resname %s and altloc '%s')" % \
        (self.chain_id, self.resid, self.resname, self.altloc)
    elif self.segid:
      return "(segid %s and resid '%s' and resname %s and altloc '%s')" % \
          (self.segid, self.resid, self.resname, self.altloc)
    else:
      raise Sorry("no chain_id or segid found for nqh flip atom selection")

  def id_str(self, ignore_altloc=False):
    if self.chain_id:
      base = "%2s%4s%1s" % (self.chain_id, self.resseq, self.icode)
    elif self.segid:
      base = "%4s%4s%1s" % (self.segid, self.resseq, self.icode)
    else:
      raise Sorry("no chain_id or segid found for nqh flip id_str")
    if (not ignore_altloc):
      base += "%1s" % self.altloc
    else :
      base += " "
    base += "%3s" % self.resname
    if (self.segid is not None):
      base += " segid='%4s'" % self.segid
    return base
  #end segid compatibility

class nqh_flips(validation):
  """
  N/Q/H sidechain flips identified by Reduce.
  """
  gui_list_headers = ["Chain", "Residue"]
  gui_formats = ["%s", "%s"]
  wx_column_widths = [75,220]
  def __init__(self, pdb_hierarchy):
    re_flip = re.compile(":FLIP")
    validation.__init__(self)
    in_lines = pdb_hierarchy.as_pdb_string()
    reduce_out = run_reduce_with_timeout(
        parameters="-BUILD -",
        stdin_lines=in_lines)
    check_and_report_reduce_failure(
        fb_object=reduce_out,
        input_lines=in_lines,
        output_fname="reduce_fail.pdb")
    from mmtbx.validation import utils
    use_segids = utils.use_segids_in_place_of_chainids(
      hierarchy=pdb_hierarchy)
    for line in reduce_out.stdout_lines:
    #chain format (2-char chain)
    #USER  MOD Set 1.1: B  49 GLN     :FLIP  amide:sc=    -2.7! C(o=-5.8!,f=-1.3!)
    #segid format (4-char segid)
    #USER  MOD Set 1.1:B     49 GLN     :FLIP  amide:sc=    -2.7! C(o=-5.8!,f=-1.3!)
      if re_flip.search(line):
        resid = line.split(":")[1]
        #reduce has slightly different outputs using chains versus segid
        if len(resid) == 15: #chain
          chain_id = resid[0:2].strip()
          segid = None
          if (len(chain_id) == 0):
            chain_id = ' '
          resid_less_chain = resid[2:]
        elif len(resid) == 17: #segid
          #self.results = []
          #return
          chain_id = None
          segid = resid[0:4].strip()
          #chain_id = resid[0:4].strip()
          resid_less_chain = resid[4:]
        else:
          raise Sorry("unexpected length of residue identifier in reduce USER MODs.")
        resname = resid_less_chain[5:8]

        assert (resname in ["ASN", "GLN", "HIS"])
        flip = nqh_flip(
          chain_id=chain_id,
          segid=segid,
          resseq=resid_less_chain[0:4].strip(),
          icode= resid_less_chain[4:5],
          altloc=resid_less_chain[12:13],
          resname=resname,
          outlier=True)
        flip.set_coordinates_from_hierarchy(pdb_hierarchy)
        self.results.append(flip)
        self.n_outliers += 1

  def show(self, out=sys.stdout, prefix=""):
    if (self.n_outliers == 0):
      print(prefix+"No backwards Asn/Gln/His sidechains found.", file=out)
    else :
      for flip in self.results :
        print(prefix+flip.as_string(), file=out)
