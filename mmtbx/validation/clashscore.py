
"""
All-atom contact analysis.  Requires Reduce and Probe (installed separately).
"""

from __future__ import division
from mmtbx.validation import validation, atoms, atom_info, residue
from libtbx.utils import Sorry
from libtbx import easy_run
import libtbx.load_env
import iotbx.pdb
import re
import sys

class clash (atoms) :
  __clash_attr__ = [
    "overlap",
    "probe_type",
    "max_b_factor",
  ]
  __slots__ = atoms.__slots__ + __clash_attr__

  @staticmethod
  def header () :
    return "%-20s %-20s  %7s" % ("Atom 1", "Atom 2", "Overlap")

  def id_str (self, spacer=" ") :
    return "%s%s%s" % (self.atoms_info[0].id_str(), spacer,
      self.atoms_info[1].id_str())

  def id_str_no_atom_name (self) :
    return "%s %s" % (self.atoms_info[0].id_str()[0:11],
      self.atoms_info[1].id_str()[0:11])

  def format_old (self) :
    return "%s :%.3f" % (self.id_str(), abs(self.overlap))

  def as_string (self) :
    return "%-20s %-20s  %7.3f" % (self.atoms_info[0].id_str(),
      self.atoms_info[1].id_str(), abs(self.overlap))

  def as_table_row_phenix (self) :
    return [ self.atoms_info[0].id_str(), self.atoms_info[1].id_str(),
             abs(self.overlap) ]

  def __cmp__ (self, other) : # sort in descending order
    return cmp(self.overlap, other.overlap)

class clashscore(validation):
  __slots__ = validation.__slots__ + [
    "clashscore",
    "clashscore_b_cutoff",
    "clash_dict",
    "clash_dict_b_cutoff",
    "list_dict",
    "b_factor_cutoff",
    "probe_file",
    "probe_clashscore_manager"
  ]
  program_description = "Analyze clashscore for protein model"
  gui_list_headers = ["Atom 1", "Atom 2", "Overlap"]
  gui_formats = ["%s", "%s", ".3f"]
  wx_column_widths = [150, 150, 150] #actually set in GUI's Molprobity/Core.py

  def get_result_class (self) : return clash

  def __init__ (self,
      pdb_hierarchy,
      keep_hydrogens=True,
      nuclear=False,
      force_unique_chain_ids=False,
      time_limit=120,
      b_factor_cutoff=None,
      save_modified_hierarchy=False,
      verbose=False,
      do_flips=False,
      out=sys.stdout) :
    validation.__init__(self)
    self.b_factor_cutoff = b_factor_cutoff
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
        print "\nUsing electron cloud x-H distances and vdW radii"
      else:
        print "\nUsing nuclear cloud x-H distances and vdW radii"
    import iotbx.pdb.hierarchy
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
        largest_occupancy=occ_max,
        b_factor_cutoff=b_factor_cutoff,
        use_segids=use_segids,
        verbose=verbose)
      if (save_modified_hierarchy) :
        self.pdb_hierarchy = iotbx.pdb.hierarchy.input(
          pdb_string=self.probe_clashscore_manager.h_pdb_string).hierarchy
      self.clash_dict[model.id] = self.probe_clashscore_manager.clashscore
      self.clash_dict_b_cutoff[model.id] = self.probe_clashscore_manager.\
                                           clashscore_b_cutoff
      self.list_dict[model.id] = self.probe_clashscore_manager.bad_clashes
      if (n_models == 1) or (self.clashscore is None) :
        self.results = self.probe_clashscore_manager.bad_clashes
        self.n_outliers = len(self.results)
        self.clashscore = self.probe_clashscore_manager.clashscore
        self.clashscore_b_cutoff = self.probe_clashscore_manager.\
                                   clashscore_b_cutoff

  def get_clashscore(self):
    return self.clashscore

  def get_clashscore_b_cutoff(self):
    return self.clashscore_b_cutoff

  def show_old_output (self, out=sys.stdout, verbose=False) :
    self.print_clashlist_old(out=out)
    self.show_summary(out=out)

  def show_summary (self, out=sys.stdout, prefix="") :
    if self.clashscore is None:
      raise Sorry("PROBE output is empty. Model is not compatible with PROBE.")
    elif (len(self.clash_dict) == 1) :
      print >> out, prefix + "clashscore = %.2f" % self.clash_dict['']
      if self.clash_dict_b_cutoff[''] is not None:
        print >> out, "clashscore (B factor cutoff = %d) = %f" % \
          (self.b_factor_cutoff,
           self.clash_dict_b_cutoff[''])
    else:
      for k in sorted(self.clash_dict.keys()) :
        print >> out, prefix + "MODEL %s clashscore = %.2f" % (k,
          self.clash_dict[k])
        if self.clash_dict_b_cutoff[k] is not None:
          print >> out, "MODEL%s clashscore (B factor cutoff = %d) = %f" % \
            (k, self.b_factor_cutoff, self.clash_dict_b_cutoff[k])

  def print_clashlist_old (self, out=sys.stdout):
    for k in self.list_dict.keys():
      if k == '':
        print >> out, "Bad Clashes >= 0.4 Angstrom:"
        for result in self.list_dict[k] :
          print >> out, result.format_old()
      else:
        print >> out, "Bad Clashes >= 0.4 Angstrom MODEL%s" % k
        for result in self.list_dict[k] :
          print >> out, result.format_old()

  def show (self, out=sys.stdout, prefix="", outliers_only=None, verbose=None) :
    if (len(self.clash_dict) == 1) :
      for result in self.list_dict[''] :
        print >> out, prefix + str(result)
    else :
      for k in self.list_dict.keys():
        for result in self.list_dict[k] :
          print >> out, prefix + str(result)
    self.show_summary(out=out, prefix=prefix)

  def as_coot_data (self) :
    data = []
    for result in self.results :
      if result.is_outlier() :
        data.append((result.atoms_info[0].id_str(),
          result.atoms_info[1].id_str(), result.overlap, result.xyz))
    return data

class probe_line_info_storage(object):
  def __init__(self, line):
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
  def is_similar(self, other):
    assert isinstance(other, probe_line_info_storage)
    return (self.srcAtom == other.srcAtom and self.targAtom == other.targAtom)
  def as_clash_obj(self, use_segids):
    atom1 = decode_atom_string(self.srcAtom,  use_segids)
    atom2 = decode_atom_string(self.targAtom, use_segids)
    if (cmp(self.srcAtom, self.targAtom) < 0):
      atoms = [ atom1, atom2 ]
    else:
      atoms = [ atom2, atom1 ]
    clash_obj = clash(
      atoms_info=atoms,
      overlap=self.gap,
      probe_type=self.type,
      outlier=abs(self.gap) > 0.4,
      max_b_factor=max(self.sBval, self.tBval),
      xyz=(self.x,self.y,self.z))
    return clash_obj

class probe_clashscore_manager(object):
  def __init__(self,
               h_pdb_string,
               nuclear=False,
               largest_occupancy=10,
               b_factor_cutoff=None,
               use_segids=False,
               verbose=False):
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
    """
    assert libtbx.env.has_module(name="probe")

    self.b_factor_cutoff = b_factor_cutoff
    self.use_segids=use_segids
    ogt = 10
    blt = self.b_factor_cutoff
    if largest_occupancy < ogt:
      ogt = largest_occupancy

    self.probe_atom_b_factor = None
    if not nuclear:
      self.probe_txt = \
        'phenix.probe -u -q -mc -het -once -NOVDWOUT "ogt%d not water" "ogt%d" -' % \
          (ogt, ogt)
      #The -NOVDWOUT probe run above is faster for clashscore to parse,
      # the full_probe_txt version below is for printing to file for coot usage
      self.full_probe_txt = \
        'phenix.probe -u -q -mc -het -once "ogt%d not water" "ogt%d" -' % \
          (ogt, ogt)
      self.probe_atom_txt = \
        'phenix.probe -q -mc -het -dumpatominfo "ogt%d not water" -' % ogt
      if blt is not None:
        self.probe_atom_b_factor = \
          'phenix.probe -q -mc -het -dumpatominfo "blt%d ogt%d not water" -' % \
            (blt, ogt)
    else: #use nuclear distances
      self.probe_txt = \
        'phenix.probe -u -q -mc -het -once -NOVDWOUT -nuclear' +\
          ' "ogt%d not water" "ogt%d" -' % (ogt, ogt)
      self.full_probe_txt = \
        'phenix.probe -u -q -mc -het -once -nuclear' +\
          ' "ogt%d not water" "ogt%d" -' % (ogt, ogt)
      self.probe_atom_txt = \
        'phenix.probe -q -mc -het -dumpatominfo -nuclear' +\
          ' "ogt%d not water" -' % ogt
      if blt is not None:
        self.probe_atom_b_factor = \
          'phenix.probe -q -mc -het -dumpatominfo -nuclear' +\
            ' "blt%d ogt%d not water" -' % (blt, ogt)

    if verbose:
      print "\nUsing input model H/D atoms...\n"
    self.h_pdb_string = h_pdb_string
    self.run_probe_clashscore(self.h_pdb_string)

  def put_group_into_dict(self, line_info, clash_hash, hbond_hash):
    key = line_info.targAtom+line_info.srcAtom
    if (cmp(line_info.srcAtom,line_info.targAtom) < 0):
      key = line_info.srcAtom+line_info.targAtom

    if (line_info.type == "so" or line_info.type == "bo"):
      if (line_info.gap <= -0.4):
        if (key in clash_hash) :
          if (line_info.gap < clash_hash[key].gap):
            clash_hash[key] = line_info
        else :
          clash_hash[key] = line_info
    elif (line_info.type == "hb"):
      if (key in hbond_hash) :
        if (line_info.gap < hbond_hash[key].gap):
          hbond_hash[key] = line_info
      else :
        hbond_hash[key] = line_info

  def filter_dicts(self, new_clash_hash, new_hbond_hash):
    temp = []
    for k,v in new_clash_hash.iteritems():
      if k not in new_hbond_hash:
        temp.append(v.as_clash_obj(self.use_segids))
    return temp

  def process_raw_probe_output(self, probe_unformatted):
    new_clash_hash = {}
    new_hbond_hash = {}
    previous_line = None
    for line in probe_unformatted:
      processed=False
      try:
        line_storage = probe_line_info_storage(line)
      except KeyboardInterrupt: raise
      except ValueError:
        continue # something else (different from expected) got into output

      if previous_line is not None:
        if line_storage.is_similar(previous_line):
          # modify previous line to store this one if needed
          previous_line.gap = min(previous_line.gap, line_storage.gap)
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

  def run_probe_clashscore(self, pdb_string):
    probe_out = easy_run.fully_buffered(self.probe_txt,
      stdin_lines=pdb_string)
    if (probe_out.return_code != 0) :
      raise RuntimeError("Probe crashed - dumping stderr:\n%s" %
        "\n".join(probe_out.stderr_lines))
    probe_unformatted = probe_out.stdout_lines
    printable_probe_out = easy_run.fully_buffered(self.full_probe_txt,
                                                  stdin_lines=pdb_string)
    self.probe_unformatted = "\n".join(printable_probe_out.stdout_lines)

    temp = self.process_raw_probe_output(probe_unformatted)

    self.n_clashes = len(temp)
    if self.b_factor_cutoff is not None:
      clashes_b_cutoff = 0
      for clash_obj in temp:
        if clash_obj.max_b_factor < self.b_factor_cutoff:
          clashes_b_cutoff += 1
      self.n_clashes_b_cutoff = clashes_b_cutoff
    used = []
    self.bad_clashes = []
    for clash_obj in sorted(temp) :
      test_key = clash_obj.id_str_no_atom_name()
      test_key = clash_obj.id_str()
      if test_key not in used:
        used.append(test_key)
        self.bad_clashes.append(clash_obj)
    probe_info = easy_run.fully_buffered(self.probe_atom_txt,
      stdin_lines=pdb_string) #.raise_if_errors().stdout_lines
    err = probe_info.format_errors_if_any()
    if err is not None and err.find("No atom data in input.")>-1:
      self.clashscore = None
      self.clashscore_b_cutoff = None
      return
    #if (len(probe_info) == 0) :
    #  raise RuntimeError("Empty PROBE output.")
    self.n_atoms = 0
    for line in probe_info.stdout_lines:
      processed=False
      try:
        dump, n_atoms = line.split(":")
      except KeyboardInterrupt: raise
      except ValueError:
        pass # something else (different from expected) got into output
    self.n_atoms = int(n_atoms)
    self.natoms_b_cutoff = None
    if self.probe_atom_b_factor is not None:
      probe_info_b_factor = easy_run.fully_buffered(self.probe_atom_b_factor,
        stdin_lines=pdb_string).raise_if_errors().stdout_lines
      for line in probe_info_b_factor :
        dump_b, natoms_b_cutoff = line.split(":")
      self.natoms_b_cutoff = int(natoms_b_cutoff)
    if self.n_atoms == 0:
      clashscore = 0.0
    else:
      clashscore = (self.n_clashes * 1000) / self.n_atoms
    self.clashscore = clashscore
    clashscore_b_cutoff = None
    if self.natoms_b_cutoff is not None and self.natoms_b_cutoff == 0:
      clashscore_b_cutoff = 0.0
    elif self.natoms_b_cutoff is not None:
      clashscore_b_cutoff = \
        (self.n_clashes_b_cutoff*1000) / self.natoms_b_cutoff
    self.clashscore_b_cutoff = clashscore_b_cutoff

def decode_atom_string (atom_str, use_segids=False) :
  # Example:
  # ' A  49 LEU HD11B'
  if (not use_segids) or (len(atom_str) == 16) :
    return atom_info(
      chain_id=atom_str[0:2],
      resseq=atom_str[2:6],
      icode=atom_str[6],
      resname=atom_str[7:10],
      altloc=atom_str[15],
      name=atom_str[11:15])
  else:
    return atom_info(
      chain_id=atom_str[0:4],
      resseq=atom_str[4:8],
      icode=atom_str[8],
      resname=atom_str[9:12],
      altloc=atom_str[17],
      name=atom_str[13:17])

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
        print >> log,"\nNo H/D atoms detected - forcing hydrogen addition!\n"
      keep_hydrogens = False
  import libtbx.load_env
  has_reduce = libtbx.env.has_module(name="reduce")
  # add hydrogen if needed
  if has_reduce and (not keep_hydrogens):
    # set reduce running parameters
    build = "phenix.reduce -oh -his -flip -keep -allalt -limit{}"
    if not do_flips : build += " -pen9999"
    if nuclear:
      build += " -nuc -"
    else:
      build += " -"
    build = build.format(time_limit)
    trim = "phenix.reduce -quiet -trim -"
    stdin_lines = r.as_pdb_string(cryst_sym)
    clean_out = easy_run.fully_buffered(trim,stdin_lines=stdin_lines)
    if (clean_out.return_code != 0) :
      msg_str = "Reduce crashed with command '%s' - dumping stderr:\n%s"
      raise Sorry(msg_str % (trim, "\n".join(clean_out.stderr_lines)))
    build_out = easy_run.fully_buffered(build,stdin_lines=clean_out.stdout_lines)
    if (build_out.return_code != 0) :
      msg_str = "Reduce crashed with command '%s' - dumping stderr:\n%s"
      raise Sorry(msg_str % (build, "\n".join(build_out.stderr_lines)))
    reduce_str = '\n'.join(build_out.stdout_lines)
    return reduce_str,True
  else:
    if not has_reduce:
      msg = 'phenix.reduce could not be detected on your system.\n'
      msg += 'Cannot add hydrogen to PDB file'
      print >> log,msg
    return r.as_pdb_string(cryst_sym),False

#-----------------------------------------------------------------------
# this isn't really enough code to justify a separate module...
#
class nqh_flip (residue) :
  """
  Backwards Asn/Gln/His sidechain, identified by Reduce's hydrogen-bond
  network optimization.
  """
  def as_string (self) :
    return self.id_str()

  def as_table_row_phenix (self) :
    return [ self.chain_id, "%1s%s %s" % (self.altloc,self.resname,self.resid) ]

class nqh_flips (validation) :
  """
  N/Q/H sidechain flips identified by Reduce.
  """
  gui_list_headers = ["Chain", "Residue"]
  gui_formats = ["%s", "%s"]
  wx_column_widths = [75,220]
  def __init__ (self, pdb_hierarchy) :
    re_flip = re.compile(":FLIP")
    validation.__init__(self)
    reduce_out = easy_run.fully_buffered("phenix.reduce -BUILD -",
      stdin_lines=pdb_hierarchy.as_pdb_string())
    for line in reduce_out.stdout_lines:
    #USER  MOD Set 1.1: B  49 GLN     :FLIP  amide:sc=    -2.7! C(o=-5.8!,f=-1.3!)
      if re_flip.search(line) :
        resid = line.split(":")[1]
        chain_id = resid[0:2].strip()
        if (len(chain_id) == 0):
          chain_id = ' '
        resname = resid[7:10]
        assert (resname in ["ASN", "GLN", "HIS"])
        flip = nqh_flip(
          chain_id=chain_id,
          resseq=resid[2:6].strip(),
          icode=resid[6:7],
          altloc=resid[14:15],
          resname=resname,
          outlier=True)
        flip.set_coordinates_from_hierarchy(pdb_hierarchy)
        self.results.append(flip)
        self.n_outliers += 1

  def show (self, out=sys.stdout, prefix="") :
    if (self.n_outliers == 0) :
      print >> out, prefix+"No backwards Asn/Gln/His sidechains found."
    else :
      for flip in self.results :
        print >> out, prefix+flip.as_string()
