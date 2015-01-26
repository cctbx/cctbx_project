
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
    return "%s :%.3f" % (self.id_str(), self.overlap)

  def as_string (self) :
    return "%-20s %-20s  %7.3f" % (self.atoms_info[0].id_str(),
      self.atoms_info[1].id_str(), self.overlap)

  def as_table_row_phenix (self) :
    return [ self.atoms_info[0].id_str(), self.atoms_info[1].id_str(),
             self.overlap ]

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
  ]
  program_description = "Analyze clashscore for protein model"
  gui_list_headers = ["Atom 1", "Atom 2", "Overlap"]
  gui_formats = ["%s", "%s", ".3f"]
  wx_column_widths = [200] * 3

  def get_result_class (self) : return clash

  def __init__ (self,
      pdb_hierarchy,
      keep_hydrogens=True,
      nuclear=False,
      force_unique_chain_ids=False,
      time_limit=120,
      b_factor_cutoff=None,
      save_probe_unformatted_file=None,
      save_modified_hierarchy=False,
      verbose=False,
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
        log=out)
      r = iotbx.pdb.hierarchy.root()
      mdc = model.detached_copy()
      r.append_model(mdc)
      occ_max = flex.max(r.atoms().extract_occ())
      pcm = probe_clashscore_manager(
        h_pdb_string=input_str,
        nuclear=nuclear,
        largest_occupancy=occ_max,
        b_factor_cutoff=b_factor_cutoff,
        use_segids=use_segids,
        verbose=verbose)
      if (save_modified_hierarchy) :
        self.pdb_hierarchy = iotbx.pdb.hierarchy.input(
          pdb_string=pcm.h_pdb_string).hierarchy
      self.clash_dict[model.id] = pcm.clashscore
      self.clash_dict_b_cutoff[model.id] = pcm.clashscore_b_cutoff
      self.list_dict[model.id] = pcm.bad_clashes
      if (n_models == 1) or (self.clashscore is None) :
        self.results = pcm.bad_clashes
        self.n_outliers = len(self.results)
        self.clashscore = pcm.clashscore
        self.clashscore_b_cutoff = pcm.clashscore_b_cutoff
      if (save_probe_unformatted_file is not None) and (n_models == 1) :
        open(save_probe_unformatted_file, "w").write(pcm.probe_unformatted)
        self.probe_file = save_probe_unformatted_file

  def get_clashscore(self):
    return self.clashscore

  def get_clashscore_b_cutoff(self):
    return self.clashscore_b_cutoff

  def show_old_output (self, out=sys.stdout, verbose=False) :
    self.print_clashlist_old(out=out)
    self.show_summary(out=out)

  def show_summary (self, out=sys.stdout, prefix="") :
    if (len(self.clash_dict) == 1) :
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

  def run_probe_clashscore(self, pdb_string):
    clash_hash = {}
    hbond_hash = {}
    probe_out = easy_run.fully_buffered(self.probe_txt,
      stdin_lines=pdb_string)
    if (probe_out.return_code != 0) :
      raise RuntimeError("Probe crashed - dumping stderr:\n%s" %
        "\n".join(probe_out.stderr_lines))
    probe_unformatted = probe_out.stdout_lines
    self.probe_unformatted = "\n".join(probe_unformatted)
    for line in probe_unformatted:
      processed=False
      try:
        name, pat, type, srcAtom, targAtom, min_gap, gap, \
        kissEdge2BullsEye, dot2BE, dot2SC, spike, score, stype, \
        ttype, x, y, z, sBval, tBval = line.split(":")
        processed=True
      except KeyboardInterrupt: raise
      except ValueError:
        pass # something else (different from expected) got into output
      if(processed):
        atom1 = decode_atom_string(srcAtom, self.use_segids)
        atom2 = decode_atom_string(targAtom, self.use_segids)
        if (cmp(srcAtom,targAtom) < 0):
          atoms = [ atom1, atom2 ]
        else:
          atoms = [ atom2, atom1 ]
        gap = float(gap)
        x, y, z = float(x), float(y), float(z)
        clash_obj = clash(
          atoms_info=atoms,
          overlap=gap,
          probe_type=type,
          outlier=abs(gap) > 0.4,
          max_b_factor=max(float(sBval), float(tBval)),
          xyz=(x,y,z))
        key = clash_obj
        if (type == "so" or type == "bo"):
          if (gap <= -0.4):
            if (key in clash_hash) :
              if (gap < clash_hash[key].overlap):
                clash_hash[key] = clash_obj
            else :
              clash_hash[key] = clash_obj
        elif (type == "hb"):
          if (key in hbond_hash) :
            if (gap < hbond_hash[key].overlap):
              hbond_hash[key] = clash_obj
          else :
            hbond_hash[key] = clash_obj
    #sort the output
    temp = []
    for k in clash_hash.keys():
      if not k in hbond_hash:
        temp.append(clash_hash[k])
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
      if test_key not in used:
        used.append(test_key)
        self.bad_clashes.append(clash_obj)
    probe_info = easy_run.fully_buffered(self.probe_atom_txt,
      stdin_lines=pdb_string).raise_if_errors().stdout_lines
    if (len(probe_info) == 0) :
      raise RuntimeError("Empty PROBE output.")
    self.n_atoms = 0
    for line in probe_info:
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
    h_count = elements.count(' H') + elements.count(' D')
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
    build = "phenix.reduce -oh -his -flip -pen9999 -keep -allalt -limit{}"
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
    return [ self.chain_id, "%s %s" % (self.resname, self.resid) ]

class nqh_flips (validation) :
  """
  N/Q/H sidechain flips identified by Reduce.
  """
  gui_list_headers = ["Chain", "Residue"]
  gui_formats = ["%s", "%s"]
  wx_column_widths = [100,220]
  def __init__ (self, pdb_hierarchy) :
    re_flip = re.compile(":FLIP")
    validation.__init__(self)
    reduce_out = easy_run.fully_buffered("phenix.reduce -BUILD -",
      stdin_lines=pdb_hierarchy.as_pdb_string())
    for line in reduce_out.stderr_lines :
    #orientation 4: A  68 HIS     :FLIP no HD1: bump=-0.607, HB=0.998, total=0.390
      if re_flip.search(line) :
        resname = line[22:25]
        assert (resname in ["ASN", "GLN", "HIS"])
        flip = nqh_flip(
          chain_id=line[16],
          resseq=line[17:21].strip(),
          icode=line[21],
          altloc=line[29],
          resname=resname)
        flip.set_coordinates_from_hierarchy(pdb_hierarchy)
        self.results.append(flip)
        self.n_outliers += 1

  def show (self, out=sys.stdout, prefix="") :
    if (self.n_outliers == 0) :
      print >> out, prefix+"No backwards Asn/Gln/His sidechains found."
    else :
      for flip in self.results :
        print >> out, prefix+flip.as_string()
