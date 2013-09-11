from __future__ import division
import libtbx.load_env
import sys, os, string
from libtbx.utils import Usage
from mmtbx import utils

try:
  from iotbx import pdb
except ImportError, e:
  print "iotbx not loaded"
  sys.exit()

import iotbx.phil
from libtbx import easy_run

def get_master_phil():
  return iotbx.phil.parse(
    input_string="""
      clashscore {
        pdb = None
        .type = path
        .help = '''Enter a PDB file name'''

        verbose = False
        .type = bool
        .help = '''Verbose'''

        keep_hydrogens = True
        .type = bool
        .help = '''Keep hydrogens in input file'''

        nuclear = False
        .type = bool
        .help = '''Use nuclear hydrogen positions'''

        time_limit = 120
        .type = int
        .help = '''Time limit (sec) for Reduce optimization'''

        b_factor_cutoff = None
        .type = int
        .help = '''B factor cutoff for use with MolProbity'''
  }

    """)

class clashscore(object):
  #flag routines-----------------------------------------------------------------------------------
  def usage(self):
    return """
phenix.clashscore file.pdb [params.eff] [options ...]

Options:

  pdb=input_file        input PDB file
  keep_hydrogens=True   keep input hydrogen files (otherwise regenerate)
  nuclear=False         use nuclear x-H distances and vdW radii
  verbose=True          verbose text output

Example:

  phenix.clashscore pdb=1ubq.pdb keep_hydrogens=True

"""
  def get_summary_and_header(self,command_name):
    header="\n"
    header+="\n#                       "+str(command_name)
    header+="\n#"
    header+="\n# Analyze clashscore for protein model"
    header+="\n# type phenix."+str(command_name)+": --help for help\n"

    summary= "phenix.%s [options] mypdb.pdb" % command_name
    return summary,header
  #------------------------------------------------------------------------------------------------

  #{{{ run
  def run(self, args, out=sys.stdout, quiet=False):
    if (len(args) == 0 or "--help" in args or "--h" in args or "-h" in args):
      raise Usage(self.usage())
    master_phil = get_master_phil()
    import iotbx.utils
    input_objects = iotbx.utils.process_command_line_inputs(
      args=args,
      master_phil=master_phil,
      input_types=("pdb",))
    work_phil = master_phil.fetch(sources=input_objects["phil"])
    work_params = work_phil.extract()
    if len(input_objects["pdb"]) != 1:
      summary, header = self.get_summary_and_header("clashscore")
      raise Usage(summary)
    file_obj = input_objects["pdb"][0]
    filename = file_obj.file_name

    command_name = "clashscore"
    summary,header=self.get_summary_and_header(command_name)
    if not quiet: print >>out, header

    #TO DO: make this a useful help section
    #if help or (params and params.clashscore.verbose):
    #  print >> out, summary
    #  pass
      #print "Values of all params:"
      #master_params.format(python_object=params).show(out=out)

    self.params=work_params # makes params available to whole class

    log=out
    if (log is None): log = sys.stdout
    if self.params.clashscore.verbose :
      print >> log, 'filename', filename
    if filename and os.path.exists(filename):
      pdb_io = pdb.input(filename)
      pass
    else:
      print >> log, "Please enter a file name"
      return
    keep_hydrogens = self.params.clashscore.keep_hydrogens
    nuclear = self.params.clashscore.nuclear
    time_limit = self.params.clashscore.time_limit
    clashscore, bad_clashes = self.analyze_clashes(pdb_io=pdb_io,
        keep_hydrogens=keep_hydrogens,
        nuclear=nuclear,
        time_limit=time_limit,
        b_factor_cutoff=self.params.clashscore.b_factor_cutoff,
        verbose=self.params.clashscore.verbose)
    if not quiet :
      self.print_clashlist(out)
      self.print_clashscore(out)
    return clashscore
  #}}}

  #{{{ analyze_clashes
  def analyze_clashes(
        self,
        pdb_io=None,
        hierarchy=None,
        keep_hydrogens=True,
        nuclear=False,
        force_unique_chain_ids=False,
        time_limit=120,
        b_factor_cutoff=None,
        verbose=False) :
    if (not libtbx.env.has_module(name="probe")):
      print "Probe could not be detected on your system.  Please make sure Probe is in your path."
      print "Probe is available at http://kinemage.biochem.duke.edu/"
      sys.exit()
    assert [pdb_io, hierarchy].count(None) == 1
    if(pdb_io is not None):
      hierarchy = pdb_io.construct_hierarchy()

    self.clash_dict = {}
    self.clash_dict_b_cutoff = {}
    self.list_dict = {}

    h_count = 0

    if verbose:
      if not nuclear:
        print "\nUsing electron cloud x-H distances and vdW radii"
      else:
        print "\nUsing nuclear cloud x-H distances and vdW radii"

    for i,m in enumerate(hierarchy.models()):
      r = iotbx.pdb.hierarchy.root()
      mdc = m.detached_copy()
      r.append_model(mdc)
      tmp_r = r
      # removed old style SEGID handling for compatibility with Probe
      # 130622 - JJH
      #bare_chains = \
      #  utils.find_bare_chains_with_segids(pdb_hierarchy=r)
      #if bare_chains:
      #  tmp_r = r.deep_copy()
      #  tmp_r.atoms().reset_i_seq()
      #  seg_dict = utils.seg_id_to_chain_id(pdb_hierarchy=tmp_r)
      #  rename_txt = utils.assign_chain_ids(pdb_hierarchy=tmp_r,
      #                                      seg_dict=seg_dict)
      #else:
      #  tmp_r = r
      #duplicate_chain_ids = \
      #  utils.check_for_duplicate_chain_ids(pdb_hierarchy=tmp_r)
      #if duplicate_chain_ids:
      #  utils.force_unique_chain_ids(pdb_hierarchy=tmp_r)
      if keep_hydrogens:
        elements = tmp_r.atoms().extract_element()
        h_count = elements.count(' H') + elements.count(' D')
        if h_count > 0:
          has_hd = True
        else:
          has_hd = False
        # if no hydrogens present, force addition for clashscore
        # calculation
        if not has_hd:
          if verbose:
            print "\nNo H/D atoms detected - forcing hydrogen addition!\n"
          keep_hydrogens = False
      input_str = tmp_r.as_pdb_string()
      largest_occupancy = self.get_largest_occupancy(atoms=tmp_r.atoms())
      pcm = probe_clashscore_manager(pdb_string=input_str,
                                     keep_hydrogens=keep_hydrogens,
                                     nuclear=nuclear,
                                     time_limit=time_limit,
                                     largest_occupancy=largest_occupancy,
                                     b_factor_cutoff=b_factor_cutoff,
                                     verbose=verbose)
      self.pdb_hierarchy = pdb.hierarchy.\
        input(pdb_string=pcm.h_pdb_string).hierarchy
      self.clash_dict[m.id]=pcm.clashscore
      self.clash_dict_b_cutoff[m.id]=pcm.clashscore_b_cutoff
      self.list_dict[m.id]=pcm.bad_clashes
    self.clashscore=pcm.clashscore
    self.clashscore_b_cutoff=pcm.clashscore_b_cutoff
    self.bad_clashes=pcm.bad_clashes
    self.probe_unformatted =pcm.probe_unformatted
    return self.clash_dict, self.list_dict
  #}}}

  #{{{ get_functions
  def get_clashscore(self):
    return self.clashscore

  def print_clashscore(self, out=sys.stdout):
    if out is None :
      out = sys.stdout
    for k in sorted(self.clash_dict.keys()):
      if k is '':
        print >> out, "clashscore = %f" % self.clash_dict[k]
        if self.clash_dict_b_cutoff[k] is not None:
          print >> out, "clashscore (B factor cutoff = %d) = %f" % \
            (self.params.clashscore.b_factor_cutoff,
             self.clash_dict_b_cutoff[k])
      else:
        print >> out, "MODEL%s clashscore = %f" % (k, self.clash_dict[k])
        if self.clash_dict_b_cutoff[k] is not None:
          print >> out, "MODEL%s clashscore (B factor cutoff = %d) = %f" % \
            (k,
             self.params.clashscore.b_factor_cutoff,
             self.clash_dict_b_cutoff[k])

  def print_clashlist(self, out=sys.stdout):
    for k in sorted(self.list_dict.keys()):
      if k is '':
        print >> out, "Bad Clashes >= 0.4 Angstrom:"
        print >> out, self.list_dict[k]
      else:
        print >> out, "Bad Clashes >= 0.4 Angstrom MODEL%s" % k
        print >> out, self.list_dict[k]

  def get_largest_occupancy(self, atoms):
    largest_occupancy = 0.0
    for atom in atoms:
      if atom.occ > largest_occupancy:
        largest_occupancy = float(atom.occ)
    return int(largest_occupancy * 100)

class probe_clashscore_manager(object):
  def __init__(self,
               pdb_string,
               keep_hydrogens=True,
               nuclear=False,
               time_limit=120,
               largest_occupancy=10,
               b_factor_cutoff=None,
               verbose=False):
    assert (libtbx.env.has_module(name="reduce") and
            libtbx.env.has_module(name="probe"))

    self.b_factor_cutoff = b_factor_cutoff
    ogt = 10
    blt = self.b_factor_cutoff
    if largest_occupancy < ogt:
      ogt = largest_occupancy

    self.trim = "phenix.reduce -quiet -trim -"
    self.probe_atom_b_factor = None
    if not nuclear:
      self.build = "phenix.reduce -oh -his -flip -pen9999" +\
                   " -keep -allalt -limit%d -" % time_limit
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
      self.build = "phenix.reduce -oh -his -flip -pen9999" +\
                   " -keep -allalt -limit%d -nuc -" % time_limit
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

    if not keep_hydrogens:
      h_pdb_string = self.run_reduce(pdb_string)
    else:
      if verbose:
        print "\nUsing input model H/D atoms...\n"
      h_pdb_string = pdb_string
    self.h_pdb_string = h_pdb_string
    self.run_probe_clashscore(self.h_pdb_string)

  def run_reduce(self, pdb_string):
    clean_out = easy_run.fully_buffered(self.trim,
                  stdin_lines=pdb_string)
    if (clean_out.return_code != 0) :
      raise RuntimeError("Reduce crashed with command '%s' - dumping stderr:\n%s"
        % (self.trim, "\n".join(clean_out.stderr_lines)))
    build_out = easy_run.fully_buffered(self.build,
                  stdin_lines=clean_out.stdout_lines)
    if (build_out.return_code != 0) :
      raise RuntimeError("Reduce crashed with command '%s' - dumping stderr:\n%s"
        % (self.build, "\n".join(build_out.stderr_lines)))
    reduce_str = string.join(build_out.stdout_lines, '\n')
    return reduce_str

  #def update_clashscore(self, pdb_string):
  #  self.run_probe_clashscore(pdb_string)

  def run_probe_clashscore(self, pdb_string):
    clash_hash={}
    hbond_hash={}
    b_factor_hash={}
    clashscore = None
    probe_out = easy_run.fully_buffered(self.probe_txt,
      stdin_lines=pdb_string)
    if (probe_out.return_code != 0) :
      raise RuntimeError("Probe crashed - dumping stderr:\n%s" %
        "\n".join(probe_out.stderr_lines))
    probe_unformatted = probe_out.stdout_lines
    self.probe_unformatted = probe_unformatted
    for line in probe_unformatted:
      name, pat, type, srcAtom, targAtom, min_gap, gap, \
      kissEdge2BullsEye, dot2BE, dot2SC, spike, score, stype, \
      ttype, x, y, z, sBval, tBval = line.split(":")
      if (cmp(srcAtom,targAtom) < 0):
        key = srcAtom+targAtom
      else:
        key = targAtom+srcAtom
      if (type == "so" or type == "bo"):
        if (float(gap) <= -0.4):
          if clash_hash.get(key) is not None:
            if (float(gap) < clash_hash[key]):
              clash_hash[key] = float(gap)
              b_factor_hash[key] = max(float(sBval), float(tBval))
          else:
            clash_hash[key] = float(gap)
            b_factor_hash[key] = max(float(sBval), float(tBval))
      elif (type == "hb"):
        if hbond_hash.get(key) is not None:
          if (float(gap) < hbond_hash[key]):
            hbond_hash[key] = float(gap)
        else:
          hbond_hash[key] = float(gap)
    clashes = len(clash_hash)

    for k in clash_hash.keys():
      if k in hbond_hash:
        clashes=clashes-1
        clash_hash[k]="Hbonded"
    bad_clashes = ''

    if self.b_factor_cutoff is not None:
      clashes_b_cutoff = 0
      for k in clash_hash.keys():
        if clash_hash[k] == "Hbonded":
          continue
        if b_factor_hash[k] < self.b_factor_cutoff:
          clashes_b_cutoff += 1

    #sort the output
    temp = []
    for k in clash_hash.keys():
      if not k in hbond_hash:
        temp.append(k)
    def get_clash(k):
      return clash_hash[k]
    temp_sorted = sorted(temp, key=get_clash)
    used = []
    for k in temp_sorted:
      test_key = k[0:11]+k[16:27]
      if test_key not in used:
        bad_clashes += k+':'+str(clash_hash[k])+'\n'
        used.append(test_key)
    probe_info = easy_run.fully_buffered(self.probe_atom_txt,
      stdin_lines=pdb_string).raise_if_errors().stdout_lines
    if (len(probe_info) == 0) :
      raise RuntimeError("Empty PROBE output.")
    natoms = 0
    for line in probe_info :
      dump, natoms = line.split(":")
    natoms_b_cutoff = None
    if self.probe_atom_b_factor is not None:
      probe_info_b_factor = easy_run.fully_buffered(self.probe_atom_b_factor,
        stdin_lines=pdb_string).raise_if_errors().stdout_lines
      for line in probe_info_b_factor :
        dump_b, natoms_b_cutoff = line.split(":")

    if int(natoms) == 0:
      clashscore = 0.0
    else:
      clashscore = (clashes*1000)/float(natoms)
    clashscore_b_cutoff = None
    if natoms_b_cutoff is not None and int(natoms_b_cutoff) == 0:
      clashscore_b_cutoff = 0.0
    elif natoms_b_cutoff is not None:
      clashscore_b_cutoff = (clashes_b_cutoff*1000)/float(natoms_b_cutoff)
    self.clashscore = clashscore
    self.clashscore_b_cutoff = clashscore_b_cutoff
    self.bad_clashes = bad_clashes
