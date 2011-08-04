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

        changes = False
        .type = bool
        .help = '''Print list of changes'''

        version = False
        .type = bool
        .help = '''Print version'''

        verbose = False
        .type = bool
        .help = '''Verbose'''

        keep_hydrogens = False
        .type = bool
        .help = '''Keep hydrogens in input file'''
  }

    """)

class clashscore(object):
  #flag routines-----------------------------------------------------------------------------------
  def usage(self):
    return """
phenix.clashscore file.pdb [params.eff] [options ...]

Options:

  pdb=input_file        input PDB file
  keep_hydrogens=False  keep input hydrogen files (otherwise regenerate)
  verbose=True          verbose text output

Example:

  phenix.clashscore pdb=1ubq.pdb keep_hydrogens=True

"""
  def changes(self):
    print "\nversion 0.10 - Development Version\n"
  def version(self):
    print "\nversion 0.10 - Copyright 2009, Jeffrey J. Headd and Mike Word\n"
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

    if self.params.clashscore.changes:
      self.changes()
      return

    if self.params.clashscore.version:
      self.version()
      return

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
    clashscore, bad_clashes = self.analyze_clashes(pdb_io=pdb_io,
        keep_hydrogens=self.params.clashscore.keep_hydrogens)
    if not quiet :
      self.print_clashlist(out)
      self.print_clashscore(out)
    return clashscore
  #}}}

  #{{{ analyze_clashes
  def analyze_clashes(self, pdb_io=None, hierarchy=None, keep_hydrogens=False) :
    if (not libtbx.env.has_module(name="probe")):
      print "Probe could not be detected on your system.  Please make sure Probe is in your path."
      print "Probe is available at http://kinemage.biochem.duke.edu/"
      sys.exit()
    assert [pdb_io, hierarchy].count(None) == 1
    if(pdb_io is not None):
      hierarchy = pdb_io.construct_hierarchy()

    self.clashscore = []
    self.bad_clashes_list = []
    self.clash_dict = {}
    self.list_dict = {}

    trim = "phenix.reduce -quiet -trim -"
    build = "phenix.reduce -oh -his -flip -pen9999 -keep -allalt -"
    probe = 'phenix.probe -u -q -mc -het -once "ogt33 not water" "ogt33" -'
    probe_atom = 'phenix.probe -q -mc -het -dumpatominfo "ogt33 not water" -'

    for i,m in enumerate(hierarchy.models()):
      r = iotbx.pdb.hierarchy.root()
      mdc = m.detached_copy()
      r.append_model(mdc)
      bare_chains = \
        utils.find_bare_chains_with_segids(pdb_hierarchy=r)
      if bare_chains:
        tmp_r = r.deep_copy()
        tmp_r.atoms().reset_i_seq()
        seg_dict = utils.seg_id_to_chain_id(pdb_hierarchy=tmp_r)
        rename_txt = utils.assign_chain_ids(pdb_hierarchy=tmp_r,
                                            seg_dict=seg_dict)
      else:
        tmp_r = r
      if keep_hydrogens is False:
        clean_out = easy_run.fully_buffered(trim,
                                    stdin_lines=tmp_r.as_pdb_string())
        build_out = easy_run.fully_buffered(build,
                                    stdin_lines=clean_out.stdout_lines)
        input_str = string.join(build_out.stdout_lines, '\n')
      else:
        input_str = tmp_r.as_pdb_string()
      self.pdb_hierarchy = pdb.hierarchy.input(pdb_string=input_str).hierarchy
      clash_hash={}
      hbond_hash={}
      probe_unformatted = easy_run.fully_buffered(probe,
                         stdin_lines=input_str).stdout_lines
      for line in probe_unformatted:
        name, pat, type, srcAtom, targAtom, min_gap, gap, kissEdge2BullsEye, dot2BE, \
        dot2SC, spike, score, stype, ttype, x, y, z, sbVal, tBval = line.split(":")
        if (cmp(srcAtom,targAtom) < 0):
          key = srcAtom+targAtom
        else:
          key = targAtom+srcAtom
        if (type == "so" or type == "bo"):
          if (float(gap) <= -0.4):
            try:
              if (float(gap) < clash_hash[key]):
                clash_hash[key] = float(gap)
            except Exception:
              clash_hash[key] = float(gap)
        elif (type == "hb"):
          try:
            if (float(gap) < hbond_hash[key]):
              hbond_hash[key] = float(gap)
          except Exception:
            hbond_hash[key] = float(gap)
      clashes = len(clash_hash)

      for k in clash_hash.keys():
        if k in hbond_hash:
          clashes=clashes-1
          clash_hash[k]="Hbonded"
      bad_clashes = ''

      for k in clash_hash.keys():
        if not k in hbond_hash:
          bad_clashes += k+':'+str(clash_hash[k])+'\n'
      probe_info = easy_run.fully_buffered(probe_atom,
        stdin_lines=input_str).raise_if_errors().stdout_lines
      for line in probe_info :
        dump, natoms = line.split(":")

      if int(natoms) == 0:
        clashscore = 0.0
      else:
        clashscore = (clashes*1000)/float(natoms)
      self.clashscore.append(clashscore)
      self.bad_clashes_list.append(bad_clashes)
      self.clash_dict[m.id]=clashscore
      self.list_dict[m.id]=bad_clashes
    self.clashscore=clashscore
    self.bad_clashes=bad_clashes
    self.probe_unformatted = probe_unformatted
    return self.clash_dict, self.list_dict
  #}}}

  #{{{ get_functions
  def get_clashscore(self):
    return self.clashscore

  def get_bad_clashes_list(self):
    return self.bad_clashes_list

  def print_clashscore(self, out=sys.stdout):
    if out is None :
      out = sys.stdout
    for k in self.clash_dict.keys():
      if k is '':
        print >> out, "clashscore = %f" % self.clash_dict[k]
      else:
        print >> out, "MODEL%s clashscore = %f" % (k, self.clash_dict[k])

  def print_clashlist(self, out=sys.stdout):
    for k in self.list_dict.keys():
      if k is '':
        print >> out, "Bad Clashes >= 0.4 Angstrom:"
        print >> out, self.list_dict[k]
      else:
        print >> out, "Bad Clashes >= 0.4 Angstrom MODEL%s" % k
        print >> out, self.list_dict[k]
