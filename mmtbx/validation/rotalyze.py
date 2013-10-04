from __future__ import division
#(jEdit options) :folding=explicit:collapseFolds=1:

import sys, os
from iotbx import pdb
from mmtbx.rotamer.sidechain_angles import SidechainAngles
from mmtbx.rotamer import rotamer_eval
from mmtbx.rotamer.rotamer_eval import RotamerID
from mmtbx.validation import utils
import iotbx.phil
from libtbx.utils import Usage, Sorry

def get_master_phil():
  return iotbx.phil.parse(
    input_string="""
    rotalyze {
      pdb = None
        .type = path
        .help = '''Enter a PDB file name'''

      outliers_only = False
      .type = bool
      .help = '''Only print rotamer outliers'''

      verbose = True
      .type = bool
      .help = '''Verbose'''

      show_errors = False
      .type = bool
      .help = '''Print out errors'''
}
""")

header = """residue:occupancy:score%:chi1:chi2:chi3:chi4:rotamer"""

class rotalyze(object):

  #{{{ flag routines
  #flag routines-----------------------------------------------------------------------------------
  def usage(self):
    return """
phenix.rotalyze file.pdb [params.eff] [options ...]

Options:

  pdb=input_file        input PDB file
  outliers_only=False   only print outliers
  verbose=False         verbose text output

Example:

  phenix.rotalyze pdb=1ubq.pdb outliers_only=True

"""
  #------------------------------------------------------------------------------------------------
  #}}}

  #{{{ get_summary_and_header
  def get_summary_and_header(self,command_name):
    header="\n"
    header+="\n#                       "+str(command_name)
    header+="\n#"
    header+="\n# Analyze protein sidechain rotamers"
    header+="\n# type phenix."+str(command_name)+": --help for help"

    summary= "phenix.%s mypdb.pdb" % command_name
    return summary,header
  #}}}

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
      summary, header = self.get_summary_and_header("rotalyze")
      raise Usage(summary)
    file_obj = input_objects["pdb"][0]
    filename = file_obj.file_name

    command_name = "rotalyze"
    summary,header=self.get_summary_and_header(command_name)
    if not quiet: print >>out, header

    #TO DO:  rework this help routine to functional
    #if help or (params and params.rotalyze.verbose):
    #  pass
      # XXX: this is really distracting!
      #print "Values of all params:"
      #master_params.format(python_object=params).show(out=out)

    self.params=work_params # makes params available to whole class

    log=out
    if (log is None): log = sys.stdout
    if filename and os.path.exists(filename):
      pdb_io = pdb.input(filename)
    else:
      print "Please enter a file name"
      return None, None

    output_text, output_list = self.analyze_pdb(
      pdb_io=pdb_io,
      outliers_only=self.params.rotalyze.outliers_only,
      show_errors=self.params.rotalyze.show_errors)
    out_count, out_percent = self.get_outliers_count_and_fraction()
    if self.params.rotalyze.verbose:
      print >> out, "residue:occupancy:score%:chi1:chi2:chi3:chi4:rotamer"
      print >> out, output_text
      print >> out, 'SUMMARY: %.2f%% outliers (Goal: %s)' % (out_percent*100,
        self.get_outliers_goal())
    todo_list = self.coot_todo(output_list)
    self.out_count = out_count
    self.out_percent = out_percent * 100.0
    return output_list, todo_list
  #}}}

  #{{{ analyze_pdb
  def analyze_pdb(self, pdb_io=None, hierarchy=None, outliers_only=False,
                  show_errors = False, out=sys.stdout):
    self.sa = SidechainAngles(show_errors)
    self.r = rotamer_eval.RotamerEval()
    if(pdb_io is not None):
      hierarchy = pdb_io.construct_hierarchy()
    use_segids = utils.use_segids_in_place_of_chainids(
                   hierarchy=hierarchy)
    analysis = ""
    output_list = []
    self.current_rotamers = {}
    self.numoutliers = 0
    self.numtotal = 0
    self.rot_id = rotamer_eval.RotamerID() # loads in the rotamer names
    for model in hierarchy.models():
      for chain in model.chains():
        if use_segids:
          chain_id = utils.get_segid_as_chainid(chain=chain)
        else:
          chain_id = chain.id
        for rg in chain.residue_groups():
          all_dict = self.construct_complete_sidechain(rg)
          #print all_dict
          for atom_group in rg.atom_groups() :
            resname = atom_group.resname
            coords = self.get_center(atom_group)
            atom_dict = all_dict.get(atom_group.altloc)
            res_key = self.get_residue_key(atom_group=atom_group)
            try:
              chis = self.sa.measureChiAngles(
                       atom_group,
                       atom_dict)#.get(conformer.altloc))
            except AttributeError:
              if show_errors:
                res_info = "%s%5s %s" % (chain_id, rg.resid(),
                  atom_group.altloc+resname)
                print >> out, '%s is missing some sidechain atoms' % res_info
                output_list.append([chain_id, rg.resid(),
                  atom_group.altloc+resname, -1, None, None, None, None,
                  "INCOMPLETE", coords])
              continue
            if (chis is not None):
              if None in chis:
                continue
              value = self.r.evaluate(resname.lower().strip(), chis)
              if value != None:
                occupancy = get_occupancy(atom_group)
                self.numtotal += 1
                s = '%s%5s %s:%.2f:%.1f' % \
                  (chain_id,
                   rg.resid(),
                   atom_group.altloc+resname.strip(),
                   occupancy,
                   value*100)
                res_out_list = [chain_id,
                                rg.resid(),
                                atom_group.altloc+resname,
                                value*100]
                wrap_chis = \
                  self.rot_id.wrap_chis(resname.strip(), chis, symmetry=False)
                sym_chis = wrap_chis[:]
                sym_chis = \
                  self.rot_id.wrap_sym(resname.strip(), sym_chis)
                for i in range(4):
                  s += ':'
                  if i < len(wrap_chis):
                    s += '%.1f' % (wrap_chis[i])
                    res_out_list.append(sym_chis[i])
                  else :
                    res_out_list.append(None)
                s += ':'
                if value < 0.01:
                  self.numoutliers += 1
                  s += "OUTLIER\n"
                  res_out_list.append("OUTLIER")
                  res_out_list.append(coords)
                  if outliers_only:
                    analysis += s
                    output_list.append(res_out_list)
                  self.current_rotamers[res_key] = "OUTLIER"
                else:
                  s += self.rot_id.identify(resname, wrap_chis) + "\n"
                  res_out_list.append(self.rot_id.identify(resname, wrap_chis))
                  res_out_list.append(coords)
                  self.current_rotamers[res_key] = \
                    self.rot_id.identify(resname, wrap_chis)
                  if self.current_rotamers[res_key] == "":
                    self.current_rotamers[res_key] = "OUTLIER"
                if not outliers_only:
                  analysis += s
                  output_list.append(res_out_list)
    return analysis.rstrip(), output_list
  #}}}

  def split_rotamer_names(self, rotamer):
    split_rotamer = []
    multi = ""
    if rotamer in ['OUTLIER', 'Cg_exo', 'Cg_endo']:
      split_rotamer.append(rotamer)
      return split_rotamer
    for i, c in enumerate(rotamer):
      if c in ['t', 'p', 'm']:
        split_rotamer.append(c)
      elif (c in ['-', '?']) or (c >= '0' and c<='9'):
        multi += c
    if len(multi) > 0:
      split_rotamer.append(multi)
    return split_rotamer

  def get_residue_key(self, atom_group):
    altloc = atom_group.altloc
    if altloc == "":
      altloc = " "
    key = None
    for atom in atom_group.atoms():
      cur_label = atom.pdb_label_columns()+atom.segid
      cur_altloc = cur_label[4:5]
      if key is None:
        if altloc == cur_altloc:
          key = cur_label[4:]
      else:
        if altloc == cur_altloc:
          if (key != cur_label[4:]) :
            raise Sorry("""\
Incompatible identifiers for one or more atoms in a residue:
%s
This is usually caused by atoms with a different segid from the rest of the
residue.  You can use phenix.pdbtools or phenix.pdb_editor to reset the
segid.""" % atom.format_atom_record())
          assert key == cur_label[4:]
    return key

  #{{{ evaluate_residue
  def evaluate_residue(
        self,
        residue_group,
        sa,
        r,
        all_dict,
        sites_cart=None):
    is_outlier = False
    for ag in residue_group.atom_groups():
      atom_dict = all_dict.get(ag.altloc)
      try:
        chis = sa.measureChiAngles(ag, atom_dict)
        value = r.evaluate(ag.resname.lower().strip(), chis, sites_cart)
      except Exception:
        #print ag.resname.lower()+residue_group.resseq+" is missing some sidechain atoms"
        value = None;
        is_outlier = None;
        return is_outlier, value
      if (value is None):
        is_outlier = False
        return is_outlier, value
      elif (value < 0.01):
        is_outlier = True
        return is_outlier, value
      else:
        return is_outlier, value
  #}}}

  #{{{ evalute_rotamer
  def evaluate_rotamer(
        self,
        atom_group,
        all_dict,
        sites_cart=None):
    atom_dict = all_dict.get(atom_group.altloc)
    resname = atom_group.resname
    try:
      chis = self.sa.measureChiAngles(atom_group, atom_dict, sites_cart)
      value = self.r.evaluate(
                atom_group.resname.lower().strip(),
                chis)
    except Exception:
      return None, None, None
    wrap_chis = \
      self.rot_id.wrap_chis(resname.strip(), chis, symmetry=False)
    rotamer_name = self.rot_id.identify(resname.strip(), wrap_chis)
    if (value is None):
      return None, None, None
    elif (value < 0.01):
      return 'OUTLIER', chis, value
    else:
      return rotamer_name, chis, value

  #{{{ get functions
  def get_outliers_count_and_fraction(self):
    if (self.numtotal != 0):
      fraction = float(self.numoutliers)/self.numtotal
      assert fraction <= 1.0
      return self.numoutliers, fraction
    return 0, 0.

  def get_outliers_goal(self):
    return "< 1%"

  def get_center(self, residue):
    for atom in residue.atoms():
      if atom.name == " CA ":
        return atom.xyz
    return None
  #}}}

  #{{{ construct_complete_sidechain
  def construct_complete_sidechain(self, residue_group):
    if (residue_group is not None):
      complete_dict = {}
      atom_dict = {}
      for ag in residue_group.atom_groups():
        for atom in ag.atoms():
          #if atom.name not in atom_dict:
          #handle hydrogen/deuterium swaps
          if atom_dict.get(atom.name) == None:
            if atom_dict.get(atom.name.replace("H","D",1)) != None:
              del(atom_dict[atom.name.replace("H","D",1)])
            elif atom_dict.get(atom.name.replace("D","H",1)) != None:
              del(atom_dict[atom.name.replace("D","H",1)])
          atom_dict[atom.name] = atom
        clone_dict = {}
        clone_dict.update(atom_dict)
        complete_dict[ag.altloc] = clone_dict
      if len(complete_dict) > 0:
        return complete_dict
    return {}
  #}}}


  #{{{ coot_todo
  def coot_todo(self, output_list):

    return ""
  #}}}

# XXX does this need to be smarter?
def get_occupancy (atom_group) :
  max_partial_occ = 0.
  for atom in atom_group.atoms() :
    if (atom.occ > max_partial_occ) and (atom.occ < 1) :
      max_partial_occ = atom.occ
  if (max_partial_occ == 0.) :
    return max([ atom.occ for atom in atom_group.atoms() ])
  else :
    return max_partial_occ
