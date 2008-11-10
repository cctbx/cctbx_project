#(jEdit options) :folding=explicit:collapseFolds=1:
# LIBTBX_SET_DISPATCHER_NAME phenix.rotalyze

import libtbx.load_env
import sys, os, getopt
#try:
from iotbx import pdb
  #from iotbx.pdb.hybrid_36 import hy36decode
#except ImportError, e:
#  print "iotbx not loaded"
#  sys.exit()

from cctbx import geometry_restraints
from mmtbx.rotamer.sidechain_angles import SidechainAngles
from mmtbx.rotamer import rotamer_eval
from mmtbx.rotamer.rotamer_eval import RotamerID
import iotbx
from phenix.autosol.UserMethods import GeneralMethods

#{{{ master_params (for using phenix built-in parameter parsing)
master_params="""
    rotalyze {
      pdb = None
        .type = path
        .help = '''Enter a PDB file name'''

      outliers_only = False
      .type = bool
      .help = '''Only print rotamer outliers'''

      changes = False
      .type = bool
      .help = '''Print change_log for rotalyze script'''

      version = False
      .type = bool
      .help = '''Print version'''
      
      verbose = True
      .type = bool
      .help = '''Verbose'''
      
      show_errors = False
      .type = bool
      .help = '''Print out errors'''
}
""" 
#}}}

class rotalyze(GeneralMethods):
  
  #{{{ flag routines
  #flag routines-----------------------------------------------------------------------------------
  def changes(self):
    print """
    version 0.02 081022 - Updated for use with GUI
    version 0.1  080326 - First version
    """
  def version(self):
    print "\nversion 0.02 081022 - created 2007, Vincent Chen\n"
  #------------------------------------------------------------------------------------------------
  #}}}
  
  #{{{ get_summary_and_header
  def get_summary_and_header(self,command_name):
    header="\n"
    header+="\n#                       "+str(command_name)
    header+="\n#"
    header+="\n# Analyze protein sidechain rotamers"
    header+="\n# type phenix."+str(command_name)+": --help for help"

    summary= "usage: phenix.%s mypdb.pdb" % command_name
    return summary,header
  #}}}
  
  #{{{ special_cases
  def special_cases(self, args):
    # special cases for input files so user doesn't need to specify:
    new_args=[]
    for arg in args:
      # special cases for input files so user doesn't need to specify:
      if (os.path.isfile(arg)):
        if arg[-3:]=='pdb':
          arg_use='pdb='+arg
        else:
          arg_use=arg
      else:
        arg_use=arg
      new_args.append(arg_use)
    return new_args
  #}}}
    
  #{{{ run
  def run(self, args, out=sys.stdout, quiet=False):
    args=self.special_cases(args)
    from mmtbx.command_line.rotalyze import master_params
    master_params = iotbx.phil.parse(master_params, process_includes=True)
    args=self.get_keyword_table(args,out=out)       # set self.keyword_table

    command_name = "rotalyze"
    summary,header=self.get_summary_and_header(command_name)
    if not quiet: print>>out, header
  
    master_params,params,changed_params,help=self.get_params(
      command_name,master_params,args,out=out)
      
    if help or (params and params.rotalyze.verbose):
      print "Values of all params:" 
      master_params.format(python_object=params).show(out=out)
  
    if help or params is None: return None, None

    #print dir(params)
    #all_vars = dir(params.rotalyze)
    #for var in all_vars:
    #  if not "__" in var: print var,getattr(params.rotalyze,var)
    self.params=params # makes params available to whole class

    if self.params.rotalyze.changes:
      self.changes()
      return None, None

    if self.params.rotalyze.version:
      self.version()
      return None, None

    log=out
    if (log is None): log = sys.stdout
    filename = self.params.rotalyze.pdb
    #print 'filename', filename
    if filename and os.path.exists(filename):
      pdb_io = pdb.input(filename)
    else:
      print "Please enter a file name"
      return None, None
      
    output_text, output_list = self.analyze_pdb(pdb_io, self.params.rotalyze.outliers_only, self.params.rotalyze.show_errors)
    if self.params.rotalyze.verbose:
      print output_text
    todo_list = self.coot_todo(output_list)
    return output_list, todo_list
  #}}}
  
  #{{{ analyze_pdb
  def analyze_pdb(self, pdb_io, outliers_only=False, show_errors = False):
    sa = SidechainAngles()
    hierarchy      = pdb_io.construct_hierarchy()
    analysis = ""
    output_list = []
    rot_id = rotamer_eval.RotamerID() # loads in the rotamer names
    #print rot_id.names
    for model in hierarchy.models():
      for chain in model.chains():
        for i_conformer, conformer in enumerate(chain.conformers()):
          #print str(conformer.residues()[-3].name)
          residues = conformer.residues()
          for i in range(len(residues)):
            residue = residues[i]
            try:
              coords = self.get_center(residue)
              chis = sa.measureChiAngles(residue)
              r = rotamer_eval.RotamerEval()
              if (chis is not None):
                value = r.evaluate(residue.resname.lower().strip(), chis)
                #print chain.id + str(residue.seq) + residue.name.strip() + str(value)
                if value != None:
                  s = '%s%4s %s:%.1f' % (chain.id,residue.resseq,residue.resname.strip(),value*100)
                  res_out_list = [chain.id,residue.resseq,residue.resname,value*100]
                  #wrap_chis = []
                  #for i in range(len(chis)):
                  #  wrap_chis.append(chis[i] % 360)
                  #  if wrap_chis[i] < 0: wrap_chis[i] += 360
                  wrap_chis = rot_id.wrap_chis(residue.resname.strip(), chis)
                  for i in range(4):
                    s += ':'
                    if i < len(wrap_chis):
                      s += '%.1f' % (wrap_chis[i])
                      res_out_list.append(wrap_chis[i])
                  s += ':'
                  if value < 0.01:
                    s += "OUTLIER\n"
                    res_out_list.append("OUTLIER")
                    res_out_list.append(coords)
                    if outliers_only: 
                      analysis += s
                      output_list.append(res_out_list)
                  else:
                    s += rot_id.identify(residue.resname, wrap_chis) + "\n"
                    res_out_list.append(rot_id.identify(residue.resname, wrap_chis))
                    res_out_list.append(coords)
                  if not outliers_only: 
                    analysis += s
                    output_list.append(res_out_list)
            except AttributeError:
              if show_errors: print '%s%4s %s is missing some sidechain atoms' % (chain.id, residue.resseq, residue.resname)
    return analysis.rstrip(), output_list
  #}}}
  
  #{{{ get_center
  def get_center(self, residue):
    coords = None
    
    for atom in residue.atoms():
      if (atom.name == " CA "): 
        coords = atom.xyz
    return coords
  #}}}
  
  #{{{ coot_todo
  def coot_todo(self, output_list):
    
    return ""
  #}}}
  
  if __name__ == "__main__":
    from mmtbx.command_line.rotalyze import rotalyze
    r = rotalyze()
    output_list, coot_todo_list = r.run(sys.argv[1:])
    #for entry in output_list:
    #  print entry
