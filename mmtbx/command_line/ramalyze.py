# (jEdit options) :folding=explicit:collapseFolds=1:
# LIBTBX_SET_DISPATCHER_NAME phenix.ramalyze

coot_script_header = """
def molprobity_fascinating_clusters_things_gui(window_name, sorting_option, cluster_list):

    ncluster_max = 75

    # a callback function
    def callback_recentre(widget, x, y, z):
        set_rotation_centre(x, y, z)

    # utility function
    def add_feature_buttons(feature_list, cluster_vbox):
        frame = gtk.Frame("Ramachandran Outliers")
        vbox = gtk.VBox(False, 0)
        cluster_vbox.pack_start(frame, False, False, 2)
        frame.add(vbox)

        # add buttons to vbox for each feature
        #
        for feature in feature_list:
            # print "feature: ", feature
            button = gtk.Button(feature[0])
            button.connect("clicked",
                           callback_recentre,
                           feature[4],
                           feature[5],
                           feature[6])
            vbox.pack_start(button, False, False, 1)

    # main body
    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    scrolled_win = gtk.ScrolledWindow()
    outside_vbox = gtk.VBox(False, 2)
    inside_vbox = gtk.VBox(False, 0)

    print "Maximum number of clusters displayed:  ", ncluster_max

    window.set_default_size(300, 200)
    window.set_title(window_name)
    inside_vbox.set_border_width(2)
    window.add(outside_vbox)
    outside_vbox.pack_start(scrolled_win, True, True, 0) # expand fill padding
    scrolled_win.add_with_viewport(inside_vbox)
    scrolled_win.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)

    count = 0

    for cluster_info in cluster_list:

        if (count == ncluster_max):
            break
        else:
            frame = gtk.Frame()
            vbox = gtk.VBox(False, 2)

            frame.set_border_width(6)
            frame.add(vbox)
            inside_vbox.pack_start(frame, False, False, 10)

            # now we have a list of individual features:
            features = cluster_info[0]
            if (len(features) > 0):
                add_feature_buttons(features, vbox)

    outside_vbox.set_border_width(2)
    ok_button = gtk.Button("  Close  ")
    outside_vbox.pack_end(ok_button, False, False, 0)
    ok_button.connect("clicked", lambda x: window.destroy())
    window.show_all()

molprobity_fascinating_clusters_things_gui(
    "MolProbity Multi-Chart",
    [],
    [[
      [
"""
import libtbx.load_env
import sys, os, getopt
try:
  from iotbx import pdb
  #from iotbx.pdb.hybrid_36 import hy36decode
except ImportError, e:
  print "iotbx not loaded"
  sys.exit()

#from mmtbx import monomer_library
#from mmtbx.monomer_library import server
#import mmtbx.monomer_library.pdb_interpretation
from mmtbx.rotamer.n_dim_table import NDimTable
from mmtbx.rotamer import rotamer_eval
from mmtbx.rotamer import ramachandran_eval
from cctbx import geometry_restraints
import iotbx
from phenix.autosol.UserMethods import GeneralMethods

#{{{ master_params (for using phenix built-in parameter parsing)
master_params="""
    ramalyze {
      pdb = None
        .type = path
        .help = '''Enter a PDB file name'''

      outliers_only = False
      .type = bool
      .help = '''Only print Ramachandran outliers'''

      changes = False
      .type = bool
      .help = '''Print change_log for ramalyze script'''

      version = False
      .type = bool
      .help = '''Print version'''
      
      verbose = True
      .type = bool
      .help = '''Verbose'''
}
""" 
#}}}
  
class ramalyze(GeneralMethods):
  
  #{{{ flag routines
  #flag routines-----------------------------------------------------------------------------------
  def changes(self):
    print """
    version 0.02 081022 - Updated for use with GUI
    version 0.1 080326 - First version
    """
  def version(self):
    print "\nversion 0.02 081022 - created 2007, Vincent Chen\n"
  def get_summary_and_header(self,command_name):
    header="\n"
    header+="\n#                       "+str(command_name)
    header+="\n#"
    header+="\n# Analyze protein backbone ramachandran"
    header+="\n# type phenix."+str(command_name)+": --help for help"

    summary= "usage: phenix.%s mypdb.pdb" % command_name
    return summary,header


  #------------------------------------------------------------------------------------------------
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
    from mmtbx.command_line.ramalyze import master_params
    master_params = iotbx.phil.parse(master_params, process_includes=True)
    args=self.get_keyword_table(args,out=out)       # set self.keyword_table

    command_name = "ramalyze"
    summary,header=self.get_summary_and_header(command_name)
    if not quiet: print>>out, header
  
    master_params,params,changed_params,help=self.get_params(
      command_name,master_params,args,out=out)
    if help or (params and params.ramalyze.verbose):
      print "Values of all params:" 
      master_params.format(python_object=params).show(out=out)
  
    if help or params is None: return None, None

    #print dir(params)
    #all_vars = dir(params.ramalyze)
    #for var in all_vars:
    #  if not "__" in var: print var,getattr(params.ramalyze,var)
    self.params=params # makes params available to whole class

    if self.params.ramalyze.changes:
      self.changes()
      return None, None

    if self.params.ramalyze.version:
      self.version()
      return None, None

    log=out
    if (log is None): log = sys.stdout
    filename = self.params.ramalyze.pdb
    #print 'filename', filename
    if filename and os.path.exists(filename):
      pdb_io = pdb.input(filename)
    else:
      print "Please enter a file name"
      return None, None
    
    output_text, output_list = self.analyze_pdb(pdb_io, self.params.ramalyze.outliers_only)
    if self.params.ramalyze.verbose:
      print output_text
    todo_list = self.coot_todo(output_list)
    return output_list, todo_list
  #}}}
  
  #{{{ analyze_pdb
  def analyze_pdb(self, pdb_io, outliers_only=None):
    hierarchy      = pdb_io.construct_hierarchy()
    analysis = ""
    output_list = []
    for model in hierarchy.models():
      for chain in model.chains():
        prevRes, prevC, resN, resCA, resO = None, None, None, None, None;
        for i_conformer, conformer in enumerate(chain.conformers()):
          #print str(conformer.residues()[-3].name)
          residues = conformer.residues()
          for i in range(len(residues)):
            residue = residues[i]
            phi = self.get_phi(residues, i)
            psi = self.get_psi(residues, i)
            coords = self.get_center(residue)
            resType = None
            if (residue.resname[0:3] == "GLY"):
              resType = "glycine"
            elif (residue.resname[0:3] == "PRO"):
              resType = "proline"
            elif (self.isPrePro(residues, i)):
              resType = "prepro"
            else:
              resType = "general"
  
            r = ramachandran_eval.RamachandranEval()
            value = 0
            if (phi is not None and psi is not None):
              value = r.evaluate(resType, [phi, psi])
              ramaType = self.evaluateScore(resType, value)
              #print params_old["outliersonly"]
              if (not outliers_only or self.isOutlier(resType, value)):
                analysis += '%s%4s %s:%.2f:%.2f:%.2f:%s:%s\n' % (chain.id,residue.resseq,residue.resname,value*100,phi,psi,ramaType,resType.capitalize())

                output_list.append([chain.id,residue.resseq,residue.resname,value*100,phi,psi,ramaType,resType.capitalize(),coords])
                
#print str(residue.seq).rjust(4,' ') + " " + residue.name + ":" + str(value) + ":" + str(phi) + ":" + str(psi) + ":OUTLIER:" + resType.capitalize()
                #print type(residue.seq)
    #print self.analysis_list
    return analysis.rstrip(), output_list
  #}}}
  
  def coot_todo(self, output_list):
    #print coot_script_header
    text=coot_script_header
    for chain_id,resnum,resname,rama_value,phi,psi,ramaType,resType,coords in output_list:
       button='       ["Ramachandran Outlier at %s%s %s (%.2f)", 0, 1, 0, %f, %f, %f],\n' %(chain_id, resnum, resname, rama_value, coords[0], coords[1], coords[2])
       #print button 
       text+=button
    text+="      ]\n"
    text+="     ]\n"
    text+="    ])\n"
    #print text
    return text

  #{{{ get_phi
  def get_phi(self, residues, i):
    prevC, resN, resCA, resC = None, None, None, None;
    if (i > 0 and i < len(residues)):
      prev = residues[i-1]
      if (prev is not None):
        for atom in prev.atoms():
          if (atom.name == " C  "): prevC = atom
      res = residues[i]
      if (res is not None):
        for atom in res.atoms():
          if (atom.name == " N  "): resN = atom
          if (atom.name == " CA "): resCA = atom
          if (atom.name == " C  "): resC = atom
      if (prevC is not None and resN is not None and resCA is not None and resC is not None):
        b = geometry_restraints.bond(
          sites=[prevC.xyz,resN.xyz],
          distance_ideal=1,
          weight=1)
        # check to see if residues are actually bonded.
        if (b.distance_model > 4): return None
        d = geometry_restraints.dihedral(
          sites=[prevC.xyz,resN.xyz,resCA.xyz,resC.xyz],
          angle_ideal=-40,
          weight=1)
        return d.angle_model
  #}}}
  
  #{{{ get_psi
  def get_psi(self, residues, i):
    resN, resCA, resC, nextN = None, None, None, None;
    if (i >= 0 and i < len(residues) - 1):
      next = residues[i+1]
      if (next is not None):
        for atom in next.atoms():
          if (atom.name == " N  "): nextN = atom
      res = residues[i]
      if (res is not None):
        for atom in res.atoms():
          if (atom.name == " N  "): resN = atom
          if (atom.name == " CA "): resCA = atom
          if (atom.name == " C  "): resC = atom
      if (nextN is not None and resN is not None and resCA is not None and resC is not None):
        b = geometry_restraints.bond(
          sites=[resC.xyz,nextN.xyz],
          distance_ideal=1,
          weight=1)
        if (b.distance_model > 4): return None
        d = geometry_restraints.dihedral(
          sites=[resN.xyz,resCA.xyz,resC.xyz,nextN.xyz],
          angle_ideal=-40,
          weight=1)
        return d.angle_model
  #}}}
 
  #{{{ get_center
  def get_center(self, residue):
    coords = None
    
    for atom in residue.atoms():
      if (atom.name == " CA "): 
        coords = atom.xyz
    return coords
  #}}}

  #{{{ isPrePro
  def isPrePro(self, residues, i):
    if (i < 0 or i >= len(residues) - 1): return False
    else:
      next = residues[i+1]
      if (next.resname[0:3] == "PRO"): return True
    return False
  #}}}
  
  #{{{ isOutlier
  def isOutlier(self, resType, value):
    if (resType == "general"):
      if (value < 0.0005): return True
      else: return False
    else:
      if (value < 0.002): return True
      else: return False
  #}}}
  
  #{{{ evaluateScore
  def evaluateScore(self, resType, value):
    if (value >= 0.02): return "Favored"
    if (resType == "general"):
      if (value >= 0.0005): return "Allowed"
      else: return "OUTLIER"
    else:
      if (value >= 0.0020): return "Allowed"
      else: return "OUTLIER"
  #}}}
  
  if __name__ == "__main__":
    #params_old = {}
    #params_old["outliersonly"]=False
    #args = parse_cmdline(params_old)
    from mmtbx.command_line.ramalyze import ramalyze
    r = ramalyze()
    output_list, coot_todo_list = r.run(sys.argv[1:])
    #for entry in output_list:
    #  print entry
