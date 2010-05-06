# (jEdit options) :folding=explicit:collapseFolds=1:

#{{{ coot_script_header
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
#}}}

import libtbx.load_env
import sys, os, getopt
try:
  from iotbx import pdb
except ImportError, e:
  print "iotbx not loaded"
  sys.exit()

from mmtbx.rotamer.n_dim_table import NDimTable
from mmtbx.rotamer import rotamer_eval
from mmtbx.rotamer import ramachandran_eval
from cctbx import geometry_restraints
import iotbx.phil

def get_master_phil():
  return iotbx.phil.parse(
    input_string="""
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
""")

class ramalyze():

  #{{{ flag routines
  #flag routines-----------------------------------------------------------------------------------
  def changes(self):
    print """
    version 0.05 090605 - Added summary and description to text output.
                          New functions to obtain percentages and goals.
    version 0.04 081216 - Added fraction output to get functions
    version 0.03 081204 - Fixes for better altconf behaviour
    version 0.02 081022 - Updated for use with GUI
    version 0.1 080326 - First version
    """
  def version(self):
    print "\nversion 0.05 090605 - created 2007, Vincent Chen\n"
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

  #{{{ run
  def run(self, args, out=sys.stdout, quiet=False):
    master_phil = get_master_phil()
    import iotbx.utils
    input_objects = iotbx.utils.process_command_line_inputs(
      args=args,
      master_phil=master_phil,
      input_types=("pdb",))
    work_phil = master_phil.fetch(sources=input_objects["phil"])
    work_params = work_phil.extract()
    assert len(input_objects["pdb"]) == 1
    file_obj = input_objects["pdb"][0]
    filename = file_obj.file_name

    command_name = "ramalyze"
    summary,header=self.get_summary_and_header(command_name)
    if not quiet: print >>out, header

    #TO DO:  make this a working help section
    #if help or (params and params.ramalyze.verbose):
    #  pass
      # XXX: disabled for GUI
      #print "Values of all params:"
      #master_params.format(python_object=params).show(out=out)

    self.params=work_params # makes params available to whole class

    if self.params.ramalyze.changes:
      self.changes()
      return None, None

    if self.params.ramalyze.version:
      self.version()
      return None, None

    log=out
    if (log is None): log = sys.stdout
    if filename and os.path.exists(filename):
      pdb_io = pdb.input(filename)
    else:
      print "Please enter a file name"
      return None, None

    output_text, output_list = self.analyze_pdb(pdb_io,
      outliers_only=self.params.ramalyze.outliers_only)
    out_count, out_percent = self.get_outliers_count_and_fraction()
    fav_count, fav_percent = self.get_favored_count_and_fraction()
    if self.params.ramalyze.verbose:
      print >> out, "residue:score%:phi:psi:evaluation:type"
      print >> out, output_text
      print >> out, 'SUMMARY: %.2f%% outliers (Goal: %s)' % \
        (out_percent*100, self.get_outliers_goal())
      print >> out, 'SUMMARY: %.2f%% favored (Goal: %s)' % \
        (fav_percent*100, self.get_favored_goal())
    todo_list = self.coot_todo(output_list)
    self.out_percent = out_percent * 100.0
    self.fav_percent = fav_percent * 100.0
    return output_list, todo_list
  #}}}

  #{{{ analyze_pdb
  def analyze_pdb(self, pdb_io=None, hierarchy=None, outliers_only=None):
    assert [pdb_io, hierarchy].count(None) == 1
    if(pdb_io is not None):
      hierarchy = pdb_io.construct_hierarchy()
    analysis = ""
    output_list = []
    self.numoutliers = 0
    self.numallowed = 0
    self.numfavored = 0
    self.numgen = 0
    self.numgly = 0
    self.numpro = 0
    self.numprepro = 0
    self.numtotal = 0
    for model in hierarchy.models():
      for chain in model.chains():
        #prevRes, prevC, resN, resCA, resO = None, None, None, None, None;
        #help(chain)
        residues = list(chain.residue_groups())
        #help(residues)
        for i, residue_group in enumerate(residues):
            #The reason I pass lists of atom_groups to get_phi and get_psi is to deal with the
            #particular issue where some residues have an A alt conf that needs some atoms from
            #a "" alt conf to get calculated correctly.
            #See 1jxt.pdb for examples.  This way I can search both the alt conf atoms and
            #the "" atoms if necessary.
            prev_rezes, next_rezes, prev_atom_list, next_atom_list, atom_list = \
              None, None, None, None, None
            if (i > 0):
              prev_rezes = self.construct_complete_residues(residues[i-1])
            rezes = self.construct_complete_residues(residues[i])
            if (i < len(residues)-1):
              next_rezes = self.construct_complete_residues(residues[i+1])
            #for alt_conf in sorted(rezes.keys()):
            for atom_group in residue_group.atom_groups():
              alt_conf = atom_group.altloc
              if rezes is not None:
                atom_list = rezes.get(alt_conf)
              if prev_rezes is not None:
                prev_atom_list = prev_rezes.get(alt_conf)
                if (prev_atom_list is None):
                  prev_keys = sorted(prev_rezes.keys())
                  prev_atom_list = prev_rezes.get(prev_keys[0])
                  #print prev_atom_list
              if next_rezes is not None:
                next_atom_list = next_rezes.get(alt_conf)
                if (next_atom_list is None):
                  next_keys = sorted(next_rezes.keys())
                  next_atom_list = next_rezes.get(next_keys[0])
              phi = self.get_phi(prev_atom_list, atom_list)
              psi = self.get_psi(atom_list, next_atom_list)
              coords = self.get_center(atom_group)

              if (phi is not None and psi is not None):
                resType = None
                self.numtotal += 1
                if (atom_group.resname[0:3] == "GLY"):
                  resType = "glycine"
                  self.numgly += 1
                elif (atom_group.resname[0:3] == "PRO"):
                  resType = "proline"
                  self.numpro += 1
                elif (self.isPrePro(residues, i)):
                  resType = "prepro"
                  self.numprepro += 1
                else:
                  resType = "general"
                  self.numgen += 1

                r = ramachandran_eval.RamachandranEval()
                #value = 0

                value = r.evaluate(resType, [phi, psi])
                ramaType = self.evaluateScore(resType, value)
                #print params_old["outliersonly"]
                if (not outliers_only or self.isOutlier(resType, value)):
                  analysis += '%s%4s %s%s:%.2f:%.2f:%.2f:%s:%s\n' % \
                    (chain.id,
                     residue_group.resseq,atom_group.altloc,
                     atom_group.resname,
                     value*100,
                     phi,
                     psi,
                     ramaType,
                     resType.capitalize())

                  output_list.append([chain.id,
                                      residue_group.resseq,
                                      atom_group.altloc+atom_group.resname,
                                      value*100,
                                      phi,
                                      psi,
                                      ramaType,
                                      resType.capitalize(),
                                      coords])

    return analysis.rstrip(), output_list
  #}}}

  #{{{ coot_todo
  def coot_todo(self, output_list):
    #print coot_script_header
    text=coot_script_header
    for chain_id,resnum,resname,rama_value,phi,psi,ramaType,resType,coords in output_list:
       if (coords is not None):
         button='       ["Ramachandran Outlier at %s%s %s (%.2f)", 0, 1, 0, %f, %f, %f],\n' % \
           (chain_id, resnum, resname, rama_value, coords[0], coords[1], coords[2])
         #print button
         text+=button
    text+="      ]\n"
    text+="     ]\n"
    text+="    ])\n"
    #print text
    return text
  #}}}

  #{{{ get_matching_atom_group
  def get_matching_atom_group(self, residue_group, altloc):
    match = None
    if (residue_group != None):
      for ag in residue_group.atom_groups():
        if (ag.altloc == "" and match == None): match = ag
        if (ag.altloc == altloc): match = ag
    return match
  #}}}

  #{{{ get_phi
  def get_phi(self, prev_atoms, atoms):
    prevC, resN, resCA, resC = None, None, None, None;
    if (prev_atoms is not None):
      for atom in prev_atoms:
        if (atom.name == " C  "): prevC = atom
    if (atoms is not None):
      for atom in atoms:
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
  def get_psi(self, atoms, next_atoms):
    resN, resCA, resC, nextN = None, None, None, None
    if (next_atoms is not None):
      for atom in next_atoms:
        if (atom.name == " N  "): nextN = atom
    if (atoms is not None):
      for atom in atoms:
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

  #{{{ construct_complete_residues
  def construct_complete_residues(self, res_group):
    if (res_group is not None):
      complete_dict = {}
      nit, ca, co, oxy = None, None, None, None
      for ag in res_group.atom_groups():
        changed = False
        for atom in ag.atoms():
          if (atom.name == " N  "): nit = atom
          if (atom.name == " CA "): ca = atom
          if (atom.name == " C  "): co = atom
          if (atom.name == " O  "): oxy = atom
          if (atom.name == " N  " or
              atom.name == " CA " or
              atom.name == " C  " or
              atom.name == " O  "):
            changed = True
        if (nit is not None and ca is not None and co is not None and oxy is not None and changed):
          # complete residue backbone found
          complete_dict[ag.altloc] = [nit, ca, co, oxy]
      if len(complete_dict) > 0:
        return complete_dict
    return None
  #}}}

  #{{{ get_center
  def get_center(self, ag):
    coords = None

    for atom in ag.atoms():
      if (atom.name == " CA "):
        coords = atom.xyz
    return coords
  #}}}

  #{{{ isPrePro
  def isPrePro(self, residues, i):
    if (i < 0 or i >= len(residues) - 1): return False
    else:
      next = residues[i+1]
      for ag in next.atom_groups():
        if (ag.resname[0:3] == "PRO"): return True
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
    if (value >= 0.02):
      self.numfavored += 1
      return "Favored"
    if (resType == "general"):
      if (value >= 0.0005):
        self.numallowed += 1
        return "Allowed"
      else:
        self.numoutliers += 1
        return "OUTLIER"
    else:
      if (value >= 0.0020):
        self.numallowed += 1
        return "Allowed"
      else:
        self.numoutliers += 1
        return "OUTLIER"
  #}}}

  #{{{ get functions
  def get_outliers_count_and_fraction(self):
    if (self.numtotal != 0):
      fraction = float(self.numoutliers)/self.numtotal
      assert fraction <= 1.0
      return self.numoutliers, fraction
    return 0, 0.

  def get_outliers_goal(self):
    return "< 0.2%"

  def get_allowed_count_and_fraction(self):
    if (self.numtotal != 0):
      fraction = float(self.numallowed)/self.numtotal
      assert fraction <= 1.0
      return self.numallowed, fraction
    return 0, 0.

  def get_allowed_goal(self):
    return "> 99.8%"

  def get_favored_count_and_fraction(self):
    if (self.numtotal != 0):
      fraction = float(self.numfavored)/self.numtotal
      assert fraction <= 1.0
      return self.numfavored, fraction
    return 0, 0.

  def get_favored_goal(self):
    return "> 98%"

  def get_general_count_and_fraction(self):
    if (self.numtotal != 0):
      fraction = (float(self.numgen)/self.numtotal)
      assert fraction <= 1.0
      return self.numgen, fraction
    return 0, 0.

  def get_gly_count_and_fraction(self):
    if (self.numtotal != 0):
      fraction = (float(self.numgly)/self.numtotal)
      assert fraction <= 1.0
      return self.numgly, fraction
    return 0, 0.

  def get_pro_count_and_fraction(self):
    if (self.numtotal != 0):
      fraction = (float(self.numpro)/self.numtotal)
      assert fraction <= 1.0
      return self.numpro, fraction
    return 0, 0.

  def get_prepro_count_and_fraction(self):
    if (self.numtotal != 0):
      fraction = (float(self.numprepro)/self.numtotal)
      assert fraction <= 1.0
      return self.numprepro, fraction
    return 0, 0.

  def get_phi_psi_residues_count(self):
    # n.b. this function returns the number of residues that have a valid phi/psi pair.
    return self.numtotal
  #}}}
