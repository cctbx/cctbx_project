from __future__ import division
import iotbx.phil
from scitbx.matrix import col, dihedral_angle, rotate_point_around_axis
from libtbx.utils import Usage
import os
import sys

def get_master_phil():
  return iotbx.phil.parse(
    input_string="""
      cbetadev {
        pdb = None
        .type = path
        .help = '''Enter a PDB file name'''

        outliers_only = False
        .type = bool
        .help = '''Only print outliers'''

        changes = False
        .type = bool
        .help = '''Print list of changes'''

        version = False
        .type = bool
        .help = '''Print version'''

        verbose = True
        .type = bool
        .help = '''Verbose'''
  }

    """)

class cbetadev(object):
  #flag routines-----------------------------------------------------------------------------------
  def usage(self):
    return """
phenix.cbetadev file.pdb [params.eff] [options ...]

Options:

  pdb=input_file        input PDB file
  outliers_only=False   only print outliers
  verbose=True          verbose text output

Example:

  phenix.cbetadev pdb=1ubq.pdb outliers_only=True

"""
  def changes(self):
    print "\nversion 0.10 080504 Initial version\nversion 0.12 090604 Added Summary line\n"
  def version(self):
    print "\nversion 0.11 090604 - Copyright 2008,2009 Jeffrey J. Headd\n"
  def get_summary_and_header(self,command_name):
    header="\n"
    header+="\n#                       "+str(command_name)
    header+="\n#"
    header+="\n# Analyze protein sidechain C-beta deviation"
    header+="\n# type phenix."+str(command_name)+": --help for help"

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
      summary, header = self.get_summary_and_header("cbetadev")
      raise Usage(summary)
    file_obj = input_objects["pdb"][0]
    filename = file_obj.file_name

    command_name = "cbetadev"
    summary,header=self.get_summary_and_header(command_name)
    if not quiet: print >>out, header

    #TO DO: make useful help section
    #if help : #or (params and params.cbetadev.verbose):
    #  print >> out, summary
      #print "Values of all params:"
      #master_params.format(python_object=params).show(out=out)

    self.params=work_params # makes params available to whole class

    if self.params.cbetadev.changes:
      self.changes()
      return

    if self.params.cbetadev.version:
      self.version()
      return

    log=out
    if (log is None): log = sys.stdout
    if self.params.cbetadev.verbose :
      print >> out, 'filename', filename
    if filename and os.path.exists(filename):
      pdb_io = iotbx.pdb.input(filename)
    else:
      print "Please enter a file name"
      return
    output_text, summary_line, output_list = self.analyze_pdb(
                                               filename=filename,
                                               pdb_io=pdb_io,
                                               outliers_only=self.params.cbetadev.outliers_only)
    if not quiet :
      print >> out, output_text
      print >> out, summary_line
    return output_list
  #}}}

  #{{{ analyze_pdb
  def analyze_pdb(self, filename=None, pdb_io=None, hierarchy=None,
                  outliers_only=False):
    relevant_atom_names = {
      " CA ": None, " N  ": None, " C  ": None, " CB ": None} # FUTURE: set
    analysis = 'pdb:alt:res:chainID:resnum:dev:dihedralNABB:Occ:ALT:\n'
    output_list = []
    cbetadev_ctr = 0
    self.num_outliers = 0
    self.expected_outliers = 0
    self.outliers = ''
    self.summary = ''
    self.beta_ideal = {}
    assert [pdb_io, hierarchy].count(None) == 1
    if(pdb_io is not None):
      hierarchy = pdb_io.construct_hierarchy()
    for model in hierarchy.models():
      for chain in model.chains():
        for rg in chain.residue_groups():
          for i_cf,cf in enumerate(rg.conformers()):
            for i_residue,residue in enumerate(cf.residues()):
              if (residue.resname == "GLY") :
                continue
              is_first = (i_cf == 0)
              is_alt_conf = False
              relevant_atoms = {}
              for atom in residue.atoms():
                if (atom.name in relevant_atom_names):
                  relevant_atoms[atom.name] = atom
                  if (len(atom.parent().altloc) != 0):
                    is_alt_conf = True
              if ((is_first or is_alt_conf) and len(relevant_atoms) == 4):
                result = calculate_ideal_and_deviation(
                  relevant_atoms=relevant_atoms,
                  resname=residue.resname)
                dev = result.deviation
                dihedralNABB = result.dihedral
                betaxyz = result.ideal
                if (dev is None) : continue
                if(dev >=0.25 or outliers_only==False):
                  if(dev >=0.25):
                    cbetadev_ctr+=1
                    self.num_outliers+=1
                  PDBfileStr = ""
                  if(filename is not None):
                    PDBfileStr = os.path.basename(filename)[:-4]
                  if (is_alt_conf):
                    altchar = cf.altloc.lower()
                  else:
                    altchar = " "
                  res=residue.resname.lower()
                  sub=chain.id
                  if(len(sub)==1):
                    sub=" "+sub
                  resnum=residue.resid()
                  resCB = relevant_atoms[" CB "]
                  occ = resCB.occ
                  analysis += '%s :%s:%s:%s:%5s:%7.3f:%7.2f:%7.2f:%s:\n' % \
                    (PDBfileStr,altchar,res,sub,resnum,dev,dihedralNABB,occ,
                     altchar)
                  if (dev >= 0.25):
                    self.outliers += \
                      '%s :%s:%s:%s:%5s:%7.3f:%7.2f:%7.2f:%s:\n' % \
                      (PDBfileStr,altchar,res,sub,resnum,dev,dihedralNABB,occ,
                       altchar)
                  key = altchar+res+sub+resnum
                  self.beta_ideal[key] = betaxyz
                  output_list.append([PDBfileStr,
                                      altchar,
                                      res,
                                      sub,
                                      resnum[0:4],
                                      resnum[4:],
                                      dev,
                                      dihedralNABB,
                                      occ,
                                      altchar,
                                      resCB.xyz])
    summary = 'SUMMARY: %d C-beta Deviation >= 0.25 Angstrom (Goal: 0)' % \
      (cbetadev_ctr)
    self.summary = 'SUMMARY: %d C-beta Deviation >= 0.25 Angstrom (Goal: 0)' %\
      (cbetadev_ctr)
    self.analysis = analysis.rstrip()
    return self.analysis, summary, output_list

  def get_outlier_count(self):
    return self.num_outliers

  def get_expected_count(self):
    return self.expected_outliers

  def get_outliers(self):
    return self.outliers

  def get_summary(self):
    return self.summary

  def get_beta_ideal(self):
    return self.beta_ideal

class calculate_ideal_and_deviation (object) :
  __slots__ = ["deviation", "ideal", "dihedral"]
  def __init__ (self, relevant_atoms, resname) :
    assert (resname != "GLY")
    self.deviation = None
    self.ideal = None
    self.dihedral = None
    resCA = relevant_atoms[" CA "]
    resN  = relevant_atoms[" N  "]
    resC  = relevant_atoms[" C  "]
    resCB = relevant_atoms[" CB "]
    dist, angleCAB, dihedralNCAB, angleNAB, dihedralCNAB, angleideal= \
      idealized_calpha_angles(resname)
    betaNCAB = construct_fourth(resN,
                                resCA,
                                resC,
                                dist,
                                angleCAB,
                                dihedralNCAB,
                                method="NCAB")
    betaCNAB = construct_fourth(resN,
                                resCA,
                                resC,
                                dist,
                                angleNAB,
                                dihedralCNAB,
                                method="CNAB")
    if (not None in [betaNCAB, betaCNAB]) :
      betaxyz = (col(betaNCAB) + col(betaCNAB)) / 2
      betadist = abs(col(resCA.xyz) - betaxyz)
      if betadist != 0:
        if(betadist != dist):
          distTemp = betaxyz - col(resCA.xyz)
          betaxyz = col(resCA.xyz) + distTemp * dist/betadist
        self.deviation = abs(col(resCB.xyz) - betaxyz)
        self.dihedral = dihedral_angle(
          sites=[resN.xyz,resCA.xyz,betaxyz.elems,resCB.xyz], deg=True)
        self.ideal = betaxyz.elems

def idealized_calpha_angles(resname):
  if (resname == "ALA"):
    dist = 1.536
    angleCAB = 110.1
    dihedralNCAB = 122.9
    angleNAB = 110.6
    dihedralCNAB = -122.6
    angleideal = 111.2
  elif (resname == "PRO"):
    dist = 1.530
    angleCAB = 112.2
    dihedralNCAB = 115.1
    angleNAB = 103.0
    dihedralCNAB = -120.7
    angleideal = 111.8
  elif (resname in ["VAL", "THR", "ILE"]) :
    dist = 1.540
    angleCAB = 109.1
    dihedralNCAB = 123.4
    angleNAB = 111.5
    dihedralCNAB = -122.0
    angleideal = 111.2
  elif (resname == "GLY"):
    dist = 1.10
    angleCAB = 109.3
    dihedralNCAB = 121.6
    angleNAB = 109.3
    dihedralCNAB = -121.6
    angleideal = 112.5
  else:
    dist = 1.530
    angleCAB = 110.1
    dihedralNCAB = 122.8
    angleNAB = 110.5
    dihedralCNAB = -122.6
    angleideal = 111.2
  return dist, angleCAB, dihedralNCAB, angleNAB, dihedralCNAB, angleideal

def construct_fourth(resN,resCA,resC,dist,angle,dihedral,method="NCAB"):
  if (not None in [resN, resCA, resC]) :
    if (method == "NCAB"):
      res0 = resN
      res1 = resC
      res2 = resCA
    elif (method == "CNAB"):
      res0 = resC
      res1 = resN
      res2 = resCA
    a = col(res2.xyz) - col(res1.xyz)
    b = col(res0.xyz) - col(res1.xyz)
    c = a.cross(b)
    cmag = abs(c)
    if(cmag > 0.000001):
      c *= dist/cmag
    c += col(res2.xyz)
    d = c
    angledhdrl = dihedral - 90
    a = col(res1.xyz)
    b = col(res2.xyz)
    # XXX is there an equivalent method for 'col'?
    newD = col(rotate_point_around_axis(
      axis_point_1=res1.xyz,
      axis_point_2=res2.xyz,
      point=d.elems,
      angle=angledhdrl,
      deg=True))
    a = newD - col(res2.xyz)
    b = col(res1.xyz) - col(res2.xyz)
    c = a.cross(b)
    cmag = abs(c)
    if(cmag > 0.000001):
      c *= dist/cmag
    angledhdrl = 90 - angle;
    a = col(res2.xyz)
    c += a
    b = c
    return rotate_point_around_axis(
      axis_point_1=a.elems,
      axis_point_2=b.elems,
      point=newD.elems,
      angle=angledhdrl,
      deg=True)
  return None

def extract_atoms_from_residue_group (residue_group) :
  """
  Given a residue_group object, which may or may not have multiple
  conformations, extract the relevant atoms for each conformer, taking into
  account any atoms shared between conformers.  This is implemented
  separately from the main validation routine, which accesses the hierarchy
  object via the chain->conformer->residue API.  Returns a list of hashes,
  each suitable for calling calculate_ideal_and_deviation.
  """
  atom_groups = residue_group.atom_groups()
  if (len(atom_groups) == 1) :
    relevant_atoms = {}
    for atom in atom_groups[0].atoms() :
      relevant_atoms[atom.name] = atom
    return [ relevant_atoms ]
  else :
    all_relevant_atoms = []
    expected_names = [" CA ", " N  ", " CB ", " C  "]
    main_conf = {}
    for atom_group in atom_groups :
      if (atom_group.altloc.strip() == '') :
        for atom in atom_group.atoms() :
          if (atom.name in expected_names) :
            main_conf[atom.name] = atom
      else :
        relevant_atoms = {}
        for atom in atom_group.atoms() :
          if (atom.name in expected_names) :
            relevant_atoms[atom.name] = atom
        if (len(relevant_atoms) == 0) : continue
        for atom_name in main_conf.keys() :
          if (not atom_name in relevant_atoms) :
            relevant_atoms[atom_name] = main_conf[atom_name]
        if (len(relevant_atoms) != 0) :
          all_relevant_atoms.append(relevant_atoms)
    if (len(main_conf) == 4) :
      all_relevant_atoms.insert(0, main_conf)
    return all_relevant_atoms
