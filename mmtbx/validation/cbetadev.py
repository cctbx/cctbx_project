import sys, os, math
from cctbx import geometry_restraints
import iotbx.phil
from libtbx.utils import Usage

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
  #{{{ local_help
  #flag routines-----------------------------------------------------------------------------------
  def local_help(self):
    version()
    print """USAGE:  phenix.cbetadev file.pdb

  FLAGS:
    -h    Print this help message
    -v    Display version information
    -c    Print a change log
    -o    Only print outliers (deviation >= 0.25)
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
              is_first = (i_cf == 0)
              is_alt_conf = False
              relevant_atoms = {}
              for atom in residue.atoms():
                if (atom.name in relevant_atom_names):
                  relevant_atoms[atom.name] = atom
                  if (len(atom.parent().altloc) != 0):
                    is_alt_conf = True
              if ((is_first or is_alt_conf) and len(relevant_atoms) == 4):
                resCA = relevant_atoms[" CA "]
                resN  = relevant_atoms[" N  "]
                resC  = relevant_atoms[" C  "]
                resCB = relevant_atoms[" CB "]
                dist, angleCAB, dihedralNCAB, angleNAB, dihedralCNAB, angleideal= \
                  self.idealized_calpha_angles(residue=residue)
                betaNCAB = self.construct_fourth(resN,
                                                 resCA,
                                                 resC,
                                                 dist,
                                                 angleCAB,
                                                 dihedralNCAB,
                                                 method="NCAB")
                betaCNAB = self.construct_fourth(resN,
                                                 resCA,
                                                 resC,
                                                 dist,
                                                 angleNAB,
                                                 dihedralCNAB,
                                                 method="CNAB")
                betaxyz = [(betaNCAB[0]+betaCNAB[0])/2,
                           (betaNCAB[1]+betaCNAB[1])/2,
                           (betaNCAB[2]+betaCNAB[2])/2]
                betadist = self.distance(resCA.xyz,betaxyz)
                if betadist == 0:
                  continue
                if(betadist != dist):
                  distTemp = [(betaxyz[0]-resCA.xyz[0]),
                              (betaxyz[1]-resCA.xyz[1]),
                              (betaxyz[2]-resCA.xyz[2])]
                  betaxyz = [(resCA.xyz[0]+distTemp[0]*dist/betadist),
                             (resCA.xyz[1]+distTemp[1]*dist/betadist),
                             (resCA.xyz[2]+distTemp[2]*dist/betadist)]
                if(residue.resname != "GLY"):
                  dev = self.distance(resCB.xyz,betaxyz)
                  if(dev >=0.25 or outliers_only==False):
                    if(dev >=0.25):
                      cbetadev_ctr+=1
                      self.num_outliers+=1
                    d = geometry_restraints.dihedral(sites=[resN.xyz,resCA.xyz,betaxyz,resCB.xyz],
                                                     angle_ideal=0,
                                                     weight=1)
                    dihedralNABB = d.angle_model

                    #internal version of dihedral calculation used to test difference
                    #dihedralTemp = dihedral4pt(resN.xyz,resCA.xyz,betaxyz,resCB.xyz)
                    #dihedralNABB = dihedralTemp

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
                    occ = resCB.occ
                    analysis += '%s :%s:%s:%s:%5s:%7.3f:%7.2f:%7.2f:%s:\n' % \
                      (PDBfileStr,altchar,res,sub,resnum,dev,dihedralNABB,occ,altchar)
                    if (dev >= 0.25):
                      self.outliers += '%s :%s:%s:%s:%5s:%7.3f:%7.2f:%7.2f:%s:\n' % \
                        (PDBfileStr,altchar,res,sub,resnum,dev,dihedralNABB,occ,altchar)
                    key = altchar+res+sub+resnum
                    self.beta_ideal[key] = betaxyz
                    output_list.append([PDBfileStr,
                                        altchar,
                                        res,
                                        sub,
                                        resnum,
                                        dev,
                                        dihedralNABB,
                                        occ,
                                        altchar,
                                        resCB.xyz])
    summary = 'SUMMARY: %d C-beta Deviation >= 0.25 Angstrom (Goal: 0)' % (cbetadev_ctr)
    self.summary = 'SUMMARY: %d C-beta Deviation >= 0.25 Angstrom (Goal: 0)' % (cbetadev_ctr)
    return analysis.rstrip(), summary, output_list
  #}}}

  #{{{ idealized_calpha_angles
  def idealized_calpha_angles(self, residue):
    if(residue.resname == "ALA"):
      dist = 1.536
      angleCAB = 110.1
      dihedralNCAB = 122.9
      angleNAB = 110.6
      dihedralCNAB = -122.6
      angleideal = 111.2
    elif(residue.resname == "PRO"):
      dist = 1.530
      angleCAB = 112.2
      dihedralNCAB = 115.1
      angleNAB = 103.0
      dihedralCNAB = -120.7
      angleideal = 111.8
    elif(residue.resname == "VAL") or (residue.resname == "THR") or (residue.resname == "ILE"):
      dist = 1.540
      angleCAB = 109.1
      dihedralNCAB = 123.4
      angleNAB = 111.5
      dihedralCNAB = -122.0
      angleideal = 111.2
    elif(residue.resname == "GLY"):
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
  #}}}

  #{{{ construct_fourth
  def construct_fourth(self, resN,resCA,resC,dist,angle,dihedral,method="NCAB"):
    if (resN is not None and resCA is not None and resC is not None):
      if (method is "NCAB"):
        res0 = resN
        res1 = resC
        res2 = resCA
      elif (method is "CNAB"):
        res0 = resC
        res1 = resN
        res2 = resCA
      a = [(res2.xyz[0] - res1.xyz[0]),(res2.xyz[1] - res1.xyz[1]),(res2.xyz[2] - res1.xyz[2])]
      b = [(res0.xyz[0] - res1.xyz[0]),(res0.xyz[1] - res1.xyz[1]),(res0.xyz[2] - res1.xyz[2])]
      c = self.cross(a,b)
      cmag = self.mag(c)
      if(cmag > 0.000001):
        c=[(c[0]*dist/cmag),(c[1]*dist/cmag),(c[2]*dist/cmag)]
      c=[(c[0] + res2.xyz[0]),(c[1] + res2.xyz[1]),(c[2] + res2.xyz[2])]
      d = c
      angledhdrl = dihedral - 90
      a = [res1.xyz[0],res1.xyz[1],res1.xyz[2]]
      b = [res2.xyz[0],res2.xyz[1],res2.xyz[2]]
      newD = self.doaxisrot(d,angledhdrl,a,b)
      a = [(newD[0]-res2.xyz[0]),(newD[1]-res2.xyz[1]),(newD[2]-res2.xyz[2])]
      b = [(res1.xyz[0]-res2.xyz[0]),(res1.xyz[1]-res2.xyz[1]),(res1.xyz[2]-res2.xyz[2])]
      c = self.cross(a,b)
      cmag = self.mag(c)
      if(cmag > 0.000001):
        c=[(c[0]*dist/cmag),(c[1]*dist/cmag),(c[2]*dist/cmag)]
      angledhdrl = 90 - angle;
      a = [res2.xyz[0],res2.xyz[1],res2.xyz[2]]
      c = [(c[0]+res2.xyz[0]),(c[1]+res2.xyz[1]),(c[2]+res2.xyz[2])]
      b = c
      return self.doaxisrot(newD,angledhdrl,a,b)
  #}}}

  #{{{ doaxisrot
  def doaxisrot(self, d,theta,res1,res2):
    LOK = True
    xx = res2[0] - res1[0]
    yy = res2[1] - res1[1]
    zz = res2[2] - res1[2]
    if(xx**2 + yy**2 + zz**2 > 0.000001):
      cosn1 = xx/(xx**2 + yy**2 + zz**2)**.5
      cosn2 = yy/(xx**2 + yy**2 + zz**2)**.5
      cosn3 = zz/(xx**2 + yy**2 + zz**2)**.5
    else:
      cosn1 = xx
      cosn2 = yy
      cosn3 = zz
    while(theta > 360):
      theta = theta - 360
    while(theta < -360):
      theta = theta + 360
    if(theta > 180):
      theta = theta - 360
    if(theta < -180):
      theta = theta + 360
    if(theta > -.0001 and theta < 0.0001):
      LOK = False
    if(LOK):
      costheta = math.cos(theta*3.14159*2/360)
      sintheta = math.sin(theta*3.14159*2/360)
      a1 = [(cosn1**2 + (1-cosn1**2)*costheta),
            (cosn1*cosn2*(1-costheta)+cosn3*sintheta),
            (cosn1*cosn3*(1-costheta) - cosn2*sintheta)]
      a2 = [(cosn1*cosn2*(1-costheta) - cosn3*sintheta),
            (cosn2**2 + (1-cosn2**2)*costheta),
            (cosn2*cosn3*(1-costheta) + cosn1*sintheta)]
      a3 = [(cosn1*cosn3*(1-costheta)+cosn2*sintheta),
            (cosn2*cosn3*(1-costheta)-cosn1*sintheta),
            (cosn3**2+(1-cosn3**2)*costheta)]
      f1 = [(d[0]-res2[0]),(d[1]-res2[1]),(d[2]-res2[2])]
      f2 = [(f1[0]*a1[0] + f1[1]*a2[0] + f1[2]*a3[0]),
            (f1[0]*a1[1] + f1[1]*a2[1] + f1[2]*a3[1]),
            (f1[0]*a1[2]+f1[1]*a2[2]+f1[2]*a3[2])]
      return [(f2[0]+res2[0]),(f2[1]+res2[1]),(f2[2]+res2[2])]
  #}}}

  #{{{ dihedral4pt
  def dihedral4pt(self, p1,p2,p3,p4):
    a = [(p1[0]-p2[0]),(p1[1]-p2[1]),(p1[2]-p2[2])]
    b = [(p3[0]-p2[0]),(p3[1]-p2[1]),(p3[2]-p2[2])]
    d = cross(a,b)
    b = [(p2[0]-p3[0]),(p2[1]-p3[1]),(p2[2]-p3[2])]
    c = [(p4[0]-p3[0]),(p4[1]-p3[1]),(p4[2]-p3[2])]
    e = cross(b,c)
    dot = dotProduct(d,e)
    dmag = mag(d)
    emag = mag(e)
    if(dmag*emag < 0.0001):
      angle = 0.0
    else:
      angle = math.acos( dot/(dmag*emag) )
    anglehdrl = angle*360.0/(2*math.pi)
    f = cross(d,b)
    dot = dotProduct(f,e)
    fmag = mag(f)
    if(fmag*emag < 0.0001):
      angle = 0.0
    else:
      angle = math.acos( dot/(fmag*emag) )
    angle = angle*360.0/(2*math.pi)
    if(angle > 90.0):
      anglehdrl = -anglehdrl
    return anglehdrl
  #}}}

  #{{{ mag
  def mag(self, v):
    return (v[0]**2 + v[1]**2 + v[2]**2)**0.5
  #}}}

  #{{{ cross
  def cross(self, v1,v2):
    newcoords = [0,0,0]
    newcoords[0] = v1[1]*v2[2]-v1[2]*v2[1]
    newcoords[1] = v1[2]*v2[0]-v1[0]*v2[2]
    newcoords[2] = v1[0]*v2[1]-v1[1]*v2[0]
    return newcoords
  #}}}

  #{{{ dotProduct
  def dotProduct(self, v1,v2):
    return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
  #}}}

  #{{{ distance
  def distance(self, a,b):
    return ((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**0.5
  #}}}

  #{{{ get_functions
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
