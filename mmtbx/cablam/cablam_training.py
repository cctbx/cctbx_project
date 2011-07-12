# (jEdit options) :folding=explicit:collapseFolds=1:
#This module contains the training/exploration components of cablam
#It can be run stand-alone with many commandline options
#It is intended for use in determining contours, motif fingerprints, etc for
#  the annotation portion
#It is probably not intended for direct use in phenix, but is included as a
#  useful tool for understanding the cablam system

import os, sys
import argparse
#argparse is a commandline parsing module in Python 2.7+
#  New phenix distributions come with phenix.python 2.7+
#Replace with optparse for compatibility between Python 2.3 and 2.6
#Commandline parsing may be be changed to phil at some point
from iotbx import pdb  #contains the very useful hierarchy
from cctbx import geometry_restraints  #contains dihedral and angle calculators
import cablam_res, cablam_math

#{{{ prunerestype function
#Deletes all members of a given residue type from a dictionary of residues
#"Residue type" is determined from pdb.hierarchy's ag.resname, the upshot being
#  that non-residues like "HOH" waters and het groups can also be pruned to
#  improve performance if their three-letter codes are known.
#-------------------------------------------------------------------------------
def prunerestype(resdata, restype):
  reslist = resdata.keys()
  for residue in reslist:
    for alt in resdata[residue].alts:
      if resdata[residue].alts[alt]['resname'].strip() == restype:
        resdata[residue].removelinks()
        trash = resdata.pop(residue)
        break
#-------------------------------------------------------------------------------
#}}}

#{{{ stripB function
#Deletes all residues containing any atom of interest with atom.b > bmax from
#  a dictionary of residues, so that the uncertainty in these atoms cannot
#  contaminate later calculations.
#Important for training, not for annotation
#-------------------------------------------------------------------------------
def stripB(resdata, bmax):
  reslist = resdata.keys()
  for residue in reslist:
    deleted = False
    for alt in resdata[residue].alts:
      if deleted:
        break
      for atom in resdata[residue].atomb[alt]:
        if resdata[residue].atomb[alt][atom] > bmax:
          resdata[residue].removelinks()
          trash = resdata.pop(residue)
          deleted = True
          break
#-------------------------------------------------------------------------------
#}}}

#{{{ add dssp info to residue records
#Given a dictionary of residues and a pre-matched dssp filename, matches each
#  dssp line to a residue and stores a one-letter dssp code in the residue.
#Dssp codes are stored in residue.dssp, which is '' by default.
#Note: dssp leaves some residues unlabeled (if they have no identifiable
#  structure). For purposes of cablam, unlabeled residues are given a dssp code
#  of 'X'.
#-------------------------------------------------------------------------------
def dsspinfo(protein, dsspfilename):
  try:
    dsspfile = open(dsspfilename)
  except KeyError: #Placeholder for real exception catching
    sys.stderr.write("How on Earth did you cause a KeyError!?")

  dssp_dict = {}
  pastheader = False
  for dsspline in dsspfile:
    if pastheader:
      dssp_resnum = dsspline[6:10].strip()
      dssp_chain = dsspline[11:12]
      dssp_id = dssp_chain + dssp_resnum
      #note: dssp_id does not handle insertion or models codes yet. Does DSSP?
      dssp_code = dsspline[16:17]

      if dssp_code != ' ':
        dssp_dict[dssp_id] = dssp_code
      else:
        dssp_dict[dssp_id] = 'X'
    elif dsspline.startswith('  #'):
      pastheader = True
    else:
      pass
  dsspfile.close()

  for resid in protein:
    residue = protein[resid]
    shortresid = residue.chain+str(residue.resnum)

    try:
      residue.dssp = dssp_dict[shortresid]
    except KeyError:
      continue
#-------------------------------------------------------------------------------
#}}}

#{{{ headers to print at start of csv or kin file
#Simply a collection of headers and other fomatting used in printing output
#-------------------------------------------------------------------------------
def csv_header(kinorder, writeto=sys.stdout):
  writeto.write(',pdb:model:chain:resname:resnum:ins,')
  writeto.write(','.join(kinorder)+'\n')

def kin_header(kinorder,kinranges, writeto=sys.stdout):
  writeto.write('@text\n')
  for arg in sys.argv:
    writeto.write(arg + '  ')
  writeto.write('\n\n@kinemage\n')
  writeto.write('@dimensions {' + '} {'.join(kinorder)+'}\n')
  writeto.write('@dimminmax '+ ' '.join(kinranges)+'\n')

  writeto.write("""
@group {Rama Frame}
@dotlist {center} color= yellow off
0 0 0
@vectorlist {frame_xy} color= yellow
P -180 -180 0
180 -180 0
180 180 0
-180 180 0
-180 -180 0
@vectorlist {frame_xz} color= yellow
P -180 0 -180
180 0 -180
180 0 180
-180 0 180
-180 0 -180
@vectorlist {frame_yz} color= yellow
P 0 -180 -180
0 180 -180
0 180 180
0 -180 180
0 -180 -180

@group {points}
""")
  writeto.write(
    '@dotlist {points} nobutton dimension='+str(len(kinorder))+'\n')

def dssp_masters(dsspcodes, writeto=sys.stdout):
  writeto.write('\n')
  for dsspcode in dsspcodes:
    writeto.write('@master {'+dsspcode+'}\n')
  writeto.write('\n')
#-------------------------------------------------------------------------------
#}}}

#{{{ skipcheck function
#Residue types can be skipped during output without pruning their influence
#  entirely  This function handles checks for skipping, and returns boolean True
#  if the residue should be skipped.
#Additional functionality is expected in this function over time. More complex
#  sequence-sensitive selection is probable as I expand my training needs.
#As with pruning, important in training, less so in annotation.
#-------------------------------------------------------------------------------
def skipcheck(residue, skiplist, inclist):
  if skiplist:      #if there's anything to skip...
    doskip = False  #...the default state is include
  elif inclist:     #if there's nothing to skip but thing to include...
    doskip = True   #...the default state is skip
  else:
    return False    #if skip and include are empty, return default 'include'

  for skip in skiplist:
    currentres = residue
    if skip.startswith('_') and skip.endswith('_'):
      sys.stderr.write('\n\
        Invalid --skip flag argument: '+skip+ ' has\'_\' on both sides\n\n')
      sys.exit()

    #Underscores are used in the commandline call to indicate position relative
    #  to a residue of interest. For example, '_PRO' refers to pre-proline, and
    #  'GLY__' (two underscores) refers to post-post-glycine.  These loops
    #  manage the underscores
    while skip.startswith('_'):
      if currentres.nextres:
        currentres = currentres.nextres
        skip = skip[1:]
      else:
        return True  #cannot determine inclusion, so exclude
    while skip.endswith('_'):
      if currentres.prevres:
        currentres = currentres.prevres
        skip = skip[:-1]
      else:
        return True  #cannot determine inclusion, so exclude

    if currentres.firstalt('CA') is not None:
      resname = currentres.alts[currentres.firstalt('CA')]['resname']
    else:
      return True
    if resname == skip.upper():
      doskip = True

  for inc in inclist:
    currentres = residue
    if inc.startswith('_') and inc.endswith('_'):
      sys.stderr.write('\n\
        Invalid --skip flag argument: '+inc+ ' has\'_\' on both sides\n\n')
      sys.exit()

    while inc.startswith('_'):
      if currentres.nextres:
        currentres = currentres.nextres
        inc = inc[1:]
      else:
        return True
    while inc.endswith('_'):
      if currentres.prevres:
        currentres = currentres.prevres
        inc = inc[:-1]
      else:
        return True #cannot determine inclusion, so exclude

    if currentres.firstalt('CA') is not None:
      resname = currentres.alts[currentres.firstalt('CA')]['resname']
    else:
      return True   #cannot determine inclusion, so exclude
    if resname == inc.upper():
      doskip = False

  return doskip


#-------------------------------------------------------------------------------
#}}}

#{{{ check for continuous dsspcode
#Given a residue i, this function checks the near-sequence residues:
#  i-2, i-1, i, i+1, i+2 If these residues all have the same (specfied) dssp
#  code, it returns True.  Used in a skipcheck-like manner.
#This is useful in training, where removing the influence of, say helix caps,
#  can provide a cleaner picture of regular helix behavior.
#-------------------------------------------------------------------------------
def continuousdssp(residue, dsspcode):
  try:
    if residue.dssp == dsspcode:
      if residue.prevres.dssp == dsspcode:
        if residue.prevres.prevres.dssp == dsspcode:
          if residue.nextres.dssp == dsspcode:
            if residue.nextres.nextres.dssp == dsspcode:
              return True
  except AttributeError:
    #catches 'NoneType' object has no attribute 'dssp' for missing residues
    #  Might catch unexpected errors, but since that would only make the check
    #  more conservative, I think I can live with it
    return False
  return False
#-------------------------------------------------------------------------------
#}}}

#{{{ print functions
#These are wrappers for printing to various formats. The print methods that they
#  wrap are in the main residue class found in cablam_res.
#-------------------------------------------------------------------------------
#prints residues in .kin format
def kin_print(protein, kinorder, skiplist=[], inclist=[], writeto=sys.stdout):
  reslist = protein.keys()
  reslist.sort()
  for resid in reslist:
    if skipcheck(protein[resid], skiplist, inclist):
      pass
    else:
      protein[resid].printtokin(kinorder, writeto)

#prints residues with selected dssp codes in .kin format
def kin_print_dssp(protein, kinorder, dsspcodes, skiplist=[], inclist=[],
  needcontinuous=False, writeto=sys.stdout):

  reslist = protein.keys()
  reslist.sort()
  for code in dsspcodes:
    writeto.write('\n@dotlist {'+pdbid+' '+code+'} dimension='+
      str(len(kinorder))+' master= {'+code+'}\n')
    for resid in reslist:
      if protein[resid].dssp != code: continue
      if needcontinuous and not continuousdssp(protein[resid], code): continue
      if skipcheck(protein[resid], skiplist, inclist): continue
      #only prints of the above checks pass
      protein[resid].printtokin(kinorder, writeto)

#prints residues with selected dssp codes in comma-separated format
def csv_print_dssp(protein, kinorder, dsspcodes, skiplist=[], inclist=[],
  needcontinuous=False, writeto=sys.stdout):

  reslist = protein.keys()
  reslist.sort()
  for code in dsspcodes:
    for resid in reslist:
      if protein[resid].dssp != code: continue
      if needcontinuous and not continuousdssp(protein[resid], code): continue
      if skipcheck(protein[resid], skiplist, inclist): continue
      #only prints if the above checks pass
      protein[resid].printtocsv(kinorder, writeto)

#Prints residues in comma-separated format, suitable for contouring and other
#  analysis
#This is currently the default behavior of cablam_training.  This output format
#  is used to generate percentile and probability contours for cablam_annote
#  using the programs Silk and kin2Dcont/kin3Dcont from the Richardson Lab.
def csv_print(protein, kinorder, skiplist=[], inclist=[], writeto=sys.stdout):
  reslist = protein.keys()
  reslist.sort()
  for resid in reslist:
    if skipcheck(protein[resid], skiplist, inclist):
      pass
    else:
      protein[resid].printtocsv(kinorder, writeto)
#-------------------------------------------------------------------------------
#}}}

#{{{ run
#The run function is currently rather messy. (Indeed, all of
#  cablam_training is a bit messy, as it's really a development tool, not a
#  general-use program.) Hopefully, everything needed for general use (structure
#  annotation) has been packaged in other modules for easy access. Good luck.
def run():
  parser = argparse.ArgumentParser()
  parser.add_argument('file_or_dir',
    help='the path to the file or directory to be operated on')
  parser.add_argument('-s','--separatefiles',action='store_true',
    help='Generate a separate, auto-named output file put each input file')
  parser.add_argument('-k','--kin', action='store_true',
    help='Print output in .kin format rather than comma-separated .csv format')

  parser.add_argument('--cad', action='store_true',
    help='2 CA pseudodihedrals')
  parser.add_argument('--caa', action='store_true',
    help='3 CA pseudo angles')
  parser.add_argument('--cod', action='store_true',
    help='2 CO pseudodihedrals')
  parser.add_argument('--rama', action='store_true',
    help='2 Ramachandran dihedrals: phi, psi')
  parser.add_argument('--exrama', action='store_true',
    help='4 Ramachandran dihedrals: psi-1, phi, psi, phi+1')
  parser.add_argument('--tau', action='store_true',
    help='1 backbone angle: tau (defined by N-CA-C)')
  parser.add_argument('-a','--allmeasures', action='store_true',
    help='Shortcut for \"all of the above\"')

  parser.add_argument('-b','--bmax',
    help='Set a max b factor, residues containing a backbone atom with higher b will be pruned, rocommended: -b=30')
  parser.add_argument('--dssp', dest='dssplist',
    help='List of DSSP codes, without spaces.  e.g. --dssp=HGI Use H,G,I,T,E,B,S X for unlabeled, just A for all, lowercase okay but poor form')
  parser.add_argument('--dsspath',
    help='Path to dssp dir or file for reference.')
  parser.add_argument('--continuous', action='store_true',
    help='Only residues i with the same DSSP code from i-2 to i+2 will print')

  parser.add_argument('--prune',
    help='List of restypes to be pruned, separated by commas, no spaces eg PRO')
  parser.add_argument('--skip',
    help='List of restypes to be skipped during printing, separated by commas')
  parser.add_argument('--include',
    help='List of restypes to be printed, all others will be skipped')

  args = parser.parse_args()

  if os.path.isdir(args.file_or_dir):
    fileset = os.listdir(args.file_or_dir)
    dirpath = args.file_or_dir
  elif os.path.isfile(args.file_or_dir):
    fileset = [args.file_or_dir]
    dirpath = None
  else:
    sys.stderr.write("Could not identify valid target file or dir.\n")
    ## print the help section
    sys.exit()

  targetatoms = ["CA","O","C","N"]

  #{{{ measurement selection
  #This section manages the user's orders for calculations
  #Note: The 'kin' in kinorder and kin ranges is a misnomer
  #-----------------------------------------------------------------------------
  if args.allmeasures:
    args.cad = True
    args.caa = True
    args.cod = True
    args.exrama = True
    args.tau = True

  kinorder, kinranges = [],[]
  if args.cad:
    kinorder.append('CA_d_in'),  kinranges.append('-180 180')
    kinorder.append('CA_d_out'), kinranges.append('-180 180')
  else:
    pass

  if args.cod:
    kinorder.append('CO_d_in'),  kinranges.append('-180 180')
    kinorder.append('CO_d_out'), kinranges.append('-180 180')
  else:
    pass

  if args.rama or args.exrama:
    if args.exrama:
      kinorder.append('psi-1'), kinranges.append('-180 180')
      kinorder.append('phi'),   kinranges.append('-180 180')
      kinorder.append('psi'),   kinranges.append('-180 180')
      kinorder.append('phi+1'), kinranges.append('-180 180')
    else:
      kinorder.append('phi'), kinranges.append('-180 180')
      kinorder.append('psi'), kinranges.append('-180 180')
  else:
    pass

  if args.caa:
    kinorder.append('CA_a_in'),  kinranges.append('0 180')
    kinorder.append('CA_a'),     kinranges.append('0 180')
    kinorder.append('CA_a_out'), kinranges.append('0 180')
  else:
    pass

  if args.tau:
    kinorder.append('tau'), kinranges.append('0 180')
  else:
    pass

  #The following lines record the order and values for kinorder and kinranges
  #  for sake of reference
  #kinorder =  ['CA_d_in', 'CA_d_out','CO_d_in', 'CO_d_out',
  #  'psi-1',   'phi',     'psi',     'phi+1',   'tau']
  #kinranges = ['-180 180','-180 180','-180 180','-180 180',
  #  '-180 180','-180 180','-180 180','-180 180','0 180']
  #-----------------------------------------------------------------------------
  #}}}

  if args.dssplist:
    dssprequest = [code.upper() for code in args.dssplist]
    realdssp = ('H','G','I','T','E','B','S','X')
    #Only real dssp codes will be accepted and passed for future use. User
    #  requests for invalid codes will be silently ignored.
    if 'A' in dssprequest:
      dsspcodes = ['H','G','I','T','E','B','S','X']
    else:
      dsspcodes = []
      for code in realdssp:
        if code in dssprequest:
          dsspcodes.append(code)
        else:
          pass

  prunelist = []
  if args.prune:
    prunelist = args.prune.split(',')
    prunelist = [res.upper() for res in prunelist] #Ha ha! List comprehension!

  skiplist = []
  inclist = []
  if args.skip:
    skiplist = args.skip.split(',')
  if args.include:
    inclist = args.include.split(',')


  if args.separatefiles:
    pass
  else:
    if args.kin:
      kin_header(kinorder,kinranges)
    else:
      csv_header(kinorder)


  for filename in fileset:
    if dirpath: #must add the path if using the listed contents of a dir
      filename = os.path.join(dirpath,filename)
    else:
      pass

    pdbid = os.path.basename(filename)
    ##if filename.endswith('.pdb'):
    ##  pdb_io = pdb.input(filename)
    ##  pdbid = os.path.basename(filename.rstrip('.pdb'))
    ##else:
    ##  sys.stderr.write(
    ##    filename + " does not end with \'.pdb\' and may not be a pdb file\n")
    ##  continue

    if args.dssplist:
      if args.dsspath:
        if os.path.isdir(args.dsspath):
          dsspfileset = os.listdir(args.dsspath)
          dsspdirpath = args.dsspath
        elif os.path.isfile(args.dsspath):
          dsspfileset = [ os.path.basename(args.dsspath) ]
          dsspdirpath =   os.path.dirname(args.dsspath)
        else:
          sys.stderr.write("could not identify valid dssp file or dir.\n")
          sys.exit()
      else:
        #Later, an empty --dsspath should cause the program to look in the same
        #  folder as the .pdb for a matching .dssp
        sys.stderr.write("Must provide --dsspath in commandline call.\n")
        sys.exit()
    else:
      pass

    sys.stderr.write(pdbid+'\n')
    pdb_io = pdb.input(filename)
    hierarchy = pdb_io.construct_hierarchy()
    resdata = cablam_res.construct_linked_residues(hierarchy,targetatoms,pdbid)
    if not resdata: #skips further processing of files not readable by hierarchy
      continue

    prunerestype(resdata, 'HOH')
    for restype in prunelist:
      prunerestype(resdata, restype)

    if args.bmax:
      stripB(resdata,float(args.bmax))
    else:
      pass

    if args.dssplist and args.dsspath:
      if os.path.isdir(args.dsspath): #Need to search for matching .dssp
        for dsspfilename in dsspfileset:
          dsspid = dsspfilename.rstrip('.dssp').upper()
          if dsspid == pdbid.upper():
            #dsspfile = open(os.path.join(dsspdirpath,dsspfilename))
            sys.stderr.write('Matched '+dsspfilename+' to '+pdbid+'\n')
            dsspinfo(resdata, os.path.join(dsspdirpath,dsspfilename))
            break
          else:
            continue
        else: #this else executes if the 'for' loop never hits a 'break'
          sys.stderr.write('No matching .dssp file found for '+pdbid+' in '+
            dsspdirpath+'\n')
      else: #If a file is specified, that specification overrides matching
        if os.path.isdir(args.file_or_dir):
          sys.stderr.write('Warning: a single .dssp has been specified for a dir of .pdb\'s.  All .pdb\'s will attempt to use this .dssp')
        else:
          pass
        #dsspfile = open(args.dsspath)
        dsspinfo(resdata, args.dsspath)
      #Now that we have an open dssp file, we can check it
    else:
      pass


    #{{{ calculation calls
    #---------------------------------------------------------------------------
    if args.cad and args.caa:
      cablam_math.CApseudos(resdata, dodihedrals = True, doangles = True)
    elif args.cad:
      cablam_math.CApseudos(resdata, dodihedrals = True, doangles = False)
    elif args.caa:
      cablam_math.CApseudos(resdata, dodihedrals = False, doangles = True)
    else: #no CA-based calculations
      pass

    if args.cod:
      cablam_math.COpseudodihedrals(resdata)
    else:
      pass

    if args.rama or args.exrama:
      cablam_math.phipsi(resdata)
    else:
      pass

    if args.tau:
      cablam_math.taucalc(resdata)
    else:
      pass
    #---------------------------------------------------------------------------
    #}}}

    if args.dssplist:
      if args.kin:
        if args.separatefiles:
          outfile = open(pdbid+'_r2b2.kin','w')
          kin_header(kinorder,kinranges,writeto=outfile)
          dssp_masters(dsspcodes,writeto=outfile)
          kin_print_dssp(resdata, kinorder, dsspcodes,
            skiplist, inclist, args.continuous, writeto=outfile)
          outfile.close()
        else:
          dssp_masters(dsspcodes)
          kin_print_dssp(resdata, kinorder, dsspcodes,
            skiplist,inclist,args.continuous)
      else:
        if args.separatefiles:
          outfile = open(pdbid_+'_r2b2.csv','w')
          csv_header(kinorder,writeto=outfile)
          csv_print_dssp(resdata, kinorder, dsspcodes,
            skiplist, inclist,args.continuous,writeto=outfile)
          outfile.close()
        else:
          csv_print_dssp(resdata,kinorder,dsspcodes,
            skiplist,inclist,args.continuous)

    else: #i.e. not args.dssplist
      if args.kin:
        if args.separatefiles:
          outfile = open(pdbid+'_r2b2.kin','w')
          kin_header(kinorder,kinranges,writeto=outfile)
          kin_print(resdata, kinorder, skiplist, inclist, writeto=outfile)
          outfile.close()
        else:
          kin_print(resdata,kinorder,skiplist,inclist)
      else:
        if args.separatefiles:
          outfile = open(pdbid+'_r2b2.csv','w')
          csv_header(kinorder,writeto=outfile)
          csv_print(resdata, kinorder, skiplist, inclist, writeto=outfile)
          outfile.close()
        else:
          csv_print(resdata,kinorder,skiplist,inclist)
#-------------------------------------------------------------------------------
#}}}

#{{{ "__main__"
#-------------------------------------------------------------------------------
if __name__ == "__main__":
  #__main__ needs to generate a hierarchy
  run()
#-------------------------------------------------------------------------------
#}}}
