# (jEdit options) :folding=explicit:collapseFolds=1:
#This module contains the training/exploration components of cablam
#It can be run stand-alone with many commandline options
#It is intended for use in determining contours, motif fingerprints, etc for
#  the annotation portion
#It is probably not intended for direct use in phenix, but is included as a
#  useful tool for understanding the cablam system
#The May 2012 update reflects a substantial change in method. from the previous
#  version.  DSSP has been abandoned in favor of directly assessing hydrogen
#  bonding patterns determined by Probe.  Hydrogen bonding pattern definitions
#  are stored in fingerprints.py, but will likely be assembled into a separate
#  library in the future.

import os, sys
import argparse
#argparse is a commandline parsing module in Python 2.7+
#  New phenix distributions come with phenix.python 2.7+
#Replace with optparse for compatibility between Python 2.3 and 2.6
#Argument parsing should be be changed to phil at some point
from iotbx import pdb  #contains the very useful hierarchy
from mmtbx.cablam import cablam_res #contains a data structure derived from
#  hierarchy, but more suited to cablam's needs - specifically it can hold
#  geometric and probe measures and can look forward and backward in sequence
from mmtbx.cablam import cablam_math #contains geometric measure calculators
from mmtbx.cablam import fingerprints #contains motif definitions
#  Storage for motif definitions subject to change

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
#Will need to make distinction between main- and side-chain eventually
#-------------------------------------------------------------------------------
def stripB(resdata, bmax):
  reslist = resdata.keys()
  for resid in reslist:
    deleted = False
    for alt in resdata[resid].alts:
      if deleted:
        break
      for atom in resdata[resid].atomb[alt]:
        if resdata[resid].atomb[alt][atom] > bmax:
          resdata[resid].removelinks()
          trash = resdata.pop(resid)
          deleted = True
          break
#-------------------------------------------------------------------------------
#}}}

#{{{ prune alts function
#Deletes all residues that have alternate conformations at one or more atoms
#  from a dictionary of residues, so that uncertainty in these atoms or in their
#  relations with other atoms in the structure cannot contaminate later
#  calculations.
#A function for robustly handling and choosing among alternates is eventually
#  forthcoming, but will be separate.
#-------------------------------------------------------------------------------
def prune_alts(resdata):
  reslist = resdata.keys()
  for resid in reslist:
    residue = resdata[resid]
    if len(residue.alts) > 1:
      residue.removelinks()
      trash = resdata.pop(resid)
#-------------------------------------------------------------------------------
#}}}

#{{{ skipcheck function
#Residue types can be skipped during output without pruning their influence
#  entirely. This function handles checks for skipping, and returns boolean True
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

#{{{ add probe data function
#Adds mainchina-mainchain hydrogen bonding information from 'unformated' Probe
#  output to a dictionary of residues.
#At the moment, reliant on precomputed .probe files, will gain run-time Probe
#May gain other contact relationship info, by mc-mc H-bonds are most important
#-------------------------------------------------------------------------------
def add_probe_data(resdata, open_probe_file):
  reskeys = resdata.keys()
  for line in open_probe_file:
    #Probe Unformatted Output:
    #name:pat:type:srcAtom:targAtom:min-gap:gap:spX:spY:spZ:spikeLen:score:stype:ttype:x:y:z:sBval:tBval
    #for condensed output we have:
    #name:pat:type:srcAtom:targAtom:*dotcount*:min-gap:gap:spX:spY:spZ:spikeLen:score:stype:ttype:x:y:z:sBval:tBval
    ###'name' is set by the user on the command line
    ###'pat' is one of 1->1, 1->2, or 2->1; where 1 is src and 2 is targ.
    ###'type' is one of wc, cc, so, bo, hb (wide/close contact, small/bad overlap, h-bond).
    ###'srcAtom' and 'targAtom' follow the pattern CNNNNITTT AAAAL, where C is chain, N is number, I is insertion code, T is residue type, A is atom name, and L is alternate conformation flag.
    ###'*dotcount*' is condensed-output-only, and gives the number of dots in the contact
    ###'min-gap' is the distance between atoms, minus their van der Waals radii; i.e., the distance of closest approach for their vdW surfaces. gap is the distance between vdW surfaces at the current dot. Negative values indicate overlap (clashes or H-bonds).
    ###'x','y','z' is a point on the vdW surface; 'spX','spY','spZ' is tip of spike, if any (same as x,y,z for contacts)
    ###'score' is "this dot's contribution to the [Probe] score" (scaled already? YES)
    ###'stype' and 'ttype' are heavy-atom element name (C, N, O, etc)

    if not line.strip(): continue #averts an IndexError problem with empty lines
    bnana = line.split(':')

    name = bnana[0]
    pattern = bnana[1]
    interactiontype = bnana[2]
    if not interactiontype == 'hb': continue #skip non-h-bonds

    srcAtom = bnana[3]
    srcChain =    srcAtom[0:2].strip()
    srcNum =      int(srcAtom[2:6].strip())
    srcIns =      srcAtom[6:7]#.strip()
    srcResname =  srcAtom[7:10].strip()
    if srcResname == 'HOH': continue #skip waters
    srcAtomname = srcAtom[11:15].strip()
    srcAlt =      srcAtom[15:16].strip()

    trgAtom = bnana[4]
    #if srcAtom == oldsrcAtom and trgAtom == oldtrgAtom: continue
    # ^^ Probe produces a line for every dot, so have to check for relationship repeats
    #going to count dots per bond as a measure of strength instead
    trgChain =    trgAtom[0:2].strip()
    trgNum =      int(trgAtom[2:6].strip())
    trgIns =      trgAtom[6:7]#.strip()
    trgResname =  trgAtom[7:10].strip()
    if trgResname == 'HOH': continue #skip waters
    trgAtomname = trgAtom[11:15].strip()
    trgAlt =      trgAtom[15:16].strip()

    dotcount = bnana[5]

    #reference: resid_string = ' '.join([modelid, chainid, '%04i' % resnum, icode])
    if (srcAtomname == 'O' and trgAtomname == 'H') or (srcAtomname == 'H' and trgAtomname == 'O'):
      srcResid_string = ' '.join(['', srcChain, '%04i' % srcNum, srcIns])
      trgResid_string = ' '.join(['', trgChain, '%04i' % trgNum, trgIns])
      #probe does not run itself on more than the first model, so the modelid
      #  bit is a problem.  It's probably either '' or '1'
      if (srcResid_string in reskeys) and (trgResid_string in reskeys):
        srcResidue = resdata[srcResid_string]
        trgResidue = resdata[trgResid_string]
        #bondjump = str(int(trgNum) - int(srcNum))
        #Probably don't need bondjump here anymore.
        #  Probably will make a function in cablam_res to find it.

        if srcAtomname == 'H' and trgResidue not in srcResidue.probeH:
          srcResidue.probeH.append(trgResidue)
        elif srcAtomname == 'O' and trgResidue not in srcResidue.probeO:
          srcResidue.probeO.append(trgResidue)
        else:
          pass
#-------------------------------------------------------------------------------
#}}}

#{{{ Output function collection
#A collection of headers, formatting, and printint functions used in output
#Default output is to stdout, but anything with a .write can be passed to the
#  'writeto=' argument of most functions.  Functions that lack a 'writeto='
#  generate or find uniquely named files in the working dir for their output.
#Print methods called by these functions are generally from cablam_res.py
#-------------------------------------------------------------------------------

#kin_frame is a 3-dimensional frame for dihedral-space (-180 to 180) kinemages
def kin_frame(writeto=sys.stdout):
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

""")

#{{{ --- CSV printing ---
#-------------------------------------------------------------------------------
#csv_header writes column names for the top of a .csv
#It starts with a comma for a reason, but I don't remember what it is
def csv_header(kinorder, doconnections=False, writeto=sys.stdout):
  writeto.write(',pdb:model:chain:resnum:ins:resname,')
  writeto.write(','.join(kinorder))
  if doconnections:
    writeto.write(',prevres,nextres')
  writeto.write('\n')

#Prints residues in comma-separated format, suitable for contouring and other
#  analysis
#This is currently the default behavior of cablam_training.  This output format
#  is used to generate percentile and probability contours for cablam_annote
#  using the programs Silk and kin2Dcont/kin3Dcont from the Richardson Lab.
def csv_print(protein, kinorder, skiplist=[], inclist=[],
  doconnections=False, writeto=sys.stdout):
  reslist = protein.keys()
  reslist.sort()
  for resid in reslist:
    if skipcheck(protein[resid], skiplist, inclist):
      pass
    else:
      protein[resid].printtocsv(kinorder, doconnections, writeto)
#-------------------------------------------------------------------------------
#}}}

#{{{ --- Generic KIN printing ---
#-------------------------------------------------------------------------------
#kin_header writes the start of a kinemage file
#@text provides self-documentation of the commandline used to generate the .kin
#@dimensions and @dimminmax allow the .kin to handle high-dimensional data
def kin_header(kinorder,kinranges, writeto=sys.stdout):
  writeto.write('@text\n')
  for arg in sys.argv:
    writeto.write(arg + '  ')
  writeto.write('\n\n@kinemage\n')
  writeto.write('@dimensions {' + '} {'.join(kinorder)+'}\n')
  writeto.write('@dimminmax '+ ' '.join(kinranges)+'\n')
  kin_frame(writeto=writeto)
  writeto.write('@group {points}\n')
  writeto.write(
    '@dotlist {points} nobutton dimension='+str(len(kinorder))+'\n')

#prints residues in .kin format
#Uses skipcheck() to select residues to print (default includes all)
def kin_print(protein, kinorder, skiplist=[], inclist=[], writeto=sys.stdout):
  reslist = protein.keys()
  reslist.sort()
  for resid in reslist:
    if skipcheck(protein[resid], skiplist, inclist):
      pass
    else:
      protein[resid].printtokin(kinorder, writeto)
#-------------------------------------------------------------------------------
#}}}

#{{{ --- Default PROBE printing ---
#-------------------------------------------------------------------------------
#Creates files and prints headers in them for generic probe output
#One .kin for each unique label in each motif. This can produce a lot of files.
def kin_print_probe_header(full_label_list, kinorder, kinranges):
  outfiles = {}
  for label in full_label_list:
    outfiles[label] = open(label+'.kin','a')
    outfiles[label].write('\n@kinemage\n')
    outfiles[label].write('@dimensions {' + '} {'.join(kinorder)+'}\n')
    outfiles[label].write('@dimminmax '+ ' '.join(kinranges)+'\n')
    kin_frame(writeto=outfiles[label])
    outfiles[label].write(
      '@group {'+label+'} dominant animate\n@dotlist {'+label+
      '} dimension='+str(len(kinorder))+'\n')
  return outfiles

#For producing distributions of points in cablam space
#Generic output is one point (many dimensions) for each residue that matches a
#  motif definition/fingerprint.
def kin_print_probe(resdata, kinorder, outfiles, skiplist=[], inclist=[]):
  reskeys = resdata.keys()
  reskeys.sort()
  for resid in reskeys:
    residue = resdata[resid]
    if skipcheck(residue, skiplist, inclist): continue
    for label in residue.motifs:
      residue.printtokin(kinorder, writeto=outfiles[label])
#-------------------------------------------------------------------------------
#}}}

#{{{ --- PROBE ANNOTE printing ---
#-------------------------------------------------------------------------------
#For annotating an existing .kin file with balls at CA's participating in
#  motifs of interest.
#Produces one .kin per input file.
#Does not require a header as such.
def kin_print_probe_annote(resdata, motif_list, writeto=sys.stdout):
  reskeys = resdata.keys()
  reskeys.sort()
  for motif_name in motif_list:
    writeto.write('@group {'+motif_name+'}\n')
    for label in fingerprints.fingerprints[motif_name].labellist:
      writeto.write('@balllist {'+label+'}\n')
      for resid in reskeys:
        residue = resdata[resid]
        if label in residue.motifs:
          firstalt = residue.firstalt('CA')
          #try:
          CAxyz = residue.atomxyz[firstalt]['CA']
          pointid = residue.pdbid +' '+ residue.chain +' '+ str(residue.resnum)+' '+ label
          writeto.write("{ "+pointid+" } "+str(CAxyz[0])+" "+str(CAxyz[1])+" "+str(CAxyz[2])+"\n")
#-------------------------------------------------------------------------------
#}}}

#{{{ --- PROBE BY INSTANCE printing ---
#-------------------------------------------------------------------------------
#Creates files and prints headers in them for instance output
#One .kin for each motif. This can produce several files.
def kin_print_by_instance_header(motif_list, kinorder, kinranges):
  outfiles = {}
  for motif_name in motif_list:
    outfiles[motif_name] = open(motif_name+'_instances.kin', 'w')
    outfiles[motif_name].write('@text\n')
    for arg in sys.argv:
      outfiles[motif_name].write(arg + '  ')
    outfiles[motif_name].write('\n@kinemage\n')
    outfiles[motif_name].write('@dimensions {' + '} {'.join(kinorder)+'}\n')
    outfiles[motif_name].write('@dimminmax '+ ' '.join(kinranges)+'\n')
    kin_frame(writeto=outfiles[motif_name])
  return outfiles

#What this means is: each instance of a full motif, printed as a vector list so
#  the path through cablam space can be followed
def kin_print_by_instance(resdata, motif_list, kinorder, outfiles):
  reskeys = resdata.keys()
  reskeys.sort()
  for motif_name in motif_list:
    motif = fingerprints.fingerprints[motif_name]
    for resid in reskeys:
      residue = resdata[resid]
      if motif.labellist[0] in residue.motifs:
        #Each instance is a group so that they are animatable in the kinemage
        outfiles[motif_name].write(
          '@group {'+residue.pdbid.rstrip('.pdb')+' '+str(residue.resnum)+
          '} dominant animate\n@vectorlist {'+motif_name+
          '} dimension='+str(len(kinorder))+'\n')
        curres = residue
        for label in motif.labellist:
          outline = ['{'+curres.id_with_resname(sep=' ')+'_'+label+'}']
          for order in kinorder:
            try:
              outline.append(str(curres.measures[order]))
            except KeyError:
              sys.stderr.write('KeyError in kin print by instance')
              break
          else:
            outfiles[motif_name].write(' '.join(outline)+'\n')
          curres = curres.nextres
#-------------------------------------------------------------------------------
#}}}

#-------------------------------------------------------------------------------
#}}}

#{{{ run
#The run function is currently rather messy. (Indeed, all of
#  cablam_training is a bit messy, as it's really a development tool, not a
#  general-use program.) Hopefully, everything needed for general use (structure
#  annotation) has been packaged in other modules for easy access. Good luck.
def run():
  #{{{ arguments
  #-----------------------------------------------------------------------------
  parser = argparse.ArgumentParser()
  parser.add_argument('file_or_dir',
    help='the path to the file or directory to be operated on')
  parser.add_argument('-s','--separatefiles',action='store_true',
    help='Generate a separate, auto-named output file put each input file')
  parser.add_argument('-k','--kin', action='store_true',
    help='Print output in .kin format rather than comma-separated .csv format')
  parser.add_argument('--doconnections', action='store_true',
    help='Adds prevres and nextres columns to .csv output')
  parser.add_argument('--debug', action='store_true',
    help='Adds some text printed to stderr for debugging esp. for fingerprints')

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
  parser.add_argument('--omega', action='store_true',
    help='1 backbone dihedral: omega (defined by CA_1-C_1-N_2-CA_2)')
  parser.add_argument('-a','--allmeasures', action='store_true',
    help='Shortcut for \"all of the above\"')
  parser.add_argument('-c', '--cablam', action='store_true',
    help='Shortcut for just cablam-relevant measures CA_d_in, CA_d_out, CA_a')

  parser.add_argument('--probe',
    help='Activates hydrogen bonding analysis, --probe=motif_name1,motif_name2,... use --listmotifs to list available fingerprints')
  #Later: --probe without --probepath should run phenix.probe
  #libtbx easyrun
  parser.add_argument('--probepath',
    help='Stores path to dir of probed files')
  parser.add_argument('--listmotifs',action='store_true',
    help='print motifs/fingerprints available')
  parser.add_argument('--probemode',
    help='=kin for dotlist kins (default) =annote for ball on model, =instance for vectorlist kins')

  parser.add_argument('-b','--bmax',
    help='Set a max b factor, residues containing a backbone atom with higher b will be pruned, rocommended: -b=30')
  parser.add_argument('--prunealts', action='store_true',
    help='Removes all residues with alternate conformations in relevant atoms')
  parser.add_argument('--prune',
    help='List of restypes to be pruned, separated by commas, no spaces eg PRO')
  parser.add_argument('--skip',
    help='List of restypes to be skipped during printing, separated by commas')
  parser.add_argument('--include',
    help='List of restypes to be printed, all others will be skipped')

  args = parser.parse_args()
  #-----------------------------------------------------------------------------
  #}}}

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

  if args.listmotifs:
    sys.stdout.write('\n\n')
    for motifname in fingerprints.fingerprints:
      sys.stdout.write(motifname + '\n')
    sys.exit()

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
    args.omega = True

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

  if args.omega:
    kinorder.append('omega'), kinranges.append('-180 180')
  else:
    pass

  if args.cablam:
    if 'CA_d_in' not in kinorder:
      kinorder.append('CA_d_in'),  kinranges.append('-180 180')
    if 'CA_d_out' not in kinorder:
      kinorder.append('CA_d_out'),  kinranges.append('-180 180')
    if 'CA_a' not in kinorder:
      kinorder.append('CA_a'),  kinranges.append('0 180')
  else:
    pass

  #The following lines record the order and values for kinorder and kinranges
  #  for sake of reference
  #kinorder =  ['CA_d_in', 'CA_d_out','CO_d_in', 'CO_d_out',
  #  'psi-1',   'phi',     'psi',     'phi+1',   'tau', 'omega']
  #kinranges = ['-180 180','-180 180','-180 180','-180 180',
  #  '-180 180','-180 180','-180 180','-180 180','0 180', '-180 180']
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ setup
  #-----------------------------------------------------------------------------
  targetatoms = ["CA","O","C","N"]

  outfiles = []
  if args.probe:
    motif_list = args.probe.split(',')
    if args.probepath:
      probefilelist = os.listdir(args.probepath)
      if args.probemode == 'kin' or args.probemode == None:
        outfiles = kin_print_probe_header(fingerprints.get_all_labels(motif_list),kinorder,kinranges)
      elif args.probemode == 'instance':
        outfiles = kin_print_by_instance_header(motif_list, kinorder, kinranges)

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
      csv_header(kinorder,args.doconnections)
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ get file, start loop
  #-----------------------------------------------------------------------------
  for filename in fileset:
    if dirpath: #must add the path if using the listed contents of a dir
      filename = os.path.join(dirpath,filename)
    else:
      pass

    pdbid = os.path.basename(filename)

    sys.stderr.write(pdbid+'\n')
    pdb_io = pdb.input(filename)
    hierarchy = pdb_io.construct_hierarchy()
    resdata = cablam_res.construct_linked_residues(hierarchy,targetatoms,pdbid)
    if not resdata: #skips further processing of files not readable by hierarchy
      continue
  #-----------------------------------------------------------------------------
  #}}}

    #{{{ preprocessing
    #---------------------------------------------------------------------------
    prunerestype(resdata, 'HOH')
    for restype in prunelist:
      prunerestype(resdata, restype)

    if args.bmax:
      stripB(resdata,float(args.bmax))

    if args.prunealts:
      prune_alts(resdata)
    #---------------------------------------------------------------------------
    #}}}

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

    if args.omega:
      cablam_math.omegacalc(resdata)
    else:
      pass

    if args.cablam:
      cablam_math.cablam_measures(resdata)
    else:
      pass
    #---------------------------------------------------------------------------
    #}}}

    #{{{ probe stuff
    #---------------------------------------------------------------------------
    #need the run phenix.probe
    if args.probe and args.probepath:
      probefilename = pdbid.rstrip('.pdb') + '.probe'
      if probefilename in probefilelist:
        probefilepath = os.path.join(args.probepath,probefilename)
        open_probe_file = open(probefilepath)
        add_probe_data(resdata,open_probe_file)
        open_probe_file.close()
      else:
        continue

    if args.probe:
      #motif_list = args.probe.split(',') (this was actually done earlier)
      for motif_name in motif_list:
        fingerprints.annote_motif_protein(resdata,motif_name)
    #---------------------------------------------------------------------------
    #}}}

    #{{{ output
    #---------------------------------------------------------------------------
    #--probemode=kin for dotlist kins, this is the default
    #--probemode=annote for balls drawn at CA positions on the model
    #--probemode=instance for kins where each veclist is one instance of motif
    if args.probe and args.probepath:
      if args.probemode == 'kin' or args.probemode == None:
        kin_print_probe(resdata, kinorder, outfiles, skiplist, inclist)
      elif args.probemode == 'annote':
        outfile = open(pdbid+'cablam_motifs.kin','w')
        kin_print_probe_annote(resdata, motif_list, writeto=outfile)
        outfile.close()
      elif args.probemode == 'instance':
        kin_print_by_instance(resdata, motif_list, kinorder, outfiles)
      else:
        sys.stderr.write('\n\nUnrecognized probemode request\n\n')
        sys.exit()
      #add if args.kin once things basically work
      outfile = sys.stdout
      #need printer from probe version

    #default printing, with no arguments, is to .csv, one line per residue
    #--separatefiles writes a separate file for each input file to working dir
    #--kin prints kinemage file, dotlist, one point per residue
    #--doconnections adds connectivity information to csv output
    else:
      if args.kin:
        if args.separatefiles:
          outfile = open(pdbid+'_cablam.kin','w')
          kin_header(kinorder,kinranges,writeto=outfile)
          kin_print(resdata, kinorder, skiplist, inclist, writeto=outfile)
          outfile.close()
        else:
          kin_print(resdata,kinorder,skiplist,inclist)
      else:
        if args.separatefiles:
          outfile = open(pdbid+'_cablam.csv','w')
          csv_header(kinorder,args.doconnections,writeto=outfile)
          csv_print(resdata, kinorder, skiplist, inclist, args.doconnections, writeto=outfile)
          outfile.close()
        else:
          csv_print(resdata,kinorder,skiplist,inclist,args.doconnections)

  if outfiles:
    for filename in outfiles:
      outfiles[filename].close()
    #---------------------------------------------------------------------------
    #}}}
#-------------------------------------------------------------------------------
#}}}

#{{{ "__main__"
#-------------------------------------------------------------------------------
if __name__ == "__main__":
  #__main__ needs to generate a hierarchy
  run()
#-------------------------------------------------------------------------------
#}}}
