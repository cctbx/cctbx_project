from __future__ import absolute_import, division, print_function
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
#2012-09-05 cablam_training can now run probe if precomputed probe files are not
#  provided. Argument parsing has been updated to libtbx.phil for phenix
#  compatibility. cablam=True now yields CA_d_in, CA_d_out, and (instead of
#  CA_a) CO_d_in.  usage() help message added.
#2012-12-04: Added cis_or_trans argument for selecting cis or non-cis peptides
#  during printing. Default returns all residues.
#2013-09-17: Major fingerprints rewrite to use new fingerprints objects/methods
#  /storage.  See cablam_fingerprints.py and the fingerprints dir.
#  add_probe_data now stores full 4-character pdb-style atom names
#  New output: probe_mode=sequence will print amino acid sequence for motif
#2014_02_07: Updates to probe output methods to match changes in
#  cablam_fingerprints. Motifs without continuous sequence now supported for all
#  probe outputs
#Next: iotbx.file_reader incorporated to control input
#To do: Collect cis-peptides for analysis. Are they identifiable in cablamspace?
#  Add clash filtering. 0.4 is sufficient clash to cull, mc-mc are the important
#  contacts, at least for base cablam

import os, sys
from iotbx import pdb  #contains the very useful hierarchy
from mmtbx.cablam import cablam_res #contains a data structure derived from
#  hierarchy, but more suited to cablam's needs - specifically it can hold
#  geometric and probe measures and can look forward and backward in sequence
from mmtbx.cablam import cablam_math #contains geometric measure calculators
#from mmtbx.cablam import fingerprints #contains motif definitions
from mmtbx.cablam import cablam_fingerprints
#import cablam_fingerprints
#  Storage for motif definitions subject to change
from libtbx import easy_run
import libtbx.phil.command_line
from iotbx import file_reader
from libtbx import group_args

#{{{ phil
#-------------------------------------------------------------------------------
master_phil = libtbx.phil.parse("""
cablam_training {
  file_or_dir = None
    .type = path
    .help = '''input pdb file or dir thereof'''
  separate_files = False
    .type = bool
    .help = '''Generate a separate, auto-named output file for each input file'''
  give_kin = False
    .type = bool
    .help = '''Print output to screen in .kin format (default is comma-separated .csv format)'''
  give_connections = False
    .type = bool
    .help = '''Add prevres and nextres columns to .csv output'''
  debug = False
    .type = bool
    .help = '''Adds some text printed to stderr for debugging esp. for fingerprints'''

  all_measures = False
    .type = bool
    .help = '''Does all measures'''
  cad = False
    .type = bool
    .help = '''2 CA pseudo dihedrals'''
  caa = False
    .type = bool
    .help = '''3 CA pseudo angles'''
  cod = False
    .type = bool
    .help = '''2 CO pseudo dihedrals'''
  rama = False
    .type = bool
    .help = '''2 Ramachandran dihedrals: phi, psi'''
  exrama = False
    .type = bool
    .help = '''4 Ramachandran dihedrals: psi-1, phi, psi, phi+1'''
  tau = False
    .type = bool
    .help = '''1 backbone angle: tau (defined by N-CA-C)'''
  omega = False
    .type = bool
    .help = '''1 backbone dihedral: omega (defined by CA_1-C_1-N_2-CA_2)'''
  cablam = False
    .type = bool
    .help = '''Shortcut for just cablam-relevant measures CA_d_in, CA_d_out, CO_in'''

  probe_motifs = None
    .type = strings
    .help = '''Activates hydrogen bonding analysis, probe=motif_name1,motif_name2,... use --listmotifs to list available fingerprints'''
  probe_path = None
    .type = path
    .help = '''Stores path to dir of probed files, probe will be called for each file if this is not provided'''
  probe_mode = *kin annote instance sequence superpose
    .type = choice
    .help = '''=kin for dotlist kins (default) =annote for ball on model, =instance for vectorlist kins'''
  list_motifs = False
    .type = bool
    .help = '''print motifs/fingerprints available to screen'''

  b_max = None
    .type = float
    .help = '''Set a max b factor, residues containing a backbone atom with higher b will be pruned, recommended: -b=30'''
  prune_alts = False
    .type = bool
    .help = '''Removes all residues with alternate conformations in relevant atoms'''
  prune = None
    .type = strings
    .help = '''List of restypes to be pruned, separated by commas, no spaces eg PRO'''
  skip_types = None
    .type = strings
    .help = '''List of restypes to be skipped during printing, separated by commas'''
  include_types = None
    .type = strings
    .help = '''List of restypes to be printed, all others will be skipped'''

  cis_or_trans = *both cis trans
    .type = choice
    .help = '''selects whether cis-peptides, trans-peptides, or both will be returned'''

  fear = False
    .type = bool
    .help = '''turns on fear-to-tread analysis (this is temporary)'''

  help = False
    .type = bool
    .help = '''print help text to screen'''
}
""", process_includes=True)
#-------------------------------------------------------------------------------
#}}}

#{{{ usage notes
#-------------------------------------------------------------------------------
def usage():
  sys.stderr.write("""
phenix.cablam_training or cablam_training.py is a program intended for the
  exploration of protein structure datasets, the annotation of motifs of
  interest, and the training of reference datasets.  It was used in the
  construction of the reference contours used by cablam_validate.  It contains a
  number of features and modes and is intended primarily as a development tool
  rather than a utility for typical users.  However, anyone interested in
  exploring protein backboen geometry may find something of use here.

--------------------------------------------------------------------------------
file_or_dir=*path*
  Path to a pdb file or dir of pdb files to operate on, the only argument that
  doesn't need an explicit flag
--------------------------------------------------------------------------------

-----Basic Printing Options-----------------------------------------------------
separate_files=True/False
  Generate a separate, auto-named output file in the current dir for each input
  file, default output prints a single file to screen

give_kin=True/False
  Print output to screen in .kin format, may be combinded with separate_files,
  default output prints comma-separated .csv format

give_connections=True/False
  If set to True, adds prevres and nextres columns to .csv output

skip_types=restype1,restype2
include_types=restype3,restype4
  Together, these control which residue types are printed to screen or file.
  Default prints all residues.
  Residue types and relationships given to skip_types are excluded from printing
  If only include_types is used, only the listed restypes will be printed
  If include_types and skip_types are both used, then the types given to
  include_types will override those skipped by skip_types.

  List restypes by their 3-letter code and separated by commas withoug spaces,
  e.g. GLY,PRO,ALA,TRP
  Sequence relationships may be represented with underscores, e.g. _PRO is
  pre-proline, and GLY__ (2 underscores) is post-post-glycine

  examples:
  skip_types=PRO would print every residue except proline
  include_types=PRO,GLY would print *only* glycines and prolines
  skip_types=_PRO include_types=GLY would skip pre-prolines unless they were
  also glycines

cis_or_trans='cis' 'trans' 'both'
  Selects printing for cis-peptides or trans-peptides exclusively. The default
  is 'both' which will print all residues.  cis is defined as -60 to +60 degrees
  trans is defined as 120 to 180 and -120 to -180 degrees for the omega dihedral
  Note that selecting 'cis' or 'trans' will also stop printing for any residue
  for which omega cannot be calculated.
--------------------------------------------------------------------------------

-----Probe and Motif Search Options---------------------------------------------
This is an alternate mode which searches for hydrogen bonding patterns defined
  in fingerprints.

probe_motifs=motif_name1,motif_name2
  This flag activates hydrogen bonding pattern analysis, which will not run
  otherwise.  The flag accepts a spaceless string of comma-separated motif names
  to search for.  Use list_motifs=True to get a list of available motifs.

probe_path=*path*
  cablam_training can use precomputed probe results to speed up runs on large
  datasets. If a path to such prepared files is not provided, Reduce and Probe
  will be run on each pdb file, which may be time-consuming.

  Running:
  phenix.probe -u -condense -self -mc -NOVDWOUT -NOCLASHOUT MC filename.pdb > filename.probe
  Should produce appropriately formatted and named files for this option

probe_mode=kin/annote/instance/sequence
  These are printing options for hydrogen bond pattern analysis, which overrides
  the Basic Printing Options above.
Choose 1 of 3:
=kin returns automatically-named kinemage files, one for each unique member
  residue in each motif. The kins are high-dimensional dotlists containing the
  measures specified in the commandline (see below for options).  This is the
  default printing.
=annote returns an automatically-named kinemage file for each pdb file. These
  kins are balllists that highlight the selected motifs of interest if appended
  to existing kinemages of the structures.
=instance returns an automatically-named vectorlist kinemage file for each motif
  of interest.  Each kin is a high-dimensional vectorlist that shows the path of
  a multi-residue motif through the measures specified in the commandline
  (see below for options)
=sequence prints to screen the animo acid sequence of the motif of interest.
  Does not behave with multiple motifs.  Uses single-letter amino acid codes, if
  a residue type is unrecognized, will print 'X' followed by the 3-letter code.

list_motifs=True/False
  Prints to screen a list of all the motifs/"fingerprints" currently available
  for hydrogen bond pattern search
--------------------------------------------------------------------------------

-----Geometric Measures---------------------------------------------------------
All of these default to False, and some output modes will not function unless at
  least one of these options is turned on.  When in doubt, cablam=True and/or
  rama=True will provide relevant information.

cad=True/False
  For each residue, calculate the 2 C-alpha pseudo dihedrals

caa=True/False
  For each residue, calculate the 3 C-alpha pseudo angles

cod=True/False
  For each residue, calculate the 2 carbonyl oxygen pseudo dihedrals

rama=True/False
  For each residue, calculate Ramachandran dihedrals phi and psi

exrama=True/False
  For each residue, calculate Ramachandran dihedrals psi-1, phi, psi, phi+1

tau=True/False
  For each residue, calculate backbone angle tau, defined by N-NA-C

omega=True/False
  For each residue, calculate backbone peptide dihedral,
  defined by CA_1,C_1,N_2,CA_2

all_measures=True/False
  For each residue, calculate all of the above measures (may be overkill)

cablam=True/False
  Recommended, but not default behavior.
  For each residue calculate the measures most relevant to cablam analysis:
  CA_d_in, CA_d_out, CO_in
--------------------------------------------------------------------------------

-----Quality Control Options----------------------------------------------------
b_max=#.#
  Set a max b factor value. Residues containing a backbone atom with higher b
  will be pruned and excluded from all calculations. Note this may affect
  neighboring residues.  Strongly Recommenced: b_max=30.0

prune_alts=True/False
  Prune and excludes from calculations all residues with alternate conformations
  for backbone atoms. Note this may affect neighboring residues. Default is
  prune_alts=False, which results in only the first alternate position for each
  residue being reported on.

prune=restype1,restype2
  Prune  and exclude from calculations the selected list of residue types. Note
  this may affect neighboring residues. Restypes should be given as 3-letter
  codes, e.g. GLY,PRO, but this option does not yet support the sequence
  relationship that skip_types= and include_types= do.
--------------------------------------------------------------------------------

-----Help Options---------------------------------------------------------------
help=True/False
  Displays this help message.

list_motifs=True/False
  Prints to screen a list of all the motifs/"fingerprints" currently available
  for hydrogen bond pattern search

debug=True/False
  Activates print-to-stderr debugging notes for hydrogen bond pattern search.
  This may be valuable when trying to define a new pattern correctly and with
  proper format.
--------------------------------------------------------------------------------

Examples:
phenix.cablam_training cad=True cod=True skip_types=GLY,PRO,_PRO,ILE,VAL b_max=30.0 kin=True file_or_dir=path/pdbfilename.pdb

phenix.cablam_training cablam=True b_max=30.0 prune=GLY probe_motifs=parallel_beta,antiparallel_beta_cwc,antiparallel_beta_wcw probe_mode=kin probe_path=path/database/probefiles file_or_dir=path/database/pdbfiles

""")
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

#{{{ fails cis check function
#Allows cis or trans peptides to be skipped during printing. Passing
#  cis_or_trans='both' will print all residues. Residues without an omega value
#  will be skipped unless cis_or_trans=='both'.
#As with pruning, important in training, less so in annotation.
#-------------------------------------------------------------------------------
def fails_cis_check(residue,cis_or_trans):
  doskip = True
  if cis_or_trans == 'both':
    doskip = False
  else:
    if 'omega' not in residue.measures:
      doskip = True
    else:
      omega = residue.measures['omega']
      if cis_or_trans == 'cis' and (omega >= -30 and omega <= 30):
        doskip = False
      if cis_or_trans == 'trans' and (omega >= 150 or omega <= -150):
        doskip = False

  return doskip
#-------------------------------------------------------------------------------
#}}}

#{{{ make probe data function
#If a precomputed probe file has not been provided, this function calls probe to
#  generate appropriate data for use in add_probe_data()
#-------------------------------------------------------------------------------
def make_probe_data(hierarchy):
  trim_command = "phenix.reduce -quiet -trim -"
  build_command = "phenix.reduce -oh -his -flip -pen9999 -keep -allalt -"
  #probe_command = "phenix.probe -u -condense -self -mc -NOVDWOUT -NOCLASHOUT MC -"
  probe_command = "phenix.probe -u -condense -self -mc -NOVDWOUT -NOCLASHOUT ALL -"

  for i,m in enumerate(hierarchy.models()):
    #multi-model compatibility coming soon?
    #probe doesn't keep model data, so add_probe_data doesn't handle that
    #so this just takes the first model
    model = m
    break
  r = pdb.hierarchy.root()
  mdc = model.detached_copy()
  r.append_model(mdc)

  sys.stderr.write('  cleaning . . .\n')
  clean_out = easy_run.fully_buffered(trim_command, stdin_lines=r.as_pdb_string())
  sys.stderr.write('  reducing . . .\n')
  build_out = easy_run.fully_buffered(build_command, stdin_lines=clean_out.stdout_lines)
  #print build_out.stdout_lines
  input_str = '\n'.join(build_out.stdout_lines)
  sys.stderr.write('  probing . . .\n')
  probe_out = easy_run.fully_buffered(probe_command, stdin_lines=input_str)
  #print '\n'.join(probe_out)

  #print '\n'.join(probe_out.stdout_lines)
  return probe_out.stdout_lines
#-------------------------------------------------------------------------------
#}}}

#{{{ add probe data function
#Adds mainchina-mainchain hydrogen bonding information from 'unformated' Probe
#  output to a dictionary of residues.
#At the moment, reliant on precomputed .probe files, will gain run-time Probe
#May gain other contact relationship info, by mc-mc H-bonds are most important
#-------------------------------------------------------------------------------
def add_probe_data(resdata, open_probe_file):
  #print open_probe_file
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
    srcAtomname = srcAtom[11:15]#.strip()
    srcAlt =      srcAtom[15:16].strip()

    trgAtom = bnana[4]
    #going to count dots per bond as a measure of strength instead
    trgChain =    trgAtom[0:2].strip()
    trgNum =      int(trgAtom[2:6].strip())
    trgNumStr =   trgAtom[2:6]
    trgIns =      trgAtom[6:7]#.strip()
    trgResname =  trgAtom[7:10].strip()
    #if trgResname == 'HOH': continue #skip waters
    trgAtomname = trgAtom[11:15]#.strip()
    trgAlt =      trgAtom[15:16].strip()

    dotcount = bnana[5]
    mingap = bnana[6]

    #new model for probe storage------------------------------------------------
    # If targ is not in resdata then it is likely a water or hetgroup. However,
    # we want to have a record of the hb info. In this case 'residue' in 'record'
    # will be an object with chain, resnum, resname, and icode.
    # If src is not in resdata then we arn't interested.
    src_key = ' '.join(['', srcChain, '%04i' % srcNum, srcIns])
    if src_key not in resdata.keys() : continue
    srcResidue = resdata[src_key]
    targ_key = ' '.join(['', trgChain, '%04i' % trgNum, trgIns])
    if targ_key not in resdata.keys():
      continue
      #trgResidue = group_args(chain   = trgChain,
      #                        resnum  = trgNum,
      #                        resname = trgResname,
      #                        icode   = trgIns)
      #recordkey  =  trgResname +' '+trgChain + trgNumStr + trgIns + trgAtomname
    else:
      trgResidue = resdata[targ_key]
      recordkey = trgResidue.id_with_resname() + trgAtomname
    record = group_args(residue  = trgResidue,
                        atom     = trgAtomname,
                        dotcount = dotcount,
                        mingap   = mingap,
                        seqdist  = srcResidue.seq_dist(trgResidue))
    if srcAtomname not in srcResidue.probe.keys():
      srcResidue.probe[srcAtomname] = {}
    #####srcResidue = resdata[' '.join(['', srcChain, '%04i' % srcNum, srcIns])]
    #####trgResidue = resdata[' '.join(['', trgChain, '%04i' % trgNum, trgIns])]
    #####recordkey = trgResidue.id_with_resname() + trgAtomname
    #####record = group_args(residue=trgResidue, atom=trgAtomname, mingap=mingap, seqdist=srcResidue.seq_dist(trgResidue))
    ######print [trgResidue.id_with_resname(),trgAtomname,dotcount,srcResidue.seq_dist(trgResidue)]
    #####if srcAtomname not in srcResidue.probe.keys():
    #####  srcResidue.probe[srcAtomname] = {}
    #probe keys first by the current residue's atom, then by the target
    #  residue's id and atom, id+atom is unique enough to handle bifurcations
    srcResidue.probe[srcAtomname][recordkey] = record
    #end new model for probe storage--------------------------------------------
    #reference: resid_string = ' '.join([modelid,chainid,'%04i' % resnum,icode])
#-------------------------------------------------------------------------------
#}}}

#{{{ Output function collection
#A collection of headers, formatting, and printint functions used in output
#Default output is to stdout, but anything with a .write can be passed to the
#  'writeto=' argument of most functions.  Functions that lack a 'writeto='
#  generate or find uniquely named files in the working dir for their output.
#Print methods called by these functions are generally from cablam_res.py
#-------------------------------------------------------------------------------

#{{{ --- kin_frame ---
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
#-------------------------------------------------------------------------------
#}}}

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
  doconnections=False, cis_or_trans='both', writeto=sys.stdout):
  reslist = protein.keys()
  reslist.sort()
  for resid in reslist:
    if skipcheck(protein[resid], skiplist, inclist):
      pass
    elif fails_cis_check(protein[resid],cis_or_trans):
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
  if len(kinorder) == 0:
    sys.stderr.write('\nNo geometric measures (e.g. rama=True) specified')
    sys.stderr.write('\nExiting . . .\n')
    sys.exit()
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
def kin_print(protein, kinorder, skiplist=[], inclist=[], cis_or_trans='both',
  writeto=sys.stdout):
  if len(kinorder) == 0:
    sys.stderr.write('\nNo geometric measures (e.g. rama=True) specified')
    sys.stderr.write('\nExiting . . .\n')
    sys.exit()
  reslist = protein.keys()
  reslist.sort()
  for resid in reslist:
    if skipcheck(protein[resid], skiplist, inclist):
      pass
    elif fails_cis_check(protein[resid],cis_or_trans):
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
  if len(kinorder) == 0:
    sys.stderr.write('\nNo geometric measures (e.g. rama=True) specified')
    sys.stderr.write('\nExiting . . .\n')
    sys.exit()
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
def kin_print_probe(motif_instances, kinorder, outfiles):
  for motif_name in motif_instances:
    for instance in motif_instances[motif_name]:
      if not instance.has_all_measures(kinorder):
        sys.stderr.write(
          '  '+motif_name+' has incomplete measures, probably due to b_max\n')
        continue
      for index in instance.names:
        residue = instance.residues[index]
        name = instance.names[index]
        residue.printtokin(kinorder, writeto=outfiles[name])
#-------------------------------------------------------------------------------
#}}}

#{{{ --- PROBE ANNOTE printing ---
#-------------------------------------------------------------------------------
#For annotating an existing .kin file with balls at CA's participating in
#  motifs of interest.
#Produces one .kin per input file.
#Does not require a header as such.
def kin_print_probe_annote(motif_instances, writeto=sys.stdout):
  for motif_name in motif_instances:
    if motif_instances[motif_name]:
      writeto.write('@group {'+motif_name+'}\n')
      ref_instance = motif_instances[motif_name][0]
      indices = ref_instance.residues.keys()
      indices.sort()
      for index in indices:
        writeto.write('@balllist {'+ref_instance.names[index]+'}\n')
        for instance in motif_instances[motif_name]:
          residue = instance.residues[index]
          firstalt = residue.firstalt('CA')
          CAxyz = residue.atomxyz[firstalt]['CA']
          pointid = residue.pdbid+' '+ residue.chain +' '+ str(residue.resnum)+' '+ instance.names[index]
          writeto.write("{ "+pointid+" } "+str(CAxyz[0])+" "+str(CAxyz[1])+" "+str(CAxyz[2])+"\n")

#{{{
def old_kin_print_probe_annote(resdata, motif_list, writeto=sys.stdout):
  reskeys = resdata.keys()
  reskeys.sort()
  motifs = cablam_fingerprints.fetch_fingerprints(motif_list)
  for motif in motifs:
    writeto.write('@group {'+motif.motif_name+'}\n')
    for label in motif.residue_names.values():
      writeto.write('@balllist {'+label+'}\n')
      for resid in reskeys:
        residue = resdata[resid]
        if label in residue.motifs:
          firstalt = residue.firstalt('CA')
          #try:
          CAxyz = residue.atomxyz[firstalt]['CA']
          pointid = residue.pdbid +' '+ residue.chain +' '+ str(residue.resnum)+' '+ label
          writeto.write("{ "+pointid+" } "+str(CAxyz[0])+" "+str(CAxyz[1])+" "+str(CAxyz[2])+"\n")
#}}}
#-------------------------------------------------------------------------------
#}}}

#{{{ --- PROBE BY INSTANCE printing ---
#-------------------------------------------------------------------------------
#Creates files and prints headers in them for instance output
#One .kin for each motif. This can produce several files.
def kin_print_by_instance_header(motif_list, kinorder, kinranges):
  if len(kinorder) == 0:
    sys.stderr.write('\nNo geometric measures (e.g. rama=True) specified')
    sys.stderr.write('\nExiting . . .\n')
    sys.exit()
  outfiles = {}
  motifs = cablam_fingerprints.fetch_fingerprints(motif_list)
  for motif in motifs:
    motif_name = motif.motif_name
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
def kin_print_by_instance(motif_instances, motif_list, kinorder, outfiles):
  for motif_name in motif_instances:
    for instance in motif_instances[motif_name]:
      if not instance.has_all_measures(kinorder):
        sys.stderr.write(
          '  '+motif_name+' has incomplete measures, probably due to b_max\n')
        continue
      indices = instance.names.keys()
      indices.sort()
      #print indices
      residue = instance.residues[indices[0]]
      outfiles[motif_name].write(
        '@group {'+residue.pdbid.rstrip('.pdb')+' '+str(residue.resnum)+
        '} dominant animate\n@vectorlist {'+motif_name+
        '} dimension='+str(len(kinorder))+'\n')
      for index in indices:#instance.names:
        residue = instance.residues[index]
        name = instance.names[index]
        outline = ['{'+residue.id_with_resname()+'_'+name+'}']
        for order in kinorder:
          outline.append(str(residue.measures[order]))
        outfiles[motif_name].write(' '.join(outline)+'\n')

#print a string of 1-char resnames for each motif instance,
#  for use with WebLogo and the like.
def res_seq_by_instance(motif_instances):
  #reshash contains the standard 3char to 1char mappings, followed by an ever-
  #  growing list of non-standard animo acids
  reshash = {'GLY':'G','ALA':'A','VAL':'V','ILE':'I','LEU':'L','PHE':'F',
  'TRP':'W','MET':'M','GLU':'E','GLN':'Q','ASP':'D','ASN':'N','SER':'S',
  'THR':'T','TYR':'Y','HIS':'H','LYS':'K','PRO':'P','CYS':'C','ARG':'R',
  'MSE':'M','SME':'M','CSO':'C','OCS':'C','CSX':'C','CME':'C','YCM':'C',
  'MLY':'K'}
  for motif_name in motif_instances:
    for instance in motif_instances[motif_name]:
      indices = instance.residues.keys()
      indices.sort()
      seq_string = []
      for index in indices:
        resname = instance.residues[index].id_with_resname()[0:3]
        if resname in reshash:
          code = reshash[resname]
        else:
          #non-standard amino acids not already handled can be found in the
          #  output by searching for 'X'
          code = 'X'+resname
        seq_string.append(code)
      seq_string.append('\n')
      sys.stdout.write(''.join(seq_string))
#-------------------------------------------------------------------------------
#}}}

#{{{ --- PROBE superposition ---
#-------------------------------------------------------------------------------
#First step: excise the relevant bits of each pdb file
def trim_motifs(motif_instances, filename, superpose_refs):
  pwd = os.getcwd()
  for motif_name in motif_instances:
    if os.path.isdir(motif_name): pass
    else: os.mkdir(motif_name)
    os.chdir(motif_name)

    instance_num = 0
    for instance in motif_instances[motif_name]:
      instance_num += 1
      outputfile = os.path.basename(filename) + "_" + str(instance_num) + ".pdb"
      resnums = []
      for residue in instance.residues.values():
        resnum = str(residue.resnum)
        resnums.append(resnum)
      selection = "resseq "+ " or resseq ".join(resnums)
      command = 'phenix.pdbtools stop_for_unknowns=False modify.keep=\"'+selection+'\" '+filename + " output.file_name=" + outputfile
      #output.file_name=*****
      #sys.stderr.write(command)
      runthis = easy_run.fully_buffered(command)
      if motif_name not in superpose_refs:
        superpose_refs[motif_name] = {"motif":instance,"filename":outputfile}
      else:
        sys.stderr.write("trying to superpose\n")
        ref = superpose_refs[motif_name]
        #phenix.superpose_pdbs fixed.pdb moving.pdb selection_fixed="name CA" selection_moving="name CA"
        command = "phenix.superpose_pdbs "+ ref["filename"] + " " + outputfile + " selection_default_fixed="+ref["motif"].superpose_thus +" selection_default_moving="+instance.superpose_thus + " output.file_name=" + outputfile
        sys.stderr.write(command)
        sys.stderr.write("\n")
        runthis = easy_run.fully_buffered(command)
    os.chdir(pwd)
  return superpose_refs
#-------------------------------------------------------------------------------
#}}}

#-------------------------------------------------------------------------------
#}}}

#{{{ run
#The run function is currently rather messy. (Indeed, all of
#  cablam_training is a bit messy, as it's really a development tool, not a
#  general-use program.) Hopefully, everything needed for general use (structure
#  annotation) has been packaged in other modules for easy access. Good luck.

def run(args):
  #{{{ phil parsing
  #-----------------------------------------------------------------------------
  interpreter = libtbx.phil.command_line.argument_interpreter(master_phil=master_phil)
  sources = []
  for arg in args:
    if os.path.isfile(arg):
      input_file = file_reader.any_file(arg)
      if (input_file.file_type == "pdb"):
        sources.append(interpreter.process(arg="file_or_dir=\"%s\"" % arg))
      elif (input_file.file_type == "phil"):
        sources.append(input_file.file_object)
    elif os.path.isdir(arg):
      sources.append(interpreter.process(arg="file_or_dir=\"%s\"" % arg))
    else:
      arg_phil = interpreter.process(arg=arg)
      sources.append(arg_phil)
  work_phil = master_phil.fetch(sources=sources)
  work_params = work_phil.extract()
  params = work_params.cablam_training
  #catch missing file or dir later?
  #if not work_params.cablam_training.file_or_dir:
  #  usage()
  #  sys.exit()
  params = work_params.cablam_training
  #-----------------------------------------------------------------------------
  #}}} end phil parsing

  if params.help:
    usage()
    sys.exit()

  if params.list_motifs:
    sys.stdout.write('\n')
    fileset = os.listdir(libtbx.env.find_in_repositories(
      "cctbx_project/mmtbx/cablam/fingerprints"))
    for filename in fileset:
      if filename.endswith(".pickle"):
        motifname = os.path.splitext(os.path.basename(filename))[0]
        sys.stdout.write(motifname + '\n')
    sys.exit()

  if not params.file_or_dir:
    usage()
    sys.exit()
  if os.path.isdir(params.file_or_dir):
    fileset = os.listdir(params.file_or_dir)
    dirpath = params.file_or_dir
  elif os.path.isfile(params.file_or_dir):
    fileset = [params.file_or_dir]
    dirpath = None
  else:
    sys.stderr.write("Could not identify valid target file or dir.\n")
    usage()
    sys.exit()

  #{{{ measurement selection
  #This section manages the user's orders for calculations
  #Note: The 'kin' in kinorder and kin ranges is a misnomer
  #-----------------------------------------------------------------------------
  if params.all_measures:
    params.cad = True
    params.caa = True
    params.cod = True
    params.exrama = True
    params.tau = True
    params.omega = True

  kinorder, kinranges = [],[]
  if params.cad:
    kinorder.append('CA_d_in'),  kinranges.append('-180 180')
    kinorder.append('CA_d_out'), kinranges.append('-180 180')
  else:
    pass

  if params.cod:
    kinorder.append('CO_d_in'),  kinranges.append('-180 180')
    kinorder.append('CO_d_out'), kinranges.append('-180 180')
  else:
    pass

  if params.caa:
    kinorder.append('CA_a_in'),  kinranges.append('0 180')
    kinorder.append('CA_a'),     kinranges.append('0 180')
    kinorder.append('CA_a_out'), kinranges.append('0 180')
  else:
    pass

  if params.cablam:
    if 'CA_d_in' not in kinorder:
      kinorder.append('CA_d_in'),  kinranges.append('-180 180')
    if 'CA_d_out' not in kinorder:
      kinorder.append('CA_d_out'),  kinranges.append('-180 180')
    if 'CO_d_in' not in kinorder:
      kinorder.append('CO_d_in'),  kinranges.append('-180 180')
    if 'CA_a' not in kinorder:
      kinorder.append('CA_a'), kinranges.append('0, 180')
  else:
    pass

  if params.rama or params.exrama:
    if params.exrama:
      kinorder.append('psi-1'), kinranges.append('-180 180')
      kinorder.append('phi'),   kinranges.append('-180 180')
      kinorder.append('psi'),   kinranges.append('-180 180')
      kinorder.append('phi+1'), kinranges.append('-180 180')
    else:
      kinorder.append('phi'), kinranges.append('-180 180')
      kinorder.append('psi'), kinranges.append('-180 180')
  else:
    pass

  if params.tau:
    kinorder.append('tau'), kinranges.append('0 180')
  else:
    pass

  if params.omega:
    kinorder.append('omega'), kinranges.append('-180 180')
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
  superpose_refs = {}

  outfiles = {}
  if params.probe_motifs:
    motif_list = params.probe_motifs[0].split(',')
    if params.probe_path:
      probefilelist = os.listdir(params.probe_path)
    if params.probe_mode == 'kin':# or params.probe_mode == None:
      outfiles = kin_print_probe_header(cablam_fingerprints.get_all_labels(motif_list),kinorder,kinranges)
    elif params.probe_mode == 'instance':
      outfiles = kin_print_by_instance_header(motif_list, kinorder, kinranges)

  prunelist = []
  if params.prune:
    prunelist = params.prune[0].split(',')
    prunelist = [res.upper() for res in prunelist] #Ha ha! List comprehension!

  skiplist = []
  inclist = []
  if params.skip_types:
    skiplist = params.skip_types[0].split(',')
  if params.include_types:
    inclist = params.include_types[0].split(',')

  if params.separate_files:
    pass
  else:
    if params.give_kin:
      kin_header(kinorder,kinranges)
    elif params.probe_motifs:
      pass
    else:
      csv_header(kinorder,params.give_connections)
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ get file, start loop
  #-----------------------------------------------------------------------------
  for filename in fileset:
    #if not filename.endswith('.pdb'):
    #  continue
    if dirpath: #must add the path if using the listed contents of a dir
      filename = os.path.join(dirpath,filename)
    else:
      pass

    pdbid = os.path.basename(filename)

    if not os.path.isfile(filename): continue
    pdb_in = file_reader.any_file(filename)
    if pdb_in.file_type != "pdb":
      sys.stderr.write(filename +" not id'd as readable file\n")
      continue
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
    cablam_res.prunerestype(resdata, 'HOH')
    for restype in prunelist:
      cablam_res.prunerestype(resdata, restype)

    if params.b_max:
      stripB(resdata,params.b_max)

    if params.prune_alts:
      prune_alts(resdata)
    #---------------------------------------------------------------------------
    #}}}

    #{{{ calculation calls
    #---------------------------------------------------------------------------
    if params.cad and params.caa:
      cablam_math.CApseudos(resdata, dodihedrals = True, doangles = True)
    elif params.cad:
      cablam_math.CApseudos(resdata, dodihedrals = True, doangles = False)
    elif params.caa:
      cablam_math.CApseudos(resdata, dodihedrals = False, doangles = True)
    else: #no CA-based calculations
      pass

    if params.cod:
      cablam_math.COpseudodihedrals(resdata)
    else:
      pass

    if params.rama or params.exrama:
      cablam_math.phipsi(resdata)
    else:
      pass

    if params.tau:
      cablam_math.taucalc(resdata)
    else:
      pass

    if params.omega or params.cis_or_trans != 'both':
      cablam_math.omegacalc(resdata)
    else:
      pass

    if params.cablam:
      cablam_math.cablam_measures(resdata)
    else:
      pass
    #---------------------------------------------------------------------------
    #}}}

    #{{{ probe stuff
    #---------------------------------------------------------------------------
    #need the run phenix.probe
    if params.probe_motifs and params.probe_path:
      probefilename = pdbid.rstrip('.pdb') + '.probe'
      if probefilename in probefilelist:
        probefilepath = os.path.join(params.probe_path,probefilename)
        open_probe_file = open(probefilepath)
        add_probe_data(resdata,open_probe_file)
        open_probe_file.close()
      else:
        continue
    elif params.probe_motifs:
      add_probe_data(resdata,make_probe_data(hierarchy))

    if params.probe_motifs:
      found_motifs = cablam_fingerprints.check_protein(resdata, motif_list)
      #found_motifs is a dictionary. The keys are motif names.
      #  The values are lists of cablam_fingerprints.motif_instance objects.
    #---------------------------------------------------------------------------
    #}}}

    #{{{ output
    #---------------------------------------------------------------------------
    #--probemode=kin for dotlist kins, this is the default
    #--probemode=annote for balls drawn at CA positions on the model
    #--probemode=instance for kins where each veclist is one instance of motif
    if params.probe_motifs:# and args.probepath:
      if params.probe_mode == 'kin':# or params.probe_mode == None:
        kin_print_probe(found_motifs, kinorder, outfiles)
      elif params.probe_mode == 'annote':
        outfile = open(pdbid+'cablam_motifs.kin','w')
        #kin_print_probe_annote(resdata, motif_list, writeto=outfile)
        kin_print_probe_annote(found_motifs, writeto=outfile)
        outfile.close()
      elif params.probe_mode == 'instance':
        #kin_print_by_instance(resdata, motif_list, kinorder, outfiles)
        kin_print_by_instance(found_motifs, motif_list, kinorder, outfiles)
      elif params.probe_mode == 'sequence':
        res_seq_by_instance(found_motifs)
        #res_seq_by_instance(resdata, motif_list)
      elif params.probe_mode == 'superpose':
        #trim_motifs(resdata, filename, motif_list)
        superpose_refs = trim_motifs(found_motifs, filename,superpose_refs)
        #superpose_motifs(motif_list)
      else:
        sys.stderr.write('\n\nUnrecognized probemode request\n\n')
        sys.exit()
      #add if args.kin once things basically work
      outfile = sys.stdout
      #need printer from probe version
      #Not sure what the stray outfile=sys.stdout is doing here anymore

    #default printing, with no arguments, is to .csv, one line per residue
    #--separatefiles writes a separate file for each input file to working dir
    #--kin prints kinemage file, dotlist, one point per residue
    #--doconnections adds connectivity information to csv output
    else:
      if params.give_kin:
        if params.separate_files:
          outfile = open(pdbid+'_cablam.kin','w')
          kin_header(kinorder,kinranges,writeto=outfile)
          kin_print(resdata, kinorder, skiplist, inclist, params.cis_or_trans, writeto=outfile)
          outfile.close()
        else:
          kin_print(resdata,kinorder,skiplist,inclist,params.cis_or_trans)
      else:
        if params.separate_files:
          outfile = open(pdbid+'_cablam.csv','w')
          csv_header(kinorder,params.give_connections,writeto=outfile)
          csv_print(resdata, kinorder, skiplist, inclist, params.give_connections,params.cis_or_trans, writeto=outfile)
          outfile.close()
        else:
          csv_print(resdata,kinorder,skiplist,inclist,params.give_connections,params.cis_or_trans,)

  if outfiles:
    for filename in outfiles:
      outfiles[filename].close()
    #---------------------------------------------------------------------------
    #}}}
#-------------------------------------------------------------------------------
#}}}
