import iotbx.phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
from mmtbx import amino_acid_codes
from mmtbx import alignment
import sys, os
from iotbx import pdb
from cctbx.array_family import flex
from scitbx.math import superpose, matrix

master_params = iotbx.phil.parse("""\
  supper{
    fixed=None
    .type=str
    moving=None
    .type=str
    moved=moved.pdb
    .type=str
    alignment_type=local *global
    .type=choice
  }
""")

def move_and_write(pdb_inp, r, t, outfile):
  for ser, lbl, atm in zip( pdb_inp.atom_serial_number_strings(),
                            pdb_inp.input_atom_labels_list(),
                            pdb_inp.atoms() ):
    xyz = r*matrix.col(atm.xyz) + t
    tmp = pdb.format_atom_record(
      serial=float(ser),
      name=lbl.name(),
      altLoc=lbl.altloc(),
      resName=lbl.resname(),
      chainID=lbl.chain(),
      site=xyz,
      occupancy=atm.occ,
      tempFactor=atm.b
      )
    print >> outfile, tmp


def extract_sequence_and_sites(pdb_input):

  seq = ""
  sites = flex.vec3_double()
  use_sites = flex.bool()

  model = pdb_input.construct_hierarchy().models()[0]
  for chain in model.chains() :
    selected_residues = chain.conformers()[0].residue_class_selection(
            class_name="common_amino_acid")
    residues = chain.conformers()[0].residues()
    for ires in selected_residues:
      resi = residues[ires]
      single = amino_acid_codes.one_letter_given_three_letter[
        resi.name[0:3] ]
      use=False
      xyz = (0,0,0)
      for atom in resi.atoms():
        if atom.name ==" CA ":
          xyz = atom.xyz
          use = True
      sites.append( xyz  )
      seq += single
      use_sites.append( use )
  return seq,sites,use_sites



def allign_sites( reference_sequence_coords,
                  aligned_sequence_coords,
                  matches,
                  reference_sites,
                  aligned_sites,
                  reference_use_flags,
                  aligned_use_flags):

  fixed_sites = flex.vec3_double()
  moving_sites = flex.vec3_double()

  assert len(reference_sequence_coords) == len(aligned_sequence_coords)
  assert len(reference_sequence_coords) == len(matches)

  n = len(reference_sequence_coords)
  #loop over the alignement
  match_id=["|","*"]
  for ii in xrange(n):
    if matches[ii] in match_id:
      assert reference_sequence_coords[ii] is not "-"
      assert aligned_sequence_coords[ii] is not "-"
      r=reference_sequence_coords[ii]
      a=aligned_sequence_coords[ii]
      if reference_use_flags[r]:
        if aligned_use_flags[a]:
          fixed_sites.append( reference_sites[r]  )
          moving_sites.append( aligned_sites[a]  )

  assert moving_sites.size() == fixed_sites.size()
  return fixed_sites, moving_sites



def run(args):

  phil_objects = []
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_params=master_params,
    home_scope="supper")

  for arg in args:
    command_line_params = None
    arg_is_processed = False
    if (os.path.isfile(arg)):
      try:
        command_line_params = iotbx.phil.parse(file_name=arg)
        if command_line_params is not None:
          phil_objects.append(command_line_params)
          arg_is_processed = True
      except KeyboardInterrupt: raise
      except : pass
    else:
      try:
        command_line_params = argument_interpreter.process(arg=arg)
        if command_line_params is not None:
          phil_objects.append(command_line_params)
          arg_is_processed = True
      except KeyboardInterrupt: raise
      except : pass

    if not arg_is_processed:
      raise Sorry("Unknown file or keyword: %s" % arg)

  scope_instance = master_params.fetch(sources=phil_objects)
  params = scope_instance.extract()
  print "#                       SUPPER                                 "
  print "#"
  print "#A lightweight sequence based structure superposition tool"
  print "#"
  print "#"
  print "#Parameters used:"
  print "#phil __ON__"
  scope_instance.show()
  print "#phil __OFF__"

  # reading in a pdb file, get sequence and sites
  fixed = pdb.input(file_name=params.supper.fixed)
  fixed_seq, fixed_sites, fixed_site_flags = extract_sequence_and_sites(fixed)

  moving = pdb.input(file_name=params.supper.moving)
  moving_seq, moving_sites, moving_site_flags = extract_sequence_and_sites(moving)

  #we now have to perform the allignment,
  #make a choice between local and global allignment
  alignment_object = None
  if params.supper.alignment_type == "local":
    alignment_object = alignment.align(
      moving_seq,
      fixed_seq,
      gop=20,gep=2,
      sim=alignment.blosum50,
      style=alignment.align.LOCAL  )

  if params.supper.alignment_type == "global":
    alignment_object = alignment.align(
      moving_seq,
      fixed_seq,
      gop=20,gep=2,
      sim=alignment.blosum50,
      style=alignment.align.GLOBAL )



  all = alignment_object.extract_alignment()
  crd = alignment_object.extract_alignment_coordinates()
  mtc = alignment.matches(all,alignment.blosum50)
  equal = ( mtc ).count("|")
  similar = (mtc  ).count("*")
  tot = len(all[0])-(all[0]).count("-")
  frac_iden = 100.0*float(equal)/float(tot)
  frac_simi = 100.0*float(equal+similar)/float(tot)

  alignment.pretty_print(all,
                         mtc,
                         block_size=50,
                         n_block=1,
                         top_name="fixed",
                         bottom_name="moving",
                         short_comment="""
The alignment used in the superposition is shown below.

The identity (fraction of |'s) between the two sequences
is %4.1f%s given over the aligned length of the fixed
molecule sequence.

The similarity (fraction of |'s and *'s) between the
sequences is equal to %4.1f%s over the aligned length
of the fixed molecule sequence.
                         """%(frac_iden,"%", frac_simi,"%") )

  assert fixed_sites.size() == len( fixed_seq )
  assert moving_sites.size() == len( moving_seq )

  fixed_sites_sel, moving_sites_sel = allign_sites(
    crd[0],
    crd[1],
    alignment.matches(all,alignment.blosum50),
    fixed_sites,
    moving_sites,
    fixed_site_flags,
    moving_site_flags)

  #now we can do the superposition
  lsq = superpose.least_squares_fit(
     fixed_sites_sel, moving_sites_sel) #, moving_sites_sel)
  #here we have the rotation and translation operators
  r = lsq.r
  t = lsq.t
  #we would like to know the rmsd on the coords used for superposition
  new_sites = lsq.other_sites_best_fit()
  deltas = fixed_sites_sel - new_sites
  rmsd = deltas.rms_length()

  print "Matching residues pairs indicated by | or * in the above alignment"
  print "have been used in a least square superposition."
  print "The rsmd between the aligned C-alpha atoms is: %5.3f"%(rmsd)
  print

  # write out the pdb file
  outfile = open(params.supper.moved,'w')
  move_and_write(moving,r,t,outfile)
  print "The shifted pdb file is written to the file: %s"%(params.supper.moved)

if __name__ == "__main__":
  run(sys.argv[1:])
