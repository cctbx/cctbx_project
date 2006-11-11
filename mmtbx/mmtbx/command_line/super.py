import mmtbx.alignment
import mmtbx.amino_acid_codes
import iotbx.pdb
from cctbx.array_family import flex
from scitbx.math import superpose
from scitbx.math import matrix
import libtbx.phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
import sys, os

master_params = libtbx.phil.parse("""\
super {
  fixed = None
    .type = str
  moving = None
    .type = str
  moved = "moved.pdb"
    .type = str
  alignment_style = *local global
    .type = choice
  gap_opening_penalty = 20
    .type = float
  gap_extension_penalty = 2
    .type = float
  similarity_matrix = *blosum50 dayhoff
    .type = choice
}
""")

def extract_sequence_and_sites(pdb_input):
  seq = []
  sites = flex.vec3_double()
  use_sites = flex.bool()
  model = pdb_input.construct_hierarchy().models()[0]
  for chain in model.chains():
    selected_residues = chain.conformers()[0].residue_class_selection(
      class_name="common_amino_acid")
    residues = chain.conformers()[0].residues()
    for ires in selected_residues:
      resi = residues[ires]
      resn = resi.name[0:3]
      single = mmtbx.amino_acid_codes.one_letter_given_three_letter[resn]
      seq.append(single)
      use = False
      xyz = (0,0,0)
      for atom in resi.atoms():
        if (atom.name == " CA "):
          xyz = atom.xyz
          use = True
          break
      sites.append(xyz)
      use_sites.append(use)
  return "".join(seq), sites, use_sites

def run(args, command_name="mmtbx.super"):
  if (len(args) == 0):
    print "usage: %s fixed.pdb moving.pdb [parameter=value ...]" % command_name
    return

  print "#"
  print "#                       ", command_name
  print "#"
  print "# A lightweight sequence-based structure superposition tool."
  print "#"
  print "#"

  phil_objects = []
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_params=master_params, home_scope="super")
  fixed_pdb_file_name = None
  moving_pdb_file_name = None
  for arg in args:
    if (os.path.isfile(arg)):
      if (fixed_pdb_file_name is None): fixed_pdb_file_name = arg
      elif (moving_pdb_file_name is None): moving_pdb_file_name = arg
      else: raise Sorry("Too many file names.")
    else:
      try: command_line_params = argument_interpreter.process(arg=arg)
      except KeyboardInterrupt: raise
      except: raise Sorry("Unknown file or keyword: %s" % arg)
      else: phil_objects.append(command_line_params)

  working_params = master_params.fetch(sources=phil_objects)
  params = working_params.extract()

  def raise_missing(what):
      raise Sorry("""\
Missing file name for %(what)s structure:
  Please add
    %(what)s=file_name
  to the command line to specify the %(what)s structure.""" % vars())

  if (fixed_pdb_file_name is None):
    if (params.super.fixed is None): raise_missing("fixed")
  else:
    params.super.fixed = fixed_pdb_file_name
  if (moving_pdb_file_name is None):
    if (params.super.moving is None): raise_missing("moving")
  else:
    params.super.moving = moving_pdb_file_name

  print "#Parameters used:"
  print "#phil __ON__"
  print
  working_params = master_params.format(python_object=params)
  working_params.show()
  print
  print "#phil __OFF__"
  print

  print "Reading fixed structure:", params.super.fixed
  fixed_pdb = iotbx.pdb.input(file_name=params.super.fixed)
  print
  print "Reading moving structure:", params.super.moving
  moving_pdb = iotbx.pdb.input(file_name=params.super.moving)
  print

  fixed_seq, fixed_sites, fixed_site_flags = extract_sequence_and_sites(
    pdb_input=fixed_pdb)
  moving_seq, moving_sites, moving_site_flags = extract_sequence_and_sites(
    pdb_input=moving_pdb)

  print "Computing sequence alignment..."
  align_obj = mmtbx.alignment.align(
    seq_a=fixed_seq,
    seq_b=moving_seq,
    gap_opening_penalty=params.super.gap_opening_penalty,
    gap_extension_penalty=params.super.gap_extension_penalty,
    similarity_function=params.super.similarity_matrix,
    style=params.super.alignment_style)
  print "done."
  print

  alignment = align_obj.extract_alignment()
  matches = alignment.matches()
  equal = matches.count("|")
  similar = matches.count("*")
  total = len(alignment.a) - alignment.a.count("-")
  alignment.pretty_print(
    matches=matches,
    block_size=50,
    n_block=1,
    top_name="fixed",
    bottom_name="moving",
    comment="""\
The alignment used in the superposition is shown below.

The sequence identity (fraction of | symbols) is %4.1f%%
of the aligned length of the fixed molecule sequence.

The sequence similarity (fraction of | and * symbols) is %4.1f%%
of the aligned length of the fixed molecule sequence.
""" % (100.*equal/max(1,total), 100.*(equal+similar)/max(1,total)))

  fixed_sites_sel = flex.vec3_double()
  moving_sites_sel = flex.vec3_double()
  for ia,ib,m in zip(alignment.i_seqs_a, alignment.i_seqs_b, matches):
    if (m not in ["|", "*"]): continue
    if (fixed_site_flags[ia] and moving_site_flags[ib]):
      fixed_sites_sel.append(fixed_sites[ia])
      moving_sites_sel.append(moving_sites[ib])

  print "Performing least-squares superposition of C-alpha atom pairs:"
  print "  Number of C-alpha atoms pairs in matching residues"
  print "  indicated by | or * above:", fixed_sites_sel.size()
  if (fixed_sites_sel.size() == 0):
    raise Sorry("No matching C-alpha atoms.")
  lsq_fit = superpose.least_squares_fit(
    reference_sites=fixed_sites_sel,
    other_sites=moving_sites_sel)
  rmsd = fixed_sites_sel.rms_difference(lsq_fit.other_sites_best_fit())
  print "  RSMD between the aligned C-alpha atoms: %.3f" % rmsd
  print

  print "Writing moved pdb to file: %s" % params.super.moved
  out = open(params.super.moved, "w")
  for serial, label, atom in zip(moving_pdb.atom_serial_number_strings(),
                                 moving_pdb.input_atom_labels_list(),
                                 moving_pdb.atoms()):
    print >> out, iotbx.pdb.format_atom_record(
      record_name={False: "ATOM", True: "HETATM"}[atom.hetero],
      serial=int(serial),
      name=label.name(),
      altLoc=label.altloc(),
      resName=label.resname(),
      resSeq=label.resseq,
      chainID=label.chain(),
      iCode=label.icode(),
      site=lsq_fit.r * matrix.col(atom.xyz) + lsq_fit.t,
      occupancy=atom.occ,
      tempFactor=atom.b,
      segID=atom.segid,
      element=atom.element,
      charge=atom.charge)
  print

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
