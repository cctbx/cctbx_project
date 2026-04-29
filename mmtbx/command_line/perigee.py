"""Interaction finder"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.perigee
import os, sys
from math import sqrt

import iotbx
import iotbx.pdb

import libtbx
from libtbx import phil
import libtbx.phil.command_line
from libtbx.option_parser import OptionParser
#from libtbx.utils import Sorry
#from libtbx import runtime_utils

master_phil_string = """

perigee
  .caption = None
{
  input
  {
    pdb_file_name = None
      .type = path
      .short_caption = model
      .help = PDB filename
      .style = bold file_type:pdb
  }
  control
  {
    distance_cutoff = 3.
      .type = float
    chain = None
      .type = str
    other_chain = None
      .type = str
    residue = None
      .type = str
  }
  output
  {
    residues_only = False
      .type = bool
  }
}
"""
old = """
    file_name = None
      .type = path
      .short_caption = Output file
      .help = Defaults to current directory
      .style = bold new_file file_type:pdb
"""
master_params = master_phil_string # need for auto documentation
if 0:
  print('-'*80)
  print(master_phil_string)
  print('-'*80)
master_phil = phil.parse(master_phil_string)

def setup_parser():
  parser = OptionParser(
    prog="phenix.perigee",
    version="""
  up-to-date version
""",
    usage="""
  phenix.perigee pdb_file_name=pdb3a37.ent
""",
    )
  # Input options
  parser.add_option("",
                    "--show_defaults",
                    dest="show_defaults",
                    default=False,
                    action="store_true",
                    help="Display defaults",
                    )
  return parser

def get_chains(hierarchy, verbose=False):
  tmp = []
  for model in hierarchy.models():
    if verbose: print('model: "%s"' % model.id)
    for chain in model.chains():
      if verbose: print('chain: "%s"' % chain.id)
      if chain.id not in tmp: tmp.append(chain.id)
  return tmp

def get_het_residues(hierarchy):
  tmp = []
  for atom in hierarchy.atoms():
    if atom.hetero:
      if atom.parent().resname not in tmp:
        tmp.append(atom.parent().resname)
  return tmp

def get_distance2(atom1, atom2):
  d2 = (atom1.xyz[0]-atom2.xyz[0])**2
  d2 += (atom1.xyz[1]-atom2.xyz[1])**2
  d2 += (atom1.xyz[2]-atom2.xyz[2])**2
  return d2

def loop_over_residues_chain(hierarchy,
                             only_chain=None,
                             exclude_chain=None,
                             verbose=False,
                             ):
  #assert len(filter(None, [only_chain, exclude_chain]))==1
  for model in hierarchy.models():
    if verbose: print('model: "%s"' % model.id)
    for chain in model.chains():
      if(only_chain is not None and
         only_chain.strip()!=chain.id.strip()
         ): continue
      if(exclude_chain is not None and
         exclude_chain.strip()==chain.id.strip()
         ): continue
      if verbose: print('chain: "%s"' % chain.id)
      for resi, residue_group in enumerate(chain.residue_groups()):
        if verbose: print('  residue_group: resseq="%s" icode="%s"' % (
          residue_group.resseq, residue_group.icode))
        yield residue_group

def loop_over_residues_residue(hierarchy,
                               only_residue=None,
                               exclude_residue=None,
                               verbose=False,
                               ):
  assert len(list(filter(None, [only_residue, exclude_residue])))==1
  for model in hierarchy.models():
    if verbose: print('model: "%s"' % model.id)
    for chain in model.chains():
      if verbose: print('chain: "%s"' % chain.id)
      for resi, residue_group in enumerate(chain.residue_groups()):
        if verbose: print('  residue_group: resseq="%s" icode="%s"' % (
          residue_group.resseq, residue_group.icode))
        for atom_group in residue_group.atom_groups():
          if verbose: print('    atom_group: altloc="%s" resname="%s"' % (
            atom_group.altloc, atom_group.resname))
          if(only_residue is not None and
             only_residue.strip()!=atom_group.resname.strip()
             ): continue
          if(exclude_residue is not None and
             exclude_residue.strip()==atom_group.resname.strip()
             ): continue
          yield residue_group
          break

def get_interacting_atoms(hierarchy,
                          distance_cutoff=3.,
                          chain1=None,
                          chain2=None,
                          residue1=None,
                          residue2=None,
                          ):
  def _compare_residues(residue_1,
                        residue_2,
                        distance_cutoff=3.,
                        exclude_water=True,
                        ):
    distance_cutoff *= distance_cutoff
    atom_pairs = []
    atom1 = residue_1.atoms()[0]
    atom2 = residue_2.atoms()[0]
    d2 = get_distance2(atom1, atom2)
    if d2>distance_cutoff*100: return atom_pairs
    if exclude_water:
      if atom1.quote().find("HOH")!=-1: return []
      if atom2.quote().find("HOH")!=-1: return []
    for atom1 in residue_1.atoms():
      for atom2 in residue_2.atoms():
        d2 = get_distance2(atom1, atom2)
        if d2<distance_cutoff:
          atom_pairs.append([atom1, atom2])
    return atom_pairs
  ###################
  if chain1:
    assert residue1 is None
    assert residue2 is None
  elif residue1:
    assert chain1 is None
    assert chain2 is None

  atom_pairs = []
  if chain1:
    for i, residue_1 in enumerate(loop_over_residues_chain(
        hierarchy,
        only_chain=chain1,
        )
                                  ):
      for j, residue_2 in enumerate(loop_over_residues_chain(
          hierarchy,
          only_chain=chain2,
          exclude_chain=chain1,
          )
                                    ):
        tmp = _compare_residues(residue_1,
                                residue_2,
                                distance_cutoff=distance_cutoff,
                                )
        atom_pairs += tmp
  else:
    for i, residue_1 in enumerate(loop_over_residues_residue(
        hierarchy,
        only_residue=residue1,
        )
                                  ):
      for j, residue_2 in enumerate(loop_over_residues_residue(
          hierarchy,
          exclude_residue=residue1,
          )
                                    ):
        tmp = _compare_residues(residue_1,
                                residue_2,
                                distance_cutoff,
                                )
        atom_pairs += tmp
  return atom_pairs

def run_probe_on_pdb_string(pdb_filename):
  cmd = 'phenix.probe -u -condense -self -mc -NOVDWOUT -NOCLASHOUT ALL -'
  assert 0

def run(rargs):
  print("""

    Running interaction finder

  """)
  rargs = list(rargs)
  parser = setup_parser()
  (options, args) = parser.parse_args(args=rargs)
  if options.show_defaults:
    for line in master_phil_string.splitlines():
      if line.strip().find(".")==0: continue
      print(line)
    sys.exit()
  if len(args)==0:
    parser.print_help()
    sys.exit()
  #
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil,
    home_scope="perigee")
  #
  phils = []
  phil_args = []
  pdbs = []
  for arg in args:
    if os.path.isfile(arg):
      if iotbx.pdb.is_pdb_file(arg):
        pdbs.append(arg)
        continue
      try :
        file_phil = phil.parse(file_name=arg)
      except RuntimeError :
        pass
      else :
        phils.append(file_phil)
    else :
      phil_args.append(arg)
      phils.append(argument_interpreter.process(arg))
  working_phil = master_phil.fetch(sources=phils)
  #working_phil.show()
  working_params = working_phil.extract()
  in_scope = working_params.perigee.input
  control_scope = working_params.perigee.control
  output_scope = working_params.perigee.output

  for i, pdb in enumerate(pdbs):
    if i==0 and not in_scope.pdb_file_name:
      in_scope.pdb_file_name = pdbs[i]
  #
  #if not output_scope.file_name:
  #  output_scope.file_name = in_scope.pdb_file_name
  #  d = os.path.dirname(output_scope.file_name)
  #  output_scope.file_name = os.path.basename(output_scope.file_name)
  #  output_scope.file_name = output_scope.file_name.split(".")[0]
  #  if True:
  #    output_scope.file_name += "_loaded.pdb"
  #  else:
  #    output_scope.file_name += "_%s.pdb" % os.path.basename(
  #      in_scope.carbohydrate_file_name).split(".")[0]
  #  output_scope.file_name = os.path.join(d, output_scope.file_name)
  #
  #preamble = output_scope.file_name.split(".")[0]
  #print "    Writing effective parameters to %s.eff\n" % preamble
  working_phil.format(python_object=working_params).show()
  #f=open("%s.eff" % preamble, "wb")
  #f.write(working_phil.format(python_object=working_params).as_str())
  #f.close()

  pdb_inp = iotbx.pdb.input(in_scope.pdb_file_name,
                            #source_info="model PDB",
                            #lines=flex.split_lines(input_lines),
                            )
  hierarchy = pdb_inp.construct_hierarchy()
  #hierarchy.show()

  if control_scope.chain and control_scope.residue:
    assert 0
  elif control_scope.chain is None and control_scope.residue is None:
    if control_scope.chain is None:
      chains = get_chains(hierarchy)
      print("\n  Chains")
      for chain in chains:
        print('    "%s"' % chain)
    if control_scope.residue is None:
      print("Residues that may be of interest")
      residues = get_het_residues(hierarchy)
      for residue in residues:
        print('    "%s"' % residue)
    return

  def display_atom_pairs(atom_pairs, residues_only=False):
    def _comp_on_d2(pair1, pair2):
      d21 = get_distance2(pair1[0], pair1[1])
      d22 = get_distance2(pair2[0], pair2[1])
      if d21<d22: return -1
      return 1
    ##########
    atom_pairs.sort(_comp_on_d2)
    outl = ""
    done = []
    for atom1, atom2 in atom_pairs:
      key = None
      if residues_only:
        key = [atom1.quote()[17:27],
               atom2.quote()[17:27],
               ]
        if key in done: continue
      d2 = get_distance2(atom1, atom2)
      if residues_only:
        outl += "%s - %s : %6.1f\n" % (atom1.quote()[17:27],
                                       atom2.quote()[17:27],
                                       sqrt(d2),
                                       )
      else:
        outl += "%s - %s : %6.1f\n" % (atom1.quote(),
                                       atom2.quote(),
                                       sqrt(d2),
                                       )
      if residues_only: done.append(key)
    print(outl)

  if control_scope.chain:
    atom_pairs = get_interacting_atoms(
      hierarchy,
      distance_cutoff=control_scope.distance_cutoff,
      chain1=control_scope.chain,
      chain2=control_scope.other_chain,
      )
  elif control_scope.residue:
    atom_pairs = get_interacting_atoms(
      hierarchy,
      distance_cutoff=control_scope.distance_cutoff,
      residue1=control_scope.residue,
      )
  print('-'*80)
  if output_scope.residues_only:
    display_atom_pairs(atom_pairs, residues_only=True)
  else:
    display_atom_pairs(atom_pairs)
  print('-'*80)

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(args)

