from __future__ import division

# chain_comparison.py
# a tool to compare main-chain from two structures with or without crystal
# symmetry
#


from libtbx import adopt_init_args
import iotbx.phil
import os,sys
from phenix.autosol.build_model import best_match
from libtbx.utils import Sorry

master_phil = iotbx.phil.parse("""

  input_files {
    pdb_in = None
      .type = path
      .multiple = True
      .help = Input PDB file
      .short_caption = Input PDB file

  }
  crystal_info {
    chain_type = *PROTEIN RNA DNA
      .type = choice
      .short_caption = Chain type
      .help = Chain type.  All residues of other chain types ignored.

    use_crystal_symmetry = True
      .type = bool
      .short_caption = Use crystal symmetry in comparison
      .help = If set, use crystal symmetry to map atoms to closest positions
  }
  comparison {
    max_dist = 3.
      .type = float
      .short_caption = Maximum close distance
      .caption = Maximum distance between atoms to be considered close
  }
  control {
      verbose = False
        .type = bool
        .help = Verbose output
        .short_caption = Verbose output

      quiet = False
        .type = bool
        .help = No printed output
        .short_caption = No printed output

  }
""", process_includes=True)
master_params = master_phil

def get_params(args,out=sys.stdout):
    command_line = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil=master_phil,
      pdb_file_def="input_files.pdb_in")

    params = command_line.work.extract()
    print >>out,"\nFind similarity between two main-chains"
    master_phil.format(python_object=params).show(out=out)
    return params


def apply_atom_selection(atom_selection,hierarchy=None):
  asc=hierarchy.atom_selection_cache()
  sel = asc.selection(string = atom_selection)
  return hierarchy.deep_copy().select(sel)  # deep copy is required

def select_atom_lines(lines):
  new_lines=[]
  for line in lines:
    if line.startswith("ATOM "):
      new_lines.append(line)
  return new_lines

class dummy_info:
  def __init__(self,i=None,j=None, distance=None):
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())

  def dist(self):
    return self.distance


def get_best_match(xyz1,xyz2,crystal_symmetry=None):
  if crystal_symmetry:
    return best_match(
      xyz1,xyz2,
      crystal_symmetry=crystal_symmetry)
  else: # do it without symmetry
    (distance,i,j)=xyz1.min_distance_between_any_pair_with_id(xyz2)
    return dummy_info(i=i,j=j,distance=distance)


def run(args=None,
   chain_hierarchy=None,
   target_hierarchy=None,
   target_file=None, # model
   chain_file=None, # query
   crystal_symmetry=None,
   max_dist=None,
   quiet=None,
   verbose=None,
   use_crystal_symmetry=None,
   chain_type=None,
   out=sys.stdout):
  if not args: args=[]
  params=get_params(args,out=out)


  if verbose is None:
    verbose=params.control.verbose
  if quiet is None:
    quiet=params.control.quiet
  if chain_type is None:
    chain_type=params.crystal_info.chain_type
  if use_crystal_symmetry is None:
    use_crystal_symmetry=params.crystal_info.use_crystal_symmetry
  if max_dist is None:
    max_dist=params.comparison.max_dist

  if not target_file and len(params.input_files.pdb_in)>0:
     target_file=params.input_files.pdb_in[0]  # model
  if not chain_file and len(params.input_files.pdb_in)>1:
     chain_file=params.input_files.pdb_in[1] # query

  # get the hierarchies
  if not chain_hierarchy or not target_hierarchy:
    assert chain_file and target_file
    from phenix.autosol.get_pdb_inp import get_pdb_inp
    pdb_inp=get_pdb_inp(file_name=chain_file  )
    if not crystal_symmetry:
      crystal_symmetry=pdb_inp.crystal_symmetry_from_cryst1()
    chain_hierarchy=pdb_inp.construct_hierarchy()
    target_pdb_inp=get_pdb_inp(file_name=target_file)
    if not crystal_symmetry:
      crystal_symmetry=target_pdb_inp.crystal_symmetry_from_cryst1()
    target_hierarchy=target_pdb_inp.construct_hierarchy()
  if not quiet:
    print >>out,"Looking for chain similarity for "+\
      "%s (%d residues) in the model %s (%d residues)" %(
     chain_file,chain_hierarchy.overall_counts().n_residues,
     target_file,target_hierarchy.overall_counts().n_residues)
    if verbose:
      print >>out,"Chain type is: %s" %(chain_type)
  # get the CA residues
  if chain_type in ["RNA","DNA"]:
    atom_selection="name P"
  else:
    atom_selection="name ca and (not element Ca)"
  chain_ca=apply_atom_selection(atom_selection,chain_hierarchy)
  chain_ca_lines=select_atom_lines(chain_ca.as_pdb_string().splitlines())
  target_ca=apply_atom_selection(atom_selection,target_hierarchy)
  target_xyz_lines=select_atom_lines(target_ca.as_pdb_string().splitlines())
  chain_xyz_cart=chain_ca.atoms().extract_xyz()
  target_xyz_cart=target_ca.atoms().extract_xyz()

  # for each xyz in chain, figure out closest atom in target and dist
  best_i=None
  best_i_dd=None
  best_pair=None
  pair_list=[]
  from scitbx.array_family import flex
  chain_xyz_fract=crystal_symmetry.unit_cell().fractionalize(chain_xyz_cart)
  target_xyz_fract=crystal_symmetry.unit_cell().fractionalize(target_xyz_cart)
  far_away_match_list=[]
  far_away_match_rmsd_list=flex.double()
  if use_crystal_symmetry:
    working_crystal_symmetry=crystal_symmetry
  else:
    working_crystal_symmetry=None
  for i in xrange(chain_xyz_fract.size()):
    best_j=None
    best_dd=None
    if working_crystal_symmetry:
      info=get_best_match(
        flex.vec3_double([chain_xyz_fract[i]]),target_xyz_fract,
        crystal_symmetry=working_crystal_symmetry)
    else:
      info=get_best_match(
        flex.vec3_double([chain_xyz_cart[i]]),target_xyz_cart)

    if info and (best_dd is None or info.dist()<best_dd):
        best_dd=info.dist()
        best_j=info.j
    if best_dd > max_dist:
      far_away_match_list.append(i)
      far_away_match_rmsd_list.append(best_dd)
      if (not quiet) and verbose:
        print >>out,"%s" %(chain_ca_lines[i])
      continue
    if best_i is None or best_dd<best_i_dd:
      best_i=i
      best_i_dd=best_dd
      best_pair=[i,best_j]
    pair_list.append([i,best_j,best_dd])
  n_forward=0
  n_reverse=0
  forward_match_list=[]
  reverse_match_list=[]
  forward_match_rmsd_list=flex.double()
  reverse_match_rmsd_list=flex.double()
  unaligned_match_list=[]
  unaligned_match_rmsd_list=flex.double()
  last_i=None
  last_j=None
  for [i,j,dd],[next_i,next_j,next_dd] in zip(
      pair_list,pair_list[1:]+[[None,None,None]]):
    if i is None or j is None: continue
    found=False
    if last_i is None: # first time
      if next_i==i+1: # starting a segment
        if next_j==j+1:
          n_forward+=1
          forward_match_list.append([i,j])
          forward_match_rmsd_list.append(dd**2)
          found=True
        elif next_j==j-1:
          n_reverse+=1
          reverse_match_list.append([i,j])
          reverse_match_rmsd_list.append(dd**2)
          found=True
    else: # not the first time
      if i==last_i+1: # continuing a segment
        if j==last_j+1:
          n_forward+=1
          forward_match_list.append([i,j])
          forward_match_rmsd_list.append(dd**2)
          found=True
        elif j==last_j-1:
          n_reverse+=1
          reverse_match_list.append([i,j])
          reverse_match_rmsd_list.append(dd**2)
          found=True
    if not found:
      last_i=None
      last_j=None
      unaligned_match_list.append([i,j])
      unaligned_match_rmsd_list.append(dd**2)
    else:
      last_i=i
      last_j=j

  if n_forward==n_reverse==0:
    direction='none'
  elif n_forward>= n_reverse:
    direction='forward'
  else:
    direction='reverse'
  if (not quiet) and verbose:
    print >>out,"%s %d  %d  N: %d" %(
     direction,n_forward,n_reverse,chain_xyz_fract.size())

  if not quiet:
    if verbose:
      print >>out,"Total CA: %d  Too far to match: %d " %(
        chain_xyz_fract.size(),len(far_away_match_list))

    if forward_match_rmsd_list.size():
      print >>out,\
          "\nResidues matching in forward direction:   %4d  RMSD: %6.2f" %(
         forward_match_rmsd_list.size(),
         forward_match_rmsd_list.min_max_mean().mean**0.5)
      if verbose:
        for i,j in forward_match_list:
          print >>out,"ID:%d:%d  RESIDUES:  \n%s\n%s" %( i,j, chain_ca_lines[i],
           target_xyz_lines[j])

    if reverse_match_rmsd_list.size():
      print >>out,\
         "Residues matching in reverse direction:   %4d  RMSD: %6.2f" %(
         reverse_match_rmsd_list.size(),
         reverse_match_rmsd_list.min_max_mean().mean**0.5)
      if verbose:
        for i,j in reverse_match_list:
          print >>out,"ID:%d:%d  RESIDUES:  \n%s\n%s" %(
           i,j, chain_ca_lines[i],
           target_xyz_lines[j])

    if unaligned_match_rmsd_list.size():
      print >>out,\
        "Residues near but not matching one-to-one:%4d  RMSD: %6.2f" %(
         unaligned_match_rmsd_list.size(),
         unaligned_match_rmsd_list.min_max_mean().mean**0.5)
      if verbose:
        for i,j in unaligned_match_list:
          print >>out,"ID:%d:%d  RESIDUES:  \n%s\n%s" %(i,j, chain_ca_lines[i],
            target_xyz_lines[j])

    if far_away_match_rmsd_list.size():
      print >>out,\
        "Residues far from target:                 %4d  RMSD: %6.2f" %(
         far_away_match_rmsd_list.size(),
         far_away_match_rmsd_list.min_max_mean().mean**0.5)
      if verbose:
        for i in far_away_match_list:
          print >>out,"ID:%d  RESIDUES:  \n%s" %(i,chain_ca_lines[i])


  return n_forward,n_reverse,len(pair_list)

if __name__=="__main__":
  args=sys.argv[1:]
  run(args=args,out=sys.stdout)
