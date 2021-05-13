##################################################################################
# This is a test program to validate that the Python wrapping of Probe worked.
#

import sys
import mmtbx_probe_ext as probe
from iotbx.map_model_manager import map_model_manager
from iotbx.data_manager import DataManager
import mmtbx

def RunProbeTests(inFileName):

  #========================================================================
  # Make sure we can get at the DotSphere objects and their methods
  print('Making DotSphere from cache and getting its dots')
  cache = probe.DotSphereCache(10)
  sphere1 = cache.get_sphere(1)
  dots = sphere1.dots()
  print(' Found',len(dots),'dots')
  print('First dot is at',dots[0][0],dots[0][1],dots[0][2])

  #========================================================================
  # Make sure we can fill in an ExtraAtomInfoList and pass it to scoring
  # Generate an example data model with a small molecule in it
  print('Generating model')
  if inFileName is not None and len(inFileName) > 0:
    # Read a model from a file using the DataManager
    dm = DataManager()
    dm.process_model_file(inFileName)
    model = dm.get_model(inFileName)
  else:
    # Generate a small-molecule model using the map model manager
    mmm=map_model_manager()         #   get an initialized instance of the map_model_manager
    mmm.generate_map()              #   get a model from a generated small library model and calculate a map for it
    model = mmm.model()             #   get the model

  # Fill in an ExtraAtomInfoList with an entry for each atom in the hierarchy.
  # We first find the largest i_seq sequence number in the model and reserve that
  # many entries so we will always be able to fill in the entry for an atom.
  print('Filling in extra atom information needed for probe score')
  atoms = model.get_atoms()
  maxI = atoms[0].i_seq
  for a in atoms:
    if a.i_seq > maxI:
      maxI = a.i_seq
  extra = []
  for i in range(maxI+1):
    extra.append(probe.ExtraAtomInfo())

  # Get the bonding information we'll need to exclude our bonded neighbors
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  model.set_pdb_interpretation_params(params = p)
  model.process_input_model(make_restraints=True) # make restraints
  geometry = model.get_restraints_manager().geometry
  sites_cart = model.get_sites_cart() # cartesian coordinates
  bond_proxies_simple, asu = \
      geometry.get_all_bond_proxies(sites_cart = sites_cart)

  # Make a mapping from sequence number to atom so that we can quickly look this
  # up when finding bonded neighbors
  # @todo There may be a faster way to find bonded neighbors that does not require
  # this.
  atomDict = {}
  for a in atoms:
    atomDict[a.i_seq] = a

  # Make a dictionary for each atom listing all of its bonded neighbors
  bondedNeighbors = {}
  for a in atoms:
    bondedNeighbors[a] = []
  for bp in bond_proxies_simple:
    bondedNeighbors[atomDict[bp.i_seqs[0]]].append(atomDict[bp.i_seqs[1]])
    bondedNeighbors[atomDict[bp.i_seqs[1]]].append(atomDict[bp.i_seqs[0]])

  # Traverse the hierarchy and look up the extra data to be filled in.
  # Get a list of all the atoms in the chain while we're at it
  atoms = []
  mon_lib_srv = model.get_mon_lib_srv()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  ph = model.get_hierarchy()
  for m in ph.models():
    for chain in m.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          md, ani = mon_lib_srv.get_comp_comp_id_and_atom_name_interpretation(
                residue_name=ag.resname, atom_names=ag.atoms().extract_name())
          atom_dict = md.atom_dict()

          for a in ag.atoms():
            atoms.append(a)
            te = atom_dict[a.name.strip()].type_energy
            extra[a.i_seq].vdwRadius = ener_lib.lib_atom[te].vdw_radius
            hb_type = ener_lib.lib_atom[te].hb_type
            if hb_type == "A":
              extra[a.i_seq].isAcceptor = True
            if hb_type == "D":
              extra[a.i_seq].isDonor = True

  # Construct a SpatialQuery and fill in the atoms.  Ensure that we can make a
  # query within 1000 Angstroms of the origin.
  sq = probe.SpatialQuery(atoms)
  nb = sq.neighbors((0,0,0), 0, 1000)
  print('Found this many atoms within 1000A of the origin:', len(nb))

  # Construct a DotScorer object.
  # Find the radius of each atom in the structure and construct dot spheres for
  # them. Find the atoms that are bonded to them and add them to an excluded list.
  # Then compute the score for each of them and report the summed score over the
  # whole molecule the way that Reduce will.
  ds = probe.DotScorer(extra)
  total = 0
  badBumpTotal = 0
  for a in atoms:
    rad = extra[a.i_seq].vdwRadius
    sphere = cache.get_sphere(rad)

    # Excluded atoms that are bonded to me or to one of my neightbors.
    # It has the side effect of excluding myself if I have any neighbors.
    # Construct as a set to avoid duplicates.
    exclude = set()
    for n in bondedNeighbors[a]:
      exclude.add(n)
      for n2 in bondedNeighbors[n]:
        exclude.add(n2)
    exclude = list(exclude)

    dots = sphere.dots()
    res = ds.score_dots(a, 1.0, sq, rad*3, 0.25, exclude, sphere.dots(), sphere.density(), False)
    total += res.totalScore()
    if res.hasBadBump:
      badBumpTotal += 1
  print('Summed probe score for molecule =',total,'with',badBumpTotal,'bad bumps')

  # @todo

  #========================================================================
  # Call the test functions for all of our files.

  print('Testing DotSphere objects')
  ret = probe.DotSpheres_test()
  assert (len(ret) == 0)

  print('Testing SpatialQuery objects')
  ret = probe.SpatialQuery_test()
  assert (len(ret) == 0)

  print('Testing Scoring objects')
  ret = probe.Scoring_test()
  assert (len(ret) == 0)

  return ret

if __name__ == '__main__':

  #==============================================================
  # Parse command-line arguments.  The 0th argument is the name
  # of the script. There can be the name of a PDB file to read.
  realParams = 0
  fileName = ""
  for i in range(1,len(sys.argv)):
    fileName = sys.argv[i]

  ret = RunProbeTests(fileName)
  if len(ret) == 0:
    print('Success!')

  assert (len(ret) == 0)
