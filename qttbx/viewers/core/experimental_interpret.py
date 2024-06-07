from collections import defaultdict

import pandas as pd
import numpy as np
import gemmi

monlib_path = "/Applications/ccp4-8.0/lib/data/monomers/"

def extract_restraint_dataframes_with_gemmi(model):
  st = gemmi.read_pdb_string(model.model_as_pdb())


  resnames = st[0].get_all_residue_names()
  monlib = gemmi.read_monomer_lib(monlib_path, resnames)
  top = gemmi.prepare_topology(st,monlib)

  bond_data = {
  "atom_id_1":[bond.atoms[0].name for bond in top.bonds],
  "atom_id_2":[bond.atoms[1].name for bond in top.bonds],
  "id_1":[bond.atoms[0].serial for bond in top.bonds],
  "id_2":[bond.atoms[1].serial for bond in top.bonds],
  "value_dist":[bond.restr.value for bond in top.bonds],
  "value_dist_esd":[bond.restr.esd for bond in top.bonds],
  }
  bonds = pd.DataFrame(bond_data)
  bonds["i_seq_1"] = bonds["id_1"]-1
  bonds["i_seq_2"] = bonds["id_2"]-1
  bonds["weight"] = 1/bonds["value_dist_esd"]**2

  angle_data = {
  "atom_id_1":[angle.atoms[0].name for angle in top.angles],
  "atom_id_2":[angle.atoms[1].name for angle in top.angles],
  "atom_id_3" :[angle.atoms[2].name for angle in top.angles],
  "id_1":[angle.atoms[0].serial for angle in top.angles],
  "id_2":[angle.atoms[1].serial for angle in top.angles],
  "id_3" :[angle.atoms[2].serial for angle in top.angles],
  "value_angle":[angle.restr.value for angle in top.angles],
  "value_angle_esd":[angle.restr.esd for angle in top.angles],
  }
  angles = pd.DataFrame(angle_data)
  angles["i_seq_1"] = angles["id_1"]-1
  angles["i_seq_2"] = angles["id_2"]-1
  angles["i_seq_3"] = angles["id_3"]-1
  angles["weight"] = 1/angles["value_angle_esd"]**2

  torsion_data = {
  "atom_id_1":[torsion.atoms[0].name for torsion in top.torsions],
  "atom_id_2":[torsion.atoms[1].name for torsion in top.torsions],
  "atom_id_3" :[torsion.atoms[2].name for torsion in top.torsions],
  "atom_id_4" :[torsion.atoms[3].name for torsion in top.torsions],
  "id_1":[torsion.atoms[0].serial for torsion in top.torsions],
  "id_2":[torsion.atoms[1].serial for torsion in top.torsions],
  "id_3" :[torsion.atoms[2].serial for torsion in top.torsions],
  "id_4" :[torsion.atoms[2].serial for torsion in top.torsions],
  "value_angle":[torsion.restr.value for torsion in top.torsions],
  "value_angle_esd":[torsion.restr.esd for torsion in top.torsions],
  "value_angle_period":[torsion.restr.period for torsion in top.torsions],

  }
  torsions = pd.DataFrame(torsion_data)
  torsions["i_seq_1"] = torsions["id_1"]-1
  torsions["i_seq_2"] = torsions["id_2"]-1
  torsions["i_seq_3"] = torsions["id_3"]-1
  torsions["i_seq_4"] = torsions["id_4"]-1
  torsions["weight"] = 1/torsions["value_angle_esd"]**2
  torsions["periodicity"] = torsions["value_angle_period"]
  torsions["alt_angle_ideals"] = [None]*len(torsions)
  torsions["origin_id"] = [None]*len(torsions)


  plane_data = {
  "atom_id_1":[plane.atoms[0].name for plane in top.planes],
  "atom_id_2":[plane.atoms[1].name for plane in top.planes],
  "atom_id_3" :[plane.atoms[2].name for plane in top.planes],
  "atom_id_4" :[plane.atoms[3].name for plane in top.planes],
  "id_1":[plane.atoms[0].serial for plane in top.planes],
  "id_2":[plane.atoms[1].serial for plane in top.planes],
  "id_3" :[plane.atoms[2].serial for plane in top.planes],
  "id_4" :[plane.atoms[2].serial for plane in top.planes],
  "plane_id":[plane.restr.esd for plane in top.planes],
  "value_dist_esd":[plane.restr.esd for plane in top.planes],

  }
  planes = pd.DataFrame(plane_data)
  planes["i_seq_1"] = planes["id_1"]-1
  planes["i_seq_2"] = planes["id_2"]-1
  planes["i_seq_3"] = planes["id_3"]-1
  planes["i_seq_4"] = planes["id_4"]-1
  planes["weight_1"] = 1/planes["value_dist_esd"]**2
  planes["weight_2"] = 1/planes["value_dist_esd"]**2
  planes["weight_3"] = 1/planes["value_dist_esd"]**2
  planes["weight_4"] = 1/planes["value_dist_esd"]**2

  return bonds, angles, torsions, planes

# check
def getallattrs(obj,return_dict=False):
  attrs = []
  for attr in dir(obj):
    if not attr.startswith("_") and not callable(getattr(obj,attr)):
      attrs.append(attr)
  if not return_dict:
    return attrs
  d = {attr:getattr(obj,attr) for attr in attrs}
  return d

def check_proxy(proxy_left,proxy_right,assertion=True,proxy_type="any",ignore_origin_id=True,weight_tol=1):
  failed = False
  attrs_left = getallattrs(proxy_left)
  attrs_right = getallattrs(proxy_right)
  assert attrs_left==attrs_right, f"Attrs mismatch between proxies, left={attrs_left} and right={attrs_right}"
  for attr in attrs_left:
    value_left = getattr(proxy_left,attr)
    value_right = getattr(proxy_right,attr)
    if proxy_type=="plane" and attr in ["weights","i_seqs"]:
      value_left = list(value_left)
      value_right = list(value_right)
    if attr == "origin_id" and ignore_origin_id:
      continue
    if value_left!=value_right:
      if assertion:
        match = value_left==value_right
        if not match and proxy_type == "plane":
          diffs = np.array([abs(int(a)-int(b)) for a,b in zip(value_left,value_right)])
          value_left = diffs
          value_right = diffs
          match = diffs.max()<weight_tol
        if not match:
          assert False, f"Value mismatch for proxy type {proxy_type} for attr: {attr}, left={value_left} and right={value_right}"
      else:
        failed = True
        print( f"Value mismatch for attr: {attr}, left={value_left} and right={value_right}")
  return failed

def is_similar_grm(grm_left,grm_right,model):


  grm_left.pair_proxies(model.get_sites_cart())
  grm_right.pair_proxies(model.get_sites_cart())

  # Bonds simple
  for proxy_left,proxy_right in zip(grm_left.get_all_bond_proxies()[0],
                                    grm_right.get_all_bond_proxies()[0]):
    check_proxy(proxy_left,proxy_right)

  # Angles
  for proxy_left,proxy_right in zip(grm_left.get_all_angle_proxies(),
                                    grm_right.get_all_angle_proxies()):
    check_proxy(proxy_left,proxy_right)

  # Torsions
  for i,(proxy_left,proxy_right) in enumerate(zip(grm_left.dihedral_proxies,
                                    grm_right.dihedral_proxies)):

    check_proxy(proxy_left,proxy_right,assertion=True)

  # Planes
  for i,(proxy_left,proxy_right) in enumerate(zip(grm_left.planarity_proxies,
                                    grm_right.planarity_proxies)):
    check_proxy(proxy_left,proxy_right,assertion=False,proxy_type="plane")

  # nonbonded
  for i,(proxy_left,proxy_right) in enumerate(zip(grm_left.pair_proxies().nonbonded_proxies.simple,
                                    grm_right.pair_proxies().nonbonded_proxies.simple)):
    check_proxy(proxy_left,proxy_right,assertion=True)



def create_atom_df_from_atoms(atoms):


  # non vectorized-by-atom attrs
  func_mapper = {

                  "chain_id":lambda atom: atom.parent().parent().parent().id.strip(),
                  "resseq":lambda atom: atom.parent().parent().resseq_as_int(),
                  "resname":lambda atom: atom.parent().resname.strip(),
                  "element":lambda atom: atom.element,
                  "altloc": lambda atom: atom.parent().altloc,
                  "model_id":lambda atom: atom.parent().parent().parent().parent().id,
                  "icode":lambda atom: atom.parent().parent().icode.strip(),
                  "charge":lambda atom: atom.charge_as_int(),
  }


  data = defaultdict(list)
  for atom in atoms:
    for key,func in func_mapper.items():
      data[key].append(func(atom))



  # vectorized-by-atom attrs non-numeric
  data["group_PDB"] = np.full(len(atoms),"ATOM",dtype='<U6')
  hetatm_mask = np.array(atoms.extract_hetero())
  if len(hetatm_mask)>0:
    data["group_PDB"][hetatm_mask] = "HETATM"
  data["id"] = np.array(atoms.extract_serial())
  data["name"] = list(atoms.extract_name())
  data["segid"] = list(atoms.extract_segid())



  # Make atom sites dataframe
  sites = pd.DataFrame(data,index=list(range(len(atoms)))).astype("string")


  # Strip strings
  for column in sites.columns:
    if sites[column].dtype == "string":
      sites[column] = sites[column].str.strip()

  # vectorized floats
  xyz = atoms.extract_xyz().as_numpy_array()
  sites["x"] = xyz[:,0]
  sites["y"] = xyz[:,1]
  sites["z"] = xyz[:,2]
  sites["b"] = atoms.extract_b().as_numpy_array()
  sites["occ"] = atoms.extract_occ().as_numpy_array()


  # Replacement of string blanks
  replacements = {'': pd.NA,
                  "?":pd.NA,
                  ".":pd.NA,
                  }

  sites.replace(replacements,inplace=True)

  return sites

def unpack_func_mapper(func_mapper,iterable,silent=True,blank=pd.NA):
  data = defaultdict(list)
  for obj in iterable:
    for key,func in func_mapper.items():
      try:
        data[key].append(func(obj))
      except:
        if not silent:
          raise
        else:
          data[key].append(blank)
  return data


def extract_processed_model_to_dataframes(model,do_nonbonded=False):

  sites = create_atom_df_from_atoms(model.get_atoms())

  # Bonds
  grm = model.restraints_manager.geometry
  grm.pair_proxies(model.get_sites_cart())
  simple_bonds, asu = grm.get_all_bond_proxies()
  func_mapper = {
    "i_seq_1":lambda proxy: proxy.i_seqs[0],
    "i_seq_2": lambda proxy: proxy.i_seqs[1],
    "value_dist":lambda proxy: proxy.distance_ideal,
    "weight":lambda proxy: proxy.weight,
  }
  bonds = pd.DataFrame(unpack_func_mapper(func_mapper,simple_bonds))
  bonds["atom_id_1"] = sites["name"].values[bonds["i_seq_1"].values]
  bonds["atom_id_2"] = sites["name"].values[bonds["i_seq_2"].values]


  # Angles
  func_mapper = {
    "i_seq_1":lambda proxy: proxy.i_seqs[0],
    "i_seq_2": lambda proxy: proxy.i_seqs[1],
    "i_seq_3": lambda proxy: proxy.i_seqs[2],
    "value_angle":lambda proxy: proxy.angle_ideal,
    "weight":lambda proxy: proxy.weight,
  }
  angles = pd.DataFrame(unpack_func_mapper(func_mapper,grm.get_all_angle_proxies()))
  angles["atom_id_1"] = sites["name"].values[angles["i_seq_1"].values]
  angles["atom_id_2"] = sites["name"].values[angles["i_seq_2"].values]
  angles["atom_id_3"] = sites["name"].values[angles["i_seq_3"].values]

  # Torsions
  func_mapper = {
    "i_seq_1":lambda proxy: proxy.i_seqs[0],
    "i_seq_2": lambda proxy: proxy.i_seqs[1],
    "i_seq_3": lambda proxy: proxy.i_seqs[2],
    "i_seq_4": lambda proxy: proxy.i_seqs[3],
    "value_angle":lambda proxy: proxy.angle_ideal,
    "weight":lambda proxy: proxy.weight,
    "periodicity":lambda proxy: proxy.periodicity,
    "alt_angle_ideals": lambda proxy: proxy.alt_angle_ideals,
    "origin_id":lambda proxy: proxy.origin_id
  }
  torsions = pd.DataFrame(unpack_func_mapper(func_mapper,grm.dihedral_proxies))
  torsions["atom_id_1"] = sites["name"].values[torsions["i_seq_1"].values]
  torsions["atom_id_2"] = sites["name"].values[torsions["i_seq_2"].values]
  torsions["atom_id_3"] = sites["name"].values[torsions["i_seq_3"].values]
  torsions["atom_id_4"] = sites["name"].values[torsions["i_seq_4"].values]


  # Planes
  max_plane_size = max([len(proxy.i_seqs) for proxy in grm.planarity_proxies])
  func_mapper = {f"i_seq_{i}": lambda proxy, i=i: proxy.i_seqs[i] for i in range(max_plane_size)}
  func_mapper.update({f"weight_{i}": lambda proxy, i=i: proxy.weights[i] for i in range(max_plane_size)})

  planes = pd.DataFrame(unpack_func_mapper(func_mapper,grm.planarity_proxies,silent=True))
  for i in range(max_plane_size):
    seq_column = f"i_seq_{i}"
    atom_id_column = f"atom_id_{i}"

    seq_values = planes[seq_column].values
    non_na_mask = ~pd.isna(seq_values)

    if non_na_mask.any():
      valid_indices = seq_values[non_na_mask].astype(int)

      # Use .iloc for positional indexing on non-NA values
      planes.loc[non_na_mask, atom_id_column] = sites["name"].iloc[valid_indices].values
    else:
      # Handle the case where all are NA, if necessary
      planes[atom_id_column] = pd.NA

  # Nonbonded
  nonbonds = None
  if do_nonbonded:
    nb_proxies = grm.pair_proxies().nonbonded_proxies.simple
    func_mapper = {
      "i_seq_1":lambda proxy: proxy.i_seqs[0],
      "i_seq_2": lambda proxy: proxy.i_seqs[1],
      "vdw_distance":lambda proxy: proxy.vdw_distance,
    }
    nonbonds = pd.DataFrame(unpack_func_mapper(func_mapper,nb_proxies))
    nonbonds["atom_id_1"] = sites["name"].values[nonbonds["i_seq_1"].values]
    nonbonds["atom_id_2"] = sites["name"].values[nonbonds["i_seq_2"].values]


  return bonds,angles,torsions,planes, nonbonds

class DummySourceInfo:
  def __init__(self,n_expected_atoms=0):
    self._n_expected_atoms = n_expected_atoms
  def labels(self):
    return "dummy_source_info_labels"
  def n_expected_atoms(self):
    return self._n_expected_atoms