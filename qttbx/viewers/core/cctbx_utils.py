from pathlib import Path

import pandas as pd
from cctbx.array_family import flex

from .pandas_utils import df_to_cif_lines



def convert_iseq_to_i(df):
  """
  Takes a dataframe and looks for i,j,k,l _seq columns
  If found, then i,j,k,l columns will be added
  """


  df.rename(columns={"i_seq":"i",
                     "j_seq":"j",
                     "k_seq":"k",
                     "l_seq":"l",
                    },inplace=True)

  return df

def cell_convert_cctbx(cell):
  # convert a cell with cctbx stuff to python
  if isinstance(cell,(flex.size_t,flex.double)):
    return list(cell)
  else:
    return cell

def obj_to_dict(obj):
    # Returns an objects non-function and non-hidden attributes as a dictionary
    d = {}
    for attr in dir(obj):
        if not attr.startswith('__'):  # Ignore built-in attributes
            value = getattr(obj, attr)
            if not callable(value):# Ignore functions
                d[attr]=value
    return d

def get_restraint_df(model,geometry_type="bond"):
    """
    Access the cctbx model's geometry restraints as pandas dataframes
    """
    assert geometry_type in ["bond","nonbonded","angle","dihedral","chirality","planarity","parallelity"]
    grm = model.get_restraints_manager()
    assert grm is not None, "Process model to build restraints manager first."
    grm = grm.geometry

    flags=None
    sites_cart=model.get_sites_cart()

    if geometry_type == "bond":
        pair_proxies = grm.pair_proxies(flags=flags, sites_cart=sites_cart)
        keys = ["i_seq", "j_seq","labels", "distance_ideal", "distance_model", "slack", "delta", "sigma", "weight", "residual", "sym_op_j", "rt_mx"]
        simple,asu = pair_proxies.bond_proxies.get_sorted("delta",model.get_sites_cart())
        records = [dict(zip(keys,info)) for info in simple]


    elif geometry_type=="nonbonded":
        pair_proxies = grm.pair_proxies(flags=flags, sites_cart=sites_cart)
        keys = ["labels", "i_seq", "j_seq", "delta", "vdw_distance", "sym_op_j", "rt_mx"]
        simple,asu = pair_proxies.nonbonded_proxies.get_sorted("delta",model.get_sites_cart())
        records = [dict(zip(keys,info)) for info in simple if info[-1]==None]
        # check for not simple, for some reason simple,asu not working like others

    elif geometry_type == "angle":
        keys = ["i_seqs", "angle_ideal", "angle_model", "delta", "sigma", "weight", "residual"]
        simple,asu =  grm.angle_proxies.get_sorted("delta",model.get_sites_cart())
        records = [dict(zip(keys,info)) for info in simple]

    elif geometry_type == "dihedral":
        keys = ["i_seqs", "angle_ideal", "angle_model", "delta", "period", "sigma", "weight", "residual"]
        simple,asu =  grm.dihedral_proxies.get_sorted("delta",model.get_sites_cart())
        records = [dict(zip(keys,info)) for info in simple]

    elif geometry_type == "chirality":
        keys = ["i_seqs","both_signs","ideal","model","delta","sigma","weight","residual"]
        simple,asu =  grm.chirality_proxies.get_sorted("delta",model.get_sites_cart())
        records = [dict(zip(keys,info)) for info in simple]

    # TODO: This code does not put actual model values in for planarity and parallelity
    elif geometry_type == "planarity":
        proxies = grm.planarity_proxies
        records = []
        for proxy in proxies:
            record = obj_to_dict(proxy)
            records.append(record)
    elif geometry_type == "parallelity":
        proxies = grm.parallelity_proxies
        records = []
        for proxy in proxies:
            record = obj_to_dict(proxy)
            records.append(record)


    df = pd.DataFrame.from_records(records)
    df = df.applymap(cell_convert_cctbx)

    # move i_seqs columns to list
    if any([("_seq" in column and column != "i_seqs") for column in df.columns]):

      indices = "ijklm"
      to_merge = []
      for i in indices:
        name = f"{i}_seq"
        if name in df.columns:
          to_merge.append(name)
      if len(to_merge)>0:
        df['i_seqs'] = df.apply(lambda row: [int(row[name]) for name in to_merge], axis=1)
        df.drop(to_merge, axis=1, inplace=True)

    # add labels if none
    if 'labels' not in df.columns:
      if 'i_seqs' in df.columns:
        df['labels'] = df['i_seqs'].apply(lambda lst: [str(x) for x in lst])




    # # move i_seqs list to column if not too many
    # if "i_seqs" in df.columns:
    #   max_len = df["i_seqs"].apply(len).max()
    #   if max_len<5:
    #     indices = "ijkl"
    #     df[[indices[e] for e in range(max_len)]] = pd.DataFrame(df['i_seqs'].tolist(), index=df.index)
    #     del df["i_seqs"]

    #df = convert_iseq_to_i(df)
    return df

def get_restraint_dfs_from_model(model,lazy_build=False):
  # convert geometry restraints manager to dataframes
  # write dataframes to cif
  if not lazy_build:
    assert model.get_restraints_manager() is not None, "Must process and build restraints first"
  else:
    if model.get_restraints_manager() is None:
      model.process(make_restraints=True)

  # these are the key,values for supported restraints
  restraint_dfs = {
      geometry_type:get_restraint_df(model,geometry_type=geometry_type)
    for geometry_type in ["bond","nonbonded","angle","dihedral","chirality","planarity","parallelity"]
    }

  return restraint_dfs


def write_phenix_restraints_cif(model,filename,data_name="data_restraints"):
  outfile = Path(filename)

  restraint_dfs = get_restraint_dfs_from_model(model)
  cif_lines_all = [data_name]
  for key,df in restraint_dfs.items():

    cif_lines = df_to_cif_lines(df,column_prefix=key,integer_padding=8)
    cif_lines_all+=cif_lines

    with outfile.open("w") as fh:
      fh.write("\n".join(cif_lines_all))


def model_to_df(model):
  atoms = model.get_atoms()
  attrs_composition = [
                        #"model_id", # model
                        "asym_id", # chain
                        "seq_id", # resid
                        "comp_id", # resname
                        "atom_id", # name
                        "type_symbol", # element
                        "alt_id"] # alt_loc




  data = {
                        "id":[atom.i_seq+1 for atom in atoms],
                        "asym_id":[atom.parent().parent().parent().id for atom in atoms],
                        "seq_id":[atom.parent().parent().resseq_as_int() for atom in atoms],
                        "comp_id":[list(atom.parent().parent().unique_resnames()) for atom in atoms],
                        "atom_id":[atom.name.strip() for atom in atoms],
                        "type_symbol":[atom.element for atom in atoms],
                        "alt_id":["." if atom.parent().altloc == "" else atom.parent().altloc for atom in atoms],
                        }
  # values non-composition
  xyz = model.get_sites_cart().as_numpy_array()
  data["x"] = xyz[:,0]
  data["y"] = xyz[:,1]
  data["z"] = xyz[:,2]


  assert len(set([len(e) for e in data["comp_id"]]))==1, "Residue groups exist with different resnames"
  data["comp_id"] = [e[0] for e in data["comp_id"]]
  df_atoms = pd.DataFrame(data,index=list(range(len(atoms))))
  df_atoms = df_atoms.astype({"id":"int",
                              "asym_id":"str",
                              "comp_id":"str",
                              "seq_id":"int",
                              "atom_id":"str",
                              "type_symbol":"str",
                              "alt_id":"str"})
  return df_atoms
