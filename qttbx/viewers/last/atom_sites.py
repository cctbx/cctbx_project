import json
import pandas as pd
import numpy as np
from collections import OrderedDict
from collections import defaultdict
from iotbx.pdb import common_residue_names_get_class
from iotbx.pdb.utils import all_label_asym_ids
from iotbx.pdb import hierarchy
from .cif_io import CifInput
from libtbx import group_args
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal

params = group_args(

# Maps internal attr name -> phenix name
attrs_map_phenix = OrderedDict([

    ("group_id","group"), # ?
    ("model_id","model"),
    ("entity_id","entity_id"), #?
    ("asym_id","chain"),
    ("seq_id","resseq"),
    ("comp_id","resname"),
    ("atom_id","atom"),
    ("id","serial"),
    ("type_symbol","element"),
    ("x", "x"),
    ("y","y"),
    ("z","z"),
    ("occupancy",'occ'),
    ("bfactor","b"),
    ('icode','icode'),
    ("charge","charge"),
    ("alt_id","altloc")
    ]),

dtypes_mmcif = {
 'auth_asym_id': 'string',
 'label_asym_id': 'string',
 'auth_seq_id': 'Int64',
 'label_seq_id': 'Int64',
 'label_comp_id': 'string',
 'type_symbol': 'string',
 'label_alt_id': 'string',
 'pdbx_PDB_model_num': 'string',
 'pdbx_PDB_ins_code': 'string',
 'Cartn_x': 'Float64',
 'Cartn_y': 'Float64',
 'Cartn_z': 'Float64',
 'B_iso_or_equiv': 'Float64',
 'occupancy': 'Float64',
 'pdbx_formal_charge': 'Int64',
 'group_PDB': 'string',
 'id': 'string',
 'label_atom_id': 'string',
 'label_entity_id': 'string'
},


# Define functions to fill missing columns.
# 'None' indicates an absolutely required initial column
fill_functions = {
 'auth_asym_id': lambda sites: sites["label_asym_id"],
 'label_asym_id': None,
 'auth_seq_id': lambda sites: sites["label_seq_id"],
 'label_seq_id': None,
 'auth_comp_id': lambda sites: sites["label_comp_id"],
 'label_comp_id': None,
 'type_symbol': None,
 'label_alt_id': None,
 'pdbx_PDB_model_num': None,
 'pdbx_PDB_ins_code': None,
 'Cartn_x': None,
 'Cartn_y': None,
 'Cartn_z': None,
 'B_iso_or_equiv': None,
 'occupancy': None,
 'pdbx_formal_charge': None,
 'group_PDB': None,
 'id': None,
 'label_atom_id': None,
 'label_entity_id': None,
},

attr_fill_values = {
  # required
  'asym_id':"?",
  'seq_id':"?",
  'comp_id':"?",
  'atom_id':"?",
  'id':"?",
  'type_symbol':"?",

  # other
  'x':np.nan,
  'y':np.nan,
  'z':np.nan,
  'group_id':"ATOM",
  'entity_id':"",
  'model_id':"1",
  'occupancy':np.nan,
  'bfactor':np.nan,
  'icode':" ",
  'charge':"",
  'alt_id':"",
},
rounding = {
   "Cartn_x":3,
   "Cartn_y":3,
   "Cartn_Z":3,
   "B_iso_or_equiv":2,
   "occupancy":2,

},
blanks = set([""," ",".","?",np.nan,pd.NA]),

)

import pandas as pd

class AtomSites(pd.DataFrame):
  params = params


  @classmethod
  def from_cif_input(cls,cif_input,name=None):
     sites = cif_input.atom_sites
     return cls.from_mmcif_df(sites,name=name)

  @classmethod
  def from_mmcif_df(cls,df,name=None,build_hierarchy=True):
    replacements = {"?":pd.NA,
                    ".":pd.NA}
    df.replace(replacements,inplace=True)
    sites = cls(df,name=name)
    sites = AtomSites._init_new_sites(sites,build_hierarchy=build_hierarchy)
    return sites

  @classmethod
  def from_mmcif_file(cls,filename,name=None,method="cifpd",build_hierarchy=True):
    assert method in ["cifpd"] # TODO: add iotbx cif
    filename = str(filename)
    cif_input = CifInput(filename)
    return cls.from_mmcif_df(cif_input.atom_sites,
                                name=name,
                                build_hierarchy=build_hierarchy,
                                )

  @classmethod
  def from_cctbx_model(cls,model,name=None):
    atoms = model.get_atoms()
    return cls.from_cctbx_atoms(atoms,name=name)

  @classmethod
  def from_cctbx_hierarchy(cls,hierarchy,name=None):
    atoms = hierarchy.atoms()
    return cls.from_cctbx_atoms(atoms,name=name)


  @classmethod
  def from_cctbx_atoms(cls,atoms,name=None):
    df_atoms=AtomSites.create_atom_df_from_atoms(atoms)
    sites = cls(df_atoms,name=name)
    sites = AtomSites._init_new_sites(sites)
    return sites

  @staticmethod
  def create_atom_df_from_atoms(atoms):

    segid_as_auth_segid = False


    # prepare label_ values
    h = atoms[0].parent().parent().parent().parent().parent()
    label_asym_per_rg= {rg:h.get_label_asym_id(rg) for rg in h.residue_groups()}
    label_seq_per_ag = {ag:h.get_label_seq_id(ag) for ag in h.atom_groups()}



    # non vectorized-by-atom attrs
    func_mapper = {

                    "auth_asym_id":lambda atom: atom.parent().parent().parent().id.strip(),
                    "label_asym_id":lambda atom: label_asym_per_rg[atom.parent().parent()],
                    "auth_seq_id":lambda atom: atom.parent().parent().resseq_as_int(),
                    "label_seq_id":lambda atom: label_seq_per_ag[atom.parent()],
                    "label_comp_id":lambda atom: atom.parent().resname.strip(),
                    "type_symbol":lambda atom: atom.element,
                    "label_alt_id": lambda atom: atom.parent().altloc,
                    "pdbx_PDB_model_num":lambda atom: atom.parent().parent().parent().parent().id,
                    "pdbx_PDB_ins_code":lambda atom: atom.parent().parent().icode.strip(),
                    "pdbx_formal_charge":lambda atom: atom.charge_as_int(),
    }


    data = defaultdict(list)
    for atom in atoms:
      for key,func in func_mapper.items():
        data[key].append(func(atom))



    # vectorized-by-atom attrs
    N = len(atoms)
    xyz = atoms.extract_xyz().as_numpy_array()
    data["Cartn_x"] = xyz[:,0]
    data["Cartn_y"] = xyz[:,1]
    data["Cartn_z"] = xyz[:,2]
    data["B_iso_or_equiv"] = atoms.extract_b().as_numpy_array()
    data["occupancy"] = atoms.extract_occ().as_numpy_array()
    data["group_PDB"] = np.full(N,"ATOM",dtype='<U6')
    hetatm_mask = np.array(atoms.extract_hetero())
    if len(hetatm_mask)>0:
      data["group_PDB"][hetatm_mask] = "HETATM"
    data["id"] = np.array(atoms.extract_serial()).astype('string')
    data["id"] = data["id"].str.strip()
    data["label_atom_id"] = np.vectorize(str.strip)(np.array(atoms.extract_name()))
    if  segid_as_auth_segid:
      data["auth_asym_id"] = np.vectorize(str.strip)(np.array(atoms.extract_segid()))

    # Need attention
    data["label_entity_id"] = np.full(N,"?") # comes from the _entity. loop


    df_atoms = pd.DataFrame(data,index=list(range(len(atoms))))

    # strip strings
    #strip_if_string = lambda x: x.strip() if isinstance(x, str) else x
    df_atoms["type_symbol"] = df_atoms["type_symbol"].apply(str.strip)

    # Replace blanks
    replacements = {'': '.',pd.NA:"."}
    cols = ["label_alt_id"]
    for col in cols:
      df_atoms[col].replace(replacements,inplace=True)

    replacements = {'': '?',pd.NA:"?"}
    cols = ["pdbx_PDB_ins_code","pdbx_formal_charge"]
    for col in cols:
      df_atoms[col].replace(replacements,inplace=True)

    replacements = {'': pd.NA,
                    "?":pd.NA,
                    ".":pd.NA,
                   }
    cols = ["label_seq_id","auth_seq_id"]
    for col in cols:
      df_atoms[col].replace(replacements,inplace=True)

    # add model num if not present
    if (df_atoms["pdbx_PDB_model_num"]=="").all():
      df_atoms["pdbx_PDB_model_num"] = "1"
    return df_atoms


  @staticmethod
  def _init_new_sites(sites,build_hierarchy=True):
    """
    This is important. It is called by constructors to do most of what
    would usually be in an __init__ method. But it occurs AFTER __init__

    Call this when starting a new molecule.
    It avoids expensive and possibly undesired initialization upon
    successive dataframe creation.
    """
    sites._validate_and_fill_input_columns()
    sites._apply_dtypes_to_mmcif()
    sites._fill_str_NA()
    #sites._unescape_characters()
    # The unescape function doesn't seem to work, doing it individually to
    # atom_id in the build hierarchy section for now
    sites._prepare_hierarchy(sites)
    sites.detect_chain_breaks(sites)
    sites._annotate(sites)
    if build_hierarchy:
      sites._hierarchy = sites._build_hierarchy(sort=True)
    #sites.add_crystal_symmetry_if_necessary(sites)
    #sites.add_label_asym_id(sites)

    return sites

  @staticmethod
  def _prepare_hierarchy(sites):
    model_keys = ["pdbx_PDB_model_num"]
    object.__setattr__(sites, "_model_keys", model_keys)

    chain_keys = [
            "auth_asym_id",
            'chain_break_idx',

           ]
    object.__setattr__(sites, "_chain_keys", chain_keys)

    rg_keys = ([
                "auth_seq_id",
                "pdbx_PDB_ins_code",
               ]
              )
    object.__setattr__(sites, "_rg_keys", rg_keys)

    ag_keys = ([
                "auth_comp_id",
                "label_alt_id",
              ]
              )
    object.__setattr__(sites, "_ag_keys", ag_keys)

    xyz_keys = ["Cartn_x","Cartn_y","Cartn_z"]
    object.__setattr__(sites, "_xyz_keys", xyz_keys)


  def __init__(self, *args,name=None,**kwargs):
    #self.custom_attribute = kwargs.pop('custom_attribute', None)
    super().__init__(*args, **kwargs)

    object.__setattr__(self, "_mapper_applied", None)
    object.__setattr__(self, "_name", name)
    object.__setattr__(self, "_hierarchy", None)

    # Initialize hierarchy keys.
    self._prepare_hierarchy(self)

    # TODO: Move this to config
    self.comp_id_key = "auth_comp_id"
    self.asym_id_key = "auth_asym_id"


  @property
  def _constructor(self):
      # Constructor for the same type
      def _c(*args, **kwargs):
          return self.__class__(*args, **kwargs).__finalize__(self, method='inherit')
      return _c

    # def __finalize__(self, other, method=None, **kwargs):
    #     # Here, you can explicitly pass attributes to the derivative DataFrame
    #     if method == 'inherit':
    #         for name in ['comp_id_key']:
    #             setattr(self, name, getattr(other, name, None))
    #     return self



  def _validate_and_fill_input_columns(self,exceptions=[]):
    for col,func in self.params.fill_functions.items():
      if col in self:
        pass # Column needed and found
      elif func is None:
        if col not in exceptions:
          assert col in self, f"Required column is missing: {col}"
      else:
        # Try to fill
        self[col] = func(self)

  def _apply_dtypes_to_mmcif(self):
      dtypes_mmcif = self.params.dtypes_mmcif
      for col, dtype in dtypes_mmcif.items():
          if col in self.columns:
              self[col] = self[col].astype(dtype)

  def _apply_rounding_to_mmcif(self):
     mapper = self.params.rounding
     self.round(mapper,inplace=True)


  def with_names(self, mapper, inplace=True):
      if inplace:
          self.rename(columns=mapper, inplace=True)
      else:
          return self.rename(columns=mapper, inplace=False)
      return self

  @property
  def model_keys(self):
    return self._model_keys
  @property
  def chain_keys(self):
    return self.model_keys + self._chain_keys
  @property
  def rg_keys(self):
    return self._model_keys + self._chain_keys + self._rg_keys
  @property
  def ag_keys(self):
    return self._model_keys + self._chain_keys + self._rg_keys + self._ag_keys
  @property
  def hierarchy_keys(self):
    return self._model_keys + self._chain_keys + self._rg_keys + self._ag_keys
  @property
  def xyz_keys(self):
     return self._xyz_keys

  def with_mmcif_names(self, inplace=True):
      mapper = {key: value for key, value in self.params.attrs_map_mmcif.items()}
      return self.with_names(mapper, inplace=inplace)

  def with_phenix_names(self, inplace=True):
      mapper = {key: value for key, value in self.params.attrs_map_phenix.items()}
      return self.with_names(mapper, inplace=inplace)

  def with_standard_names(self, inplace=True):
      if self._mapper_applied:
          mapper = {value: key for key, value in self._mapper_applied.items()}
          return self.with_names(mapper, inplace=inplace)

  def as_DICT(self,mode="flat"):
    """
    Define the 'standard' way to get a dict

    mode='flat': returns a list of dicts for each site (row)
    """
    d = self.to_dict(orient="records")
    return d

  def as_JSON(self,mode="flat",indent=2):
    """
    Define the 'standard' way to get a json string

    mode='flat': returns a list of entries for each site (row)
    """
    d = self.as_DICT()
    j =  json.dumps(d,indent=indent)
    return j

  def as_DF_str(self):
    """
    Convert to string types.
    Creates a copy

    """
    df_str = self.astype("string").fillna(".")
    return df_str

  def as_DF_numeric(self):
    """
    Convert to prescribed types

    """
    self._apply_dtypes_to_mmcif()
    return self

  # def cleanup_columns(self):
  #   """
  #   Cleanup extra columns that have accumulated
  #   """
  #   self.drop(columns=self._cleanup_columns,inplace=True)
  #   self._cleanup_columns = []


  def get_attrs_map_phenix_to_mmcif(self):
    attrs_map_phenix_to_mmcif = {}

    for a, b in params.attrs_map_mmcif.items():
      if a in params.attrs_map_phenix:
        attrs_map_phenix_to_mmcif[b] = params.attrs_map_phenix[a]

    return {v:k for k,v in attrs_map_phenix_to_mmcif.items()}

  def _unescape_characters(sites):
    for col in ["label_atom_id","auth_atom_id"]:
      if col in sites:
        sites[col] = sites[col].str.replace("\\'", "'",regex=False)

  def _fill_str_NA(self):
    # fill with non-na values
    dtype_placeholders = {
    'string': '.',
    }
    col_name_placeholders = {
      "pdbx_PDB_ins_code":"?",
      "label_entity_id":"?",
    }

    for col in self.columns:
      if col in col_name_placeholders:
        placeholder = col_name_placeholders[col]
        self[col] = self[col].fillna(placeholder)

      else:
        dtype = str(self[col].dtype)
        if dtype in ["string","str"]:
          placeholder = dtype_placeholders['string']

          self[col] = self[col].fillna(placeholder)

  @property
  def xyz(self):

     return self[self.xyz_keys].values

  def get_sites_cart(self):
    sites_cart = flex.vec3_double(self.xyz)
    sites_cart_h = self.hierarchy.atoms().extract_xyz()
    assert approx_equal(sites_cart,sites_cart_h), (
        "Mismatch between dataframe cartesian sites and those in hierarchy"
     )
    return sites_cart_h


  @property
  def hierarchy(self):
    # if self._hierarchy is None:
    #    self._hierarchy = self._build_hierarchy()

    return self._hierarchy

  def add_i_seqs_from_index(self):
    # In theory this should always work, should never need the hierarchy method
    self.reset_index(drop=True,inplace=True)
    self["i_seq"] = self.index.values

  def _format_atom_name(self,atom_name, atom_type):
    if len(atom_name) < 4:


      # Strip and conditionally prefix atom_name based on atom_type match and length
      atom_name = atom_name.strip()
      if atom_name.upper().startswith(atom_type.upper()) and len(atom_type) == 1:
          atom_name = f" {atom_name}"

      # Ensure atom_name is exactly 4 characters, padded with spaces as needed
      atom_name = atom_name.ljust(4)
    atom_name =  atom_name.replace("\\'", "'").strip('"')
    return atom_name

  def _json_serializer(self,obj):
    if isinstance(obj, type(pd.NA)):
        return "."  # Or format however you prefer
    # Add more custom handling as necessary
    raise TypeError(f"Object of type {obj.__class__.__name__} is not JSON serializable")

  def _sort_with_hierarchy(sites,hierarchy):
    # sort according to the serial values in a hierarchy
    # NOTE: This MIGHT introduce issues with previously applied annotations...

    # resort according to hierarchy
    sorted_id = [atom.serial.strip() for atom in hierarchy.atoms()]
    sites['id'] = pd.Categorical(sites['id'], categories=sorted_id, ordered=True)
    # Sort the DataFrame by the 'id' column
    sites.sort_values('id',inplace=True)
    sites["id"] = sorted_id
    sites["id"] = sites["id"].astype("string")
    sites.reset_index(drop=True, inplace=True)

  def sort_ligand_atom_names(sites):
    # TODO: Deprecated, use sort_with_hierarchy
    df = sites
    sort_col = "label_atom_id"
    check_col = "residue_class"
    check_list = ["ligand","other"]
    group_cols = sites.ag_keys

    # Define groups to be sorted
    groups_to_sort = ["ligand","other"]  # Groups to sort
    sort_col = "label_atom_id"  # Column to sort within each sort group


    # Identify start indices for each group to be sorted
    df['Sort_Group_Idx'] = df.groupby(sites.ag_keys,sort=False).ngroup() + 1  # +1 to ensure 0 is not used as a valid group ID
    df['To_Sort'] = df['residue_class'].isin(groups_to_sort) * df['Sort_Group_Idx']

    # Sort each group individually
    for group_id in df['To_Sort'].unique():
        if group_id > 0:  # Check if it's a group we want to sort
            # Identify the start and end index of the current group
            group_indices = df.index[df['To_Sort'] == group_id].tolist()
            if len(group_indices)>1:
              start_idx, end_idx = group_indices[0], group_indices[-1] #+ 1  # +1 to make end_idx exclusive
              # Sort the subsection and update the DataFrame
              idcol = df.loc[start_idx:end_idx,"id"].values.copy()
              sorted_subsection = df.loc[start_idx:end_idx].sort_values(by=sort_col, ascending=True)
              sorted_subsection["id"] = idcol
              df.loc[start_idx:end_idx, sorted_subsection.columns] = sorted_subsection.values

    # Clean up: Remove the auxiliary columns
    df.drop(columns=['Sort_Group_Idx', 'To_Sort'], inplace=True)


    return df



  @staticmethod
  def detect_chain_breaks(sites,cleanup_columns=True):
    resname_key = "auth_comp_id"

    blanks = sites.params.blanks
    def add_residue_class(resname):
       if resname in blanks:
          return pd.NA
       else:
          return common_residue_names_get_class(resname)

    sites["residue_class"] = sites[resname_key].apply(add_residue_class)
    residue_categories = defaultdict(lambda: 'ligand')
    residue_categories.update({
      'common_amino_acid':'poly',
      'modified_amino_acid':'poly',
      'common_rna_dna':'poly',
      'modified_rna_dna':'poly',
      'common_water':'water',
    })
    sites['residue_category'] = [residue_categories[v] for v in sites["residue_class"]]
    # sites['category_change'] = sites['category'].ne(sites['category'].shift())


    # sites['from_poly'] = (sites['category'].shift() == 'poly') & sites['category'].ne('poly') # end poly chain
    # sites['from_poly_to_water'] = (sites['category'].shift() == 'poly') & sites['category'].eq('water') # end poly chain
    # sites['from_poly_to_ligand'] = (sites['category'].shift() == 'poly') & sites['category'].eq('ligand') # end poly chain
    # sites['from_water_to_non_water'] = (sites['category'].shift() == 'water') & sites['category'].ne('water')
    # sites['from_other_to_any'] = ((sites['category'].shift() != 'poly') & (sites['category'].shift() != 'water')) & sites['category'].isin(['poly', 'water', 'ligand'])

    sites['pdbx_PDB_model_num_change'] = sites['pdbx_PDB_model_num'].ne(sites['pdbx_PDB_model_num'].shift())
    sites['auth_asym_id_change'] = sites['auth_asym_id'].ne(sites['auth_asym_id'].shift())
    sites['label_asym_id_change'] = sites['label_asym_id'].ne(sites['label_asym_id'].shift())
    sites['auth_seq_id_change'] = sites['auth_seq_id'].ne(sites['auth_seq_id'].shift())

    sites['is_aa_or_rna_or_dna'] = (
        sites['residue_class'].eq('common_amino_acid') |
        sites['residue_class'].eq('modified_amino_acid') |
        sites['residue_class'].eq('common_rna_dna') |
        sites['residue_class'].eq('modified_rna_dna')
    )

    # Initialize the first row for change columns to False explicitly
    for col in ['pdbx_PDB_model_num_change', 'auth_asym_id_change', 'label_asym_id_change', 'auth_seq_id_change']:
        sites.loc[0, col] = False

    # Adjusting the chain_break calculation
    # Ensuring NA values resulting from .shift() are treated as False

    sites['chain_break_poly'] =(
         sites['label_asym_id_change'] &
         sites['is_aa_or_rna_or_dna'].shift(fill_value=False) #&
         #sites['is_aa_or_rna_or_dna']
        )

    # sites["chain_break_water"] = (
    #      sites['label_asym_id_change'] &
    #      sites['residue_category'].eq("water") &
    #      sites['residue_category'].shift(fill_value="other").eq("poly"))

    # sites["chain_break_ligand"] = (
    #      sites['label_asym_id_change'] &
    #      sites['residue_category'].eq("ligand") &
    #      sites['residue_category'].shift(fill_value="other").ne("ligand"))

    sites['chain_break'] = (sites['pdbx_PDB_model_num_change'] |
                           sites['auth_asym_id_change'] |
                           sites['chain_break_poly'] #|
                           # sites["chain_break_water"] |
                           # sites['chain_break_ligand']
                           ).astype(int)



    # Calculating chain_break_idx
    sites['chain_break_idx'] = sites['chain_break'].cumsum()

    cleanup_columns =[
      'chain_break_poly',
      "chain_break",
      #"chain_break_idx",
      'is_aa_or_rna_or_dna',
      'auth_seq_id_change',
      'label_asym_id_change',
      'auth_asym_id_change',
      'pdbx_PDB_model_num_change',
      'residue_category',
      #'residue_class',
    ]
    if cleanup_columns:
      sites.drop(columns=cleanup_columns,inplace=True)



  @staticmethod
  def add_label_asym_id(sites,cleanup_columns=True):

    assert "chain_break_idx" in sites

        # add the labels
    label_asym_ids = all_label_asym_ids()
    label_asym_ids_array = np.array(all_label_asym_ids())
    sites["label_asym_id"] = label_asym_ids_array[sites['chain_break_idx'].values]
    sites["label_asym_id"] = sites["label_asym_id"].astype("string")


    cleanup_columns = [
      "label_asym_id_check",
      "number_label_asym_id_check",
    ]
    if cleanup_columns:
      sites.drop(columns=cleanup_columns,inplace=True)



  def _build_hierarchy(sites,cleanup_columns=True,sort=True):
    """
    Build hierarchy. Assumes standard mmcif keys present.
    """

    # assemble some useful intermediate columns
    # for charge

    sites["charge_sign"] = np.sign(sites.pdbx_formal_charge.fillna(0))
    sites["charge_sign"] = sites["charge_sign"].astype('string')
    sites["charge_sign"].replace({"-1":'-',"1":"+"},inplace=True)
    sites["charge_abs"] = np.abs(sites.pdbx_formal_charge.values)
    sites["charge_abs"] = sites["charge_abs"].astype('string')
    sites["charge_string"] = sites["charge_abs"] + sites["charge_sign"]
    isna = sites["charge_string"].isna()
    sites.loc[isna,'charge_string'] = "  " # charge blank to match pdbh

    # Element preprocessing
    sites["element"] = sites["type_symbol"]
    isna = sites["element"].isna()
    sites.loc[isna,"element"] = ""
    #sites["element"] = sites["element"].apply(lambda x: x.rjust(2))

    # Serial preprocessing
    sites["serial"] = sites["id"]
    #sites["serial"] = sites["serial"].apply(lambda x: x.rjust(5))


    blanks = sites.params.blanks
    h_keys = ["model_idx","chain_idx","rg_idx","ag_idx"]

    h_dict = {}
    h = hierarchy.root()
    for i,t in enumerate(sites.itertuples()):
      m,c,rg,ag,a = t.model_idx,t.chain_idx,t.rg_idx,t.ag_idx,i
      m_id = (m,None,None,None,None)
      c_id = (m,c,None,None,None)
      rg_id = (m,c,rg,None,None)
      ag_id = (m,c,rg,ag,None)
      a_id =  (m,c,rg,ag,i)

      # Models
      if m_id not in h_dict:
        h_dict[m_id] = {}
        model = hierarchy.model(id=t.pdbx_PDB_model_num)
        h_dict[m_id]["cctbx_obj"] = model
        h.append_model(model)

      # chains
      if c_id not in h_dict:
        h_dict[c_id] = {}
        # make
        chain = hierarchy.chain(id=t.auth_asym_id)
        h_dict[c_id]["cctbx_obj"] = chain

        # attach
        model = h_dict[m_id]["cctbx_obj"]
        model.append_chain(chain)

      # rgs
      if rg_id not in h_dict:
        h_dict[rg_id] = {}
        #make
        resseq = t.auth_seq_id
        icode = t.pdbx_PDB_ins_code
        if icode in blanks:
          icode = " "
        rg = hierarchy.residue_group(resseq=str(resseq).rjust(4),icode=icode)
        h_dict[rg_id]["cctbx_obj"] = rg

        # attach
        chain = h_dict[c_id]["cctbx_obj"]
        chain.append_residue_group(rg)

      # ags
      if ag_id not in h_dict:
        h_dict[ag_id] = {"atoms":[]}

        # make ag
        altloc = t.label_alt_id
        resname = t.label_comp_id
        if altloc in blanks:
          altloc = ""
        ag = hierarchy.atom_group(altloc=altloc, resname=resname)
        h_dict[ag_id]["cctbx_obj"] = ag


        # attach
        rg = h_dict[rg_id]["cctbx_obj"]
        if altloc=="":
          rg.insert_atom_group(0, ag)
        else:
          rg.append_atom_group(ag)

      # Atoms
      assert a_id not in h_dict, "Duplicate atoms detected"
      atom = hierarchy.atom()
      atom.set_serial(t.serial)

      # store element


      # it's still not clear what to do with justification
      # atom.set_element(t.type_symbol.rjust(2))
      atom.set_element(t.element)

      # store name
      name = sites._format_atom_name(t.label_atom_id,t.type_symbol)
      atom.set_name(name)

      atom.set_xyz((t.Cartn_x,t.Cartn_y,t.Cartn_z))
      atom.set_b(t.B_iso_or_equiv)
      atom.set_occ(t.occupancy)
      atom.set_segid("    ")

      # set pdb group
      if t.group_PDB == "HETATM":
        atom.hetero = True

      #if t.charge_string is not pd.NA:
      #print("adding charge string")
      atom.set_charge(t.charge_string)

      ag = h_dict[ag_id]["cctbx_obj"]

      ag.append_atom(atom)


    # Sort
    h.sort_atoms_in_place()
    h.reset_atom_i_seqs()
    # DON'T reset serial, it's needed to map i_seqs in hierarchy to the dataframe

    if cleanup_columns:
      cols = [
         "charge_sign",
         "charge_abs",
         "charge_string",
      ]
      sites.drop(columns=cols,inplace=True)

    if sort:
       sites._sort_with_hierarchy(h)
    return h



  @staticmethod
  def _annotate(sites):
    # Annotations

    # Identify hierarchy membership
    sites["ag_idx"] = sites.groupby(sites.ag_keys,as_index=False,sort=False).ngroup()
    sites["rg_idx"] = sites.groupby(sites.rg_keys,as_index=False,sort=False).ngroup()
    sites["chain_idx"] = sites.groupby(sites.chain_keys,as_index=False,sort=False).ngroup()
    sites["model_idx"] = sites.groupby(sites.model_keys,as_index=False,sort=False).ngroup()

    # Annotate first and last residues in chain

    # Identify chain breaks
    sites['residue_change_shift_back'] = sites['rg_idx'] != sites['rg_idx'].shift(-1)
    sites['residue_change_shift_forward'] = sites['residue_change_shift_back'].shift(1).fillna(True)
    sites['chain_change_shift_back'] = sites['chain_idx'] != sites['chain_idx'].shift(-1)
    sites['chain_change_shift_forward'] = sites['chain_change_shift_back'].shift(1).fillna(True)

    sites['chain_first_atom'] =  ( sites['residue_change_shift_forward'] &
                                sites['chain_change_shift_forward']
                                )
    sites['chain_last_atom'] =  (sites['residue_change_shift_back'] &
                                sites['chain_change_shift_back']
                                )


    # propagate up
    sites['chain_first_residue'] = sites.groupby('rg_idx')['chain_first_atom'].transform('max')
    sites['chain_last_residue'] = sites.groupby('rg_idx')['chain_last_atom'].transform('max')


    # N-terminal or C-terminal residues would be

    # sites['chain_first_residue'] & sites['residue_class'].isin(["common_amino_acid","modified_amino_acid"])
    # sites['chain_last_residue'] & sites['residue_class'].isin(["common_amino_acid","modified_amino_acid"])

    # N-terminal or C-terminal atoms would be
    # (sites['chain_first_residue'] &
    #  sites['residue_class'].isin(["common_amino_acid","modified_amino_acid"]) &
    #  sites['atom_id'].eq("N"))
    #
    # (sites['chain_last_residue'] &
    #  sites['residue_class'].isin(["common_amino_acid","modified_amino_acid"]) &
    #  sites['atom_id'].eq("C"))

    # cleanup
    cleanup_cols = [
      'residue_change_shift_back',
      'residue_change_shift_forward',
      'chain_change_shift_back',
      'chain_change_shift_forward',
    ]
    sites.drop(columns=cleanup_cols,inplace=True)

    # End of "Annotate first and last residues in chain"

    # End Annotate
    return sites

def compare_hierarchies(h_left,h_right):

  print("starting hierarchical compare")
  atoms_left = list(h_left.atoms())
  atoms_right = list(h_right.atoms())
  assert len(atoms_left) == len(atoms_right)
  for model_left,model_right in zip(h_left.models(),h_right.models()):
    assert len(h_left.models()) == len(h_right.models())
    if model_left.id.strip() == model_right.id.strip(): # TODO: would be nice to not need strip
      model_left.id = model_right.id
    assert len(model_left.chains()) == len(model_right.chains())
    for chain_left,chain_right in zip(model_left.chains(),model_right.chains()):
      assert chain_left.id == chain_right.id
      for rg_left,rg_right in zip(chain_left.residue_groups(),chain_right.residue_groups()):
        assert int(rg_left.resseq)==int(rg_right.resseq)
        assert rg_left.unique_resnames()==rg_right.unique_resnames()
        assert rg_left.icode == rg_right.icode
        for ag_left,ag_right in zip(rg_left.atom_groups(),rg_right.atom_groups()):
          assert ag_left.altloc == ag_right.altloc
          for atom_left,atom_right in zip(ag_left.atoms(),ag_right.atoms()):
            assert atom_left.name.strip() == atom_right.name.strip()
            if atom_left.name != atom_right.name:
               print("Atom names differ in whitespace.")
            assert round(atom_left.b,2) == round(atom_right.b,2)
            xyz_left = [round(v,3) for v in atom_left.xyz]
            xyz_right = [round(v,3) for v in atom_right.xyz]
            assert xyz_left == xyz_right
            assert round(atom_left.occ,2) == round(atom_right.occ,2)
            assert atom_left.charge.strip() == atom_right.charge.strip()
            if atom_left.charge != atom_right.charge:
               print("Atom charges differ in whitespace.")
            assert atom_left.element.strip() == atom_right.element.strip()
            if atom_left.element != atom_right.element:
               print("Atom elements differ in whitespace.")
            #assert atom_left.serial_as_int() == atom_right.serial_as_int(), f"Failed at atom index: {atoms_left.index(atom_left)}"
            #assert atom_left.i_seq == atom_right.i_seq, f"Failed at atom index: {atoms_left.index(atom_left)}"


  # Finally, compare hierarchies
  if not h_right.is_similar_hierarchy(h_left):
     print("Hierarchies are similar but do not pass is_similar_hierarchy()")
  return True


def compare_sites(sites1, sites2,skip_columns=[]):
    sites1 = sites1.copy()
    sites2 = sites2.copy()
    df1,df2 = sites1,sites2
    columns_df1 = set(df1.columns)
    columns_df2 = set(df2.columns)

    # Find columns in df1 not in df2 and vice versa
    missing_in_df2 = columns_df1.difference(columns_df2)
    missing_in_df1 = columns_df2.difference(columns_df1)
    # Display the results
    #print("Columns in df1 missing in df2:", missing_in_df2)
    #print("Columns in df2 missing in df1:", missing_in_df1)

    # Compare shapes
    if df1.shape != df2.shape:
        return f"DataFrames have different shapes: {df1.shape} vs {df2.shape}"

    # Ensure both DataFrames have the same column order for comparison
    df1 = df1.reindex(columns=df2.columns)

    # Compare column names
    if not df1.columns.equals(df2.columns):
        return f"DataFrames have different columns."

    # Compare index
    if not df1.index.equals(df2.index):
        return f"DataFrames have different indexes."

    # Initialize a mask for differences
    diff_mask = pd.DataFrame(False, index=df1.index, columns=df1.columns)

    # Iterate over columns to compare values
    for col in df1.columns:
        if col not in skip_columns:
          values1 = df1[col].values
          values2 = df2[col].values
          assert np.all(values1==values2), f"Failed at column: {col}"

    if diff_mask.any().any():
        diff = pd.concat([df1[diff_mask], df2[diff_mask]], axis=0, keys=['df1', 'df2'])
        #return f"DataFrames differ:\n{diff}"
    else:
        return "DataFrames are identical."



class StringifyJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        # # Convert any object to its string representation
        # if isinstance(obj, (list, dict, int, float, bool)):
        #     return super().default(str(obj))
        return str(obj)
