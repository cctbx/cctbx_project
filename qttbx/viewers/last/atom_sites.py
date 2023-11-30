from collections import OrderedDict
import numpy as np
import pandas as pd
import networkx as nx
import copy

from .selection_common import CommonSelectionParser
from .selection_utils import (
  df_nodes_group_to_query,
  find_simplest_selected_nodes,
  form_simple_str_common,
  SelConverterPhenix,
  Selection,
  SelectionQuery
)
from .pandas_utils import cctbx_atoms_to_df, cif_df_to_numeric
from .cif_io import read_cif_file, find_key_in_dict


class AtomSites(pd.DataFrame):
  @classmethod
  def from_file_cif(cls,filename,method="pdbe",**kwargs):
    # read from _atom_site in cif file
    cif_dict = read_cif_file(filename,method=method)
    return cls.from_atom_site_dict(cif_dict,**kwargs)

  @classmethod
  def from_atom_site_dict(cls,atom_site_dict,**kwargs):
    counter, value = find_key_in_dict(atom_site_dict, '_atom_site')
    assert counter==1, f"Did not find exactly one '_atom_site' key in dict: {atom_site_dict}"
    df = pd.DataFrame(value)
    df = AtomSites.insert_defaults_for_missing_columns(df,AtomSites._attrs_default_values)
    df = cif_df_to_numeric(df)
    return cls(df,**kwargs)

  @classmethod
  def from_mmtbx_model(cls,model,**kwargs):
    atoms = model.get_atoms()
    return cls.from_cctbx_atoms(atoms,params={"model"},**kwargs)

  @classmethod
  def from_cctbx_atoms(cls,atoms,**kwargs):
    df = cctbx_atoms_to_df(atoms)
    df = AtomSites.insert_defaults_for_missing_columns(df,AtomSites._attrs_default_values)
    return cls(df,**kwargs)
    
  _attrs_map_default = OrderedDict([
      ("asym_id",["auth_asym_id","label_asym_id","chain"]),
      ("seq_id",["auth_seq_id","label_seq_id","resseq"]),
      ("comp_id",["auth_comp_id","label_seq_id","resname"]),
      ("atom_id",["auth_atom_id","label_atom_id","name"]),
      ("id",["id"]), # i_seq?
      ("type_symbol",["type_symbol","element"]),
      ("x", ["Cartn_x"]),
      ("y",["Cartn_y"]),
      ("z",["Cartn_z"]),
      ("occupancy",['occupancy','occ']),
      ("B_iso_or_equiv",['bfactor','b','B'])])
  
  @staticmethod
  def map_core_attrs_to_columns(df, attrs_core_map):
    attrs_core = list(attrs_core_map.keys())
    mapping = {}  # Initialize the mapping dictionary
  
    for attr in attrs_core:
      if attr in df.columns:
        # If the core attribute is already in the DataFrame columns, map it to itself
        mapping[attr] = attr
      else:
        # Check if any of the aliases are present in the DataFrame columns
        if attr in attrs_core_map:
          found = False
          for alias in attrs_core_map[attr]:
            if alias in df.columns:
              mapping[attr] = alias
              found = True
              break
          if not found:
            raise ValueError(f"For mapping: {attrs_core_map} Neither attribute '{attr}' nor its aliases are present in DataFrame columns: {list(df.columns)}")
        else:
          raise ValueError(f"For mapping: {attrs_core_map} Attribute '{attr}' and its aliases are not present in DataFrame columns: {list(df.columns)}")
    
    return mapping

  _attrs_hierarchy_core_default = [
    'asym_id',
    'seq_id',
    'comp_id',
    'atom_id',
    'id',
    'type_symbol',
    'occupancy',
    'B_iso_or_equiv'
    ]

  _attrs_default_values = OrderedDict([
      ("asym_id","?"),
      ("seq_id","?"),
      ("comp_id","?"),
      ("atom_id","?"),
      ("id","?"), # i_seq?
      ("type_symbol","?"),
      ("x",np.nan),
      ("y",np.nan),
      ("z",np.nan),
      ("occupancy",np.nan),
      ("B_iso_or_equiv",np.nan)])

  @staticmethod
  def insert_defaults_for_missing_columns(df,attrs_default_values):
    for attr,value in attrs_default_values.items():
      if attr not in df.columns:
        df[attr] = value
    return df
    
  def __init__(self,data,*args,attrs_map=None, params={},attrs_hierarchy_core=None,**kwargs):
    if attrs_hierarchy_core is None:
      attrs_hierarchy_core = self._attrs_hierarchy_core_default
    if isinstance(data, pd.DataFrame):
      # Pass the internal data block manager directly to avoid a copy
      super().__init__(data._data, *args, **kwargs)
    else:
      super().__init__(data, *args, **kwargs)
    # Custom attributes need to be defined differently to avoid pandas adding as columns
    # Use object.__setattr__ to safely set the custom attribute
    
    
    object.__setattr__(self,"params", params)
    object.__setattr__(self,"_default_str_format", "phenix")
    
    
    # First build a mapper between "core" attrs and column names ("original" attrs)
    if attrs_map is None:
      attrs_map = copy.deepcopy(self._attrs_map_default)
  
    object.__setattr__(self,"attrs_map", self.map_core_attrs_to_columns(self,attrs_map))
    object.__setattr__(self,"attrs_core_map",{v:k for k,v in self.attrs_map.items()})
    object.__setattr__(self,"attrs_hierarchy_core", attrs_hierarchy_core)
    attrs_hierarchy = [self.attrs_map[attr] for attr in self.attrs_hierarchy_core]
    object.__setattr__(self,"attrs_hierarchy", attrs_hierarchy)


    # build hierarchy graph
    object.__setattr__(self, '_G', self._create_hierarchy_graph(self.attrs_hierarchy_core))
    
    self._validate()
    
  def _validate(self):
    assert len(set(self.core["id"].values)) == len(self), "The 'id' field must be unique for each atom"
  
  @property
  def attrs_core(self):
    return list(self.attrs_core_map.keys())
  @property
  def core(self):
    # set instance dataframes
    return self[self.attrs_core_map.keys()].rename(columns=self.attrs_core_map)
  
  @property
  def model(self):
    assert "model" in self.params, "A mmtbx model was not used to initialize this instance"
    return self.params["model"]
    
  @property
  def model(self):
    assert "filename" in self.params, "A filename was not provided for this instance"
    return self.params["filename"]
    
  @property
  def G(self):
    # the networkx graph for the macromolecular hierarchy
    if self._G is None:
      self._G = self._create_hierarchy_graph(self.attrs_hierarchy_core)
    return self._G

  def _create_hierarchy_graph(self,columns):
      df = self.core
      G = nx.DiGraph()
    
      # Create root node
      root = "root"
      G.add_node(root, ids=[])
    
      # Function to add IDs to a node and its ancestors
      def add_id_to_ancestors(G, node, id_value):
        G.nodes[node]['ids'].append(id_value)
        for parent in G.predecessors(node):
          add_id_to_ancestors(G, parent, id_value)
    
      for row in df.itertuples(index=False):
        cur_node = root  # Start at the root for each row
        id_value = getattr(row, "id")
    
        # Initialize a list of "*" values with the same length as "columns"
        cur_level_columns = ["*" for _ in columns]
    
        for idx, col in enumerate(columns):
          parent_node = cur_node  # The current node becomes the parent
    
          # Fill in the next value in "cur_level_columns"
          cur_level_columns[idx] = getattr(row, col)
    
          # Create a new node for the current row and column
          cur_node = tuple(cur_level_columns.copy())  # Make a copy of the list and convert to tuple
    
          # Add the node to the graph if it doesn't exist, then connect it to its parent
          if cur_node not in G:
            G.add_node(cur_node, ids=[])
            G.add_edge(parent_node, cur_node)
    
          # Add the ID to the current node and its ancestors
          add_id_to_ancestors(G, cur_node, id_value)
      return G
    
  ###############################
  #### Starting Selection Fresh #
  ###############################
  def select(self):
    raise NotImplementedError

  def select_from_query(self,query):
    if len(query.selections)==0:
      return self.__class__(self[self.index < 0]) # empty
    else:
      return self._convert_query_to_sites(query)

  def _convert_query_to_sites(self,query):
    """
    Convert a SelectionQuery object to a 
    subset sites data frame (This class)
    """
    # get a pandas query
    return self._pandas_query_to_sites(query.pandas_query)


  

  def _pandas_query_to_sites(self,pandas_query):
    # rename columns
    self.rename(columns=self.attrs_core_map, inplace=True)
    
    # perform query
    if copy:
      # returns a copy
      df_sel = self.query(pandas_query).copy()
      # must rename both
      self.rename(columns=self.attrs_map, inplace=True)
      df_sel.rename(columns=self.attrs_map, inplace=True)
      
    else:
      # returns a subset
      df_sel = self.query(pandas_query)
      # rename back (only one necessary)
      self.rename(columns=self.attrs_map, inplace=True)
    
    # return new selected df
    return self.__class__(df_sel,attrs_map=self.attrs_map)

  def _select_query_from_str_phenix(self,str_phenix):
    converter = SelConverterPhenix()
    sel_str_common = converter.convert_phenix_to_common(str_phenix)
    return self._select_query_from_str_common(sel_str_common)

  def _select_query_from_str_common(self,str_common):
    """
    Given a selection string, go to sites -> query
    """
    sites_sel = self._select_sites_from_str_common(str_common)
    query = self._convert_sites_to_query(sites_sel)
    return query


  def _find_simplest_nodes_from_sites(self,sites_sel):
    """
    Utility function to reduce a subset sites dataframe
    together with a graph (self.G) to get the minimum number
    of nodes to fullfill the subset.
    """
    assert isinstance(sites_sel,self.__class__), f"Provide an atom sites object, not: {type(sites_sel), {print(sites_sel)}}"
    nodes = find_simplest_selected_nodes(self.G,sites_sel.core["id"])
    if nodes == ["root"]:
      nodes = [["*" for _ in sites_sel.attrs_hierarchy_core]]
      nodes = pd.DataFrame(nodes,columns=sites_sel.attrs_hierarchy_core)
    else:
      nodes = pd.DataFrame(nodes,columns=sites_sel.attrs_hierarchy_core)
    return nodes

  def _convert_sites_to_query(self,sites_sel):
    """
    With a subset sites data frame, convert to a
    SelectionQuery object
    """
    assert isinstance(sites_sel,self.__class__), f"Provide an atom sites object, not: {type(sites_sel), {print(sites_sel)}}"
    # find the minimum selections required
    nodes = self._find_simplest_nodes_from_sites(sites_sel)

    # Group by seq range and convert to query
    query = df_nodes_group_to_query(nodes)
    return query

  def _convert_sites_to_str_common(self,sites_sel):
    str_common = form_simple_str_common(sites_sel)
    return str_common

  def _select_sites_from_str_common(self,str_common):
    """
    This is the main path from a Phenix/common string to an
    atom sites dataframe, Going through an the ast.
    """
    sel_str_common = str_common
    # remove 'sel' statements
    if sel_str_common.startswith("select" ):
      sel_str_common = sel_str_common[7:]
    elif sel_str_common.startswith("sel "):
      sel_str_common = sel_str_common[4:]

    # parse to ast
    parser = CommonSelectionParser(sel_str_common, debug=False)
    parser.parse()
    
    # interpret ast as pandas selection query
    query = parser.to_pandas_query()
    #print("Pandas query:")
    #print(query)    
    return self._pandas_query_to_sites(query)

  # ## Public selection methods
  # #############
  # def select(self,arg,input_format='phenix',return_format="sites"):
  #   """
  #   Select atoms with various input/output types
  #   """
  #   return_type = return_format
  #   str_format = input_format
  #   assert str_format in ["phenix","common",'json'], f"Str format specified: {str_format} not recognized"
  #   assert return_type in ["phenix",'common','sites','query','json'], f"Provide valid return type, not {return_type}"

  #   # get the sites first
  #   sites = self._select(arg,str_format=str_format)
  #   if return_type == 'sites':
  #     return sites

  #   # possibly convert to query object
  #   elif return_type in ['query','json','phenix','common']:

  #     # 'simplify' it to a more compact form than per-atom
  #     df_sel_simple = self._simplify_selected_df(sites)

  #     if return_type in ['query','json']:
  #       # convert to a selection query object
  #       query = self._simple_df_to_query(df_sel_simple)
  #       if return_type == 'query':
  #         return query
  #       elif return_type == 'json':
  #         return query.to_json()
  #     elif return_type in ['phenix','common']:

  #       str_simple = form_simple_str_common(df_sel_simple)
  #       if return_type == "common":
  #         return str_simple
  #       elif return_type == 'phenix':
  #         converter = SelConverterPhenix()
  #         str_simple_phenix = converter.convert_common_to_phenix(str_simple)
  #         return str_simple_phenix

  # def _select(self,arg,str_format=None):
  #   """
  #   Select but it only returns one type, another instance of self.__class__
  #   """
  #   if isinstance(arg,str):
  #     if str_format is None:
  #       str_format = self._default_str_format
  #     assert str_format in ["phenix","common",'json'], f"Str format specified: {str_format} not recognized"
  #     if str_format == "common":
  #       return self._select_from_common_str(arg)
  #     elif str_format == 'phenix':
  #       return self._select_from_phenix_str(arg)
  #     elif str_format == 'json':
  #       query = SelectionQuery.from_json(arg)
  #       raise NotImplementedError
  #       #query = self.create_query_string(records,self.columns)
  #       #df_sel = self.query(query)
  #       #return self.__class__(df_sel,attrs_map=self.attrs_map)
  #   elif isinstance(arg,pd.DataFrame):
  #     if not isinstance(arg,self.__class__):
  #       arg = self.__class__(arg,attrs_map=self.attrs_map)
  #     return arg
  
  # def _select_from_query_obj(self,query):
  #   assert isinstance(query,SelectionQuery), "Provide SelectionQuery object"


  # #############
  # ## Private selection methods
  # #############

  # def _simple_df_to_query(self,df):
  #   records = df.to_dict(orient='records')

  #   # Create Selection objects and store them in a list
  #   selections = [Selection.from_df_record(record) for record in records]
  #   query = SelectionQuery(selections=selections)
  #   return query

                         
  # def _select_from_phenix_str(self,sel_str_phenix,copy=False):
  #   # convert to common selection string (very similar)
  #   converter = SelConverterPhenix()
  #   sel_str_common = converter.convert_phenix_to_common(sel_str_phenix)
  #   return self._select_from_common_str(sel_str_common,copy=False)
    
  # def _select_from_common_str(self,sel_str_common,copy=False):

  #   # remove 'sel' statements
  #   if sel_str_common.startswith("select" ):
  #     sel_str_common = sel_str_common[7:]
  #   elif sel_str_common.startswith("sel "):
  #     sel_str_common = sel_str_common[4:]

  #   # parse to ast
  #   parser = CommonSelectionParser(sel_str_common, debug=False)
  #   parser.parse()
    
  #   # interpret ast as pandas selection query
  #   query = parser.to_pandas_query()
    
  #   # rename columns
  #   self.rename(columns=self.attrs_core_map, inplace=True)
    
  #   # perform query
  #   if copy:
  #     # returns a copy
  #     df_sel = self.query(query).copy()
  #     # must rename both
  #     self.rename(columns=self.attrs_map, inplace=True)
  #     df_sel.rename(columns=self.attrs_map, inplace=True)
      
  #   else:
  #     # returns a subset
  #     df_sel = self.query(query)
  #     # rename back (only one necessary)
  #     self.rename(columns=self.attrs_map, inplace=True)
    
  #   # return new selected df
  #   return self.__class__(df_sel,attrs_map=self.attrs_map)

  # def _simplify_selected_df(self,df_sel):
  #   """
  #   Two major transformations on a selected dataframe
  #   1. Find highest graph level that includes only the selected atoms (simplify)
  #   2. Merge seq_id ranges (group)

  #   This output is usually what you want to transform into string selections
  #   """
  #   # Simplify selected atoms to highest level nodes
  #   simplest_nodes = find_simplest_selected_nodes(self.G,df_sel.core["id"])

  #   # Merge residue ranges
  #   df_simple_sel = group_seq_range(simplest_nodes,columns=self.attrs_hierarchy_core)

  #   return df_simple_sel

  # def _simplify_selection_str_common(self,sel_str_common):
  #   # select an atom df
  #   df_sel = self._select_from_common_str(sel_str_common)

  #   # simplify the df
  #   df_simple_sel = self._simplify_selection(df_sel) #TODO: defined?
    
  #   # make a final additive statement of rows in final_df
  #   sel_str_simple = form_simple_str_common(df_simple_sel)
  #   return sel_str_simple

  # def _simplify_selection_str_phenix(self,sel_str_phenix):

  #   converter = SelConverterPhenix()
  #   sel_str_common = converter.convert_phenix_to_common(sel_str_phenix)
  #   sel_str_simple_common = self.simplify_selection_str_common(sel_str_common)
  #   converter = SelConverterPhenix()
  #   sel_str_simple_phenix = converter.convert_common_to_phenix(sel_str_simple_common)
  #   return sel_str_simple_phenix

  # def _select_from_selection(self,selection):
  #   if not isinstance(selection,(np.ndarray,list)):
  #     # assume cctbx object
  #     selection = selection.as_numpy_array()
    
  #   # check type
  #   e = selection[0]
  #   if not isinstance(e,bool):
  #     selection_bool = np.full(len(self),False)
  #     selection_bool[selection] = True
  #     selection = selection_bool
    
  #   # select atoms using phenix selection
  #   df_sel = self.core[selection]
  #   return df_sel

  

  # def _select_simple_from_selection(self,selection):
  #   df_sel = self._select_from_selection(selection)
    
  #   # Simplify selected atoms to highest level nodes
  #   simplest_nodes = find_simplest_selected_nodes(self.G,df_sel.index.values)

  #   # Merge residue ranges
  #   df_simple_sel = group_seq_range(simplest_nodes,columns=self.attrs_hierarchy_core)
  #   return df_simple_sel
  
  # def _str_common_from_selection(self,selection):
  #   df_simple = self._select_simple_from_selection(selection)
  #   str_simple = form_simple_str_common(df_simple)
  #   return str_simple
    
  # def _str_phenix_from_selection(self,selection):

  #   str_simple_common = self._str_common_from_selection(selection)
  #   converter = SelConverterPhenix()
  #   str_simple_phenix = converter.convert_common_to_phenix(str_simple_common)
  #   return str_simple_phenix

  # def create_query_string(self,subset_records, columns):
  #   # selection to df query string.
  #   # TODO: Need to unify all selection language (again)
  #   query_parts = []
  #   for record in subset_records:
  #     conditions = []
  #     for col in columns:
  #       if col in record:
  #         conditions.append(f"{col} == {repr(record[col])}")
  #       elif f"{col}_start" in record and f"{col}_stop" in record:
  #         start_val = repr(record[f"{col}_start"])
  #         stop_val = repr(record[f"{col}_stop"])
  #         conditions.append(f"({col} >= {start_val} and {col} < {stop_val})")

  #     if conditions:  # Skip if there are no conditions for this record
  #       query_part = " and ".join(conditions)
  #       query_parts.append(f"({query_part})")

  #   query_string = " or ".join(query_parts)
  #   return query_string
    




    