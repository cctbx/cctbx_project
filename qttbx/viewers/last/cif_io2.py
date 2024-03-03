from io import StringIO
import warnings
from collections import UserDict, defaultdict
import re
from pathlib import Path
import sys

from iotbx import cif
from iotbx.cif.builders import crystal_symmetry_builder
from iotbx.cif.model import block
import pandas as pd
from pandas.errors import ParserWarning

from libtbx.utils import null_out
from .python_utils import find_key_path, get_value_by_path
from .cif_io import convert_dataframes_to_iotbx_cif, convert_dict_to_dataframes



class CifInput(UserDict):

  # Define explicit dtypes, what isn't here become 'string'
  dtypes = defaultdict(lambda: 'string')
  dtypes.update({
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
  })

  def __init__(self,*args,**kwargs):
    """
    Initialize with regular dictionary args,
      or a filename of a cif file
    """

    # Check if initialized with a filename
    self.filename = None
    if len(args) == 1 and isinstance(args[0], (Path,str)):
      self.filename = args[0]
      d = self._read_mmcif_file_to_dataframes(self.filename)
      args = args[1:]

    super().__init__(*args[1:],**kwargs)
    if self.filename is not None:
      self.update(d)

  def read_mmcif_file_to_dataframes(self,file):
    d = self._read_mmcif_file_to_dataframes(file)
    self.update(d)

  def _read_mmcif_file_to_dataframes(self,file):
    with open(file, 'r') as fh:
      lines = fh.readlines()

    self.parser= CifParser(lines)
    self.parser.run()
    dataframes = self.parser.process_to_pandas()
    return dataframes

  @property
  def atom_sites(self):
    # Search rather than assume location
    atom_site_path = find_key_path(self,"_atom_site")
    if atom_site_path is None:
      return None
    df =  get_value_by_path(self,atom_site_path)
    return df





  def extract_crystal_symmetry(self):
    """
    Extract the crystal symmetry as an object from a nested cif dict
    """

    # search the nested dict and retrieve value
    key_path = find_key_path(self,"_symmetry")
    d_symmetry = get_value_by_path(self,key_path)

    # value is dataframe, convert to dict
    d_cs =d_symmetry.to_dict(orient="records")[0]

    # join cif keys with dot
    d_cs = {"_symmetry."+key:value for key,value in d_cs.items()}

    # build the cctbx cif block
    b = block()
    for k,v in d_cs.items():
      b[k] = v

    # get crystal symmetry
    cs_builder = crystal_symmetry_builder(b)
    cs = cs_builder.crystal_symmetry
    return cs

  def find_dataframe(self,name):
    key_path = find_key_path(self,name)
    if key_path is None:
      # try with prefix
      name = "_"+name
      key_path = find_key_path(self,name)

    if key_path is None:
      return None
    d = get_value_by_path(self,key_path)
    return d


class CifParser:
  def __init__(self,lines):
    self.lines = lines
    self.i = 0
    self.columns = []
    self.data = []
    self.results = {}
    self.loop_parsing = False
    self.block_parsing = False
    self.current_prefix = None
    self.ready_to_eat_loop = False
    self.ready_to_eat_block = False
    self.ready_to_eat_key_value = False
    self.loop_marker = "loop_"
    self.last_result_key = None
    self.current_result = {}

  @property
  def line(self):
    return self.lines[self.i]

  def run(self):
    while self.i < len(self.lines)-1:
      self.read_next_line()
      self.end_line()
    self.finish_file()

  def clear(self):
    self.columns = []
    self.data = []
    self.loop_parsing = False
    self.block_parsing = False
    self.ready_to_eat_loop = False
    self.ready_to_eat_block = False
    self.ready_to_eat_key_value = False

  def eat(self,is_loop=False,is_block=False):

    if is_loop or is_block:
      key = self.current_prefix
      # if key is None:
      #   key = self.columns[0].split(".")[0]
      assert key is not None
      # if self.i ==3:
      #   import pdb
      #   pdb.set_trace()
      self.current_result["columns"] = self.columns
      self.current_result["data"] = self.data
      self.current_result["is_loop"] = is_loop
      self.current_result["is_block"] = is_block
      self.current_result["is_single_value"] = False
      self.current_result["prefix"]  = key
      self.current_result["stop"] = self.i

      if is_loop: # get the last line of a loop
        self.data.append(self.i)
      self.results[key] = self.current_result
    else:
      # simple key value pair
      assert len(self.columns)==1 and len(self.data)==1
      self.current_result["is_loop"] = False
      self.current_result["is_block"] = False
      self.current_result["is_single_value"] = True
      key = self.columns[0]
      value = self.data[0]
      self.current_result["key"] = key
      self.current_result["value"] = value
      assert key is not None
      self.results[key] = self.current_result

    # Finish eating
    self.last_result_key = key
    self.current_result = {"start":self.i+1}
    self.clear()


  def end_line(self):
    if self.ready_to_eat_key_value:
      self.eat(is_loop=False,is_block=False)
    elif self.ready_to_eat_loop:
      self.eat(is_loop=True,is_block=False)
    elif self.ready_to_eat_block:
      self.eat(is_loop=False,is_block=True)

    self.i += 1
    #print("\n\n\n")

  def finish_file(self):
    if self.loop_parsing:
      self.ready_to_eat_loop = True
    elif self.block_parsing:
      self.ready_to_eat_block = True
    elif len(self.columns)>0:
      self.ready_to_eat_key_value = True
    self.end_line()

  def process_block_line(self,line):
    #print("Processing block line: ",self.i)
    #print(line)
    #print()
    if line.startswith("_"):
      prefix,key,value = self.process_single_line_data(line,log=False)
      if prefix is not None:
        key = ".".join([prefix,key])
      self.current_prefix = prefix
      self.columns.append(key)
      self.data.append(value)
    else:
      self.current_result["is_strange"] = True
      if len(self.data)>0:
        if self.data[-1] == "":
          self.data[-1] = line # newline value without comment



  def start_block(self,prefix,key,value):
    #print("Starting block with (prefix,key):",prefix,key)
    self.current_prefix = prefix
    self.columns.append(key)
    self.data.append(value)
    self.block_parsing = True
    #print("Start block on line:",self.i)

  def start_loop(self):
    #print("Start loop on line:",self.i)
    self.current_prefix = None
    self.loop_parsing = True

  def process_single_line_data(self,line,log=False):
    if log:
      #print("Processing single value line: ",self.i)
      pass
    line_split = split_with_quotes(line)
    if len(line_split)==1:
      line_split = (line_split[0],"")
    assert len(line_split)==2, f"Failed to split key,value pair into two parts: {line} was split to: {line_split}"
    # check prefix
    prefix = None
    key,value = line_split
    prefix_split = key.split(".")
    if len(prefix_split)>1:
      assert len(prefix_split)==2
      prefix,key = prefix_split

    return prefix,key,value


  def read_next_line(self):
    #print("Reading line: ",self.i)

    line = self.lines[self.i]
    stripped_line = line.strip()
    is_key_row = stripped_line.startswith('_')
    is_break_row = (len(stripped_line)==0 or
                     stripped_line.startswith("#")
                     )
    is_loop_entry_row = stripped_line.startswith(self.loop_marker)
    is_multi_line = stripped_line.startswith(";")

    # First, handle multi line
    if is_multi_line:
      self.current_result["is_strange"] = True
      #print("multi line start",self.i)
      multiline_value, last_index = handle_multiline(self.lines, self.i + 1)
      if len(self.columns)==0:

        if isinstance(self.results[self.last_result_key],dict):
          self.results[self.last_result_key]["data"][-1] = multiline_value
        else:
          assert self.last_result_key in self.results
          assert self.last_result_key is not None
          self.results[self.last_result_key] = multiline_value
      else:
        self.data[-1] = multiline_value
      self.i = last_index  # Skip processed lines
      return


    # Now consider block or loop parsing...
    # Handle the case of no obvious separator

    is_no_sep_break = False
    if self.loop_parsing or self.block_parsing:
      if not is_break_row and (is_key_row or is_loop_entry_row):
        if is_loop_entry_row:
          is_no_sep_break = True
        elif self.current_prefix is not None and not stripped_line.startswith(self.current_prefix):
          #print("setting sep break",self.current_prefix,self.line)
          is_no_sep_break = True



    # Handle blocks
    if self.block_parsing:
      if is_break_row or is_no_sep_break:
        if is_no_sep_break:
          #self.i-=1
          pass
        self.ready_to_eat_block = True
        return
      else:
        # Not breaking
        self.process_block_line(stripped_line)
        return



    # Handle loops
    if self.loop_parsing:
      if is_break_row or is_no_sep_break:
        if is_no_sep_break:
          self.i-=1
        self.ready_to_eat_loop = True
        return
      elif is_key_row:
        #print("Processing loop column line:",self.i)
        self.columns.append(stripped_line)
        if self.current_prefix is None:
          prefix,columns,data = self.process_single_line_data(stripped_line,log=False)
          self.current_prefix = prefix
        return

      else:
        # assume loop data line
        #print("Processing loop data line:",self.i)
        self.data.append(self.i)
        return


    # Check for new things
    assert not self.loop_parsing and not self.block_parsing
    if is_loop_entry_row:
      # start loop
      self.start_loop()
      return

    elif is_key_row:
      prefix,key,value = self.process_single_line_data(stripped_line)
      if prefix is not None:
        self.start_block(prefix,key,value)
        return
      else:
        # simple key value
        self.columns.append(key)
        self.data.append(value)
        self.ready_to_eat_key_value = True
        return

    # Should never reach here
    self.current_result["is_strange"] = True


  def process_to_pandas(self):
    """
    Process the output of parse_mmcif_file. This function
    separates all the pandas functionality from the file io
    """
    lines = self.lines
    output_dfs = {}
    for i,(df_name,result) in enumerate(self.results.items()):
      if not isinstance(result,dict):
        output_dfs[df_name] = result
        continue

      # inline iotbx function alternative
      def do_iotbx_cif_parsing(string):

        reader = cif.reader(file_object=StringIO(string))
        model = reader.model()
        d = convert_iotbx_cif_to_dict(model)
        d = d[f"loop_{i}"]
        assert len(d.keys())==1 and list(d.keys())[0] == loop_name, f"Expected loop to have 1 stem key: {loop_name}, got: {d.keys()}"
        d = d[loop_name]
        df = convert_dict_to_dataframes(d)
        return df

      if "is_strange" in result:
        try:
          start = result["start"]
          stop = result["stop"]
          strange_lines = lines[start:stop]
          string =  "\n".join(strange_lines)
          df = do_iotbx_cif_parsing(string)
          output_dfs[loop_name] = df
          continue
        except:
          continue

      # determine column stem name
      columns = result["columns"]
      has_prefix = "prefix" in result
      if has_prefix:
        prefix = result["prefix"]
        col_stem = prefix # alias
        col_names = [col.replace(prefix,"") for col in columns]
      # Process loop
      if result["is_loop"]:
        loop_name = df_name
        start = result["data"][0]
        stop = result["data"][-1]
        columns = result["columns"]


        # prepare file io object with the loop string
        loop_data_string = "".join(lines[start:stop])
        loop_data = StringIO(loop_data_string)


        # set dtypes for each column
        if has_prefix:
          dtype_keys = col_names
        else:
          dtype_keys = columns
        dtypes = {col:CifInput.dtypes[col] for col in dtype_keys}


        with warnings.catch_warnings(record=True) as w:
          warnings.simplefilter("always", ParserWarning)
          try:
            loop_df = pd.read_csv(loop_data,
                            sep='\s+',
                            engine='c',
                            quotechar="'",
                            dtype=dtypes,
                            quoting=0,
                            names=columns,
                            index_col=False,
                            low_memory=False,
                            keep_default_na=False
                            )


            if has_prefix:
              loop_df.rename(columns={col:col.replace(col_stem+".","") for col in loop_df.columns},inplace=True)
              output_dfs[loop_name] = loop_df
            else:
              for col in loop_df.columns:
                output_dfs[col] = loop_df[col].values
          except ParserWarning:
            print(f"ParserWarning caught: Mismatch between column names and data. Attempting iotbx parsing method for loop: {loop_name}")
            loop_string =  f"data_loop_{i}\n" + "\n".join(["loop_"]+columns)+"\n"+loop_data_string
            df = do_iotbx_cif_parsing(loop_string)
            output_dfs[loop_name] = df
          except Exception as e:
              print(f"An unexpected error occurred: {e}")
              print(f"Failed with: {loop_name}")
      elif result["is_block"]:
        assert has_prefix, "The term 'block' here requires a _prefix."
        # Process data block
        columns, data = result["columns"], result["data"]
        if len(data)==len(columns):
          data = [data]
        assert len(data)==1, f"Expected a single row for key,value cif block, got: {data}"
        assert len(columns) == len(data[0]), f"Failed to interpret key,value cif block for columns: {columns} and data: {data}"

        if has_prefix:
          # transpose
          data = [list(row) for row in zip(*data)]
          block_df_data = {col:d for col,d in zip(columns,data)}

          block_df = pd.DataFrame(block_df_data)
          block_df.rename(columns={col:col.replace(df_name+".","") for col in block_df.columns},inplace=True)
          output_dfs[df_name] = block_df
      elif result["is_single"]:
        # Single data line, just add to the top level dict
        for col,value in zip(columns,data[0]):
          output_dfs[col] = value

    # Postprocess
    for df_name,df in output_dfs.items():
      for col in df.select_dtypes(include=[object,"string"]):
        df[col] = df[col].str.replace("\\'", "'",regex=False)
        output_dfs[df_name] = df

    return output_dfs




# compile re and inline split function
QUOTES_REGEX = re.compile(r'"[^"]*"|\'[^\']*\'|\S+')
def split_with_quotes(line):
  parts = QUOTES_REGEX.findall(line)
  return [part[1:-1] if (part[0]=="'" or part[0]=='"') else part for part in parts]

# inline multiline handler
def handle_multiline(lines, current_index):
  multiline_content = [";"]
  for j in range(current_index, len(lines)):
      line = lines[j].strip()
      if line == ";":  # End of multiline value
          multiline_content.append(";")
          break
      multiline_content.append(line)
  return "\n".join(multiline_content), j  # Return the concatenated multiline value and the last index
