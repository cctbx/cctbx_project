"""
A companion to the Selection data type. Used to convert between different selection syntax
"""
import json
import re

from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, Sorry

#Parameters for parsing and translating Phenix-type syntax (chain, resname, resseq, etc)
# to the mmcif-type syntax used in Molstar (asym_id, comp_id, seq_id, etc)

# A compositional label, should be able to uniquely identify an atom. Analogous to id_str
compositional_labels_justify_map = { # What to include, and how much space for each attribute
    "asym_id":4,
    "comp_id":5,
    "seq_id":5,
    "atom_id":4,
    "alt_id": 2,
}

core_map_to_mmcif = {
# Maps a 'core' attribute key to a default 'real' attribute key
'asym_id':"auth_asym_id",
'seq_id':"auth_seq_id",
'comp_id':"label_comp_id",
'atom_id':"label_atom_id",
'type_symbol':"type_symbol",
'alt_id':'label_alt_id',
'id':'id',
}
core_map_to_core = {v:k for k,v in core_map_to_mmcif.items()}

# Maps mmcif attribute name to phenix name
attrs_map_to_phenix = {
'auth_asym_id': 'chain',
'label_asym_id': 'chain',
'auth_seq_id': 'resseq',
'label_seq_id': 'resseq',
'label_comp_id': 'resname',
'type_symbol': 'element',
'label_alt_id': 'altloc',
'pdbx_PDB_model_num': 'model',
'pdbx_PDB_ins_code': 'icode',
'Cartn_x': 'x',
'Cartn_y': 'y',
'Cartn_z': 'z',
'B_iso_or_equiv': 'bfactor',
'occupancy': 'occupancy',
'pdbx_formal_charge': 'charge',
'group_PDB': 'group',
'id': 'id',
'label_atom_id': 'name',
'label_entity_id': 'segid'
}
# The reverse, 
attrs_map_to_mmcif = {}
for k,v in attrs_map_to_phenix.items():
    core_key = k.replace("auth_","").replace("label_","")
    if core_key in core_map_to_mmcif:
        mmcif_key = core_map_to_mmcif[core_key]
        attrs_map_to_mmcif[v] = mmcif_key # set the core:real map
        attrs_map_to_mmcif["auth_"+core_key] = mmcif_key # alias
        attrs_map_to_mmcif["label_"+core_key] = mmcif_key # alias


    else:
        attrs_map_to_mmcif[v] = k

#Create list of all keywords for parsing
keywords_all = list((set(list(attrs_map_to_mmcif.keys())) | 
                     set(list(attrs_map_to_mmcif.values())) | 
                     set(list(core_map_to_mmcif.keys())) | 
                     set(list(core_map_to_mmcif.values()))))

# Convert logic operators between pandas query and phenix selection
logic_map_to_pandas ={ 
    'or': '|',
    'and': '&',
    'not': '~',
    "and":"&",
    "or":"|",
    "not":"~",
    " ": "==",
    ">=": ">=",
    ">": ">",
    "<=": "<=",
    "<": "<",
    "=":"==",
    "==":"==",
}
# The reverse
logic_map_to_phenix = {v:k for k,v in logic_map_to_pandas.items()}

logic_map_molstar = {
    "==": "eq",
    ">=": "gre",
    ">": "gr",
    "<=": "lte",
    "<": "lt",
    "!=": "neq",
}



class PhenixParser:
  """
  This class parses text selections from varied sources (phenix, molstar, chimera) and builds
    a common abstract syntax tree (AST). From that tree, any other selection statement can be made. 

  So there is some duplication of functionality between this class and iotbx.pdb.atom_selection. 
    The rationale for that is:
      1. There is behavior in iotbx.pdb.atom_selection that is not convertable to other selection languages. 
          This parser aims to be a common parser where only convertable expressions are processed.
      2. Adding conversion functionality directly to iotbx.pdb.atom_selection would 
          require extensive modifications to the existing, well-tested code.
  """
  @classmethod
  def from_phenix_string(cls,phenix_string,debug=False,verify=True):
    return cls(phenix_string,debug=debug,verify=verify)

  @staticmethod
  def remove_whitespace_around_colon(s):
    return re.sub(r'\s*:\s*', ':', s)
  
  @staticmethod
  def remove_top_level_parentheses(expr):
    start = 0
    end = len(expr) - 1
    while start <= end and expr[start] == '(' and expr[end] == ')':
      count = 0
      remove = True
      for i in range(start + 1, end):
        if expr[i] == '(':
          count += 1
        elif expr[i] == ')':
          count -= 1
        if count < 0:
          remove = False
          break
      if remove:
        start += 1
        end -= 1
      else:
        break
    return expr[start:end+1]

  @staticmethod
  def is_compatible_string(text,model):
    """
    This pre-processing function filters out expressions that 
      are known to be unsupported. 

    Importantly, 'resid' and 'through' are substituted for 'resseq' and colon ':'
      because the corner cases for those are not easily convertable to other selection languages.
      If the substitution does not result in identical selections, it returns a failure. 
    """
    output = text
    fail_reason = None
    failed= False
    if model:
      try:
        sel = model.selection(text)
      except:
        failed = True
        fail_reason = f"invalid phenix selection"
        raise

    if not failed and "*" in text:
      fail_reason = "wildcards not currently supported"
      failed = True
    if not failed and "\\" in text:
      fail_reason = "backslashes not currently supported"
      failed = True
    if not failed and "segid" in text:
      fail_reason = "segid not currently supported"
      failed = True
    if not failed and "icode" in text:
      fail_reason = "insertion code (icode) not currently supported"
      failed = True
    if not failed and "within" in text:
      fail_reason = "within syntax not currently supported"
      failed = True
    # Test through and resid substitution
    if not failed and model:
      s = text
      s2 = s.replace("through",":").replace("resid","resseq")
      model2 = model.select(model.selection(s2))
      xyz1 = model.get_hierarchy().atoms().extract_xyz()
      xyz2 = model2.get_hierarchy().atoms().extract_xyz()
      if not approx_equal(xyz1,xyz1,eps=3):
        fail_reason = """
        through and/or resid could not be substituted for colon and/or resseq. 
        This is a known limitation of the viewer. 
        """
        failed = True
      else:
        output = s2
    passed = not failed
    return passed, fail_reason, output

  def __init__(self, input_str, debug=False, verify=True):
    self.debug = debug
    self.verify = verify  # Do round trip verifications
    # preprocessing
    input_str = self.remove_whitespace_around_colon(input_str)
    input_str = self.remove_top_level_parentheses(input_str)

    self.input_str = input_str
    
    self.keywords = keywords_all

    if self.debug:
      print(f"String to parse: {self.input_str}")
    self.tokens, self.regular_str = self.lexer(self.input_str)
    if self.debug:
      print("\nTokenize:")
      print(self.tokens)
      print()
      print(f"Regularized string: {self.regular_str}")

    # instance vars
    self.ast = None
    self.current_token = None
    self.index = 0
    self._has_parsed = False
    self._tree = None
    
  def lexer(self, input_str):
    keywords = "|".join(self.keywords)
    token_specification = [
      ("ALIAS", r'\b(all|protein|water)\b'),  # Specific keywords (aliases)
      ('ID_QUOTED', r"'[^']*'"),  # Matches fully quoted strings, including empty ones like ''
      ('KEYWORD', rf'\b({keywords})\b'),  # Matches specific keywords defined in self.keywords
      ('LOGICAL', r'\b(and|AND|&|or|OR|not|NOT|~|\|)\b'),  # Ensure logical operators are correctly identified
      ('OPERATOR', r'[<>]=?|==|=|!='),
      ('RANGE', r'\b\d+:\d+\b'),
      ('FLOAT', r'\b\d+\.\d+\b'),
      ('INT', r'\b\d+\b'),  # Match integers before PRIME_ID to avoid misclassification
      ('PRIME_ID', r"[A-Za-z_][A-Za-z0-9_]*'?[A-Za-z0-9_]*"),  # Matches identifiers possibly ending in a prime
      ('ID', r'[A-Za-z_][A-Za-z0-9_]*'),  # Matches regular identifiers
      ('OPEN_PAREN', r'\('),
      ('CLOSE_PAREN', r'\)'),
      ('WILDCARD', r'\*'),
      ('COLON', r':'),
      ('WHITESPACE', r'\s+'),
    ]


    tok_regex = '|'.join('(?P<%s>%s)' % pair for pair in token_specification)
    tokens = []
    for mo in re.finditer(tok_regex, input_str):
      type = mo.lastgroup
      value = mo.group(type)
      token = {"type": type, "value": value}
      tokens.append(token)
      if self.debug:
        print("TOKEN: ",token)

    
    
    if self.verify:
      reconstructed = "".join(token["value"] for token in tokens)
      try:
        assert input_str == reconstructed, f"Input str and reconstructed str do not match\n{input_str}\n{reconstructed}"
      except:
        print(input_str)
        print(tokens)
        print(reconstructed)
        raise
    tokens = [token for token in tokens if token["type"] != "WHITESPACE"]


    regular_str = " ".join(token["value"] for token in tokens)
    return tokens, regular_str

  def consume(self):
    if self.index < len(self.tokens):
      self.current_token = self.tokens[self.index]
      self.index += 1
      if self.debug:
        print(f"Consumed: {self.current_token}")
    else:
      self.current_token = None

  def regularize_ast(self,ast):
    conversion_map = {k:v.upper() for k,v in logic_map_to_phenix.items()}
    def update_values(nested_dict, conversion_dict):
      def recurse(d):
        for key, value in d.items():
          if isinstance(value, dict):
            recurse(value)
          elif key == "value" and value in conversion_dict:
            d[key] = conversion_dict[value]
          elif key == "type" and value in conversion_dict:
            d[key] = conversion_dict[value]
      
      recurse(nested_dict)
      return nested_dict
    return update_values(ast,conversion_map)

  def parse(self):
    if self.debug:
      print("\nParsing...")
    self.consume()
    ast = self.parse_expression()
    ast = self.regularize_ast(ast)
    self.ast = ast

    self._has_parsed = True

  def parse_expression(self):
    node = self.parse_logical_core()
    while self.current_token and self.current_token["type"] == 'LOGICAL':
      logical_type = self.current_token["value"]
      self.consume()
      right_node = self.parse_logical_core()
      if node['type'] == logical_type.upper():
        node['children'].append(right_node)
      else:
        node = {
          'type': logical_type.upper(),
          'children': [node, right_node]
        }
    return node

  def parse_logical_core(self):
    if self.current_token["type"] == 'OPEN_PAREN':
      open_paren_node = {'type': 'OPEN_PAREN', 'value': self.current_token['value']}
      self.consume()
      expr_node = self.parse_expression()
      if self.current_token["type"] != 'CLOSE_PAREN':
        raise Exception("Expected closing parenthesis")
      close_paren_node = {'type': 'CLOSE_PAREN', 'value': self.current_token['value']}
      self.consume()
      return {
        'type': 'PAREN_GROUP',
        'children': [open_paren_node, expr_node, close_paren_node]
      }
    elif self.current_token["type"] == 'LOGICAL' and self.current_token["value"].lower() == 'not':
      self.consume()
      node = self.parse_logical_core()
      return {
        'type': 'NOT',
        'children': [node]
      }
    elif self.current_token["type"] == 'ALIAS':
      return self.parse_alias()
    else:
      return self.parse_keyword_or_value()

  def parse_alias(self):
    token = self.current_token
    self.consume()
    return {
      'type': 'ALIAS',
      'value': token['value'],
      'children': []
    }

  def parse_keyword_or_value(self):
    if self.current_token["type"] == 'KEYWORD':
      return self.parse_keyword()
    else:
      return self.parse_value()

  def parse_keyword(self):
    token = self.current_token
    self.consume()
    if self.current_token and self.current_token["type"] == 'OPERATOR':
      operator = self.current_token
      self.consume()
      if self.current_token["type"] in ["FLOAT", "INT", "ID_QUOTED", "RANGE"]:
        value = self.current_token
        self.consume()
        # if operator["value"] == "=":
        #   operator["value"] = "==" # phenix uses " " and "=" to represent "=="
        return {
          'type': 'COMPARISON',
          'children': [
            {'type': 'KEYWORD', 'value': token},
            {'type': 'OPERATOR', 'value': operator},
            {'type': 'VALUE', 'value': value}
          ]
        }
      else:
        raise Exception("Expected integer, float, or quoted ID after operator")
    elif self.current_token and self.current_token["type"] in ["FLOAT", "INT", "ID_QUOTED", "ID", "PRIME_ID", "RANGE"]:
      value = self.current_token
      self.consume()
      return {
        'type': 'COMPARISON',
        'children': [
          {'type': 'KEYWORD', 'value': token},
          {'type': 'OPERATOR', 'value': {'type': 'OPERATOR', 'value': '=='}},
          {'type': 'VALUE', 'value': value}
        ]
      }
    else:
      raise Exception(f"Expected operator or value after keyword, got: {self.current_token}")

  def parse_value(self):
    token = self.current_token
    if token["type"] in ["ID", "PRIME_ID","ID_QUOTED", "FLOAT", "INT", "RANGE", "WILDCARD"]:
      self.consume()
      return {'type': 'VALUE', 'value': token}
    else:
      raise Exception(f"Expected value token, got {token['type']}")

  @property
  def tree(self):
    if self._tree is None:
      if not self._has_parsed:
        self.parse()
      root = parse_ast(self.ast)
      tree = SelectionTree(root)
      self._tree 
    
    return tree
  @property
  def phenix_string(self):
    return self.tree.phenix_string
  @property
  def pandas_string(self):
    return self.tree.pandas_string
  @property
  def molstar_syntax(self):
    s = self.tree.molstar_syntax
    return s

#########################
### Molstar Generation
#########################

class SelectionTree:
  """
  Container to hold nodes of the AST
  """
  def __init__(self,root):
    self.root = root


  @property
  def phenix_string(self):
    return self.root.phenix_string()

  @property
  def pandas_string(self):
    return self.root.pandas_string()

  @property
  def molstar_syntax(self):
    return self.root.molstar_syntax()

class Node:
  """
  The base class of the AST. Each node should be able to express
    itself in any of the supported selection languages. 
  """
  indent_padding = 2

  ms_operator_map = logic_map_molstar

  @staticmethod
  def ms_test_text(keyword):
    if keyword == "type_symbol":
      return "MS.struct.atomProperty.core.elementSymbol()"
    elif keyword == "id":
      return "MS.struct.atomProperty.macromolecular.id()"
    else:
      return f"MS.ammp('{keyword}')"

  @staticmethod
  def unquote_string(s):
    return s.replace("'","").replace('"','')

  @staticmethod
  def format_value(value, value_type):
    
    if value_type == "INT":
      value = Node.unquote_string(value)  # remove quotes
      return int(value)
    elif value_type == "FLOAT":
      value = Node.unquote_string(value)  # remove quotes
      return float(value)
    else:
      return value.upper() #f"'{value.upper()}'"  # re-quote

  @staticmethod
  def format_keyword(keyword):
    if keyword in attrs_map_to_mmcif:
      return attrs_map_to_mmcif[keyword]
    return keyword

  def __init__(self, node_type):
    self.type = node_type
    self.children = []

  def add_child(self, child):
    if child:  # Ensure that only valid nodes are added
      self.children.append(child)

  def __repr__(self):
    return f"{self.type}({self.children})"

  def molstar_syntax(self, level=0):
    raise NotImplementedError("Subclasses should implement this method")

  def pandas_string(self, level=0):
    raise NotImplementedError("Subclasses should implement this method")

  def phenix_string(self,level=0):
    raise NotImplementedError("Subclasses should implement this method")

  def chimerax_string(self,level=0):
    raise NotImplementedError("Subclasses should implement this method")



class Alias(Node):
  alias_map_pd = {
    "all":"index>=0",
    #"protein":"(" + "|".join([f"label_comp_id == '{aa}'" for aa in aas]) + ")"
  }
  alias_map_ms = {
    "all":"MS.struct.generator.all()",
    #"protein":"MS.core.rel.eq([MS.ammp('entityType'), 'protein'])"

  }
  def __init__(self, value):
    super().__init__('ALIAS')
    self.value = value

  def pandas_string(self, level=0):
    indent = ' ' * (level * Node.indent_padding)
    return f"{indent}{self.value}"

  def pandas_string(self, level=0):
    indent = ' ' * (level * Node.indent_padding)
    return f"{indent}{self.alias_map_pd[self.value]}"

  def molstar_syntax(self, level=0):
    indent = ' ' * (level * Node.indent_padding)
    return f"{indent}{self.alias_map_ms[self.value]}"

  def __repr__(self):
    return f"{self.type}{self.value}"


class Comparison(Node):

  def __init__(self, keyword_token, operator_token, value_token):
    super().__init__('COMPARISON')
    self.keyword_token = keyword_token
    self.operator_token = operator_token
    self.value_token = value_token

  @property
  def keyword(self):
    keyword = self.keyword_token["value"]
    keyword = attrs_map_to_mmcif[keyword] 
    return keyword

  @property
  def keyword_phenix(self):
    keyword = self.keyword_token["value"]
    # Convert to mmcif and back to phenix for standarization
    keyword = attrs_map_to_mmcif[keyword] 
    keyword = attrs_map_to_phenix[keyword]
    return keyword


  @property
  def operator(self):
    operator =  self.operator_token["value"]
    operator = logic_map_to_pandas[operator]
    return operator

  @property
  def operator_phenix(self):
    operator =  self.operator_token["value"]
    operator = logic_map_to_pandas[operator]
    if operator == "==":
      if (self.value_token["type"] == "FLOAT" or
          self.keyword in ["B_iso_or_equiv","occupancy"]
      ):
        return "="
      else:
        return ""
    else:
      return logic_map_to_phenix[operator]

  @property
  def value(self):
    return self.value_token["value"]

  @property
  def value_phenix(self):
    return self.value

  @property
  def value_molstar(self):
    value = self.value
    if "seq_id" in self.keyword: # convert to int for seq
      if self.value_token["type"] != "RANGE":
        if isinstance(self.value,str):
          value = int(self.unquote_string(self.value))
    return value
  

  def phenix_string(self):
    return f"( {self.keyword_phenix} {self.operator_phenix} {self.value_phenix} )"
    

  def pandas_string(self):
    if "seq_id" in self.keyword: # convert to int for seq
      if self.value_type != "RANGE":
        if isinstance(self.value,str):
          self.value = int(self.unquote_string(self.value))
        else:
          self.value = self.value
    if self.value_type == "RANGE":
      value = self.unquote_string(self.value)
      low, high = value.split(":")
      return f"( {low} <= {self.keyword} <= {high} )"
    else:
      return f"( {self.keyword} {self.operator} {self.value} )"

  def molstar_syntax(self, level=0):
    if self.keyword == "id": # convert to int for id for molstar
      if isinstance(self.value,str):
        self.value = int(self.unquote_string(self.value))
      else:
        self.value = self.value
    indent = ' ' * (level * Node.indent_padding)
    if self.value_token["type"] == "RANGE":
      value = self.unquote_string(self.value)
      low, high = value.split(":")
      return f"""
{indent}  MS.core.logic.and([
{indent}    MS.core.rel.gre([{self.ms_test_text(self.keyword)}, {low}]),
{indent}    MS.core.rel.lt([{self.ms_test_text(self.keyword)}, {high}])
{indent}  ])
"""
    else:
      return f"{indent}MS.core.rel.{self.ms_operator_map[self.operator]}([{self.ms_test_text(self.keyword)}, {self.value}])"

  def __repr__(self):
    return f"{self.type}({self.keyword} {self.operator} {self.value})"



class And(Node):
  def __init__(self):
    super().__init__('AND')

  def phenix_string(self):
    children_syntax = ' and '.join(child.phenix_string() for child in self.children)
    return f"{children_syntax}"

  def pandas_string(self):
    children_syntax = ' & '.join(child.pandas_string() for child in self.children)
    return f"{children_syntax}"

  def molstar_syntax(self, level=0):
    indent = ' ' * (level * Node.indent_padding)
    children_syntax = ',\n'.join(child.molstar_syntax(level + 1) for child in self.children)
    return f"{indent}MS.core.logic.and([\n{children_syntax}\n{indent}])"


class Or(Node):
  def __init__(self):
    super().__init__('OR')

  def phenix_string(self):
    children_syntax = ' or '.join(child.phenix_string() for child in self.children)
    return f"{children_syntax}"

  def pandas_string(self):
    children_syntax = ' | '.join(child.pandas_string() for child in self.children)
    return f"{children_syntax}"

  def molstar_syntax(self, level=0):
    indent = ' ' * (level * Node.indent_padding)
    children_syntax = ',\n'.join(child.molstar_syntax(level + 1) for child in self.children)
    return f"{indent}MS.core.logic.or([\n{children_syntax}\n{indent}])"


class Not(Node):
  def __init__(self):
    super().__init__('NOT')

  def phenix_string(self):
    children_syntax = ' '.join(child.phenix_string() for child in self.children)
    return f"not {children_syntax}"

  def pandas_string(self):
    children_syntax = ' '.join(child.pandas_string() for child in self.children)
    return f"~{children_syntax}"


  def molstar_syntax(self, level=0):
    indent = ' ' * (level * Node.indent_padding)
    children_syntax = ',\n'.join(child.molstar_syntax(level + 1) for child in self.children)
    return f"{indent}MS.core.logic.not([\n{children_syntax}\n{indent}])"


class ParenGroup(Node):
  def __init__(self):
    super().__init__('PAREN_GROUP')

  def phenix_string(self):
    children_syntax = ' '.join(child.phenix_string() for child in self.children)
    return f"({children_syntax})"

  def pandas_string(self):
    children_syntax = ' '.join(child.pandas_string() for child in self.children)
    return f"({children_syntax})"

  def molstar_syntax(self, level=0):
    indent = ' ' * (level * Node.indent_padding)
    children_syntax = ',\n'.join(child.molstar_syntax(level + 1) for child in self.children)
    #return f"{indent}MS.struct.modifier.union([\n{children_syntax}\n{indent}])"
    return f"{indent}{children_syntax}{indent}"

def parse_ast(ast):
  if ast['type'] in ['AND', 'OR', 'NOT', 'PAREN_GROUP']:
    if ast['type'] == 'AND':
      node = And()
    elif ast['type'] == 'OR':
      node = Or()
    elif ast['type'] == 'NOT':
      node = Not()
    elif ast['type'] == 'PAREN_GROUP':
      node = ParenGroup()

    for child in ast['children']:
      node.add_child(parse_ast(child))
    return node

  elif ast['type'] == 'COMPARISON':
    keyword_token = ast['children'][0]['value']
    operator_token = ast['children'][1]['value']
    value_token = ast['children'][2]['value']
    return Comparison(keyword_token,operator_token,value_token)
  elif ast["type"] == "ALIAS":
    value = ast["value"]
    return Alias(value)
  return None  # Ensure that invalid nodes are not processed