

import json
import re

from .parameters import attrs_map_to_phenix, blanks, core_map_to_mmcif, logic_map_to_phenix
from .parameters import params

class PhenixParser:
  """
  Parse a 'common' selection string. ie, Phenix logical syntax with mmcif-like keywords.
  It can build an AST then translate that to other outputs.
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

  def __init__(self, input_str, debug=False, verify=True):
    self.debug = debug
    self.verify = verify  # Do round trip verifications
    # preprocessing
    input_str = self.remove_whitespace_around_colon(input_str)
    input_str = self.remove_top_level_parentheses(input_str)

    self.input_str = input_str
    
    self.keywords = params.keywords_all

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
      ("ALIAS", r'\b(all|protein|water)\b'),
      ('ID_QUOTED', r"'[^']*'"),
      ('KEYWORD', rf'\b({keywords})\b'),
      ('OPEN_PAREN', r'\('),
      ('CLOSE_PAREN', r'\)'),
      ('LOGICAL', r'(and|AND|&|or|OR|not|NOT|~|\|)'),
      ('OPERATOR', r'[<>]=?|==|!='),
      ('RANGE', r'\b\d+:\d+\b'),
      ('FLOAT', r'\b\d+\.\d+\b'),
      ('INT', r'\b\d+\b'),
      ('ID', r'[A-Za-z0-9_][A-Za-z_0-9]*'),
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

    reconstructed = "".join(token["value"] for token in tokens)
    if self.verify:
      assert input_str == reconstructed, f"Input str and reconstructed str do not match\n{input_str}\n{reconstructed}"

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
    elif self.current_token and self.current_token["type"] in ["FLOAT", "INT", "ID_QUOTED", "ID", "RANGE"]:
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
      raise Exception("Expected operator or value after keyword")

  def parse_value(self):
    token = self.current_token
    if token["type"] in ["ID", "ID_QUOTED", "FLOAT", "INT", "RANGE", "WILDCARD"]:
      self.consume()
      return {'type': 'VALUE', 'value': token}
    else:
      raise Exception(f"Expected value token, got {token['type']}")



  # def _postprocess_ast(self,d):
  #   # Capitalize certain values, etc
    
  #   if isinstance(d, dict):
  #     # Check if this node is the target type
  #     if d.get('type') in ['ID',"ID_QUOTED"]:
  #       # Check if 'value' is present and is a string
  #       if 'value' in d and not isinstance(d['value'], dict):
  #         # Convert the 'value' string to uppercase
  #         value = SelConverterPhenix.unquote_string(d['value'])
  #         d['value'] = f"'{d['value'].upper()}'"
  #     elif d.get('type') == "INT":
  #       if 'value' in d and not isinstance(d['value'], dict):
  #         d['value'] = int(d['value'])
  #     elif d.get('type') == "FLOAT":
  #       if 'value' in d and not isinstance(d['value'], dict):
  #         d['value'] = float(d['value'])
  #     # Recurse into each dictionary value
  #     for key, value in d.items():
  #       if isinstance(value, dict):
  #         self._postprocess_ast(value)
  #       elif isinstance(value, list):
  #         for item in value:
  #           if isinstance(item, dict):
  #             self._postprocess_ast(item)

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
  indent_padding = 2

  ms_operator_map = params.logic_map_molstar

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
    value = Node.unquote_string(value)  # remove quotes
    if value_type == "INT":
      return int(value)
    elif value_type == "FLOAT":
      return float(value)
    else:
      return f"'{value.upper()}'"  # re-quote

  @staticmethod
  def format_keyword(keyword):
    if keyword in params.attrs_map_to_mmcif:
      return params.attrs_map_to_mmcif[keyword]
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

  def __init__(self, keyword, operator, value, value_type='STRING'):
    super().__init__('COMPARISON')
    self.keyword = self.format_keyword(keyword)
    self.operator = operator
    self.value_type = value_type
    self.value = self.format_value(value, value_type)
    
  @staticmethod
  def unquote_string(s):
    return s.replace("'","").replace('"','')
    
  def phenix_string(self):

    if "atom_id" in self.keyword: # remove quotes
      if " " not in self.value:
        self.value = self.unquote_string(self.value)
    if "seq_id" in self.keyword: # convert to int
      if isinstance(self.value,str):
        self.value = int(self.unquote_string(self.value))
      else:
        self.value = self.value
    keyword = attrs_map_to_phenix[self.keyword]
    # Implicit operator
    return f"( {keyword} {self.value} )"

  def pandas_string(self):
    if "seq_id" in self.keyword: # convert to int for seq
      if isinstance(self.value,str):
        self.value = int(self.unquote_string(self.value))
      else:
        self.value = self.value
    if self.value_type == "RANGE":
      value = self.unquote_string(self.value)
      low, high = sorted(value.split(":"))
      return f"( {low} <= {keyword} <= {high} )"
    else:
      return f"( {self.keyword} {self.operator} {self.value} )"

  def molstar_syntax(self, level=0):
    if self.keyword == "id": # convert to int for id for molstar
      if isinstance(self.value,str):
        self.value = int(self.unquote_string(self.value))
      else:
        self.value = self.value
    indent = ' ' * (level * Node.indent_padding)
    if self.value_type == "RANGE":
      value = self.unquote_string(self.value)
      low, high = sorted(value.split(":"))
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
    keyword = ast['children'][0]['value']['value']
    operator = ast['children'][1]['value']['value']
    value = ast['children'][2]['value']['value']
    value_type = ast['children'][2]['value']['type']
    return Comparison(keyword, operator, value,value_type)
  elif ast["type"] == "ALIAS":
    value = ast["value"]
    return Alias(value)
  return None  # Ensure that invalid nodes are not processed