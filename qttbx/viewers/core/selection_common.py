import re
#from .selection_utils import SelConverterPhenix

import re
from .parameters import params

import re

class PhenixParser:
  """
  Parse a 'common' selection string. ie, Phenix logical syntax with mmcif-like keywords.
  It can build an AST then translate that to other outputs.
  """
  @classmethod
  def from_phenix_string(cls,phenix_string,debug=False,verify=True):
    converter = SelConverterPhenix()
    pandas_str = converter.convert_phenix_to_pandas(phenix_string)
    return cls(pandas_str,debug=debug,verify=verify)



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
    
  def lexer(self,input_str):

    keywords = "|".join(self.keywords)
    token_specification = [
      ('ID_QUOTED', r"'[^']*'"),  # recognize strings enclosed in single quotes
      ('KEYWORD', rf'\b({keywords})\b'),  # recognize keywords
      ('OPEN_PAREN', r'\('),
      ('CLOSE_PAREN', r'\)'),
      ('LOGICAL', r'\b(and|AND|&|or|OR|re.escape(|)|not|NOT|~)\b'),  # recognize logic
      ('OPERATOR', r'[<>]=?|==|!='),
      ('RANGE', r'\b\d+:\d+\b'),  # recognize ranges like "10:20" as single token
      ('FLOAT', r'\b\d+\.\d+\b'),  # recognize floating-point numbers
      ('INT', r'\b\d+\b'),  # recognize integers
      ('ID', r'[A-Za-z0-9_][A-Za-z_0-9]*'),  # recognize identifiers starting with digits
      ('WILDCARD', r'\*'),  # recognize *
      ('COLON', r':'),
      ('WHITESPACE', r'\s+'),  # recognize whitespace
    ]

    tok_regex = '|'.join('(?P<%s>%s)' % pair for pair in token_specification)
    tokens = []
    for mo in re.finditer(tok_regex, input_str):
      type = mo.lastgroup
      value = mo.group(type)
      token = {"type": type, "value": value}
      tokens.append(token)

    reconstructed = ""
    for token in tokens:
      reconstructed += token["value"]
    assert input_str == reconstructed, f"Input str and reconstructed str do not match\n{input_str}\n{reconstructed}"

    for token in tokens:
      if token["type"] != "WHITESPACE":
        assert len(token["value"]) == len(token["value"].strip()), f"Token value contains whitespace: {token}"

    tokens = [token for token in tokens if token["type"] != "WHITESPACE"]

    regular_str = ""
    for i, token in enumerate(tokens):
      regular_str += " "
      regular_str += token["value"]

    return tokens, regular_str

  def consume(self):
    if self.index < len(self.tokens):
      self.current_token = self.tokens[self.index]
      self.index += 1
      if self.debug:
        print(f"Consumed: {self.current_token}")
    else:
      self.current_token = None

  def parse(self):
    if self.debug:
      print("\nParsing...")
    self.consume()
    ast = self.parse_expression()
    #ast = self._postprocess_ast(ast)
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
    elif self.current_token["type"] == 'LOGICAL' and self.current_token["value"] == 'not':
      self.consume()
      node = self.parse_logical_core()
      return {
        'type': 'NOT',
        'children': [node]
      }
    else:
      return self.parse_keyword_or_value()

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
      if self.current_token["type"] in ["FLOAT", "INT", "ID_QUOTED","RANGE"]:
        value = self.current_token
        self.consume()
        return {
          'type': 'COMPARISON',
          'children': [{'type': 'KEYWORD', 'value': token}, {'type': 'OPERATOR', 'value': operator}, {'type': 'VALUE', 'value': value}]
        }
      else:
        raise Exception("Expected integer, float or quoted ID after operator")
    elif self.current_token and self.current_token["type"] in ["FLOAT", "INT", "ID_QUOTED", "ID","RANGE"]:
      # Implicit equality
      value = self.current_token
      self.consume()
      return {
        'type': 'COMPARISON',
        'children': [{'type': 'KEYWORD', 'value': token}, {'type': 'OPERATOR', 'value': {'type': 'OPERATOR', 'value': '=='}}, {'type': 'VALUE', 'value': value}]
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
  def pandas_string(self):
    return self.tree.as_pandas_string()
  @property
  def molstar_syntax(self):
    return self.tree.as_molstar_syntax()

#########################
### Molstar Generation
#########################

class SelectionTree:

  def __init__(self,root):
    self.root = root

  def as_pandas_string(self):
    return self.root.pandas_string()


  def as_molstar_syntax(self):
    return self.root.molstar_syntax()

class Node:
  indent_padding = 2

  ms_operator_map = params.logic_map_molstar

  @staticmethod
  def ms_test_text(keyword):
    if keyword == "type_symbol":
      return "MS.struct.atomProperty.core.elementSymbol()"
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
      return f"'{value}'"  # re-quote

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

  def common_syntax(self, level=0):
    raise NotImplementedError("Subclasses should implement this method")


class Comparison(Node):

  def __init__(self, keyword, operator, value, value_type='STRING'):
    super().__init__('COMPARISON')
    self.keyword = self.format_keyword(keyword)
    self.operator = operator
    self.value_type = value_type
    self.value = self.format_value(value, value_type)

  def pandas_string(self):
    if self.value_type == "RANGE":
      value = SelConverterPhenix.unquote_string(self.value)
      low, high = sorted(value.split(":"))
      return f"( {low} <= {keyword} <= {high} )"
    else:
      return f"( {self.keyword} {self.operator} {self.value} )"

  def common_string(self):
    if self.value_type == "RANGE":
      value = SelConverterPhenix.unquote_string(self.value)
      low, high = sorted(value.split(":"))
      return f"{self.keyword} {self.operator} {low}:{high}"
    else:
      return f"{self.keyword} {self.operator} {self.value}"

  def common_syntax(self, level=0):
    indent = ' ' * (level * Node.indent_padding)
    return f"{indent}{self.common_string()}"


  def molstar_syntax(self, level=0):
    indent = ' ' * (level * Node.indent_padding)
    if self.value_type == "RANGE":
      value = SelConverterPhenix.unquote_string(self.value)
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

  def pandas_string(self):
    children_syntax = ' & '.join(child.pandas_string() for child in self.children)
    return f"{children_syntax}"

  def common_string(self):
    children_syntax = ' and '.join(child.common_string() for child in self.children)
    return f"{children_syntax}"

  def common_syntax(self, level=0):
    indent = ' ' * (level * Node.indent_padding)
    if level == 0:
      children_syntax = '\n'.join(child.common_syntax(level + 1) for child in self.children)
      return f"{children_syntax}\n"
    children_syntax = '\n'.join(child.common_syntax(level + 1) for child in self.children)
    return f"{indent}{self.children[0].common_syntax(level)} \n{indent}and\n{indent}{self.children[1].common_syntax(level + 1)}"


  def molstar_syntax(self, level=0):
    indent = ' ' * (level * Node.indent_padding)
    children_syntax = ',\n'.join(child.molstar_syntax(level + 1) for child in self.children)
    return f"{indent}MS.core.logic.and([\n{children_syntax}\n{indent}])"


class Or(Node):
  def __init__(self):
    super().__init__('OR')

  def pandas_string(self):
    children_syntax = ' | '.join(child.pandas_string() for child in self.children)
    return f"{children_syntax}"


  def common_string(self):
    children_syntax = ' or '.join(child.common_string() for child in self.children)
    return f"{children_syntax}"

  def common_syntax(self, level=0):
    indent = ' ' * (level * Node.indent_padding)
    children_syntax = '\n'.join(child.common_syntax(level + 1) for child in self.children)
    return f"{indent}(\n{children_syntax}\n{indent})"

  def molstar_syntax(self, level=0):
    indent = ' ' * (level * Node.indent_padding)
    children_syntax = ',\n'.join(child.molstar_syntax(level + 1) for child in self.children)
    return f"{indent}MS.core.logic.or([\n{children_syntax}\n{indent}])"


class Not(Node):
  def __init__(self):
    super().__init__('NOT')

  def pandas_string(self):
    children_syntax = ' '.join(child.pandas_string() for child in self.children)
    return f"~{children_syntax}"

  def common_string(self):
    children_syntax = ' '.join(child.common_string() for child in self.children)
    return f"not {children_syntax}"

  def common_syntax(self, level=0):
    indent = ' ' * (level * Node.indent_padding)
    children_syntax = '\n'.join(child.common_syntax(level + 1) for child in self.children)
    return f"{indent}not \n{children_syntax}\n{indent}"

  def molstar_syntax(self, level=0):
    indent = ' ' * (level * Node.indent_padding)
    children_syntax = ',\n'.join(child.molstar_syntax(level + 1) for child in self.children)
    return f"{indent}MS.core.logic.not([\n{children_syntax}\n{indent}])"


class ParenGroup(Node):
  def __init__(self):
    super().__init__('PAREN_GROUP')

  def pandas_string(self):
    children_syntax = ' '.join(child.pandas_string() for child in self.children)
    return f"({children_syntax})"

  def common_string(self):
    children_syntax = ' '.join(child.common_string() for child in self.children)
    return f"({children_syntax})"

  def common_syntax(self, level=0):
    indent = ' ' * (level * Node.indent_padding)
    children_syntax = '\n'.join(child.common_syntax(level + 1) for child in self.children)
    return f"{indent}(\n{children_syntax}\n{indent})"

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

  return None  # Ensure that invalid nodes are not processed


# class CommonSelectionParser:
#   """
#   Parse a 'common' selection string. ie, Phenix logical syntax with mmcif-like keywords.
#   It can build an AST then translate that to other outputs
#   """

#   @staticmethod
#   def remove_whitespace_around_colon(s):
#     return re.sub(r'\s*:\s*', ':', s)
#   @staticmethod
#   def remove_top_level_parentheses(expr):
#     start = 0
#     end = len(expr) - 1
#     while start <= end and expr[start] == '(' and expr[end] == ')':
#       count = 0
#       remove = True
#       for i in range(start + 1, end):
#         if expr[i] == '(':
#           count += 1
#         elif expr[i] == ')':
#           count -= 1
#         if count < 0:
#           remove = False
#           break
#       if remove:
#         start += 1
#         end -= 1
#       else:
#         break
#     return expr[start:end+1]

#   @classmethod
#   def from_phenix_string(cls,phenix_string,debug=False,verify=True):
#     converter = SelConverterPhenix()
#     common_str = converter.convert_phenix_to_common(phenix_string)
#     return cls(common_str,debug=debug,verify=verify)

#   def __init__(self, input_str, debug=False,verify=True):
#     self.debug = debug
#     self.verify=verify # Do round trip verifications
#     # preprocessing
#     input_str = self.remove_whitespace_around_colon(input_str)
#     input_str = self.remove_top_level_parentheses(input_str)

#     self.input_str = input_str

#     if self.debug:
#       print(f"String to parse: {self.input_str}")
#     self.tokens, self.regular_str = self.lexer(self.input_str)
#     if self.debug:
#       print("\nTokenize:")
#       print(self.tokens)
#       print()
#       print(f"Regularized string: {self.regular_str}")


#     # instance vars
#     self.ast = None
#     self.current_token = None
#     self.index = 0
#   @staticmethod
#   def lexer(input_str):
#     token_specification = [
#     ('ID_QUOTED', r"'[^']*'"),  # recognize strings enclosed in single quotes
#     ('KEYWORD', r'\b(asym_id|seq_id|comp_id|atom_id|id|type_symbol|B_iso_or_equiv|occupancy|alt_id)\b'),  # recognize keywords
#     ('OPEN_PAREN', r'\('),
#     ('CLOSE_PAREN', r'\)'),
#     ('LOGICAL', r'\b(and|or|not)\b'),  # recognize logic
#     ('OPERATOR', r'[<>]=?|==|!='),
#     ('RANGE', r'\b\d+:\d+\b'),  # recognize ranges like "10:20" as single token
#     ('FLOAT', r'\b\d+\.\d+\b'),  # recognize floating-point numbers
#     ('INT', r'\b\d+\b'),  # recognize integers
#     ('ID', r'[A-Za-z0-9_][A-Za-z_0-9]*'),  # recognize identifiers starting with digits
#     ('WILDCARD', r'\*'),  # recognize *
#     ('COLON', r':'),
#     ('WHITESPACE', r'\s+'),  # recognize whitespace
#     ]

#     tok_regex = '|'.join('(?P<%s>%s)' % pair for pair in token_specification)
#     tokens = []
#     for mo in re.finditer(tok_regex, input_str):
#       type = mo.lastgroup
#       value = mo.group(type)
#       token = {"type":type,
#                "value": value}
      
#       tokens.append(token)

#     # Regularize

#     # Test round trip
#     #if self.verify:
#     reconstructed = ""
#     for token in tokens:
#       reconstructed+=token["value"]
#     assert input_str == reconstructed, f"Input str and reconstructed str do not match\n{input_str}\n{reconstructed}"


#     # Check for whitespace in values
#     for token in tokens:
#       if token["type"] != "WHITESPACE":
#         assert len(token["value"]) == len(token["value"].strip()), f"Token value contains whitespace: {token}"

#     # Remove whitespace tokens
#     tokens = [token for token in tokens if token["type"] != "WHITESPACE"]

#     # Create standardized whitespace string.
#     regular_str = ""
#     for i,token in enumerate(tokens):
#       regular_str+=" "
#       regular_str+=token["value"]

#     return tokens, regular_str

#   def consume(self):
#     if self.index < len(self.tokens):
#       self.current_token = self.tokens[self.index]
#       self.index += 1
#       if self.debug:
#         print(f"Consumed: {self.current_token}")
#     else:
#       self.current_token = None

#   def parse(self):
#     # Entry point for parsing
#     if self.debug:
#       print("\nParsing...")
#     self.consume()
#     ast = self.parse_expression()
#     if self.debug:
#       print("\n\nAST:")
#       print(json.dumps(ast, indent=2))
#     tokens_out = self.unparse_ast(ast)
#     if self.debug:
#       self.debug_tokens(self.tokens,tokens_out)
#     assert self.tokens==tokens_out, f"Parsing failed round trip check:\n{self.input_str}\n{self.tokens}\n{tokens_out}"
#     reconstructed_regular_str = ""
#     for i,token in enumerate(tokens_out):
#       reconstructed_regular_str+=" "
#       reconstructed_regular_str+=token["value"]
    
#     if self.debug:
#       print("original vs reconstructed: ")
#       print(self.regular_str)
#       print(reconstructed_regular_str)
#     if self.verify:
#       assert self.regular_str == reconstructed_regular_str, "Lexer + Parsing failed round trip check"
#     self.ast = ast

#   def parse_expression(self):
#       # Initialize node with the first logical core expression
#       node = self.parse_logical_core()
#       # Continue building the AST with subsequent logical cores if they exist
#       while self.current_token and self.current_token["type"] == 'LOGICAL':
#         token = self.current_token
#         self.consume()
#         node = self.handle_logical_combinations(node, token)
#       return node

#   def parse_logical_core(self):
#     # If the expression starts with an open parenthesis, parse that first
#     if self.current_token["type"] == 'OPEN_PAREN':
#       self.consume()
#       node = self.parse_expression()  # Note: recursion here
#       if self.current_token["type"] != 'CLOSE_PAREN':
#         raise Exception("Expected closing parenthesis")
#       self.consume()
#       return {
#         'type': 'ParenthesizedExpression',
#         'operand': node
#       }
#     else:  # Otherwise, parse keyword or value
#       return self.parse_keyword_or_value()

#   def handle_logical_combinations(self, node, token):
#     # Handle 'and not' combination or just AND/OR
#     if token['value'] == 'and' and self.current_token and self.current_token['value'] == 'not':
#       not_token = self.current_token
#       self.consume()
#       return {
#         'type': 'BinaryLogical',
#         'left': node,
#         'operator': token,
#         'right': {
#           'type': 'UnaryLogical',
#           'operator': not_token,
#           'operand': self.parse_keyword_or_value()
#         }
#       }
#     else:
#       return {
#         'type': 'BinaryLogical',
#         'left': node,
#         'operator': token,
#         'right': self.parse_expression()
#       }

#   def parse_keyword_or_value(self):
#     # choose between keyword and value for the initial token
#     if self.current_token["type"] == 'KEYWORD':
#       return self.parse_keyword()
#     else:
#       return self.parse_value()

#   def parse_keyword(self):
#     # handle keyword and what follows it
#     token = self.current_token
#     self.consume()
#     if self.current_token["type"] == 'OPERATOR':
#       return self.handle_keyword_with_operator(token)
#     else:
#       return self.handle_keyword_with_value(token)

#   def handle_keyword_with_operator(self, keyword_token):
#     # Handle keyword followed by an operator
#     operator_token = self.current_token
#     self.consume()
#     if self.current_token["type"] in ["FLOAT","INT"]:
#       node = {
#         'type': 'KeywordComparison',
#         'keyword': keyword_token,
#         'operator': operator_token,
#         'value': self.current_token
#       }
#       self.consume()
#       return node
#     else:
#       raise Exception("Expected integer after operator")

#   def handle_keyword_with_value(self, keyword_token):
#     # Handle keyword followed by a value
#     if self.current_token["type"] in ["ID", "ID_QUOTED","FLOAT","INT", "RANGE", "WILDCARD"]:
#       node = {'type': 'Keyword', 'keyword': keyword_token, 'value': self.current_token}
#       self.consume()
#       return node
#     else:
#       raise Exception("Expected value or operator after keyword")

#   def parse_value(self):
#     # Handle value types like ID, INT, etc.
#     token = self.current_token
#     if token["type"] in ["ID","ID_QUOTED", "FLOAT","INT", "RANGE", "WILDCARD"]:
#       node = {'type': 'Value', 'value': token}
#       self.consume()
#       return node
#     else:
#       raise Exception(f"Expected value token, got {token['type']}")

#   def _postprocess_ast(self,ast):
#     # Capitalize certain values, etc
#     self._capitalize_ast(ast)

#   def _capitalize_ast(self,ast):
#     d = ast
#     if isinstance(d, dict):
#       # Check if this node is the target type
#       if d.get('type') in ['ID',"ID_QUOTED"]:
#         # Check if 'value' is present and is a string
#         if 'value' in d and isinstance(d['value'], str):
#           # Convert the 'value' string to uppercase
#           d['value'] = d['value'].upper()

#       # Recurse into each dictionary value
#       for key, value in d.items():
#         if isinstance(value, dict):
#           self._capitalize_ast(value)
#         elif isinstance(value, list):
#           for item in value:
#             if isinstance(item, dict):
#               self._capitalize_ast(item)

  # def consume(self):
  #   if self.index < len(self.tokens):
  #     self.current_token = self.tokens[self.index]
  #     self.index += 1
  #     if self.debug:
  #       print(f"Consumed: {self.current_token}")
  #   else:
  #     self.current_token = None

  # def parse_expression(self):
  #   if self.current_token["type"] == 'OPEN_PAREN':
  #     self.consume()
  #     node = self.parse_logical_core()  # Parse inside parenthesis
  #     if self.current_token["type"] != 'CLOSE_PAREN':
  #       raise Exception("Expected closing parenthesis")
  #     self.consume()
  #     # Add a ParenthesizedExpression node to capture the parenthesis
  #     return {
  #       'type': 'ParenthesizedExpression',
  #       'operand': node
  #     }
  #   else:
  #     return self.parse_logical_core()

  # def parse(self):
  #   if self.debug:
  #     print("\nParsing...")
  #   self.consume()
  #   ast = self.parse_expression()  # Kick-start with parse_expression
    # if self.debug:
    #   print("\n\nAST:")
    #   print(json.dumps(ast, indent=2))
    # tokens_out = self.unparse_ast(ast)
    # self.debug_tokens(self.tokens,tokens_out)
    # assert self.tokens==tokens_out, f"Parsing failed round trip check:\n{self.input_str}\n{self.tokens}\n{tokens_out}"
    # reconstructed_regular_str = " ".join([token["value"] for token in tokens_out])
    # if self.debug:
    #   print("original vs reconstructed: ")
    #   print(self.regular_str)
    #   print(reconstructed_regular_str)
    # assert self.regular_str == reconstructed_regular_str, "Lexer + Parsing failed round trip check"
    # return ast


  # def parse_logical(self):
  #   return self.parse_logical_core()  # Start with the core logic


  # def parse_logical_core(self):
  #   # Allow nested parentheses by calling parse_expression
  #   if self.current_token["type"] == 'OPEN_PAREN':
  #     node = self.parse_expression()
  #   else:
  #     node = self.parse_keyword()  # Otherwise start with keyword or value

  #   while self.current_token and self.current_token["type"] == 'LOGICAL':
  #     token = self.current_token
  #     self.consume()

  #     # Check for the 'and not' combination
  #     if token['value'] == 'and' and self.current_token and self.current_token['value'] == 'not':
  #       not_token = self.current_token
  #       self.consume()
  #       node = {
  #         'type': 'BinaryLogical',
  #         'left': node,
  #         'operator': token,
  #         'right': {
  #           'type': 'UnaryLogical',
  #           'operator': not_token,
  #           'operand': self.parse_keyword()
  #         }
  #       }
  #     else:
  #       # Instead of parse_keyword, call parse_expression
  #       node = {
  #         'type': 'BinaryLogical',
  #         'left': node,
  #         'operator': token,
  #         'right': self.parse_expression()
  #       }

  #   return node


  # def parse_keyword(self):
  #   token = self.current_token
  #   if token["type"] == 'KEYWORD':
  #     self.consume()

  #     if self.current_token["type"] == 'OPERATOR':
  #       operator_token = self.current_token
  #       self.consume()
  #       if self.current_token["type"] == "INT":
  #         node = {
  #           'type': 'KeywordComparison',
  #           'keyword': token,
  #           'operator': operator_token,
  #           'value': self.current_token
  #         }
  #         self.consume()
  #         return node
  #       else:
  #         raise Exception("Expected integer after operator")
  #     elif self.current_token["type"] in ["ID", "INT", "RANGE", "WILDCARD"]:
  #       node = {'type': 'Keyword', 'keyword': token, 'value': self.current_token}
  #       self.consume()
  #       return node
  #     else:
  #       raise Exception("Expected value or operator after keyword")
  #   return self.parse_value()


  # def parse_value(self):
  #   token = self.current_token
  #   if token["type"] in ["ID", "INT", "RANGE","WILDCARD"]:
  #     node = {'type': 'Value', 'value': token}
  #     self.consume()
  #     return node
  #   else:
  #     raise Exception(f"Expected value token, got {token['type']}")

  def unparse_ast(self, node):
    if node['type'] == 'BinaryLogical':
      left_tokens = self.unparse_ast(node['left'])
      right_tokens = self.unparse_ast(node['right'])
      operator_token = [{'type': 'LOGICAL', 'value': node['operator']['value']}]

      # Ensure correct parenthesization
      if node['left']['type'] == 'BinaryLogical' and node['left']['operator']['value'] != node['operator']['value']:
        left_tokens = [{'type': 'OPEN_PAREN', 'value': '('}] + left_tokens + [{'type': 'CLOSE_PAREN', 'value': ')'}]

      if node['right']['type'] == 'BinaryLogical' and node['right']['operator']['value'] != node['operator']['value']:
        right_tokens = [{'type': 'OPEN_PAREN', 'value': '('}] + right_tokens + [{'type': 'CLOSE_PAREN', 'value': ')'}]
      print(f"BinaryLogical: {left_tokens} {operator_token} {right_tokens}\n")
      return left_tokens + operator_token + right_tokens

    elif node['type'] == 'UnaryLogical':
      operator_token = [{'type': 'LOGICAL', 'value': node['operator']['value']}]
      operand_tokens = self.unparse_ast(node['operand'])
      print(f"UnaryLogical: {operator_token} {operand_tokens}\n")
      return operator_token + operand_tokens

    elif node['type'] == 'Keyword':
      keyword_token = [{'type': 'KEYWORD', 'value': node['keyword']['value']}]
      value_token = [{'type': node['value']['type'], 'value': node['value']['value']}]
      print(f"Keyword: {keyword_token} {value_token}\n")
      return keyword_token + value_token

    elif node['type'] == 'KeywordComparison':
      keyword_token = [{'type': 'KEYWORD', 'value': node['keyword']['value']}]
      operator_token = [{'type': 'OPERATOR', 'value': node['operator']['value']}]
      value_token = [{'type': node['value']['type'], 'value': node['value']['value']}]
      print(f"KeywordComparison: {keyword_token} {operator_token} {value_token}\n")
      return keyword_token + operator_token + value_token

    elif node['type'] == 'ParenthesizedExpression':
      operand_tokens = self.unparse_ast(node['operand'])
      print(f"ParenthesizedExpression: ( {operand_tokens} )\n")
      return [{'type': 'OPEN_PAREN', 'value': '('}] + operand_tokens + [{'type': 'CLOSE_PAREN', 'value': ')'}]

    else:
      raise Exception(f"Unrecognized node type: {node['type']}")



  def debug_tokens(self,tokens1, tokens2):
    print("\nOriginal Tokens:")
    for token in tokens1:
      print(token)

    unparsed_tokens = tokens2

    print("\nUnparsed Tokens:")
    for token in unparsed_tokens:
      print(token)

    assert tokens1 == unparsed_tokens, "The original and unparsed tokens do not match!"

  # def to_pandas_query(self):
  #   assert self.ast is not None, "Parse first to get ast"
  #   self._postprocess_ast(self.ast)
  #   query = self._to_pandas_query(self.ast)
    
  #   if self.debug:
  #     print(f"Pandas query:\n{query}")
  #   return query

  # def _to_pandas_query(self,ast):

  #   # handle the BinaryLogical nodes
  #   if ast['type'] == 'BinaryLogical':
  #     left = self._to_pandas_query(ast['left'])
  #     right = self._to_pandas_query(ast['right'])
  #     operator = ast['operator']['value']

  #     if operator == 'and':
  #       operator = '&'
  #     elif operator == 'or':
  #       operator = '|'

  #     return f"({left} {operator} {right})"

  #   # handle the KeywordComparison nodes
  #   elif ast['type'] == 'KeywordComparison':
  #     keyword = ast['keyword']['value']
  #     operator = ast['operator']['value']
  #     value = ast['value']['value']
  #     return f"{keyword} {operator} {value}"

  #   # handle the Keyword nodes
  #   elif ast['type'] == 'Keyword':
  #     keyword = ast['keyword']['value']
  #     value = ast['value']

  #     if value['type'] =='ID':
  #       return f"{keyword} == '{value['value']}'"
  #     elif value['type'] =='ID_QUOTED':
  #       return f"{keyword} == {value['value']}"
  #     elif value['type'] == 'RANGE':
  #       start, end = value['value'].split(':')
  #       return f"{start} <= {keyword} <= {end}"
  #     elif value['type'] == 'INT':
  #       return f"{keyword} == {value['value']}"
  #     elif value['type'] == 'WILDCARD':
  #       return f"{keyword} == {keyword}"

  #   # handle ParenthesizedExpression nodes
  #   elif ast['type'] == 'ParenthesizedExpression':
  #     operand = self._to_pandas_query(ast['operand'])
  #     return f"({operand})"

  #   # handle the UnaryLogical type node
  #   elif ast['type'] == 'UnaryLogical':
  #     operator = ast['operator']['value']
  #     operand = self._to_pandas_query(ast['operand'])

  #     if operator == 'not':
  #       return f"~({operand})"

  #   # handle Value nodes
  #   elif ast['type'] == 'Value':
  #     value = ast['value']['value']
  #     return f"{value}"

  #   # if none of the above, return an empty string
  #   else:
  #     return ''
