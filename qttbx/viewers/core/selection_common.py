import re
import json


class CommonSelectionParser:
  """
  Parse a 'common' selection string. ie, Phenix logical syntax with mmcif-like keywords.
  It can build an AST then translate that to other outputs
  """

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


  def __init__(self, input_str, debug=False,verify=True):
    self.debug = debug
    self.verify=verify # Do round trip verifications
    # preprocessing
    input_str = self.remove_whitespace_around_colon(input_str)
    input_str = self.remove_top_level_parentheses(input_str)

    self.input_str = input_str

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

  def lexer(self,input_str):
    token_specification = [
      ('ID_QUOTED', r"'[^']*'"),  # recognize strings enclosed in single quotes
      ('KEYWORD', r'\b(asym_id|seq_id|comp_id|atom_id|id|type_symbol|B_iso_or_equiv|occupancy|alt_id)\b'),  # recognize keywords
      ('OPEN_PAREN', r'\('),
      ('CLOSE_PAREN', r'\)'),
      ('LOGICAL', r'\b(and|or|not)\b'),  # recognize logic
      ('OPERATOR', r'[<>]=?|==|!='),
      ('RANGE', r'\b\d+:\d+\b'),  # recognize ranges like "10:20" as single token
      ('INT', r'\b\d+\b'),  # recognize integers
      ('ID', r'[A-Za-z0-9_][A-Za-z_0-9]*'),  # modified to recognize identifiers starting with digits
      ('WILDCARD', r'\*'),  # recognize *
      ('COLON', r':'),
      ('WHITESPACE', r'\s+'),
    ]
    tok_regex = '|'.join('(?P<%s>%s)' % pair for pair in token_specification)
    tokens = []
    for mo in re.finditer(tok_regex, input_str):
      type = mo.lastgroup
      value = mo.group(type)
      token = {"type":type,
               "value": value}
      
      tokens.append(token)

    # Regularize

    # Test round trip
    if self.verify:
      reconstructed = ""
      for token in tokens:
        reconstructed+=token["value"]
      assert input_str == reconstructed, f"Input str and reconstructed str do not match\n{input_str}\n{reconstructed}"


    # Check for whitespace in values
    for token in tokens:
      if token["type"] != "WHITESPACE":
        assert len(token["value"]) == len(token["value"].strip()), f"Token value contains whitespace: {token}"

    # Remove whitespace tokens
    tokens = [token for token in tokens if token["type"] != "WHITESPACE"]

    # Create regularized string
    regular_str = ""
    for i,token in enumerate(tokens):
      regular_str+=" "
      regular_str+=token["value"]

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
    # Entry point for parsing
    if self.debug:
      print("\nParsing...")
    self.consume()
    ast = self.parse_expression()
    if self.debug:
      print("\n\nAST:")
      print(json.dumps(ast, indent=2))
    tokens_out = self.unparse_ast(ast)
    if self.debug:
      self.debug_tokens(self.tokens,tokens_out)
    assert self.tokens==tokens_out, f"Parsing failed round trip check:\n{self.input_str}\n{self.tokens}\n{tokens_out}"
    reconstructed_regular_str = ""
    for i,token in enumerate(tokens_out):
      reconstructed_regular_str+=" "
      reconstructed_regular_str+=token["value"]
    
    if self.debug:
      print("original vs reconstructed: ")
      print(self.regular_str)
      print(reconstructed_regular_str)
    if self.verify:
      assert self.regular_str == reconstructed_regular_str, "Lexer + Parsing failed round trip check"
    self.ast = ast

  def parse_expression(self):
      # Initialize node with the first logical core expression
      node = self.parse_logical_core()
      # Continue building the AST with subsequent logical cores if they exist
      while self.current_token and self.current_token["type"] == 'LOGICAL':
        token = self.current_token
        self.consume()
        node = self.handle_logical_combinations(node, token)
      return node

  def parse_logical_core(self):
    # If the expression starts with an open parenthesis, parse that first
    if self.current_token["type"] == 'OPEN_PAREN':
      self.consume()
      node = self.parse_expression()  # Note: recursion here
      if self.current_token["type"] != 'CLOSE_PAREN':
        raise Exception("Expected closing parenthesis")
      self.consume()
      return {
        'type': 'ParenthesizedExpression',
        'operand': node
      }
    else:  # Otherwise, parse keyword or value
      return self.parse_keyword_or_value()

  def handle_logical_combinations(self, node, token):
    # Handle 'and not' combination or just AND/OR
    if token['value'] == 'and' and self.current_token and self.current_token['value'] == 'not':
      not_token = self.current_token
      self.consume()
      return {
        'type': 'BinaryLogical',
        'left': node,
        'operator': token,
        'right': {
          'type': 'UnaryLogical',
          'operator': not_token,
          'operand': self.parse_keyword_or_value()
        }
      }
    else:
      return {
        'type': 'BinaryLogical',
        'left': node,
        'operator': token,
        'right': self.parse_expression()
      }

  def parse_keyword_or_value(self):
    # choose between keyword and value for the initial token
    if self.current_token["type"] == 'KEYWORD':
      return self.parse_keyword()
    else:
      return self.parse_value()

  def parse_keyword(self):
    # handle keyword and what follows it
    token = self.current_token
    self.consume()
    if self.current_token["type"] == 'OPERATOR':
      return self.handle_keyword_with_operator(token)
    else:
      return self.handle_keyword_with_value(token)

  def handle_keyword_with_operator(self, keyword_token):
    # Handle keyword followed by an operator
    operator_token = self.current_token
    self.consume()
    if self.current_token["type"] == "INT":
      node = {
        'type': 'KeywordComparison',
        'keyword': keyword_token,
        'operator': operator_token,
        'value': self.current_token
      }
      self.consume()
      return node
    else:
      raise Exception("Expected integer after operator")

  def handle_keyword_with_value(self, keyword_token):
    # Handle keyword followed by a value
    if self.current_token["type"] in ["ID", "ID_QUOTED","INT", "RANGE", "WILDCARD"]:
      node = {'type': 'Keyword', 'keyword': keyword_token, 'value': self.current_token}
      self.consume()
      return node
    else:
      raise Exception("Expected value or operator after keyword")

  def parse_value(self):
    # Handle value types like ID, INT, etc.
    token = self.current_token
    if token["type"] in ["ID","ID_QUOTED", "INT", "RANGE", "WILDCARD"]:
      node = {'type': 'Value', 'value': token}
      self.consume()
      return node
    else:
      raise Exception(f"Expected value token, got {token['type']}")

  def _postprocess_ast(self,ast):
    # Capitalize certain values, etc
    self._capitalize_ast(ast)

  def _capitalize_ast(self,ast):
    d = ast
    if isinstance(d, dict):
      # Check if this node is the target type
      if d.get('type') in ['ID',"ID_QUOTED"]:
        # Check if 'value' is present and is a string
        if 'value' in d and isinstance(d['value'], str):
          # Convert the 'value' string to uppercase
          d['value'] = d['value'].upper()

      # Recurse into each dictionary value
      for key, value in d.items():
        if isinstance(value, dict):
          self._capitalize_ast(value)
        elif isinstance(value, list):
          for item in value:
            if isinstance(item, dict):
              self._capitalize_ast(item)


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
      if node['left']['type'] == 'BinaryLogical' and node['left']['operator']['value'] != node['operator']['value']:
        left_tokens = [{'type': 'OPEN_PAREN', 'value': '('}] + left_tokens + [{'type': 'CLOSE_PAREN', 'value': ')'}]

      if node['right']['type'] == 'BinaryLogical' and node['right']['operator']['value'] != node['operator']['value']:
        right_tokens = [{'type': 'OPEN_PAREN', 'value': '('}] + right_tokens + [{'type': 'CLOSE_PAREN', 'value': ')'}]
      return left_tokens + operator_token + right_tokens

    elif node['type'] == 'UnaryLogical':
      operator_token = [{'type': 'LOGICAL', 'value': node['operator']['value']}]
      operand_tokens = self.unparse_ast(node['operand'])
      return operator_token + operand_tokens

    elif node['type'] == 'Keyword':
      keyword_token = [{'type': 'KEYWORD', 'value': node['keyword']['value']}]
      value_token = [{'type': node['value']['type'], 'value': node['value']['value']}]
      return keyword_token + value_token

    elif node['type'] == 'KeywordComparison':
      keyword_token = [{'type': 'KEYWORD', 'value': node['keyword']['value']}]
      operator_token = [{'type': 'OPERATOR', 'value': node['operator']['value']}]
      value_token = [{'type': node['value']['type'], 'value': node['value']['value']}]
      return keyword_token + operator_token + value_token

    elif node['type'] == 'ParenthesizedExpression':
      operand_tokens = self.unparse_ast(node['operand'])
      return [{'type': 'OPEN_PAREN', 'value': '('}] + operand_tokens + [{'type': 'CLOSE_PAREN', 'value': ')'}]


    else:
      raise Exception(f"Unrecognized node type: {node['type']}")


  def debug_tokens(self,tokens1, tokens2):
    N = max([len(tokens1),len(tokens2)])
    for i in range(N):
      if i< len(tokens1):
        t1 = tokens1[i]
      else:
        t1 = ""
      if i< len(tokens2):
        t2 = tokens2[i]
      else:
        t2 = ""
      print(t1,"\t",t2)

  def to_pandas_query(self):
    assert self.ast is not None, "Parse first to get ast"
    self._postprocess_ast(self.ast)
    query = self._to_pandas_query(self.ast)
    
    if self.debug:
      print(f"Pandas query:\n{query}")
    return query

  def _to_pandas_query(self,ast):

    # handle the BinaryLogical nodes
    if ast['type'] == 'BinaryLogical':
      left = self._to_pandas_query(ast['left'])
      right = self._to_pandas_query(ast['right'])
      operator = ast['operator']['value']

      if operator == 'and':
        operator = '&'
      elif operator == 'or':
        operator = '|'

      return f"({left} {operator} {right})"

    # handle the KeywordComparison nodes
    elif ast['type'] == 'KeywordComparison':
      keyword = ast['keyword']['value']
      operator = ast['operator']['value']
      value = ast['value']['value']
      return f"{keyword} {operator} {value}"

    # handle the Keyword nodes
    elif ast['type'] == 'Keyword':
      keyword = ast['keyword']['value']
      value = ast['value']

      if value['type'] =='ID':
        return f"{keyword} == '{value['value']}'"
      elif value['type'] =='ID_QUOTED':
        return f"{keyword} == {value['value']}"
      elif value['type'] == 'RANGE':
        start, end = value['value'].split(':')
        return f"{start} <= {keyword} <= {end}"
      elif value['type'] == 'INT':
        return f"{keyword} == {value['value']}"
      elif value['type'] == 'WILDCARD':
        return f"{keyword} == {keyword}"

    # handle ParenthesizedExpression nodes
    elif ast['type'] == 'ParenthesizedExpression':
      operand = self._to_pandas_query(ast['operand'])
      return f"({operand})"

    # handle the UnaryLogical type node
    elif ast['type'] == 'UnaryLogical':
      operator = ast['operator']['value']
      operand = self._to_pandas_query(ast['operand'])

      if operator == 'not':
        return f"~({operand})"

    # handle Value nodes
    elif ast['type'] == 'Value':
      value = ast['value']['value']
      return f"{value}"

    # if none of the above, return an empty string
    else:
      return ''
