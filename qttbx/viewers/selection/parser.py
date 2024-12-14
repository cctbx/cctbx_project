import re
from qttbx.viewers.selection.ast import *
from libtbx.test_utils import approx_equal

ms_operator_map = {
  "==": "eq",
  ">=": "gre",
  ">": "gr",
  "<=": "lte",
  "<": "lt",
  "!=": "neq",
}

attrs_map_to_phenix = {
'auth_asym_id': 'chain',
'label_asym_id': 'chain',
'auth_seq_id': 'resseq',
'label_seq_id': 'resseq',
'label_comp_id': 'resname',
'auth_comp_id': 'resname',
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
'auth_atom_id': 'name',
'label_entity_id': 'segid'
}
phenix_map_to_mmcif = {
"chain":'auth_asym_id',    # for chain
"resname":'label_comp_id',   # for resname
"resseq":'auth_seq_id',     # for resseq
"name":'label_atom_id'    # for name
}
attrs_map_to_mmcif = {v:k for k,v in attrs_map_to_phenix.items()}
attrs_map_to_mmcif.update(phenix_map_to_mmcif)
attrs_map_to_mmcif.update(ms_operator_map) 

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

class Token:
  rexp = None # Define in subclass

  def __init__(self, content):
    assert isinstance(content, str), f"Received content of type: {type(content)}, expected string"
    self.content = content

  def __repr__(self):
    name = f"{self.__class__.__name__}"
    return f"{name}('{self.content}')"
      
  def __str__(self):
    return self.content

  @property
  def ms(self):
    return attrs_map_to_mmcif[str(self)]

  def phenix_string(self):
    raise NotImplementedError


class KeywordToken(Token):
  keywords = ["model","chain","resname","resid","resseq","name","bfactor","occupancy","elment","altloc","icode","segid","charge","id"]
  
  @classmethod
  @property
  def rexp(cls):
    rexp = rf'\b({"|".join(cls.keywords)})\b'
    return rexp
  
  def phenix_string(self):
    return self.content

class OperatorToken(Token):
  rexp = r'[<>]=?|==|=|!='
  
  def __init__(self,*args,is_implicit_eq=False,**kwargs):
    super().__init__(*args,**kwargs)
    self.is_implicit_eq = is_implicit_eq

  def phenix_string(self):
    if self.content == "==":
      return "" # equality implicit in phenix
    else:
      return self.content


class ValueToken(Token):
  def phenix_string(self):
    return self.content

class AliasToken(ValueToken):
  rexp = r'\b(all|protein|water)\b'

class IdQuotedToken(ValueToken):
  rexp = r"'[^']*'"

class AndToken(Token):
  rexp = r'\b(and|AND|&)\b'

class OrToken(Token):
  rexp = r'\b(or|OR|\|)\b'

class NotToken(Token):
  rexp = r'\b(not|NOT|~)\b'

class ColonRangeToken(Token):
  rexp = r'\b\d+\s*:\s*\d+\b'

class FloatToken(ValueToken):
  rexp = r'\b\d+\.\d+\b'

class IntToken(ValueToken):
  rexp = r'\b\d+\b'

class PrimeIDToken(ValueToken):
  rexp =  r"[A-Za-z_][A-Za-z0-9_]*'?[A-Za-z0-9_]*"

class IDToken(ValueToken):
  rexp = r'[A-Za-z_][A-Za-z0-9_]*'

class OpenParenToken(Token):
  rexp = r'\('

class CloseParenToken(Token):
  rexp = r'\)'

class WildcardToken(ValueToken):
  rexp = r'\*'

class ColonToken(Token):
  rexp = r':'

class ThroughToken(Token):
  rexp = r'\b(through)\b'

class WithinToken(Token):
  rexp = r'\b(within)\b'

class CommaToken(Token):
  rexp = r'\s*,\s*'

class RangeValueToken(Token):
  # Not directly detected from text, 
  # made later in the parser from components (ColonToken)
  def __init__(self,start_token,range_init_token, stop_token):
    super().__init__(range_init_token.content)
    self.start_token = start_token
    self.range_init_token = range_init_token
    self.stop_token = stop_token

class WhitespaceToken(Token):
  rexp = r'\s+'

class Lexer:
  # The order of these classes is super important, 
  # due to how the regular expressions are processed
  token_classes_default = [
    AliasToken,
    ThroughToken,
    IdQuotedToken,
    KeywordToken,
    AndToken,
    OrToken,
    NotToken,
    OperatorToken,
    ColonRangeToken,
    FloatToken,
    IntToken,
    WithinToken,
    CommaToken,
    PrimeIDToken,
    IDToken,
    OpenParenToken,
    CloseParenToken,
    WildcardToken,
    ColonToken,
    WhitespaceToken,

  ]
  def __init__(self,input_str,token_classes=None,debug=False):
    self.input_str = input_str
    self.token_classes = token_classes
    if not token_classes:
      self.token_classes = self.token_classes_default
    self.debug = debug
    self.tokens = []
    
    tok_dict = {token.__name__:token for token in self.token_classes}
    self.rexp = '|'.join(f"(?P<{token.__name__}>{token.rexp})" for token in self.token_classes)
    if self.debug:
      self._print("Regular exp: ")
      self._print(self.rexp)

    for mo in re.finditer(self.rexp, self.input_str):
      name = mo.lastgroup
      content = mo.group(name)
      token_class = tok_dict[name]
      token = token_class(content)
      self.tokens.append(token)
      if self.debug:        
        self._print("TOKEN: ",token)
    
    self._verify()
  
  def _verify(self):
    reconstructed = "".join(token.content for token in self.tokens)
    assert self.input_str == reconstructed, f"Input str and reconstructed str do not match\n{self.input_str}\n{reconstructed}"

  def _print(self,*args):
    print(*args)

class Parser:
    ignored_token_classes = [WhitespaceToken]

    @classmethod
    def from_string(cls, input_str):
        compatible, error, sel = is_compatible_string(input_str,None)
        lexer = Lexer(sel)
        tokens = lexer.tokens
        return cls(tokens)

    def __init__(self, tokens, debug=False):
        self.tokens = tokens
        self.debug = debug
        self._nodes = {}  # Node storage: nid -> node
        assert all([isinstance(t, Token) for t in tokens]), f"Types provided as tokens: {set([type(t) for t in tokens])}"
        self.input_tokens = [t for t in tokens if t.__class__ not in self.ignored_token_classes]
        self.current_index = 0
        self.root = AstRoot(None)
        self.count_node(self.root)

    def parse_range_value(self, colon_range_token):
        """Parse a range value from a ColonRangeToken"""
        # Split the range token value (e.g., "4:10" into start and end)
        start_str, end_str = str(colon_range_token).split(':')
        
        # Create tokens for start and end values
        start_token = IntToken(start_str)
        end_token = IntToken(end_str)
        
        # Create and return a range node
        range_node = AstRange(start_token, colon_range_token, end_token)
        return self.count_node(range_node)

    def parse_comparison(self):
        """Parses comparison expressions (keyword operator value), now with range support"""
        
        keyword = self.current_token()
        self.advance()
        
        # Handle the operator
        if isinstance(keyword,AliasToken):
          return self.count_node(AstAlias(keyword))
        elif isinstance(self.current_token(), (IntToken, PrimeIDToken, IdQuotedToken, ColonRangeToken)):
            operator = OperatorToken("==", is_implicit_eq=True)
        else:
            operator = self.current_token()
            self.advance()
        
        # Handle the value, checking for range tokens
        value_token = self.current_token()
        self.advance()
        
        if isinstance(value_token, ColonRangeToken):
            value = self.parse_range_value(value_token)
        else:
            value = value_token
        
        node = self.count_node(AstComparison(keyword, operator, value))
        return node

    # Rest of the Parser class methods remain the same
    def _print(self, *args):
        if self.debug:
            print(*args)
            
    def count_node(self, node):
        nid = len(self._nodes)
        node.nid = nid
        self._nodes[nid] = node
        return node
        
    def current_token(self):
        if self.current_index < len(self.input_tokens):
            return self.input_tokens[self.current_index]
        return None
        
    def peek_next(self):
        if self.current_index + 1 < len(self.input_tokens):
            return self.input_tokens[self.current_index + 1]
        return None
        
    def advance(self):
        self.current_index += 1
        return self.current_token()
        
    def parse(self):
        result = self.parse_expression()
        if not isinstance(self.current_token(),(type(None),NotToken)):
            raise SyntaxError(f"Unexpected tokens after expression: {self.current_token()}")
        return result
        
    def parse_expression(self):
        left = self.parse_term()
        
        while self.current_token() and isinstance(self.current_token(), (AndToken, OrToken)):
            operator = self.current_token()
            self.advance()
            
            if isinstance(operator, AndToken):
                # Handle AND with following NOT as a unit
                if isinstance(self.current_token(), NotToken):
                    not_operator = self.current_token()
                    self.advance()
                    right = self.parse_term()
                    not_node = self.count_node(AstNot(not_operator))
                    not_node.left = right
                    and_node = self.count_node(AstAnd(operator))
                    and_node.left = left
                    and_node.right = not_node
                    # Important: treat this whole AND-NOT group as a term
                    left = and_node
                else:
                    right = self.parse_term()
                    node = self.count_node(AstAnd(operator))
                    node.left = left
                    node.right = right
                    left = node
            elif isinstance(operator, OrToken):
                right = self.parse_expression()  # Use parse_expression to handle nested AND-NOT
                node = self.count_node(AstOr(operator))
                node.left = left
                node.right = right
                left = node
        return left
    def parse_term(self):
        token = self.current_token()
        
        if isinstance(token, OpenParenToken):
            self.advance()
            node = self.count_node(AstParenGroup(token))
            node.left = self.parse_expression()
            
            if not isinstance(self.current_token(), CloseParenToken):
                raise SyntaxError("Expected closing parenthesis")
            self.advance()
            return node
            
        elif isinstance(token, KeywordToken):
            return self.parse_comparison()
        
        elif isinstance(token, WithinToken):
            within_operator = token
            self.advance()  # Move past within token
            
            # Expect and consume opening parenthesis
            if not isinstance(self.current_token(), OpenParenToken):
                raise SyntaxError("Expected opening parenthesis after within")
            self.advance()
            
            # Get and validate radius
            radius = self.current_token()
            if not isinstance(radius, (IntToken,FloatToken)):
                raise SyntaxError("Expected numeric radius in within expression")
            self.advance()
            
            # Expect and consume comma
            if not isinstance(self.current_token(), CommaToken):
                raise SyntaxError("Expected comma after radius in within expression")
            self.advance()
            
            # Create within node with radius
            within_node = self.count_node(AstWithin(within_operator, radius))
            
            # Parse the expression after the comma
            within_node.left = self.parse_expression()
            
            # Expect and consume closing parenthesis
            if not isinstance(self.current_token(), CloseParenToken):
                raise SyntaxError("Expected closing parenthesis after within expression")
            self.advance()
            
            return within_node

        elif isinstance(token, AliasToken):
            return self.parse_comparison()
        else:
            raise SyntaxError(f"Unexpected token: {token}")


# Utils
def check_phenix_strings(a,b,model):
  selA = model.selection(a).iselection()
  selB = model.selection(b).iselection()
  assert(selA == selB).all()

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
    # Test through and resid substitution
    if not failed:
      s = text
      s2 = s.replace("through",":").replace("resid","resseq")
      output = s2
      if model:
        model2 = model.select(model.selection(s2))
        xyz1 = model.get_hierarchy().atoms().extract_xyz()
        xyz2 = model2.get_hierarchy().atoms().extract_xyz()
        if not approx_equal(xyz1,xyz1,eps=3):
          fail_reason = """
          through and/or resid could not be substituted for colon and/or resseq. 
          This is a known limitation of the viewer. 
          """
          failed = True
    passed = not failed
    return passed, fail_reason, output
