import re


pymol_operator_map = {
    "=": "", 
    "==": "",  # implicit equality
    ">=": ">=",
    ">": ">",
    "<=": "<=",
    "<": "<",
    "!=": "!="
}

pymol_keyword_map = {
  "occupancy":"q",
  "bfactor":"b",
  "charge":"fc.",
  "type_symbol": "element",
  "id": "id",
  "resname": "resn",
  "chain": "chain",
  "resseq": "resi",
  "icode": "alt",
  "altloc": "alt",
  "name": "name"
}
pymol_value_map = {} 

pymol_alias_map = {
    "all": "all",
    "protein": "polymer.protein",
    "backbone": "backbone",
    "sidechain": "sidechain",
    "water": "solvent",
    "ligand": "organic"
}


class AstNode:
  indent_padding = 2 # formatting for ms expressions

  @staticmethod
  def unquote_string(s):
    return s.replace("'","").replace('"','')
    
  @staticmethod
  def pymol_keyword(kw_token):
    if str(kw_token) in pymol_keyword_map:
      return pymol_keyword_map[str(kw_token)]
    return str(kw_token)

  @staticmethod
  def pymol_operator(op_token):
    if str(op_token) in pymol_operator_map:
      return pymol_operator_map[str(op_token)]
    return str(op_token)

  @staticmethod
  def pymol_value(val_token):
    if str(val_token) in pymol_value_map:
      return pymol_value_map[str(val_token)]
    return str(val_token)

  @staticmethod
  def pymol_alias(token):
    if str(token) in pymol_alias_map:
      return pymol_alias_map[str(token)]
    else:
      raise ValueError("Unable to convert alias keyword to pymol: ",str(token))



  def __init__(self, initiator_token, nid=None):
    self.initiator_token = initiator_token
    self.nid = nid
    self._left = None
    self._right = None
    self._parent = None

  @property
  def left(self):
    return self._left
  @left.setter
  def left(self,value):
    self._left = value
    if value:
      value.parent = self
  @property
  def right(self):
    return self._right
  @right.setter
  def right(self,value):
    self._right = value
    if value:
      value.parent = self
  @property
  def parent(self):
    return self._parent

  @parent.setter
  def parent(self, value):
    if self.parent == value:
      return 
    assert self.parent is None, f"Tried to set a parent on node: {self} to be: {value}, but is already: {self.parent}"
    self._parent = value

  def phenix_string(self):
    raise NotImplementedError

  def __repr__(self):
    return f"{self.__class__.__name__}: nid: {self.nid}"

  def add_child(self, child_node):
    print(f"\tAdding child {child_node} to {self}")
    if not self.left:
      self.left = child_node
    elif not self.right:
      self.right = child_node
    else:
      assert False, "left/right children full"
    child_node.parent = self

  @property
  def children(self):
    return [child for child in [self.left,self.right] if child]

  def get_all_children(self):
    # Start with the current object
    objects = [self]
    # Iterate over each child of the current object
    if self.children:
      for child in self.children:
        # Recursively get all objects from the child and extend the list
        if child:
          objects.extend(child.get_all_children())
    return objects

  @staticmethod
  def clean_string(string):
    # replace any sequence of whitespace characters with a single space
    return re.sub(r'\s+', ' ', string).strip()

  def print(self, level=0):
    """Recursively print the AST, indented by the node's depth in the tree."""
    indent = '  ' * level  # Indentation for the current level
    print(f"{indent}{repr(self)}")

    # Recursively print children, if any
    if hasattr(self, 'children') and self.children:
      for child in self.children:
        child.print(level + 1)  # Recursive call for each child


class AstRoot(AstNode):
  def phenix_string(self):
    return "".join([child.phenix_string() for child in self.children])

  def pymol_string(self):
    return "".join([child.pymol_string() for child in self.children])


class AstAnd(AstNode):
    def pymol_string(self):
        if not self.left or not self.right:
            return ""
        left_str = self.left.pymol_string()
        right_str = self.right.pymol_string()
        # Add parentheses around NOT conditions when combined with AND
        if isinstance(self.left, AstNot) or isinstance(self.right, AstNot):
            return f"({left_str} and {right_str})"
        return f"{left_str} and {right_str}"

    def phenix_string(self):
        if not self.left or not self.right:
            return ""
        
        left_str = self.left.phenix_string()
        right_str = self.right.phenix_string()
        return f"{left_str} and {right_str}"



class AstOr(AstNode):

    def pymol_string(self):
        if not self.left or not self.right:
            return ""
        left_str = self.left.pymol_string()
        right_str = self.right.pymol_string()
        
        # Add parentheses only if needed (i.e., if child is not already a ParenGroup)
        # and the content contains multiple terms that need grouping
        if " and " in left_str and not isinstance(self.left, AstParenGroup):
            left_str = f"({left_str})"
        if " and " in right_str and not isinstance(self.right, AstParenGroup):
            right_str = f"({right_str})"
            
        return f"{left_str} or {right_str}"


    def phenix_string(self):
        if not self.left or not self.right:
            return ""
            
        left_str = self.left.phenix_string()
        right_str = self.right.phenix_string()
        return f"{left_str} or {right_str}"

class AstNot(AstNode):
    def __init__(self, initiator_token, nid=None):
        super().__init__(initiator_token, nid)
        self._right = None  # NOT is unary, only needs left child
    
    @property
    def right(self):
        return None
    
    @right.setter
    def right(self, value):
        raise NotImplementedError("NOT node does not support right child")

    def pymol_string(self):
      if not self.left:
          return "not"
      # Don't automatically wrap in parentheses
      return f"not {self.left.pymol_string()}"
    
    def phenix_string(self):
        if not self.left:
            return "not"
        return f"not {self.left.phenix_string()}"


class AstWithin(AstNode):
    def __init__(self, initiator_token, radius_token, nid=None):
        super().__init__(initiator_token, nid)
        self.radius_token = radius_token # Int token
        self._right = None  # within is unary, only needs left child
    
    @property
    def right(self):
        return None
    
    @right.setter
    def right(self, value):
        raise NotImplementedError("Within node does not support right child")
    
    def pymol_string(self):
      if not self.left:
        return ""
      r = float(self.pymol_value(self.radius_token))
      s = f"({self.left.pymol_string()} expand {r})"
      return s

    def phenix_string(self):
        #print("Generating phenix_string for:", self)
        if not self.left:
            return "within"
        return f"within({self.radius_token.content}, {self.left.phenix_string()})"

class AstParenGroup(AstNode):

    def pymol_string(self):
        if not self.left:
            return "( )"
        # If parent is an OR node, we don't need to add extra parens
        # as the OR node will handle necessary grouping
        if hasattr(self, 'parent') and isinstance(self.parent, AstOr):
            return self.left.pymol_string().strip()
        inner_content = self.left.pymol_string().strip()
        return f"({inner_content})"


    def phenix_string(self):
        if not self.left:
            return "( )"
            
        # Only include the first child (left) as ParenGroups shouldn't have right children
        inner_content = self.left.phenix_string().strip()
        return f"({inner_content})"
  
class AstComparison(AstNode):
  
    def __init__(self, keyword, operator, value, nid=None):
        super().__init__(keyword, nid)
        self.keyword = keyword
        self.operator = operator
        self.value = value
    
    @property
    def left(self):
        return None

    @left.setter
    def left(self, value):
        raise NotImplementedError

    @property
    def right(self):
        return None
        
    @right.setter
    def right(self, value):
        raise NotImplementedError
        
    @property
    def children(self):
        return None
      
    def __repr__(self):
        kw_content = ""
        if self.keyword:
          kw_content = str(self.keyword)
        op_content = ""
        if self.operator:
          op_content = str(self.operator)
        v_content = ""
        if self.value:
          v_content = str(self.value)

        return f"{super().__repr__()}: nid: {self.nid}, keyword: {kw_content}, operator: {op_content}, value:{v_content}"

    def pymol_string(self):
        
        # Handle different value types
        keyword = self.pymol_keyword(self.keyword)
        if isinstance(self.value, AstRange):
            return self.value.pymol_string(keyword)
            
        value = self.pymol_value(self.value)
        if isinstance(value, str):
            value = value.replace("'", "").replace('"', '')
            
        # PyMOL uses implicit equality
        operator = self.pymol_operator(self.operator)
        if operator == "":
            return f"{keyword} {value}"
        return f"{keyword} {operator} {value}"


    def phenix_string(self):
        #print("Generating phenix_string for:", self)
        if self.operator.is_implicit_eq:
            return f"{self.keyword.phenix_string()} {self.value.phenix_string()}"
        return f"{self.keyword.phenix_string()} {self.operator.phenix_string()} {self.value.phenix_string()}"


    def add_child(self, child):
        self.parent.add_child(child)

class AstRange(AstNode):
    def __init__(self, start_value, colon_token, end_value, nid=None):
        super().__init__(colon_token, nid)  # Use colon_token as the initiator_token
        self.start_value = start_value
        self.colon_token = colon_token
        self.end_value = end_value
    
    @property
    def content(self):
        """Required for compatibility with printing"""
        return f"{self.start_value.content}:{self.end_value.content}"
    
    def pymol_string(self, keyword=None):
        """Generate PyMOL-compatible range selection"""
        start = self.pymol_value(self.start_value).replace("'", "").replace('"', '')
        end = self.pymol_value(self.end_value).replace("'", "").replace('"', '')
        
        if keyword:
            # PyMOL uses dash for ranges
            return f"{keyword} {start.rstrip()}-{end.lstrip()}"
        return f"{start}-{end}"

    def phenix_string(self):
        """Generate the phenix-compatible string representation"""
        ret =  f"{self.start_value.phenix_string().replace('ColonRangeToken(','')}:{self.end_value.phenix_string()}"
        ret = ret.replace("'","").replace("(","").replace(")","")
        return ret
            
    def __repr__(self):
        return f"{super().__repr__()}: range: {self.start_value.content}:{self.end_value.content}"
    
    @property
    def left(self):
        return None
        
    @left.setter
    def left(self, value):
        raise NotImplementedError("Range nodes don't support children")
        
    @property
    def right(self):
        return None
        
    @right.setter
    def right(self, value):
        raise NotImplementedError("Range nodes don't support children")
        
    @property
    def children(self):
        return None


class AstAlias(AstNode):
  alias_map_ms = {
  "all":"MS.struct.generator.all()",
  "protein":"MS.core.rel.eq([MS.ammp('entityType'), 'protein'])",
      }

  def phenix_string(self):
      return self.initiator_token.content

  def pymol_string(self):
    return self.pymol_alias(self.initiator_token)
      