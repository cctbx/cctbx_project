import re
import numpy as np

from cctbx.array_family import flex
from iotbx.pdb.atom_selection import selection_string_from_selection


__doc__ = """
Definition of terms:

Selection string: Has parenthetic logic, "(chain 'A' and ((resid 19 and (name CB or name OG )) or (resid 20 and (name N or name CG))))"

Selection expression: No parenthesis, "name N or name CG or name CD1 or name CD2"

Selection statement: A single key,value selection: "name N"
                     

"""



# Utility functions

def is_int(s):
    try:
        i = int(s)
        return True
    except:
        return False

      
      
# Simple selection statement translation functions

def translate_model(x,model_number=1):
    assert str(x)==str(model_number)
    return "#"+str(model_number)

def translate_chain(chain_id,model_number=1):
    return "#"+str(model_number)+'/'+chain_id

def translate_residue_name(residue_name,model_number=1):
    return "#"+str(model_number)+":"+residue_name.lower()




def translate_resid_colon(residue_id,model_number=1):
  split = residue_id.split(":")
  assert len(split)==2, "Residue id/number must be a single range, cannot parse: "+residue_id
  resid_start,resid_stop = split
  resid_start,resid_stop = resid_start.strip(),resid_stop.strip()
  # use single value if start == stop
  if resid_start == resid_stop:
    return "#"+str(model_number)+":"+str(resid_start)
  else:
    return "#"+str(model_number)+":"+str(resid_start)+"-"+str(resid_stop)
  

def translate_resid_through(residue_id,model_number=1):
  split = residue_id.split("through")
  assert len(split)==2, "Residue id/number must be a single range, cannot parse: "+residue_id
  resid_start,resid_stop = split
  resid_start,resid_stop = resid_start.strip(),resid_stop.strip()

  # use single value if start == stop
  if resid_start == resid_stop:
    return "#"+str(model_number)+":"+str(resid_start)
  else:
    return "#"+str(model_number)+":"+str(resid_start)+"{-}"+str(resid_stop)+"" # brackets added to denote through range
  


def translate_resseq_colon(residue_id,model_number=1):
  split = residue_id.split(":")
  assert len(split)==2, "Residue id/number must be a single range, cannot parse: "+residue_id
  resid_start,resid_stop = split
  resid_start,resid_stop = resid_start.strip(),resid_stop.strip()

  # use dash with an addition stop:stop statement
  return "((#"+str(model_number)+":"+str(resid_start)+"-"+str(resid_stop)+")|(#"+str(model_number)+":"+str(resid_stop)+"-"+str(resid_stop)+"))"


def translate_resseq_through(residue_id,model_number=1):
  assert False, "Resseq is not supported with through"

def translate_resid(residue_id,model_number=1):
    if ":" in residue_id:
       return translate_resid_colon(residue_id,model_number=model_number)

    elif "through" in residue_id:
        return translate_resid_through(residue_id,model_number=1)
    else:
        return "#"+str(model_number)+":"+str(residue_id)

def translate_resseq(residue_id,model_number=1):
    if ":" in residue_id:
       return translate_resseq_colon(residue_id,model_number=model_number)
    elif "through" in residue_id:
        return translate_resseq_through(residue_id,model_number=1)
    else:
        return "#"+str(model_number)+":"+str(residue_id)+"-"+str(residue_id)

def translate_atom_name(atom_name,model_number=1):
    return "#"+str(model_number)+"@"+atom_name.lower()
    

def translate_element_name(element_name,model_number=1):
    if len(element_name)==2:
        element_name = element_name[0].upper()+element_name[1].lower()
    else:
        assert len(element_name)==1, "Element name must be 1 or two characters, got: "+element_name
    return element_name
        
def translate_attr(key,operator,value):
    return "@@"+key+operator+str(value)


def translate_altloc(value,model_number=1):
    return "@@"+"alt_loc"+"=="+str(value)

def translate_protein(value,model_number=1):
    return "protein"

def translate_nucleotide(value,model_number=1):
    return "nucleic"

def translate_hetatm(value,model_number=1):
    return "::"+"isHet"
  
def translate_water(value,model_number=1):
    return translate_residue_name("HOH",model_number=model_number)
  
def translate_all(value,model_number=1):
    return "all"
  
def translate_not(value,model_number=1):
    return "~"
  
def translate_within(value,model_number=1):
  return "@<__START_WITHIN__"
    
translation_dict = {
    
# compositional keywords, implicit "="
"model":translate_model,
"chain":translate_chain,
"resname":translate_residue_name,
"resid":translate_resid, 
"resseq":translate_resseq, 
"resi":translate_resid, # alias fo resid
"name":translate_atom_name,
"element":translate_element_name,

# attribute keywords, explicit "=",">","<"
"bfactor":translate_attr,
"occupancy":translate_attr,
"altloc":translate_altloc,
    
# single value keywords
"protein":translate_protein,
"nucleotide":translate_nucleotide,
"hetatm":translate_hetatm,
"hetero":translate_hetatm,
"water":translate_water,
"all":translate_all,

#other
"within":translate_within,
"not":translate_not
}

def translate_phenix_expression(selection_expression,model_number=1):
  """
  Translate a Phenix selection expression to Chimera selection. A selection expression
  is defined as anything that fits in a single parenthetical level.
  
  Example:
    resname lys or resname arg and chain A
    This is a valid selection expression
    
    Not:
    (chain B and (resname lys or resname arg and chain A))
    This is not a single expression, it is a full selection string.
  
  Attributes:
    selection_expression (str): A Phenix-style selection expression
    model_number (int): The Chimera model number of the model to appear in the Chimera selection
    
    
  Returns: 
    translated_string (str): The selection expression translated to Chimera-style syntax
    
  """
  translated_expression = ""
  keywords = [" and ", " or "," not ","not "]
  attribute_operators = ["=",">","<",]
  operator_conversion={"and":"&","or":"|","not":"~","":""}
  pattern = "|".join(keywords)
  result = re.split(pattern, selection_expression)
  statements = [r.strip() for r in result]
  operators = [""]+re.findall(pattern,selection_expression)

  for operator,statement in zip(operators,statements):
      operator = operator.strip()

      # Debug
      #print(operator,',',statement)

      # check if attribute statement
      operator_attrs = [o for o in attribute_operators if o in statement]
      assert len(operator_attrs)<=1, f"Multiple attribute operators found in statement: {statement}"
      if len(operator_attrs)==1:
          # an attribute statment
          pattern = "|".join(attribute_operators)
          operators_attr_statement = re.findall(pattern,statement)
          operators_attr_statement = [o for o in operators_attr_statement if len(o.strip())>0]
          assert len(operators_attr_statement)==1, f"Multiple attribute operators found: {statement}"
          operator_attr_statement = operators_attr_statement[0]
          result = re.split(pattern,statement)
          assert len(result)==2, f"Attribute test requires a keyword and value, cannot interpret: {statement}"
          key, value= result
          key, value = key.strip(), value.strip()
          assert key in translation_dict.keys(), f"Statement keyword not recognized: {key}"
          translate_func_attr = translation_dict[key]
          translated_statement = translate_func_attr(key,operator_attr_statement,value)


      else:
          # not an attribute statement
          if len(statement.strip())==0:
              translated_expression+=operator_conversion[operator]
              continue
          statement_split = statement.split()
          if len(statement_split)==1:
              # try single value statements
              key = statement.strip()
              value = ""

          else:
              if len(statement_split)>2:
                  search = statement.split()[0].strip()
                  statement_split = re.split(r"(%s)"%search,statement)
                  statement_split = [s for s in statement_split if len(s)>0]

              assert len(statement_split)==2,f"Unable to parse selection statement: {selection_expression}"
              key,value = statement_split
              key = key.strip()
              value = value.strip("'")
              value = value.strip('"')

          assert key in translation_dict.keys(), f"Statement keyword not recognized: {key}"
          translate_func = translation_dict[key]
          translated_statement = translate_func(value,model_number=model_number)


      # update translated expression
      translated_operator = operator_conversion[operator]
      translated_expression+="".join([translated_operator,translated_statement])

  return translated_expression


def process_within_statements(translated_string,staged_within=[]):
  """
  Attributes:
    translated_string (str): The chimera string with special keywords defining the scope of withins 
    staged_within (list): Information collected about withins during translation (symbol,distance)
                        
  Returns:
    translated_string (str): The selection string translated to Chimera-style syntax
  """
    
  keywords = ["__START_WITHIN__"]
  pattern = "|".join(keywords)
  results = re.split(pattern,translated_string)
  operators = [""]+re.findall(pattern,translated_string)

  within_starts = [(result,operator) for result,operator in zip(results,operators) if len(operator)>0]
  assert len(within_starts)==len(staged_within), "Error parsing within statement"
  
  reformatted_withins = []
  for staged,(result,operator) in zip(staged_within,within_starts):
    assert "__END_WITHIN__" in result, "Error parsing within statement"
    result = result.replace("__END_WITHIN__",")")
    symbol, distance = staged


    # format for chimera
    prepend = " "
    if "&" in symbol:
      prepend = "&"
      symbol = symbol.replace("&","")
    elif "|" in symbol:
      prepend = "| "
      symbol = symbol.replace("|","")

    reformatted = prepend+result[0:-2]+" "+symbol+str(distance)+")"
    reformatted_withins.append(reformatted)
  return "".join(reformatted_withins)



def translate_phenix_selection_string(selection_string,model_number=1):
  """
  Attributes:
    selection_string (str): The full Phenix-style selection string to translate.
    model_number (int): The Chimera model number of the model to appear in the Chimera selection
                        
  Returns:
    translated_string (str): The selection string translated to Chimera-style syntax
  """
  # if there is no parenthetical logic, we don't need this function
  open_count = selection_string.count("(")
  closed_count = selection_string.count(")")
  if open_count==0 and closed_count==0:
    return translate_phenix_expression(selection_string,model_number=model_number)
  


  # preprocessing
  # remove disabled wildcards
  selection_string = selection_string.replace("\*","")
  
  start = None # the character index that starts the content
  stop = None # the index that ends the content
  translated_string = "" # the translated string we build on
  
  open_stack = 0 # this keeps track of the number of open parenthesis
  entered_within = False # This is flipped upon encountering the 'within' keyword
  within_stack_start = None # This is the state of the open stack upon encountering the 'within' keyword
  staged_within = [] # This holds the within information that needs to be incorporated (distance,chimerax symbol)
  
  
  for i, c in enumerate(selection_string): # step through the string for each character 'c'
      
      #debug
      #print(c,entered_within,open_stack,within_stack_start)
          
      last_i = max(i-1,0) # the index preceeding i
      next_i = min(i+1,len(selection_string)-1) # the index following i
      
      last_c = selection_string[last_i] # the character preceding c
      next_c = selection_string[next_i] # the character following c

      if c in ["(",")"] or i== len(selection_string)-1: # a parenthesis encountered

          # trigger the end of a content block
          if last_c not in ["(",")"]:
              # if the preceeding character was not a parenthesis, 
              # we are are likely at the end of a selection expression
              # we call it 'content' 
            
              stop = i # the index that ends the content
              
              content = selection_string[start:stop] # the content string
              #print("content:",content)
              
              # check if we are currently processing a within statement
              if entered_within ==True:
                # if so, then the content will likely be uninterpretable to the translate_expression function
                # Something like: '5, resname lys'
                # We need to remove the 5, store it for later, and pass the rest on.s

                if "," in content:
                  content_split = content.split(",")
                  assert len(content_split)==2,  "Error parsing within statement"
                  
                  # split content into within distance and the actual content
                  d,content = content_split
                  
                  # add the distance to the last staged within list
                  staged_within[-1].append(d)

              # translate the content   
              translated_content = translate_phenix_expression(content,model_number=model_number)
      
              # check if we encountered a 'within' statement
              if "__START_WITHIN__" in translated_content:
                # if the 'within' keyword appeared, the function above will
                #  embed __START_WITHIN__ in the translated content
                
                # sanity check that we are not already handling another within
                assert entered_within == False, "Error parsing within statement"
                entered_within = True # set the flag that we are processing a within
                within_stack_start = open_stack+1 # this is the parenthesis stack state when within encountered
                # if we encountered a within, we need to save the content and not put it immediatetly in the 
                # translated string, because chimera's zone grammar has the subject of the selection first
                staged_within.append([translated_content.replace("__START_WITHIN__","")]) 
                
                # we will put the placeholder into the translated string to replace later
                translated_content = "__START_WITHIN__"
                
              # append the content to the translated string
              translated_string+=translated_content
              

              
              
          
        
      # trigger the start of a content block
      if c in ["(",")"] and next_c not in ["(",")"]:
          # if the character following c is not a parenthesis, but c is,
          # we are at the start of a content block
          start = next_i
      
      # if c is parenthesis, also add it to the translated string
      if i<len(selection_string) and c in ["(",")"]:
        translated_string+=c
      
      # # debug
      # print(c,entered_within,open_stack,within_stack_start)
      # print("")
      
      # finally we need to check if we are currently processing a special within statement
      # If so, then a parenthesis encountered means that this is the end of the within statement
      if c in [")"]:
        if entered_within ==True:

          # But.... if the subject of the within statement also has parentheses, we need to make sure
          # that the parentheses stack has returned to where it was when the within began
          #print("open_stack",open_stack,"within_stack",within_stack_start)
          if open_stack==within_stack_start:
            entered_within = False
            within_stack_start = None
            translated_string+="__END_WITHIN__"
        
      # decide whether to increase to decrease the stack
      if c == "(":   
        open_stack+=1
      elif c==")":
        open_stack-=1


  if "__START_WITHIN__" in translated_string:
    translated_string = process_within_statements(translated_string,staged_within=staged_within)
  return translated_string


def convert_selection_str_to_int(model,string_selection):
  s = model.selection(string_selection)
  int_selection =  s.as_numpy_array().nonzero()[0].tolist()
  return int_selection



def convert_selection_int_to_str(model,int_selection):
  try:
    int_selection = list(int_selection)
  except:
    assert type(int_selection) in [int,float,np.int64], "Pass int_selection as iterable or single number"
    int_selection = [int(int_selection)]
  
  return selection_string_from_selection(pdb_h=model.get_hierarchy(),selection=flex.int(int_selection))