"""Methods to convert between a hierarchy object and a forward_compatible_pdb compatible string.
"""
from __future__ import absolute_import, division, print_function

'''
Rationale: Hierarchy object and mmcif representations can contain
  chain ID values with n-characters and residue names with 3 or 5
  characters.  PDB format only allows 2 chars for chain ID and 3 for
  residue names.

Approach: Convert all non-forward_compatible_pdb-compliant chain ID and residue names
  to suitable number of characters and save the conversion information
  as a conversion_info object and as RESNAM records (for residue names)
  and REMARK records (for chain ID) in PDB string representations of
  the hierarchy.

Examples of typical uses:

A. Write a forward_compatible_pdb compatible string with conversion information in REMARK
   and RESNAM records from any hierarchy (ph):
   NOTE: any kw and args for as_pdb_string() can be supplied

  from iotbx.pdb.forward_compatible_pdb_cif_conversion import hierarchy_as_forward_compatible_pdb_string
  forward_compatible_pdb_string =  hierarchy_as_forward_compatible_pdb_string(ph)

B. Read a forward_compatible_pdb compatible string (forward_compatible_pdb_string) with conversion
   information in REMARK/RESNAM records and convert to a hierarchy (
   inverse of A).
   NOTE: same function will read any mmcif string as well.

  from iotbx.pdb.forward_compatible_pdb_cif_conversion import pdb_or_mmcif_string_as_hierarchy
  ph = pdb_or_mmcif_string_as_hierarchy(forward_compatible_pdb_string).hierarchy

C. Get conversion info from any hierarchy (ph):

  from iotbx.pdb.forward_compatible_pdb_cif_conversion import forward_compatible_pdb_cif_conversion
  conversion_info = forward_compatible_pdb_cif_conversion(hierarchy = ph)

D. Get conversion info from unique chain_ids and residue names (
    unique_values_dict):

  from iotbx.pdb.forward_compatible_pdb_cif_conversion import forward_compatible_pdb_cif_conversion
  conversion_info = forward_compatible_pdb_cif_conversion(
    unique_values_dict = unique_values_dict)

D. Get conversion info as REMARK and RESNAM string

  from iotbx.pdb.forward_compatible_pdb_cif_conversion import forward_compatible_pdb_cif_conversion
  remark_hetnam_string = forward_compatible_pdb_cif_conversion(ph).conversion_as_remark_hetnam_string()

E. Convert a forward_compatible_pdb compatible hierarchy to a full hierarchy with
   conversion information in conversion_info. This approach can be
   used to (1) save conversion information from a hierarchy,
   (2) write a forward_compatible_pdb file, (3) do something with the forward_compatible_pdb file that loses
   the header information, (4) read back the forward_compatible_pdb file that does not have
   REMARK records, and (5) restore the original structure in the new
   hierarchy.

  from iotbx.pdb.forward_compatible_pdb_cif_conversion import forward_compatible_pdb_cif_conversion
  from iotbx.pdb.forward_compatible_pdb_cif_conversion import hierarchy_as_forward_compatible_pdb_string
  from iotbx.pdb.forward_compatible_pdb_cif_conversion import pdb_or_mmcif_string_as_hierarchy

  # Get conversion information
  conversion_info = forward_compatible_pdb_cif_conversion(ph)

  # Get a forward_compatible_pdb string with no remarks
  forward_compatible_pdb_string =  hierarchy_as_forward_compatible_pdb_string(ph)
  forward_compatible_pdb_string_no_remarks = remove_remarks(forward_compatible_pdb_string)

  # convert back to hierarchy (this can be a new pdb string obtained
  #  after manipulations of the model but with residue names and chain id
  #  values matching the forward_compatible_pdb_string)
  ph = pdb_or_mmcif_string_as_hierarchy(forward_compatible_pdb_string_no_remarks).hierarchy

  # Apply the conversions to obtain a full representation in ph
  conversion_info.convert_hierarchy_to_full_representation(ph)
  '''

def hierarchy_as_forward_compatible_pdb_string(ph, conversion_info = None, *args, **kw):
  '''Convert a hierarchy into a forward_compatible_pdb compatible string, with any
     conversion information written as REMARK records

    parameters:
      ph: hierarchy object
      conversion_info: optional conversion_info object specifying conversion
      args, kw: any args and kw suitable for the hierarchy
          method ph.as_pdb_string()

    returns:  string
  '''

  if not conversion_info:
    conversion_info = forward_compatible_pdb_cif_conversion(hierarchy = ph)

  if (not conversion_info.conversion_required()):
    return ph.as_pdb_string(*args, **kw)
  else:
    ph_forward_compatible_pdb = ph.deep_copy()
    conversion_info.convert_hierarchy_to_forward_compatible_pdb_representation(ph_forward_compatible_pdb)
    remark_hetnam_string = conversion_info.conversion_as_remark_hetnam_string()
    forward_compatible_pdb_string = ph_forward_compatible_pdb.as_pdb_string(*args, **kw)
    full_string = remark_hetnam_string + forward_compatible_pdb_string
    return full_string

def pdb_or_mmcif_string_as_hierarchy(pdb_or_mmcif_string,
       conversion_info = None):
  '''Convert an mmcif string or a forward_compatible_pdb compatible string into a
      hierarchy object, using any conversion information written as
      REMARK records in the forward_compatible_pdb string, or using any supplied
      conversion information.

    parameters:
      pdb_or_mmcif_string: mmcif string or a forward_compatible_pdb compatible string
      conversion_info: optional forward_compatible_pdb_cif_conversion object to apply

    returns: group_args (hierarchy, pdb_inp, crystal_symmetry, conversion_info)
  '''
  import iotbx.pdb
  from iotbx.pdb.forward_compatible_pdb_cif_conversion import forward_compatible_pdb_cif_conversion
  inp = iotbx.pdb.input(lines=pdb_or_mmcif_string, source_info=None)
  remark_hetnam_string = "\n".join(inp.remark_section())
  hetnam_string = "\n".join(inp.heterogen_section())
  remark_hetnam_string += "\n"+ hetnam_string
  crystal_symmetry = inp.crystal_symmetry()

  if (not conversion_info):
    conversion_info = forward_compatible_pdb_cif_conversion()
    conversion_info.set_conversion_tables_from_remark_hetnam_records(
      remark_hetnam_records = remark_hetnam_string.splitlines())
  assert conversion_info.is_initialized()

  # Get the hierarchy
  ph = inp.construct_hierarchy()
  from libtbx import group_args
  result = group_args(
    group_args_type = 'hierarchy and crystal_symmetry from text string',
    hierarchy = ph,
    pdb_inp = inp,
    crystal_symmetry = crystal_symmetry,
    conversion_info = conversion_info)

  # Determine if this is already in full format
  if forward_compatible_pdb_cif_conversion(ph).conversion_required(): # already set
    assert not conversion_info.conversion_required(), \
      "Cannot apply forward_compatible_pdb conversions to a hierarchy that is not forward_compatible_pdb"
  elif conversion_info.conversion_required(): # convert it
    conversion_info.convert_hierarchy_to_full_representation(ph)
  else: # nothing to convert
    pass
  return result

class forward_compatible_pdb_cif_conversion:
  ''' Class to generate and save forward_compatible_pdb representation of 5-character
    residue names and n-character chain IDs. Used to convert between
    forward_compatible_pdb and mmcif formatting.

    NOTE 1: marked as self._is_initialized when values are available
    NOTE 2: hierarchy object that has been converted to forward_compatible_pdb compatible
    will be marked with the attribute
      self._is_forward_compatible_pdb_representation=True


    To modify these tables to add another field to check:
    1. Add new field to self._keys and self._max_chars_dict
    2. Add new methods like "def _unique_chain_ids_from_hierarchy"
    3. Use these new methods in "def _set_up_conversion_table"
    4. Add code at "Modify hierarchy here to convert to forward_compatible_pdb"
    5. Add code at "Modify hierarchy here to convert from forward_compatible_pdb"
    6. Add code to regression test at iotbx/regression/tst_hierarchy_forward_compatible_pdb.py
    '''

  def __init__(self, hierarchy = None,
     unique_values_dict = None,
     residue_conversion_as_remark = True,
     residue_conversion_as_hetnam = True,
     end_residue_names_with_tilde_if_possible = True,
     ):
    ''' Identify all unique chain_ids and residue names that are not compatible
        with forward_compatible_pdb. Generate dictionary relating original names and
        compatible names and for going backwards.

    parameters:  iotbx.pdb.hierarchy object (required unless unique_values_dict
           is supplied)
        unique_values_dict:  Optional dict with unique values for each key
        residue_conversion_as_remark:   read and write conversion for residue
            name as a REMARK
        residue_conversion_as_hetnam:   read and write conversion for residue
            name as a HETNAM record
        end_residue_names_with_tilde_if_possible:  try to make 3-char residue
                                                    names as 2 chars + "~"

    returns:  None

    '''


    # Fields in hierarchy that are limited in number of characters in forward_compatible_pdb
    self._keys = ['chain_id', 'resname']
    self._max_chars_dict = {'chain_id':2, 'resname':3}
    self._end_with_tilde_dict = {'chain_id':False, 'resname':True}

    if unique_values_dict is not None:
      for key in self._keys:
        assert key in list(unique_values_dict.keys())

    self._remark_keys = ['chain_id']
    self._hetnam_keys = []
    self._residue_conversion_as_remark = residue_conversion_as_remark
    if self._residue_conversion_as_remark:
      self._remark_keys.append('resname')

    self._residue_conversion_as_hetnam = residue_conversion_as_hetnam
    if self._residue_conversion_as_hetnam:
      self._hetnam_keys.append('resname')

    self._is_initialized = False

    if hierarchy is None and unique_values_dict is None:
      self._conversion_table_info_dict = None
      self._conversion_required = None
      return

    # Set up conversion tables
    self._conversion_table_info_dict = {}

    # Flag that indicates if any conversion is necessary
    self._conversion_required = False

    for key in self._keys:
      self._set_up_conversion_table(key, hierarchy,
        unique_values_dict = unique_values_dict)

    self._is_initialized = True

  def is_initialized(self):

    '''Public method to return True if this is initialized
    parameters:  None
    returns: True if initialized
    '''
    return self._is_initialized

  def conversion_required(self):

    '''Public method to return True if conversion for forward_compatible_pdb is necessary
    parameters:  None
    returns: True if conversion is necessary
    '''
    assert self.is_initialized(), "Need to initialize"
    return self._conversion_required

  def conversion_as_remark_hetnam_string(self):
    '''Public method to return a PDB REMARK/HETNAM string representing all the
    conversions that are necessary
    '''

    assert self.is_initialized(), "Need to initialize"

    if not self.conversion_required():
      return ""  # No info needed


    from six.moves import cStringIO as StringIO
    f = StringIO()
    print(
       "REMARK 987 PDB_V3_CONVERSION  CONVERSIONS MADE FOR PDB_V3 COMPATIBILITY",
           file = f)

    # Set up conversion info that goes in REMARK records
    for key in self._remark_keys:
      info =  self._conversion_table_info_dict[key]
      if info:
        for full_text, forward_compatible_pdb_text in zip(
           info.full_representation_list,
           info.forward_compatible_pdb_representation_list,
           ):
          print(
            "REMARK 987 PDB_V3_CONVERSION  %s: %s  PDB_V3_TEXT: %s" %(
              key.upper(),
              full_text,
              forward_compatible_pdb_text),
            file = f)
    print(file = f)

    # Set conversion info that goes in HETNAM records
    """
HETNAM
Overview

This record gives the chemical name of the compound with the given hetID.

Record Format

COLUMNS       DATA  TYPE    FIELD           DEFINITION
----------------------------------------------------------------------------
 1 -  6       Record name   "HETNAM"
 9 - 10       Continuation  continuation    Allows concatenation of multiple records.
12 - 14       LString(3)    hetID           Het identifier, right-justified.
16 - 70       String        text            Chemical name.
    """
    for key in self._hetnam_keys:
      info =  self._conversion_table_info_dict[key]
      if info:
        for full_text, forward_compatible_pdb_text in zip(
           info.full_representation_list,
           info.forward_compatible_pdb_representation_list,
           ):
          print("%6s  %2s %3s %55s%10s" %(
              "HETNAM".ljust(6),
              "".ljust(2),  # continuation chars
              forward_compatible_pdb_text.ljust(3),  # 3-char version
              "PDB_V3_CONVERSION (FULL NAME IN COLS 71-80)".ljust(55),  # any text for 55 chars
              full_text.ljust(10)),  # full version
            file = f)
    print(file = f)


    return f.getvalue()


  def convert_hierarchy_to_forward_compatible_pdb_representation(self, hierarchy):

    '''Public method to convert a hierarchy in place to forward_compatible_pdb compatible
       hierarchy using information in self._conversion_table_info_dict
    parameters: hierarchy (modified in place)
    output: None

    '''

    assert self.is_initialized(), "Need to initialize"
    assert hierarchy is not None, "Need hierarchy for conversion"

    if hasattr(hierarchy, '_is_forward_compatible_pdb_representation') and (
        hierarchy._is_forward_compatible_pdb_representation):
      return # nothing to do because it was already converted

    if not self.conversion_required():
      return # nothing to do because no conversion is necessary

    # Modify hierarchy here to convert to forward_compatible_pdb

    for model in hierarchy.models():
      for chain in model.chains():
        new_id = self.get_forward_compatible_pdb_text_from_full_text(
           key = 'chain_id',
           full_text = chain.id)
        if new_id and new_id != chain.id:
          chain.id = new_id  # Modify chain ID here

        for residue_group in chain.residue_groups():
          for atom_group in residue_group.atom_groups():
            new_resname = self.get_forward_compatible_pdb_text_from_full_text('resname',
                atom_group.resname)
            if new_resname and (new_resname != atom_group.resname):
              atom_group.resname = new_resname  # Modify residue name here

    hierarchy._is_forward_compatible_pdb_representation = True

  def convert_hierarchy_to_full_representation(self, hierarchy):
    '''Public method to convert a hierarchy in place from forward_compatible_pdb compatible
       hierarchy using information in self._conversion_table_info_dict
    parameters: hierarchy (modified in place)
    output: None

    '''
    assert hierarchy is not None, "Need hierarchy for conversion"
    assert self.is_initialized(), "Need to initialize"

    if hasattr(hierarchy, '_is_forward_compatible_pdb_representation') and (
        hierarchy._is_forward_compatible_pdb_representation):
      return # nothing to do because it was already converted

    if not self.conversion_required():
      return # nothing to do because no conversion is necessary

    # Modify hierarchy here to convert from forward_compatible_pdb

    for model in hierarchy.models():
      for chain in model.chains():
        new_id = self.get_full_text_from_forward_compatible_pdb_text(
          key = 'chain_id',
          forward_compatible_pdb_text = chain.id)
        if new_id and new_id != chain.id:
          chain.id = new_id  # Modify chain_id here

        for residue_group in chain.residue_groups():
          for atom_group in residue_group.atom_groups():
            new_resname = self.get_full_text_from_forward_compatible_pdb_text(
              key = 'resname',
              forward_compatible_pdb_text = atom_group.resname)
            if new_resname and (new_resname != atom_group.resname):
              atom_group.resname = new_resname # Modify residue name here

    hierarchy._is_forward_compatible_pdb_representation = False

  def set_conversion_tables_from_remark_hetnam_records(
       self, remark_hetnam_records, add_to_existing = False):
    ''' Public method to set conversion tables based on remarks and hetnam
        records written in standard form as by this class

    parameters:
      remark_hetnam_records:  list of lines, containing REMARK
                        and HETNAM lines with information
                        conversion_as_remark_hetnam_string
      add_to_existing: do not re-initialize if already initialized
    returns: None
    '''


    if not remark_hetnam_records:
      self._is_initialized = True
      return # nothing to do

    self._conversion_required = False

    if (self._is_initialized and add_to_existing):
      pass # keep existing dicts
    else: # usual...initialize
      full_representation_list_dict = {}
      forward_compatible_pdb_representation_list_dict = {}

      for key in self._keys:
        full_representation_list_dict[key] = []
        forward_compatible_pdb_representation_list_dict[key] = []
    self._is_initialized = True

    for line in remark_hetnam_records:
      if not line: continue
      spl = line.split()
      if (spl[0] == "REMARK") and (spl[1] == "987") and (spl[2] == "PDB_V3_CONVERSION"):
        if len(spl) != 7: continue
        key = spl[3].lower()[:-1] # take off ":"
        if not key in self._remark_keys: continue
        full = spl[4]
        forward_compatible_pdb = spl[6]
      elif self._residue_conversion_as_hetnam and (spl[0] == "HETNAM"):
        key = "resname"
        if not key in self._hetnam_keys: continue
        forward_compatible_pdb = line[11:14].strip()
        full = line[69:80].strip()
        if not forward_compatible_pdb: continue
        if not full: continue
      else:
        continue

      if not full in full_representation_list_dict[key]:
        full_representation_list_dict[key].append(full)
        forward_compatible_pdb_representation_list_dict[key].append(forward_compatible_pdb)

      # there was something needing conversion
      self._conversion_required = True

    self._is_initialized = True

    if not self._conversion_required: # nothing to do
      return

    self._conversion_table_info_dict = {}
    from libtbx import group_args
    for key in self._keys:
      self._conversion_table_info_dict[key] = group_args(
        group_args_type = 'conversion tables for %s' %(key),
        full_representation_list = full_representation_list_dict[key],
        forward_compatible_pdb_representation_list =  forward_compatible_pdb_representation_list_dict[key])


  def _set_up_conversion_table(self, key, hierarchy, unique_values_dict = None):
    ''' Private method to set up conversion table from a hierarchy for
        field named by key and put it in self._conversion_table_info_dict[key].
        Also set self._conversion_required if conversion is needed.
        also set self._is_initialized'''

    if unique_values_dict is not None:
      unique_values = unique_values_dict[key]  # use supplied values
    elif key == 'chain_id':
      unique_values = self._unique_chain_ids_from_hierarchy(hierarchy)
    elif key == 'resname':
      unique_values = self._unique_resnames_from_hierarchy(hierarchy)
    else:
      raise "NotImplemented"

    end_with_tilde = self._end_with_tilde_dict[key]

    max_chars = self._max_chars_dict[key]
    allowed_ids, ids_needing_conversion = self._choose_allowed_ids(
        unique_values,
        max_chars = max_chars)

    forward_compatible_pdb_representation_list = self._get_any_forward_compatible_pdb_representation(
        ids_needing_conversion, max_chars, exclude_list = allowed_ids,
        end_with_tilde = end_with_tilde)

    if ids_needing_conversion:
      assert len(ids_needing_conversion) == len(forward_compatible_pdb_representation_list)

    from libtbx import group_args
    self._conversion_table_info_dict[key] = group_args(
      group_args_type = 'conversion tables for %s' %(key),
      full_representation_list = ids_needing_conversion,
      forward_compatible_pdb_representation_list = forward_compatible_pdb_representation_list)

    if forward_compatible_pdb_representation_list:  # there was something needing conversion
      self._conversion_required = True

  def _unique_chain_ids_from_hierarchy(self, hierarchy):
    ''' Private method to identify all unique chain IDs in a hierarchy
    parameters:  hierarchy
    returns:  list of unique chain ids

    '''
    chain_ids = []
    for model in hierarchy.models():
      for chain in model.chains():
        if (not chain.id in chain_ids):
          chain_ids.append(chain.id)
    return chain_ids

  def _unique_resnames_from_hierarchy(self, hierarchy):
    ''' Private method to identify all unique residue names in a hierarchy
    parameters:  hierarchy
    returns:  list of unique residue names

    '''
    resnames = []
    for model in hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          for atom_group in residue_group.atom_groups():
            if (not atom_group.resname in resnames):
              resnames.append(atom_group.resname)
    return resnames

  def _get_any_forward_compatible_pdb_representation(self, ids_needing_conversion,
     max_chars, exclude_list = None,
     end_with_tilde = None):
    '''Private method to try a few ways to generate unique forward_compatible_pdb
      representations for a set of strings.  Order to try:
      1. take first max_chars of each.
      2. if end_with_tilde, then generate max_chars-1 of numbers plus tilde,
      3. generate anything up to max_chars
    parameters:
      ids_needing_conversion:  list of strings to convert
      max_chars:  maximum characters in converted strings
      exclude_list: list of strings not to use as output
      end_with_tilde: try to end strings with a tilde ("~")
    returns:
      forward_compatible_pdb_representation_list: list of converted strings, same order and
        length as ids_needing_conversion
    '''

    if (not ids_needing_conversion):
       return [] # ok with nothing in it

    # Try just taking first n_chars of strings...ok if they are all unique
    forward_compatible_pdb_representation_list = self._get_forward_compatible_pdb_representation(
        ids_needing_conversion, max_chars, exclude_list = exclude_list,
        take_first_n_chars = True)
    if forward_compatible_pdb_representation_list:
      return forward_compatible_pdb_representation_list

    # Generate unique strings for all the ids needing conversion, preventing
    #   duplications of existing ids
    forward_compatible_pdb_representation_list = self._get_forward_compatible_pdb_representation(
        ids_needing_conversion, max_chars, exclude_list = exclude_list,
        end_with_tilde = end_with_tilde)
    if forward_compatible_pdb_representation_list:
      return forward_compatible_pdb_representation_list

    # Failed to get forward_compatible_pdb representation...
    from libtbx.utils import Sorry
    raise Sorry("Unable to generate forward_compatible_pdb representation of %s" %(key))

  def _get_forward_compatible_pdb_representation(self, ids, max_chars,
      exclude_list = None, take_first_n_chars = False,
      end_with_tilde = False):

    '''Private method to try and get forward_compatible_pdb representation of ids that fit in
       max_chars and do not duplicate anything in exclude_list
    parameters:
      ids:  strings to convert
      max_chars:  maximum characters in output
      exclude_list: strings to not include in output
      take_first_n_chars: just take the first max_chars if set
      end_with_tilde: try to end strings with a tilde ("~")

    returns:
      list of converted strings of same order and length as ids, if successful
      otherwise, None
    '''

    forward_compatible_pdb_representation_list = []
    for id in ids:
      if take_first_n_chars:  # Just take the first n_chars
        new_id = id[:max_chars]
        if new_id in exclude_list + forward_compatible_pdb_representation_list:
          return None # cannot do it this way
      else:  # generate a new id
        new_id = self._get_new_unique_id(id, max_chars,
           exclude_list + forward_compatible_pdb_representation_list,
           end_with_tilde = end_with_tilde)
        if not new_id:
          return None # could not do this
      forward_compatible_pdb_representation_list.append(new_id)
    return forward_compatible_pdb_representation_list

  def _get_new_unique_id(self, id, max_chars, exclude_list,
     end_with_tilde):
    ''' Private method to get a unique ID with up to max_chars that is not
    in exclude_list. Start with max_chars and work down and use reverse order
    so as to generally create codes that are unlikely for others to have used.
    Also start with numbers, then numbers and letters, then everything
    '''
    for z in (
        [False, False, True, False],
        [True, False, True, False],
        [True, True, True, False],
        [True, True, True, True],):
      include_upper, include_lower, include_numbers, include_special_chars = z

      if end_with_tilde:
        id = self._get_new_id(max_chars, exclude_list, end_with_tilde = True,
          include_upper = include_upper,
          include_lower = include_lower,
          include_numbers = include_numbers,
          include_special_chars = include_special_chars,)
        if id:
          return id

      for n_chars_inv in range(max_chars):
        n_chars = max_chars - n_chars_inv
        id = self._get_new_id(n_chars, exclude_list,
          include_upper = include_upper,
          include_lower = include_lower,
          include_numbers = include_numbers,
          include_special_chars = include_special_chars,)
        if id:
          return id

  def _get_new_id(self, n_chars, exclude_list, end_with_tilde = None,
      include_upper = True,
      include_lower = True,
      include_numbers = True,
      include_special_chars = True,
       ):
    ''' Private method to get a unique ID with exactly n_chars that is not
    in exclude_list
    '''
    from iotbx.pdb.utils import generate_n_char_string
    x = generate_n_char_string(n_chars = n_chars,
       reverse_order = (not end_with_tilde),
       end_with_tilde = end_with_tilde,
       include_upper = include_upper,
       include_lower = include_lower,
       include_numbers = include_numbers,
       include_special_chars = include_special_chars,
      )
    while 1:
      new_id = x.next()
      if (not new_id):
        return None # failed
      elif (not new_id in exclude_list):
        return new_id

  def _choose_allowed_ids(self, unique_values, max_chars):
    ''' Private method to separate unique_values into those that are and
        are not compatible with forward_compatible_pdb (i.e., have max_chars or fewer)
    '''
    allowed = []
    not_allowed = []
    for u in unique_values:
      if self._is_allowed(u, max_chars):
        allowed.append(u)
      else:
        not_allowed.append(u)
    return allowed, not_allowed

  def _is_allowed(self, u, max_chars):
    ''' Private method to identify whether the string u is or is not
        compatible with forward_compatible_pdb (i.e., has max_chars or fewer)
    '''
    if len(u) <= max_chars:
      return True
    else:
      return False


  def _get_conversion_table_info(self, key):
    ''' Private method to return conversion table info for
        specified key (e.g., chain_id, resname)
    '''

    if not key in self._keys:
      return None
    elif (not self._conversion_required):
      return None
    else:
      return self._conversion_table_info_dict[key]

  def get_full_text_from_forward_compatible_pdb_text(self, key = None, forward_compatible_pdb_text = None):
    '''Public method to return full text from forward_compatible_pdb_text based on
       conversion table

    parameters:
      key: field to convert (e.g., chain_id, resname)
      forward_compatible_pdb_text: text to convert from forward_compatible_pdb to full text
    '''

    assert key is not None
    assert forward_compatible_pdb_text is not None

    conversion_table_info = self._get_conversion_table_info(key)

    if conversion_table_info and (
        forward_compatible_pdb_text in conversion_table_info.forward_compatible_pdb_representation_list):
      index = conversion_table_info.forward_compatible_pdb_representation_list.index(
        forward_compatible_pdb_text)
      full_text = conversion_table_info.full_representation_list[index]
    else:
      full_text = forward_compatible_pdb_text

    return full_text

  def get_forward_compatible_pdb_text_from_full_text(self, key = None, full_text = None):
    '''Public method to return forward_compatible_pdb text from full text based on
       conversion table

    parameters:
      key: field to convert (e.g., chain_id, resname)
      full_text: text to convert to forward_compatible_pdb
    '''

    assert key is not None
    assert full_text is not None

    conversion_table_info = self._get_conversion_table_info(key)
    if conversion_table_info and (
        full_text in conversion_table_info.full_representation_list):
      index = conversion_table_info.full_representation_list.index(
        full_text)
      forward_compatible_pdb_text = conversion_table_info.forward_compatible_pdb_representation_list[index]
    else:
      forward_compatible_pdb_text = full_text

    # Make sure that the resulting text is allowed in forward_compatible_pdb
    assert self._is_allowed(forward_compatible_pdb_text, self._max_chars_dict[key])

    return forward_compatible_pdb_text
