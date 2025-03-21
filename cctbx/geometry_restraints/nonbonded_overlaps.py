from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
from scitbx.array_family import flex
from libtbx import easy_run
import iotbx.pdb
import math
import sys
import six
from six.moves import range, zip

class nonbonded_overlaps_results(object):
  """ Container for non-bonded overlaps results """

  def __init__(self):
    """
    - nonbonded_overlaps is number of overlaps based on the non-bonded
      proxies pairs
    - Overlapping proxies is a list containing the information on the
      overlapping atoms
      in the format:([pdb labels],i_seq,j_seq,model,vdw_distance,sym_op_j,rt_mx)

    The following non-bonded overlaps are calculated:
      nb_overlaps_due_to_sym_op: (calculated with the complete model)
      nb_overlaps_macro_molecule: (protein, DNA and RNA)
        excluding symmetry related overlaps
      nb_overlaps_all: (calculated with the complete model)

    For each of the overlaps a list of overlapping proxies (atoms) is provided:
      nb_overlaps_proxies_due_to_sym_op
      nb_overlaps_proxies_macro_molecule
      nb_overlaps_proxies_all

    The following non-bonded overlaps are calculated:
      normalized_nbo_sym: (calculated with the complete model)
      normalized_nbo_macro_molecule: (protein, DNA and RNA)
        excluding symmetry related overlaps
      normalized_nbo_all: (calculated with the complete model)

    """
    self.nb_overlaps_due_to_sym_op = 0
    self.nb_overlaps_macro_molecule = 0
    self.nb_overlaps_all = 0
    #
    self.nb_overlaps_proxies_due_to_sym_op = []
    self.nb_overlaps_proxies_macro_molecule = []
    self.nb_overlaps_proxies_all = []
    #
    self.normalized_nbo_sym = 0
    self.normalized_nbo_macro_molecule = 0
    self.normalized_nbo_all = 0

class compute(object):
  """
  Compute non-bonded overlap using geometry restraints manager.

  @author: Youval Dar (LBL 2013)
  """

  def __init__(
    self,
    nonbonded_list,
    hd_sel,
    full_connectivity_table,
    connectivity_table_2,
    sites_cart):
    """
    Arguments:
    nonbonded_list: a list with items in the following format
                    ([pdb labels],i_seq,j_seq,model,vdw_distance,sym_op_j,rt_mx)
                    i_seq,j_seq: position of residues in the pdb file
                    model: The pdb or other model distance
                    vdw_distance: Van Der Waals distance
                    sym_op_j: is this a product of a symmetry operation
                    rt_mx: Rotation matrix for symmetry operation
    hd_sel: hd_sel[i] returns True of False, indicating whether an atom i
            is a Hydrogen or not
    full_connectivity_table: full_connectivity_table[i] is a dictionary
                            containing a list of all atoms connected to atom i
    connectivity_table_2: connectivity_table[i] is a dictionary containing a
                         list of all atoms separated by three bonds from atom i
                         (1 - 4 interaction)
    sites_cart: sites_cart[i] tuple containing the x,y,z coordinates of atom i
    """
    self.nonbonded_list = nonbonded_list
    self.hd_sel = hd_sel
    self.full_connectivity_table = full_connectivity_table
    self.connectivity_table_2 = connectivity_table_2
    self.sites_cart = sites_cart
    if nonbonded_list != []:
      try:
        nbo_list = self.nb_overlaps_list()
      except TypeError as e:
        # When proxies_info_nonbonded are not available
        if e == "vec3_double' object is not callable":
          raise Sorry(e)
        else:
          nbo_list = [[],[]]
          print(e)
      except Sorry as e:
        raise Sorry(e)
      except Exception as e:
        m='Failed processing proxies_info_nonbonded in nb_overlaps_list()'
        raise Sorry(m)
      # Collect overlaps information
      self.nb_overlaps_proxies_due_to_sym_op = nbo_list[0]
      self.nb_overlaps_non_sym_overlaps = nbo_list[1]
      self.nb_overlaps_proxies_all = nbo_list[0] + nbo_list[1]
      # overlap is the of steric overlaps (>0.4A)
      self.n_atoms = len(self.sites_cart)
      nbo_sym = len(nbo_list[0])
      nbo_non_sym = len(nbo_list[1])
      nbo_all_overlaps = nbo_sym + nbo_non_sym
      #
      self.nb_overlaps_due_to_sym_op = nbo_sym
      self.nb_overlaps_non_sym = nbo_non_sym
      self.nb_overlaps_all = nbo_all_overlaps
      # compute clashscores for testing
      # number of steric overlaps (>0.4A) per 1000 atoms
      clashscore_sym = nbo_sym*1000/self.n_atoms
      clashscore_non_sym = nbo_non_sym*1000/self.n_atoms
      clashscore_all_clashes = clashscore_sym + clashscore_non_sym
      self.normalized_nbo_sym = clashscore_sym
      self.cctbx_clashscore_non_sym = clashscore_non_sym
      self.normalized_nbo_all = clashscore_all_clashes
    else:
      self.nb_overlaps_due_to_sym_op = 0
      self.nb_overlaps_non_sym = 0
      self.nb_overlaps_all = 0
      self.normalized_nbo_sym = 0
      self.cctbx_clashscore_non_sym = 0
      self.normalized_nbo_all = 0
      #
      self.nb_overlaps_proxies_due_to_sym_op = []
      self.nb_overlaps_non_sym_overlaps = []
      self.nb_overlaps_proxies_all = []

  @staticmethod
  def is_1_5_interaction(i_seq,j_seq,hd_sel,full_connectivity_table):
    """(int,int,bool array,list of lists of int) -> bool
    Check if we have 1-5 interaction between two hydrogen and heavy atom

    Args:
      i_seq,j_seq:  are the number of the atoms we are checking in the pdb table
      hd_sel: hd_sel[i] returns True of False, indicating whether an atom i is
        a Hydrogen or not
      full_connectivity_table: full_connectivity_table[i] is a dictionary
        constraining a list of all atoms connected to atom

    Returns:
      True if we have 1-5 hydrogen and heavy atom interaction
    """
    # check if we have hydrogen - heavy atom interaction
    xor = lambda a,b: (a or b) and not (a and b)
    if xor(hd_sel[i_seq],hd_sel[j_seq]):
      # starting with hydrogen will make process shorter
      if not hd_sel[i_seq]:
        i_seq,j_seq = j_seq,i_seq
      # build connection table of i_seq, 4 steps deep
      atoms_numbers = dict([(i_seq, 0)]) # i_seq is in zero distance
      used_connections = {i_seq}
      new_connections = {i_seq}
      for i in range(2,6):
        connections = set()
        for key in new_connections:
          # add all new connections for the current step
          connections = connections.union(set(full_connectivity_table[key]))
        # Remove the connection that were already used
        new_connections = connections - used_connections
        # Add the new connection to the used once
        used_connections = used_connections.union(connections)
        # Add the new atoms with their distance from key
        for new_atom in new_connections:
          atoms_numbers[new_atom] = i
      # return true if j_seq in the is 1-5 connection
      return (j_seq in atoms_numbers) and (atoms_numbers[j_seq] == 5)
    else:
      return False

  def get_nbo_keys(self,record):
    """(str) -> str,str,str

    Collect atoms keys and coordinates

    Args:
      record : an overlap record, for example:
    (['pdb=" HA  LEU A  38 "','pdb="HD23 LEU A  38 "'],523,532,1.81,2.44,'',None)

    Returns:
      nbo_vec: a unique key for an overlap. for example:
        '0.144,0.323,1.776,0.154,0.327,1.786'
      key1,key2: are the atoms keys. for example
        ' HA  LEU A  38 ::120::' or ' CB ASN A  55 ::52::sym.op.'
        "atom1 label::atom1 number::sym.op""
    """
    # get atoms keys
    record = list(record)
    key = '{0}::{1}::{2}'
    key1 = key.format(record[0][0][5:-1],record[1],record[5])
    key2 = key.format(record[0][1][5:-1],record[2],record[5])
    # get coordinates
    atomi_xyz = self.sites_cart[record[1]]
    atomj_xyz = self.sites_cart[record[2]]
    # make tupe containing both atom coordinates
    vec = atomi_xyz + atomj_xyz
    # creat a string from vec
    nbo_vec = '{0:.3f},{1:.3f},{2:.3f},{3:.3f},{4:.3f},{5:.3f}'.format(*vec)
    # make key1 always smaller than key2
    if key1 > key2:
      key1,key2 = key2,key1
    return nbo_vec,key1,key2

  @staticmethod
  def cos_vec(u,v):
    """(tuple,tuple) -> float

    Calculate the cosine used to evaluate whether the atoms
    colliding are inline or not

    Args:
      u,v: lists containing x1,y1,z1,x2,y2,z2 coordinates
    Returns:
      cos_angle: the cosine of the angle between center of the common atom
        and the mid point between the other two atoms.
    """
    # Calculate dot product
    u_dot_v = lambda u,v: (u[0]*v[0]+u[1]*v[1]+u[2]*v[2])
    # Calculate mid point between to atoms
    u_mid_v = lambda u,v: [(x+y)/2 for (x,y) in zip(u,v)]
    # find the length of a vector
    u_len = lambda u: math.sqrt(u[0]**2 + u[1]**2 + u[2]**2)
    # Move vector to the origin
    make_vec = lambda u1,u2: [(x-y) for (x,y) in zip(u1,u2)]
    v1 = v[0:3]
    v2 = v[3:6]
    u1 = u[0:3]
    u2 = u[3:6]
    # find the common atoms
    if v1 == u1:
      mid_uv = u_mid_v(u2,v2)
      vec1 = make_vec(v1,mid_uv)
      vec2 = make_vec(u2,v2)
    elif v1 == u2:
      mid_uv = u_mid_v(u1,v2)
      vec1 = make_vec(v1,mid_uv)
      vec2 = make_vec(u1,v2)
    elif v2 == u1:
      mid_uv = u_mid_v(u2,v1)
      vec1 = make_vec(v2,mid_uv)
      vec2 = make_vec(u2,v1)
    else: #v2 == u2
      mid_uv = u_mid_v(u1,v1)
      vec1 = make_vec(v2,mid_uv)
      vec2 = make_vec(u1,v1)
    # return the cosine of the angle between the two vectors
    try:
      cos_angle = abs(u_dot_v(vec1,vec2)/u_len(vec1)/u_len(vec2))
    except ZeroDivisionError:
      cos_angle = 1
    return cos_angle

  def nb_overlaps_list(self):
    """(self) -> [list,list]
    Collect information of overlapping nonbonded proxies (neighboring atoms)
    for overlap count calculation.
    - proxies (atoms) considered to overlap when:
        model_distance - van_der_waals_distance < -0.40
    - For every atom consider only worst overlap
    - Do not double count in-line overlaps.
       for example, C-H1  H2-...
       if both C and H1 are overlapping with H2 count only the worst of them
    - Exclude 1-5 interaction of Hydrogen and heavy atom

    Returns:
      Overlapping_atoms_list: [[list of overlaps due to sym operation],
                            [list of all other overlaps]]
    """
    Overlapping_atoms_list = [[],[]]
    overlap_atoms_dict = {}
    # overlap_atoms_dict[atom_key] = [all overlaps information for this atom]
    # overlap_atoms_dict[' HA  ASP A  44 '] =
    #   [ [atom1,vec1,i_seq,j_seq1,model)],[atom2, vec2,i_seq, j_seq2,model]...]
    # vec = [delta_x,delta_y,delta_z] , difference between the atoms coordinates
    overlaps_dict = {}
    # overlaps_dict[vec] = [atom1, atom2, nonbonded_list nbo_record]
    # vec uniquely define a overlap, regardless of atoms order
    for rec in self.nonbonded_list:
      i_seq = rec[1]
      j_seq = rec[2]
      model = rec[3]
      vdw = rec[4]
      symop = rec[5]
      delta = model - vdw
      # check for overlap
      if (delta < -0.40):
        # Check of 1-5 interaction
        if not self.is_1_5_interaction(
                i_seq, j_seq,self.hd_sel,self.full_connectivity_table):
          nbo_vec,key1,key2 = self.get_nbo_keys(rec)
          # nbo_key is a string a string that uniquely identify each overlap
          if key1>key2:
            nbo_key = '::'.join([key2,key1])
          else:
            nbo_key = '::'.join([key1,key2])
          if nbo_key not in overlaps_dict:
            # record new overlap
            overlaps_dict[nbo_key] = [key1,key2,rec]
            if key1 not in overlap_atoms_dict: overlap_atoms_dict[key1] = []
            if key2 not in overlap_atoms_dict: overlap_atoms_dict[key2] = []
            overlap_atoms_dict[key1].append([key2,nbo_vec,nbo_key,symop,model])
            overlap_atoms_dict[key2].append([key1,nbo_vec,nbo_key,symop,model])
    for key in overlap_atoms_dict:
      # for atoms that overlap more than once, check for inline overlaps
      if len(overlap_atoms_dict[key]) > 1:
        temp_nbo_list = []
        n_overlaps = len(overlap_atoms_dict[key])
        # iterate over all overlap combination
        for i in range(n_overlaps-1):
          for j in range(i+1,n_overlaps):
            vec_i = overlap_atoms_dict[key][i][1]
            vec_j = overlap_atoms_dict[key][j][1]
            u = [float(_v) for _v in vec_i.split(',')]
            v = [float(_v) for _v in vec_j.split(',')]
            cos_angle = 0
            # test inline only if the two atoms, overlapping with the
            # common atom, are connected
            overlapping_atom_1 = int(overlap_atoms_dict[key][i][0].split('::')[1])
            overlapping_atom_2 = int(overlap_atoms_dict[key][j][0].split('::')[1])
            if overlapping_atom_1 in self.full_connectivity_table[overlapping_atom_2]:
              if not [0.0,0.0,0.0] in [u,v]:
                # Do not calculate cos_angle when atoms are overlapping
                if overlap_atoms_dict[key][i][3] == '':
                  # ignore overlaps that are cause by symmetry operation
                  cos_angle = self.cos_vec(u,v)
            # atoms consider to be inline if cosine of
            # the angle between vectors > 0.707
            if abs(cos_angle) > 0.707 and (vec_i != vec_j):
              # for inline overlaps keep the closer two(compare models)
              if overlap_atoms_dict[key][i][4] < overlap_atoms_dict[key][j][4]:
                temp_nbo_list.append(overlap_atoms_dict[key][i])
                # remove overlap from overlaps_dict
                remove_key = overlap_atoms_dict[key][j][2]
                if remove_key in overlaps_dict: del overlaps_dict[remove_key]
              else:
                temp_nbo_list.append(overlap_atoms_dict[key][j])
                # remove overlap from overlaps_dict
                remove_key = overlap_atoms_dict[key][i][2]
                if remove_key in overlaps_dict: del overlaps_dict[remove_key]
            else:
              # overlaps are not inline, keep both
              temp_nbo_list.append(overlap_atoms_dict[key][j])
              temp_nbo_list.append(overlap_atoms_dict[key][i])
        overlap_atoms_dict[key] = temp_nbo_list
    for (key,val) in six.iteritems(overlaps_dict):
      if key.split('::')[2] != '':
        # not to symmetry operation
        Overlapping_atoms_list[0].append(val[2])
      else:
        # not due to symmetry operation
        Overlapping_atoms_list[1].append(val[2])
    return Overlapping_atoms_list

class info(object):

  def __init__(self,
    model,
#    geometry_restraints_manager,
    macro_molecule_selection,
#    sites_cart,
#    hd_sel,
#    site_labels=None,
    do_only_macro_molecule=False,
    check_for_unknown_pairs=True):

    geometry_restraints_manager = model.get_restraints_manager().geometry
    xrs = model.get_xray_structure()
    sites_cart = model.get_sites_cart()
    site_labels = xrs.scatterers().extract_labels()
    hd_sel = model.get_hd_selection()
    '''
    Construct nonbonded_overlaps_info, the non-bonded overlaps number and list

    Args:
      geometry_restraints_manager:
      macro_molecule_selection : selection array typically of corresponding to a
        selection string of "protein of dna or rna"
      sites_cart: sites_cart[i] tuple containing the x,y,z coordinates of atom i
      site_labels: a list of lables such as " HA  LEU A  38 ", for each atom
      hd_sel: hd_sel[i] retruns True of False, indicating whether an atom i is
            a Hydrogen or not

    NOTE:
    As of Dec. 2013 manipulation of scatterers can produce scatteres which
    have no lables. In that case, water interaction score will not be accurate

    NOTE:
    Be aware to the parameters:
       assume_hydrogens_all_missing=False,
       hard_minimum_nonbonded_distance=0.0

    The default values of these are True and 0.001, which will alter
    the size of the vdw bonds and the overlaps that being counted

    As of Dec. 2013 manipulation of scatterers can produce scatterers which
    have no labels. In that case, water interaction score will not be accurate

    @author: Youval Dar, LBL 12-2013
    '''
    selection_list = [flex.bool([True]*sites_cart.size())]
    results = []
    # second_grm_selection = macro_molecule_selection.count(False) > 0
    # This is 10 times faster and produces the same result
    second_grm_selection = not macro_molecule_selection.all_eq(True)
    if second_grm_selection:
      selection_list.append(macro_molecule_selection)

    for i, sel in enumerate(selection_list):
      if do_only_macro_molecule and i == 0:
        results.append(compute(
            nonbonded_list=[],
            hd_sel=None,
            full_connectivity_table=None,
            connectivity_table_2=None,
            sites_cart=None))
        continue
      grm = geometry_restraints_manager.select(sel)
      cart = sites_cart.select(sel)
      if site_labels:
        site_label = site_labels.select(sel)
      else:
        site_label = site_labels
      if (check_for_unknown_pairs and
          unknown_pairs_present(grm=grm,sites_cart=cart,site_labels=site_label)):
        msg = "nonbonded overlaps can't be calculated."
        msg += " PDB file contains unknown type pairs. Please provide cif file."
        raise Sorry(msg)
      pair_proxies = grm.pair_proxies(sites_cart=cart,site_labels=site_label)
      proxies_info_nonbonded = pair_proxies.nonbonded_proxies.get_sorted(
        by_value="delta",
        sites_cart=cart,
        site_labels=site_label)

      if proxies_info_nonbonded != None:
        nonbonded_list = proxies_info_nonbonded[0]
      else:
        nonbonded_list = []

      fsc0=grm.shell_sym_tables[0].full_simple_connectivity()
      fsc2=grm.shell_sym_tables[2].full_simple_connectivity()

      results.append(compute(
        nonbonded_list=nonbonded_list,
        hd_sel=hd_sel.select(sel),
        full_connectivity_table=fsc0,
        connectivity_table_2=fsc2,
        sites_cart=cart))

    self.result = nonbonded_overlaps_results()
    r_complete = results[0]
    # all
    self.result.nb_overlaps_all = r_complete.nb_overlaps_all
    self.result.nb_overlaps_proxies_all = r_complete.nb_overlaps_proxies_all
    # Symmetry
    self.result.nb_overlaps_due_to_sym_op = r_complete.nb_overlaps_due_to_sym_op
    self.result.nb_overlaps_proxies_due_to_sym_op = \
      r_complete.nb_overlaps_proxies_due_to_sym_op
    # CCTBX clashscore
    self.result.normalized_nbo_all = r_complete.normalized_nbo_all
    self.result.normalized_nbo_sym = \
      r_complete.normalized_nbo_sym
    # macro molecule
    if second_grm_selection:
      r_macro_mol = results[1]
      nb_overlaps = r_macro_mol.nb_overlaps_non_sym
      self.result.nb_overlaps_macro_molecule = nb_overlaps
      self.result.nb_overlaps_proxies_macro_molecule = \
        r_macro_mol.nb_overlaps_non_sym_overlaps
      clashscore = r_macro_mol.cctbx_clashscore_non_sym
      self.result.normalized_nbo_macro_molecule = clashscore
    else:
      self.result.nb_overlaps_macro_molecule = \
        r_complete.nb_overlaps_non_sym
      self.result.nb_overlaps_proxies_macro_molecule = \
        r_complete.nb_overlaps_non_sym_overlaps
      self.result.normalized_nbo_macro_molecule = \
        r_complete.cctbx_clashscore_non_sym

  def show(self, log=None, nbo_type='all',normalized_nbo=False):
    """
    Show (prints to log) nonbonded_overlaps_info on overlapping atoms

    Args:
      show_normalized_nbo=False Show non-bonded overlaps per 1000 atoms
      log : when no log is given function will print to sys.stdout
      nbo_type (str): The type of overlaps to show
        'all': Show all overlapping atoms
        'sym': Show symmetry related overlaps
        'macro': Show macro molecule overlaps (not including sym related overlaps)

    Returns:
      out_string (str): the output string that is printed to log
    """
    nb_overlaps = self.result
    if not log: log = sys.stdout
    out_list = []
    result_str = '{:<54} :{:5.2f}'
    if normalized_nbo:
      names = ['Total normalized NBO',
             'normalized NBO macro molecule (Protein, RNA, DNA)',
             'normalized NBO due to symmetry']
      scores = [nb_overlaps.normalized_nbo_all,
              nb_overlaps.normalized_nbo_macro_molecule,
              nb_overlaps.normalized_nbo_sym]
      for name,score in zip(names,scores):
        out_list.append(result_str.format(name,round(score,2)))
    names = ['Total non-bonded overlaps',
             'non-bonded overlaps macro molecule (Protein, RNA, DNA)',
             'non-bonded overlaps due to symmetry']
    scores = [nb_overlaps.nb_overlaps_all,
              nb_overlaps.nb_overlaps_macro_molecule,
              nb_overlaps.nb_overlaps_due_to_sym_op]
    for name,score in zip(names,scores):
      out_list.append(result_str.format(name,round(score,2)))
    if nbo_type == 'sym':
      nbo_proxies = nb_overlaps.nb_overlaps_proxies_due_to_sym_op
      title = 'Overlaps due to symmetry operation,'
    elif nbo_type == 'macro':
      nbo_proxies = nb_overlaps.nb_overlaps_proxies_macro_molecule
      title = 'Overlaps in macro molecule (not including sym. related overlaps),'
    else:
      nbo_proxies = nb_overlaps.nb_overlaps_proxies_all
      title = 'Overlapping atoms, complete model,'
    title += ' based on pair_proxies.nonbonded_proxies'
    out_list.append(title)
    out_list.append('='*len(title))
    labels =  ["Overlapping residues info","i_seq","j_seq","model-vdw",
               "sym overlap"]
    lbl_str = '{:^33}|{:^7}|{:^7}|{:^11}|{:<10}'
    out_str = '{:>16}|{:>16}|{:^7}|{:^7}|  {:>6.3f}   |{:^10}|'
    out_list.append(lbl_str.format(*labels))
    out_list.append('-'*73)
    argmented_counts = [0,0]
    def _adjust_count(d):
      d=(abs(d)-0.4)*10
      return d+1
    for data in nbo_proxies:
      # clean and order info for output string
      d = list(data)
      rec_list = [x.replace('pdb=','') for x in d[0]]
      rec_list = [x.replace('"','') for x in rec_list]
      rec_list.extend(d[1:3])
      rec_list.append(d[3]-d[4])
      rec_list.append('1'*bool(d[5]) + ' '*(not bool(d[5])))
      #print(rec_list)
      ptr = 0
      if rec_list[5].strip(): ptr=1
      argmented_counts[ptr] += _adjust_count(rec_list[4])
      out_list.append(out_str.format(*rec_list))
    out_string = '\n'.join(out_list)
    print(out_string, file=log)
    #print(argmented_counts)
    return out_string

def get_macro_mol_sel(pdb_processed_file,selection='protein or dna or rna'):
  """
  Get macro molecule selection from a PBD interpretation object, for non bonded
  overlaps calculation

  Args:
    pdb_processed_file : Object result from pdb_interpretation of a file
    selection (str): By default a string to select the macro molecule.

  Return:
    macro_mol_sel (flex.bool): selection array
  """
  proxies = pdb_processed_file.all_chain_proxies
  cache = proxies.pdb_hierarchy.atom_selection_cache()
  macro_mol_sel = proxies.selection(
    cache=cache,
    string=selection )
  return macro_mol_sel

def create_cif_file_using_ready_set(
        pdb_hierarchy=None,
        crystal_symmetry=None,
        file_name=None,
        log=None):
  """
  When model contains unknown pairs, create a cif file for nonbonded_overlaps
  calculation using READY_SET.

  Args:
    pdb_hierarchy : pdb hierarchy
    file_name (str): pdb file name
    log : output location
    crystal_symmetry : must provide crystal symmetry when using pdb_hierarchy

  Returns:
    (str): the cif file name that was created
  """
  import libtbx.load_env
  has_ready_set = libtbx.env.has_module(name="phenix")
  if file_name:
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    cryst_sym = pdb_inp.crystal_symmetry()
  else:
    assert crystal_symmetry
    cryst_sym = crystal_symmetry
  assert pdb_hierarchy
  if not log: log = sys.stdout
  models = pdb_hierarchy.models()
  assert len(models) == 1
  if not has_ready_set:
    msg = 'phenix.ready_set could not be detected on your system.\n'
    msg += 'Cannot process PDB file'
    print(msg, file=log)
    return [False,False]
  if not file_name:
    file_name = pdb_hierarchy.write_pdb_or_mmcif_file(
        target_format = 'pdb',
        crystal_symmetry=cryst_sym,
        target_filename = 'input_pdb_file_for_ready_set.pdb')
  cmd = "phenix.ready_set {} --silent".format(file_name)
  out = easy_run.fully_buffered(cmd)
  if (out.return_code != 0):
    msg_str = "ready_set crashed - dumping stderr:\n%s"
    raise RuntimeError(msg_str % ( "\n".join(out.stderr_lines)))
  fn_pdb = file_name.replace('.pdb','.updated.pdb')
  fn_cif = file_name.replace('.pdb','.ligands.cif')
  return [fn_cif,fn_pdb]

def unknown_pairs_present(grm,sites_cart,site_labels):
  """
  Test if PDB file contains unknown type pairs

  Args:
    grm (obj): geometry restraints manager
    sites_cart (flex.vec3): atoms sites cart (coordinates)
    site_labels: a list of lables such as " HA  LEU A  38 ", for each atom

  Return:
    (bool): True if PDB file contains unknown type pairs
  """
  pp= grm.pair_proxies(sites_cart=sites_cart,site_labels=site_labels)
  return (pp.nonbonded_proxies.n_unknown_nonbonded_type_pairs != 0)
