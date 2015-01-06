from libtbx.utils import Sorry
import math

class compute(object):
  """
  Compute clashscore using cctbx geometry restraints manager. Similar but not
  identical to MolProbity clashscore. Symmetry interactions are included.

  @author: Youval Dar (LBL 2013)
  """

  def __init__(
    self,
    nonbonded_list,
    hd_sel,
    full_connectivty_table,
    connectivty_table_2,
    sites_cart):
    """
    Arguments:
    nonbonded_list: a list with items in the following format
                    ([pdb labels], i_seq, j_seq, model, vdw_distance, sym_op_j, rt_mx)
                    i_seq,j_seq: position of residues in the pdb file
                    model: The pdb or other model distance
                    vdw_distance: Van Der Waals distance
                    sym_op_j: is this a product of a symmetry opeation
                    rt_mx: Rotation matrix for symmetry operation
    hd_sel: hd_sel[i] retruns True of False, indicating whether an atom i is a Hydrogen or not
    full_connectivty_table: full_connectivty_table[i] is a dictinary constaining a
                            list of all atoms connected to atom i
    connectivty_table_2: connectivty_table[i] is a dictinary constaining a
                         list of all atoms connected to atom i
    sites_cart: sites_cart[i] tuple containing the x,y,z coordinates of atom i
    """
    self.nonbonded_list = nonbonded_list
    self.hd_sel = hd_sel
    self.full_connectivty_table = full_connectivty_table
    self.connectivty_table_2 = connectivty_table_2
    self.sites_cart = sites_cart
    if nonbonded_list != []:
      try:
        clashlist = self.clashscore_clash_list()
      except TypeError as e:
        # When proxies_info_nonbonded are not available
        if e == "vec3_double' object is not callable":
          raise Sorry(e)
        else:
          clashlist = [[],[],[]]
          print e
      except Sorry as e:
        raise Sorry(e)
      except Exception as e:
        m='Failed processing proxies_info_nonbonded in clashscore_clash_list()'
        raise Sorry(m)
      # Collect, seperatly, clashes due to symmetry operation and those that are not
      self.nb_clash_proxies_due_to_sym_op = clashlist[0]
      self.nb_clash_proxies_solvent_solvent = clashlist[1]
      self.nb_clash_proxies_simple = clashlist[2]
      self.nb_clash_proxies_all_clashes = clashlist[0] + clashlist[1] + clashlist[2]
      # clashscore is the number of steric overlaps (>0.4A) per 1000 atoms
      n_atoms = len(self.sites_cart)
      clashscore_due_to_sym_op = len(clashlist[0])*1000/n_atoms
      clashscore_solvent_solvent = len(clashlist[1])*1000/n_atoms
      clashscore_simple = len(clashlist[2])*1000/n_atoms
      clashscore_all_clashes = clashscore_due_to_sym_op + \
        clashscore_solvent_solvent + clashscore_simple
      #
      self.nb_clashscore_due_to_sym_op = clashscore_due_to_sym_op
      self.nb_clashscore_solvent_solvent = clashscore_solvent_solvent
      self.nb_clashscore_simple = clashscore_simple
      self.nb_clashscore_all_clashes = clashscore_all_clashes
    else:
      self.nb_clashscore_due_to_sym_op = 0
      self.nb_clashscore_solvent_solvent = 0
      self.nb_clashscore_simple = 0
      self.nb_clashscore_all_clashes = 0
      #
      self.nb_clash_proxies_due_to_sym_op = []
      self.nb_clash_proxies_solvent_solvent = []
      self.nb_clash_proxies_simple = []
      self.nb_clash_proxies_all_clashes = []

  @staticmethod
  def is_1_5_interaction(i_seq,j_seq,hd_sel,full_connectivty_table):
    '''(int,int,bool array,list of lists of int) -> bool
    Check if we have 1-5 interaction between two hydrogen and heavy atom

    arguments:
    i_seq,j_seq:  are the number of the atoms we are checking in the pdb table
    hd_sel: hd_sel[i] retruns True of False, indicating whether an atom i is a Hydrogen or not
    full_connectivty_table: full_connectivty_table[i] is a dictinary constaining a
                            list of all atoms connected to atom

    retruns:
    True if we have 1-5 hydrogen and heavy atom interaction
    '''
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
          # add all new connetions for the current step
          connections = connections.union(set(full_connectivty_table[key]))
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

  def get_clash_keys(self,record):
    '''(str) -> str,str,str

    Collect atoms keys and coordinates

    Argumets:
    record : a clash record, for example:
    (['pdb=" HA  LEU A  38 "', 'pdb="HD23 LEU A  38 "'], 523, 532, 1.8108674716831155, 2.44, '', None)

    Returns:
    clash_vec: a unique key for a clash. for example: '0.144,0.323,1.776,0.154,0.327,1.786'
    key1,key2: are the atoms keys. for example ' HA  LEU A  38 ::120::' or ' CB ASN A  55 ::52::sym.op.'
               "atom1 label::atom1 number::sym.op""
    '''
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
    clash_vec = '{0:.3f},{1:.3f},{2:.3f},{3:.3f},{4:.3f},{5:.3f}'.format(*vec)
    # make key1 always smaller than key2
    if key1 > key2:
      key1,key2 = key2,key1
    return clash_vec,key1,key2

  @staticmethod
  def cos_vec(u,v):
    '''(tuple,tuple) -> float

    Calculate the cosine used to evaluate wheather the atoms
    coliding are inline or not

    Arguments:
    u,v: lists containing x1,y1,z1,x2,y2,z2 coordinates
    Returns:
    cos_angle: the cosine of the angle between center of the common atom
    and the mid point between the other two atoms.
    '''
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

  @staticmethod
  def are_both_solvent(key,connectivty_table_2):
    '''(str,dict) -> bool
    retrun true if the two atoms, which number is specified in "key" are single, double or tripple atom molecule

    Arguments:
    connectivty_table_2: connectivty_table[i] is a dictinary constaining a
                         list of all atoms connected to atom i
    key: a list "atom1 lable::atom1 number::sym.op.::atom2 label::atom2 number::sym.op."

    Return True is both atoms have no connection to atoms two bonds apart
    '''
    key = key.split('::')
    i = int(key[1])
    j = int(key[4])
    i_is_solvent = list(connectivty_table_2[i])==[]
    j_is_solvent = list(connectivty_table_2[j])==[]
    return i_is_solvent and j_is_solvent

  def clashscore_clash_list(self):
    '''(self) -> [list,list,list]
    Collect inforamtion of clashing nonbonded proxies (neighboring atoms) for clash
    score calculation.
    - proxies considered to clash when model_distance - van_der_waals_distance < -0.41
    - For every atom consider only worst clash
    - Do not double count in-line clashes.
       for example, C-H1  H2-...
       if both C and H1 are clashing with H2 count only the worst of them
    - Exclude 1-5 interaction of Hydrogen and heavy atom

    Retruns:
    clashing_atoms_list: [[list of clashes due to sym operation],
      [list of solvent-solvent clashes],[list of all other clashes]]
    '''
    clashing_atoms_list = [[],[],[]]
    clashing_atoms_dict = {}
    # clashing_atoms_dict[atom_key] = [all clashes information for this atom]
    # clashing_atoms_dict[' HA  ASP A  44 '] = [[atom1,vec1, i_seq, j_seq1,model)],
    #                                           [atom2, vec2,i_seq, j_seq2,model]...]
    # vec = [delta_x,delta_y,delta_z] , the difference between the atoms coordinates
    clashes_dict = {}
    # clashes_dict[vec] = [atom1, atom2, clash_record]
    # vec uniquly define a clash, regardless of atoms order
    for rec in self.nonbonded_list:
      i_seq = rec[1]
      j_seq = rec[2]
      model = rec[3]
      vdw = rec[4]
      symop = rec[5]
      delta = model - vdw
      # check for clash
      if (delta < -0.41):
        # Check of 1-5 interaction
        if not self.is_1_5_interaction(i_seq, j_seq,self.hd_sel,self.full_connectivty_table):
          clash_vec,key1,key2 = self.get_clash_keys(rec)
          # clash_key is a string a string that uniquly identify each clash
          if key1>key2:
            clash_key = '::'.join([key2,key1])
          else:
            clash_key = '::'.join([key1,key2])
          if clash_key not in clashes_dict:
            # record new clash
            clashes_dict[clash_key] = [key1,key2,rec]
            if key1 not in clashing_atoms_dict: clashing_atoms_dict[key1] = []
            if key2 not in clashing_atoms_dict: clashing_atoms_dict[key2] = []
            clashing_atoms_dict[key1].append([key2,clash_vec,clash_key,symop,model])
            clashing_atoms_dict[key2].append([key1,clash_vec,clash_key,symop,model])
    for key in clashing_atoms_dict:
      # for atoms that clash more than once, check for inline clashes
      if len(clashing_atoms_dict[key]) > 1:
        temp_clash_list = []
        n_clashes = len(clashing_atoms_dict[key])
        # itereate over all clash combination
        for i in range(n_clashes-1):
          for j in range(i+1,n_clashes):
            vec_i = clashing_atoms_dict[key][i][1]
            vec_j = clashing_atoms_dict[key][j][1]
            u = map(float,vec_i.split(','))
            v = map(float,vec_j.split(','))
            cos_angle = 0
            # test inline only if the two atoms, clashing with the
            # common atom, are connected
            clashing_atom_1 = int(clashing_atoms_dict[key][i][0].split('::')[1])
            clashing_atom_2 = int(clashing_atoms_dict[key][j][0].split('::')[1])
            if clashing_atom_1 in self.full_connectivty_table[clashing_atom_2]:
              if not [0.0,0.0,0.0] in [u,v]:
                # Do not calculate cos_angle when atoms are overlapping
                if clashing_atoms_dict[key][i][3] == '':
                  # ignore clashes that are cause by symmetry operation
                  cos_angle = self.cos_vec(u,v)
            # atoms consider to be inline if cosine of the angle between vectors > 0.866
            if cos_angle > 0.867 and (vec_i != vec_j):
              # for inline clashes keep the closer two(compare models)
              if clashing_atoms_dict[key][i][4] < clashing_atoms_dict[key][j][4]:
                temp_clash_list.append(clashing_atoms_dict[key][i])
                # remove clash from clashes_dict
                remove_key = clashing_atoms_dict[key][j][2]
                if remove_key in clashes_dict: del clashes_dict[remove_key]
              else:
                temp_clash_list.append(clashing_atoms_dict[key][j])
                # remove clash from clashes_dict
                remove_key = clashing_atoms_dict[key][i][2]
                if remove_key in clashes_dict: del clashes_dict[remove_key]
            else:
              # clashes are not inline, keep both
              temp_clash_list.append(clashing_atoms_dict[key][j])
              temp_clash_list.append(clashing_atoms_dict[key][i])
        clashing_atoms_dict[key] = temp_clash_list
    for (key,val) in clashes_dict.iteritems():
      if key.split('::')[2] != '':
        # not to symmetry operation
        clashing_atoms_list[0].append(val[2])
      elif self.are_both_solvent(key,self.connectivty_table_2):
        # solvent-solvent clashes
        clashing_atoms_list[1].append(val[2])
      else:
        # not due to symmetry operation or solvent-solvent, simple clashes
        clashing_atoms_list[2].append(val[2])
    return clashing_atoms_list

class info(object):
  def __init__(self,
    geometry_restraints_manager,
    sites_cart,
    hd_sel,
    site_labels=None):
    '''
    Construct nonbonded_clash_info, the non-bonded clashss lists and scores

    Arguments:
    sites_cart: sites_cart[i] tuple containing the x,y,z coordinates of atom i
    site_labels: a list of lables such as " HA  LEU A  38 ", for each atom
    hd_sel: hd_sel[i] retruns True of False, indicating whether an atom i is a Hydrogen or not

    NOTE:
    As of Dec. 2013 manipulation of scatterers can produce scatteres which
    have no lables. In that case, water interaction score will not be accurate

    NOTE:
    Be aware to the parameters:
       assume_hydrogens_all_missing=False,
       hard_minimum_nonbonded_distance=0.0

    The default values of these are True and 0.001, which will alter
    the size of the vdw bonds and the clashes that being counted

    As of Dec. 2013 manipulation of scatterers can produce scatteres whish
    have no lables. In that case, water interaction score will not be accurate

    @author: Youval Dar, LBL 12-2013
    '''
    pair_proxies = geometry_restraints_manager.pair_proxies(
          sites_cart=sites_cart,
          site_labels=site_labels)
    if pair_proxies.nonbonded_proxies.n_unknown_nonbonded_type_pairs != 0:
      msg = "nonbonded clashscore can't be calculated."
      msg += " PDB file contains unknown type pairs. Please provide cif file."
      raise Sorry(msg)
    proxies_info_nonbonded = pair_proxies.nonbonded_proxies.get_sorted(
      by_value="delta",
      sites_cart=sites_cart,
      site_labels=site_labels)

    if proxies_info_nonbonded != None:
      nonbonded_list = proxies_info_nonbonded[0]
    else:
      nonbonded_list = []

    grm=geometry_restraints_manager
    fsc0=grm.shell_sym_tables[0].full_simple_connectivity()
    fsc2=grm.shell_sym_tables[2].full_simple_connectivity()
    self.result = compute(
      nonbonded_list=nonbonded_list,
      hd_sel=hd_sel,
      full_connectivty_table=fsc0,
      connectivty_table_2=fsc2,
      sites_cart=sites_cart)

  def show(self, log=None):
    """
    Show (prints to log) CCTBX nonbonded_clash_info on clashing atoms

    Args:
      sites_cart: sites_cart[i] tuple containing the x,y,z coordinates of atom i
      site_labels: a list of lables such as " HA  LEU A  38 ", for each atom
      hd_sel: hd_sel[i] retruns True of False, indicating whether an
        atom i is a Hydrogen or not
      log : when no log is given function will print to sys.stdout

    Returns:
      out_string (str): the output string that is printed to log
    """
    nb_clash = self.result
    if not log: log = sys.stdout
    out_list = []
    all_clashes = round(nb_clash.nb_clashscore_all_clashes,2)
    out_list.append('Total non-bonded CCTBX clashscore: {}'.format(all_clashes))
    simple = round(nb_clash.nb_clashscore_simple,2)
    out_list.append(
      'Clashscore without symmetry or solvent clashes: {}'.format(simple))
    due_to_sym_op = round(nb_clash.nb_clashscore_due_to_sym_op,2)
    out_list.append('Clashscore due to symmetry: '.format(due_to_sym_op))
    solvent = round(nb_clash.nb_clashscore_solvent_solvent,2)
    out_list.append('Clashscore solvent - solvent clashes: '.format(solvent))
    title = 'All clashing atoms info table'
    title += ' based on pair_proxies.nonbonded_proxies data'
    out_list.append(title)
    out_list.append('='*len(title))
    labels =  ["pdb labels","i_seq","j_seq","model","vdw","sym_op_j","rt_mx"]
    out_str = '{:^50} | {:<5} | {:<5} | {:<6.4} | {:<6.4} | {:<10} | {:<8}'
    out_list.append(out_str.format(*labels))
    out_list.append('-'*108)
    for data in nb_clash.nb_clash_proxies_all_clashes:
      out_list.append(out_str.format(*data))
    out_string = '\n'.join(out_list)
    print >> log,out_string
    return out_string
