from libtbx.utils import Sorry
from scitbx.array_family import flex
from libtbx import easy_run
import iotbx.pdb
import string
import math
import sys

class cctbx_clashscore_results(object):
  """ Container for cctbx clashscore results """

  def __init__(self):
    """
    - Clashscore is number of clashes per 1000 atoms
    - Clashscore for All and symmetry related clashes are evaluated
      relative to the total number of atoms. Clashscore for macro molecule
      does not include symmetry related clashes and use the number of atoms
      in the macro molecule.
    - Clash proxies is a list containing the information on the clashing atoms
      in the format:([pdb labels],i_seq,j_seq,model,vdw_distance,sym_op_j,rt_mx)

    The following clash scores are calculated:
      cctbx_clashscore_due_to_sym_op: (calculated with the complete model)
      cctbx_clashscore_macro_molecule: (protein, DNA and RNA)
        excluding symmetry related clashes
      cctbx_clashscore_all: (calculated with the complete model)

    For each of the clashscores a list of clashing proxies (atoms) is provided:
      cctbx_clash_proxies_due_to_sym_op
      cctbx_clash_proxies_macro_molecule
      cctbx_clash_proxies_all

    """
    self.cctbx_clashscore_due_to_sym_op = 0
    self.cctbx_clashscore_macro_molecule = 0
    self.cctbx_clashscore_all = 0
    #
    self.cctbx_clash_proxies_due_to_sym_op = []
    self.cctbx_clash_proxies_macro_molecule = []
    self.cctbx_clash_proxies_all = []

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
                    ([pdb labels],i_seq,j_seq,model,vdw_distance,sym_op_j,rt_mx)
                    i_seq,j_seq: position of residues in the pdb file
                    model: The pdb or other model distance
                    vdw_distance: Van Der Waals distance
                    sym_op_j: is this a product of a symmetry operation
                    rt_mx: Rotation matrix for symmetry operation
    hd_sel: hd_sel[i] returns True of False, indicating whether an atom i
            is a Hydrogen or not
    full_connectivty_table: full_connectivty_table[i] is a dictionary containing
                            a list of all atoms connected to atom i
    connectivty_table_2: connectivty_table[i] is a dictionary containing a
                         list of all atoms separated by three bonds from atom i
                         (1 - 4 interaction)
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
          clashlist = [[],[]]
          print e
      except Sorry as e:
        raise Sorry(e)
      except Exception as e:
        m='Failed processing proxies_info_nonbonded in clashscore_clash_list()'
        raise Sorry(m)
      # Collect clashes information
      self.cctbx_clash_proxies_due_to_sym_op = clashlist[0]
      self.cctbx_non_sym_clashes = clashlist[1]
      self.cctbx_clash_proxies_all = clashlist[0] + clashlist[1]
      # clashscore is the number of steric overlaps (>0.4A) per 1000 atoms
      self.n_atoms = len(self.sites_cart)
      clashscore_sym = len(clashlist[0])*1000/self.n_atoms
      clashscore_non_sym = len(clashlist[1])*1000/self.n_atoms
      clashscore_all_clashes = clashscore_sym + clashscore_non_sym
      #
      self.cctbx_clashscore_due_to_sym_op = clashscore_sym
      self.cctbx_clashscore_non_sym = clashscore_non_sym
      self.cctbx_clashscore_all = clashscore_all_clashes
    else:
      self.cctbx_clashscore_due_to_sym_op = 0
      self.cctbx_clashscore_non_sym = 0
      self.cctbx_clashscore_all = 0
      #
      self.cctbx_clash_proxies_due_to_sym_op = []
      self.cctbx_non_sym_clashes = []
      self.cctbx_clash_proxies_all = []

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
    colliding are inline or not

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

  def clashscore_clash_list(self):
    '''(self) -> [list,list]
    Collect information of clashing nonbonded proxies (neighboring atoms)
    for clash score calculation.
    - proxies (atoms) considered to clash when:
        model_distance - van_der_waals_distance < -0.41
    - For every atom consider only worst clash
    - Do not double count in-line clashes.
       for example, C-H1  H2-...
       if both C and H1 are clashing with H2 count only the worst of them
    - Exclude 1-5 interaction of Hydrogen and heavy atom

    Returns:
    clashing_atoms_list: [[list of clashes due to sym operation],
                          [list of all other clashes]]
    '''
    clashing_atoms_list = [[],[]]
    clashing_atoms_dict = {}
    # clashing_atoms_dict[atom_key] = [all clashes information for this atom]
    # clashing_atoms_dict[' HA  ASP A  44 '] =
    #   [ [atom1,vec1,i_seq,j_seq1,model)],[atom2, vec2,i_seq, j_seq2,model]...]
    # vec = [delta_x,delta_y,delta_z] , difference between the atoms coordinates
    clashes_dict = {}
    # clashes_dict[vec] = [atom1, atom2, nonbonded_list clash_record]
    # vec uniquely define a clash, regardless of atoms order
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
        if not self.is_1_5_interaction(
                i_seq, j_seq,self.hd_sel,self.full_connectivty_table):
          clash_vec,key1,key2 = self.get_clash_keys(rec)
          # clash_key is a string a string that uniquely identify each clash
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
        # iterate over all clash combination
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
            # atoms consider to be inline if cosine of
            # the angle between vectors > 0.707
            if abs(cos_angle) > 0.707 and (vec_i != vec_j):
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
      else:
        # not due to symmetry operation
        clashing_atoms_list[1].append(val[2])
    return clashing_atoms_list

class info(object):

  def __init__(self,
    geometry_restraints_manager,
    macro_molecule_selection,
    sites_cart,
    hd_sel,
    site_labels=None):
    '''
    Construct nonbonded_clash_info, the non-bonded clashss lists and scores

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
    the size of the vdw bonds and the clashes that being counted

    As of Dec. 2013 manipulation of scatterers can produce scatterers which
    have no labels. In that case, water interaction score will not be accurate

    @author: Youval Dar, LBL 12-2013
    '''
    selection_list = [flex.bool([True]*sites_cart.size())]
    results = []
    second_grm_selection = macro_molecule_selection.count(False) > 0
    if second_grm_selection:
      selection_list.append(macro_molecule_selection)

    for sel in selection_list:
      grm = geometry_restraints_manager.select(sel)
      cart = sites_cart.select(sel)
      if site_labels:
        site_label = site_labels.select(sel)
      else:
        site_label = site_labels
      if unknown_pairs_present(grm=grm,sites_cart=cart,site_labels=site_label):
        msg = "nonbonded clashscore can't be calculated."
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
        full_connectivty_table=fsc0,
        connectivty_table_2=fsc2,
        sites_cart=cart))

    self.result = cctbx_clashscore_results()
    r_complete = results[0]
    # all
    self.result.cctbx_clashscore_all = r_complete.cctbx_clashscore_all
    self.result.cctbx_clash_proxies_all = r_complete.cctbx_clash_proxies_all
    # Symmetry
    self.result.cctbx_clashscore_due_to_sym_op = \
      r_complete.cctbx_clashscore_due_to_sym_op
    self.result.cctbx_clash_proxies_due_to_sym_op = \
      r_complete.cctbx_clash_proxies_due_to_sym_op
    # macro molecule
    if second_grm_selection:
      r_macro_mol = results[1]
      clashscore = r_macro_mol.cctbx_clashscore_non_sym
      self.result.cctbx_clashscore_macro_molecule = clashscore
      self.result.cctbx_clash_proxies_macro_molecule = \
        r_macro_mol.cctbx_non_sym_clashes
    else:
      self.result.cctbx_clashscore_macro_molecule = \
        r_complete.cctbx_clashscore_non_sym
      self.result.cctbx_clash_proxies_macro_molecule = \
        r_complete.cctbx_non_sym_clashes

  def show(self, log=None, clash_type='all'):
    """
    Show (prints to log) CCTBX nonbonded_clash_info on clashing atoms

    Args:
      sites_cart: sites_cart[i] tuple containing the x,y,z coordinates of atom i
      site_labels: a list of lables such as " HA  LEU A  38 ", for each atom
      hd_sel: hd_sel[i] retruns True of False, indicating whether an
        atom i is a Hydrogen or not
      log : when no log is given function will print to sys.stdout
      clash_type (str): The type of clashes to show
        'all': Show all clashing atoms
        'sym': Show symmetry related clashes
        'macro': Show macro molecule clashes (not including sym related clashes)


    Returns:
      out_string (str): the output string that is printed to log
    """
    cctbx_clash = self.result
    if not log: log = sys.stdout
    out_list = []
    result_str = '{:<50} :{:5.2f}'
    names = ['Total non-bonded CCTBX clashscore',
             'Clashscore macro molecule (Protein, RNA, DNA)',
             'Clashscore due to symmetry']
    scores = [cctbx_clash.cctbx_clashscore_all,
              cctbx_clash.cctbx_clashscore_macro_molecule,
              cctbx_clash.cctbx_clashscore_due_to_sym_op]
    for name,score in zip(names,scores):
      out_list.append(result_str.format(name,round(score,2)))
    if clash_type == 'sym':
      clash_proxies = cctbx_clash.cctbx_clash_proxies_due_to_sym_op
      title = 'Clashes due to symmetry operation,'
    elif clash_type == 'macro':
      clash_proxies = cctbx_clash.cctbx_clash_proxies_macro_molecule
      title = 'Clashes in macro molecule (not including sym. related clashes),'
    else:
      clash_proxies = cctbx_clash.cctbx_clash_proxies_all
      title = 'Clashing atoms, complete model,'
    title += ' based on pair_proxies.nonbonded_proxies'
    out_list.append(title)
    out_list.append('='*len(title))
    labels =  ["Clashing residues info","i_seq","j_seq","model-vdw","sym clash"]
    lbl_str = '{:^33}|{:^7}|{:^7}|{:^11}|{:<10}'
    out_str = '{:>16}|{:>16}|{:^7}|{:^7}|  {:>6.3f}   |{:^8}|'
    out_list.append(lbl_str.format(*labels))
    out_list.append('-'*71)
    for data in clash_proxies:
      # clean and order info for output string
      d = list(data)
      rec_list = [x.replace('pdb=','') for x in d[0]]
      rec_list = [x.replace('"','') for x in rec_list]
      rec_list.extend(d[1:3])
      rec_list.append(d[3]-d[4])
      rec_list.append('1'*bool(d[5]) + ' '*(not bool(d[5])))
      out_list.append(out_str.format(*rec_list))
    out_string = '\n'.join(out_list)
    print >> log,out_string
    return out_string

def check_and_add_hydrogen(
        pdb_hierarchy=None,
        file_name=None,
        nuclear=False,
        keep_hydrogens=True,
        verbose=False,
        model_number=0,
        n_hydrogen_cut_off=0,
        time_limit=120,
        allow_multiple_models=True,
        crystal_symmetry=None,
        log=None):
  """
  If no hydrogens present, force addition for clashscore calculation.
  Use REDUCE to add the hydrogen atoms.

  Args:
    pdb_hierarchy : pdb hierarchy
    file_name (str): pdb file name
    nuclear (bool): When True use nuclear cloud x-H distances and vdW radii,
      otherwise use electron cloud x-H distances and vdW radii
    keep_hydrogens (bool): when True, if there are hydrogen atoms, keep them
    verbose (bool): verbosity of printout
    model_number (int): the number of model to use
    time_limit (int): limit the time it takes to add hydrogen atoms
    n_hydrogen_cut_off (int): when number of hydrogen atoms < n_hydrogen_cut_off
      force keep_hydrogens tp True
    allow_multiple_models (bool): Allow models that contain more than one model
    crystal_symmetry : must provide crystal symmetry when using pdb_hierarchy

  Returns:
    (str): PDB string
    (bool): True when PDB string was updated
  """
  if file_name:
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    cryst_sym = pdb_inp.crystal_symmetry()
    pdb_hierarchy = pdb_inp.construct_hierarchy()
  elif not allow_multiple_models:
    assert crystal_symmetry
    cryst_sym = crystal_symmetry
  else:
    cryst_sym = None
  assert pdb_hierarchy
  assert model_number < len(pdb_hierarchy.models())
  if not log: log = sys.stdout
  models = pdb_hierarchy.models()
  if (len(models) > 1) and (not allow_multiple_models):
    raise Sorry("When using CCTBX clashscore, provide only a single model.")
  model = models[model_number]
  r = iotbx.pdb.hierarchy.root()
  mdc = model.detached_copy()
  r.append_model(mdc)
  if keep_hydrogens:
    elements = r.atoms().extract_element()
    h_count = elements.count(' H') + elements.count(' D')
    if h_count > n_hydrogen_cut_off:
      has_hd = True
    else:
      has_hd = False
    if not has_hd:
      if verbose:
        print >> log,"\nNo H/D atoms detected - forcing hydrogen addition!\n"
      keep_hydrogens = False
  import libtbx.load_env
  has_reduce = libtbx.env.has_module(name="reduce")
  # add hydrogen if needed
  if has_reduce and (not keep_hydrogens):
    # set reduce running parameters
    build = "phenix.reduce -oh -his -flip -pen9999 -keep -allalt -limit{}"
    if nuclear:
      build += " -nuc -"
    else:
      build += " -"
    build = build.format(time_limit)
    trim = "phenix.reduce -quiet -trim -"
    stdin_lines = r.as_pdb_string(cryst_sym)
    clean_out = easy_run.fully_buffered(trim,stdin_lines=stdin_lines)
    if (clean_out.return_code != 0) :
      msg_str = "Reduce crashed with command '%s' - dumping stderr:\n%s"
      raise RuntimeError(msg_str % (trim, "\n".join(clean_out.stderr_lines)))
    build_out = easy_run.fully_buffered(build,stdin_lines=clean_out.stdout_lines)
    if (build_out.return_code != 0) :
      msg_str = "Reduce crashed with command '%s' - dumping stderr:\n%s"
      raise RuntimeError(msg_str % (build, "\n".join(build_out.stderr_lines)))
    reduce_str = string.join(build_out.stdout_lines, '\n')
    return reduce_str,True
  else:
    if not has_reduce:
      msg = 'phenix.reduce could not be detected on your system.\n'
      msg += 'Cannot add hydrogen to PDB file'
      print >> log,msg
    return r.as_pdb_string(cryst_sym),False

def get_macro_mol_sel(pdb_processed_file):
  """
  A central place to get macro molecule selection for CCTBX clashscore
  calculation

  Args:
    pdb_processed_file : Object result from pdb_interpretation of a file

  Return:
    macro_mol_sel (flex.bool): selection array
  """
  macro_molecule_selection_str = 'protein or dna or rna'
  proxies = pdb_processed_file.all_chain_proxies
  cache = proxies.pdb_hierarchy.atom_selection_cache()
  macro_mol_sel = proxies.selection(
    cache=cache,
    string=macro_molecule_selection_str )
  return macro_mol_sel

def create_cif_file_using_ready_set(
        pdb_hierarchy=None,
        crystal_symmetry=None,
        file_name=None,
        log=None):
  """
  When model contains unknown pairs, create a cif file for clashscore
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
    print >> log,msg
    return [False,False]
  if not file_name:
    file_name = 'input_pdb_file_for_ready_set.pdb'
    open(file_name,'w').write(pdb_hierarchy.as_pdb_string(cryst_sym))
  cmd = "phenix.ready_set {} --silent".format(file_name)
  out = easy_run.fully_buffered(cmd)
  if (out.return_code != 0) :
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
