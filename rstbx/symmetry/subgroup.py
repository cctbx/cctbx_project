import pickle
from libtbx import adopt_init_args
from scitbx import matrix
from cctbx import sgtbx
#only need this for legacy table lookup; deprecate later:
from rstbx.symmetry.sgtbx_adaptor import get_patterson_group
from rstbx.command_line.lattice_symmetry import metric_subgroups as base_subgroups

group_classification = {
  1:{'system': 'triclinic','bravais': 'aP'},
  3:{'system': 'monoclinic','bravais': 'mP'},
  5:{'system': 'monoclinic','bravais': 'mC'},
  16:{'system': 'orthorhombic','bravais': 'oP'},
  21:{'system': 'orthorhombic','bravais': 'oC'},
  22:{'system': 'orthorhombic','bravais': 'oF'},
  23:{'system': 'orthorhombic','bravais': 'oI'},
  75:{'system': 'tetragonal','bravais': 'tP'}, # P4
  89:{'system': 'tetragonal','bravais': 'tP'},
  79:{'system': 'tetragonal','bravais': 'tI'}, # I4
  97:{'system': 'tetragonal','bravais': 'tI'},
  146:{'system': 'rhombohedral','bravais': 'hR'}, # R3
  155:{'system': 'rhombohedral','bravais': 'hR'},
  143:{'system': 'hexagonal','bravais': 'hP'}, # P3
  149:{'system': 'hexagonal','bravais': 'hP'}, # P312
  150:{'system': 'hexagonal','bravais': 'hP'}, # P321
  168:{'system': 'hexagonal','bravais': 'hP'}, # P6
  177:{'system': 'hexagonal','bravais': 'hP'},
  195:{'system': 'cubic','bravais': 'cP'}, # P23
  207:{'system': 'cubic','bravais': 'cP'},
  197:{'system': 'cubic','bravais': 'cI'}, # I23
  211:{'system': 'cubic','bravais': 'cI'},
  196:{'system': 'cubic','bravais': 'cF'}, # F23
  209:{'system': 'cubic','bravais': 'cF'},
}

#For LABELIT, the derive from the iotbx subgroup list such that the
#  input cell is NOT reduced to "minimum symmetry".  This allows alignment
#  of the basis with a previous data wedge.

class metric_subgroups(base_subgroups):
  def __init__(self,
        input_symmetry,
        max_delta,
        enforce_max_delta_for_generated_two_folds=True,
        bravais_types_only=True,force_minimum=False,
        best_monoclinic_beta=True):
    adopt_init_args(self,locals())
    self.parse_reference()
    self.result_groups = []
    self.cb_op_best_cell_vector =[]
    self.minimum_symmetry = self.input_symmetry
    self.cb_op_inp_minimum = sgtbx.change_of_basis_op() #identity
    if force_minimum:
      self.change_input_to_minimum_cell()
    self.derive_result_group_list()

    #as a future reindexing reference, record the cb_op to best cell
    for i,j in zip(self.result_groups,self.cb_op_best_cell_vector):
      i['cb_op_best_cell']=j

  reasonable_cutoff = 10.0 # maximum direct-space mean-square deviation of two orientations
                           # measured in Angstrom^2; larger value indicates no match
                           # later, this could be redefined as a fraction of cell size

  def change_of_basis_op_to_best_cell(self,ref_subsym):

    #plain algorithm for routine work.  Triclinic input symmetry is already
    #  in the minimum (reduced) form.  Determine the best cell for conventional
    #  monoclinic and orthorhombic crystal systems.
    subgroup_cb_op_best_cell = ref_subsym.change_of_basis_op_to_best_cell(
      best_monoclinic_beta=True)

    #fancy algorithm for reindexing work.  Triclinic symmetry is aligned
    #  with a reference orientation, but is not necessarily in reduced
    #  form.  Guarantee that the best cell
    #  conforms to the one that has already been determined for a reference
    #  dataset collected on the same cystal.  This option is only
    #  executed if the reference information is present in a pickled file.
    try:
      assert self.reference_subgroups != None #subsequent code only for reindexing cases
      import inspect
        #do not save any frame references, thus avoiding reference circles:
      current_supersym = inspect.currentframe().f_back.f_locals['supersym']
      current_bravais = group_classification[sgtbx.space_group_info(
          group=current_supersym.space_group()).type().number()]['bravais']
      current_cb_op_minimum_ref = inspect.currentframe().f_back.f_locals['cb_op_minimum_ref']
      guess1_cb_op_inp_best = subgroup_cb_op_best_cell * current_cb_op_minimum_ref *\
                              self.cb_op_inp_minimum
        #assume this frame is called in a very specific way such that four frames
        # back gives the lepage module, class character.  It's the only way to
        # pull out the current triclinic orientation.
      current_orientation = inspect.currentframe().f_back.f_back.f_back.f_back.f_locals['orientation']
      guess1_orient = current_orientation.__copy__()
      guess1_orient.change_basis(matrix.sqr(
        guess1_cb_op_inp_best.c().as_double_array()[0:9]).transpose().elems)

      for j,refitem in enumerate(self.reference_subgroups):
        if current_bravais != refitem['bravais']: continue

        # first guess.  test plain algorithm best cell for match to refitem
        dmsd = guess1_orient.direct_mean_square_difference(refitem['orient'])
        if dmsd < self.reasonable_cutoff:
          break

        #second guess.  use the cb_op_best_cell from the reference indexing solution
        guess2_cb_op_inp_best = refitem['cb_op_best_cell'] * current_cb_op_minimum_ref *\
                              self.cb_op_inp_minimum
        guess2_orient = current_orientation.__copy__()
        guess2_orient.change_basis(matrix.sqr(
          guess2_cb_op_inp_best.c().as_double_array()[0:9]).transpose().elems)
        dmsd = guess2_orient.direct_mean_square_difference(refitem['orient'])
        if dmsd < self.reasonable_cutoff:
          subgroup_cb_op_best_cell = refitem['cb_op_best_cell']
          break

    except:
      pass # do nothing; use the plain algorithm

    self.cb_op_best_cell_vector.append(subgroup_cb_op_best_cell)
    return subgroup_cb_op_best_cell

  def parse_reference(self):
    try:
      o_file = ('crystal_orientation')
      f = open(o_file,'r')
      reference_ori = pickle.load(f)
      self.reference_subgroups = pickle.load(f)
      f.close()

    except:
      self.reference_subgroups = None

class MetricSubgroup(dict):

  def import_iotbx_style(self,subgroup):
    #subgroup is a dictionary from iotbx.lattice_symmetry
    self.update(subgroup)
    return self

  def number(self):
    return self['subsym'].space_group().type().number()

  #formerly 'matrix'
  def to_reference_setting_as_double_array_transpose(self):
    return matrix.sqr(
      self['cb_op_inp_best'].c_inv().as_double_array()[0:9]).transpose().elems

  def to_reference_setting_as_double_array_inverse_transpose(self):
    return matrix.sqr(
      self['cb_op_inp_best'].c().as_double_array()[0:9]).transpose().elems

  def digest(self,add_inv=False):
    return '#'.join([
      '%.5f'%self['max_angular_difference'],
      self.pretty_cb_op(),
      self['bravais'],
      sgtbx.space_group_info(group=self['reduced_group']).type().lookup_symbol(),
      sgtbx.space_group_info(group=self['best_group']).type().lookup_symbol(),
      ' '.join([str(int(x)) for x in self['constraints']]),
      ])

  def short_digest(self,add_inv=False):
    return "%7.4f %2s (%-9s)"%(self['max_angular_difference'],      self['bravais'],
      sgtbx.space_group_info(group=self['best_group']).type().lookup_symbol(),
      )



  def reference_lookup_symbol(self):
    return sgtbx.space_group_info(group=self['best_group']).type().lookup_symbol()

  def pretty_cb_op(self):
    for x in self['cb_op_inp_best'].c_inv().as_double_array():
      assert int(x)==float(x)
    return ' '.join(['%2d'%int(x) for x in matrix.sqr(
      self['cb_op_inp_best'].c_inv().as_double_array()[0:9]
       ).transpose().elems])

  def __getitem__(self,key):
    # legacy code to support table lookup of group type
    # patterson:=the centric group built from the reference setting of Subgroup
    if key=='patterson':
      PG = get_patterson_group(self['group'])
      return sgtbx.space_group_info(group = PG)
    return dict.__getitem__(self,key)
