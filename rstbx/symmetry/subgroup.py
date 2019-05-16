from __future__ import absolute_import, division, print_function
from six.moves import cPickle as pickle
from libtbx import adopt_init_args
from cctbx.sgtbx import subgroups
from scitbx.array_family import flex
from cctbx import crystal
from scitbx import matrix
from cctbx import sgtbx
#only need this for legacy table lookup; deprecate later:
from rstbx.symmetry.sgtbx_adaptor import get_patterson_group
from cctbx.sgtbx.lattice_symmetry import metric_supergroup,metric_subgroups as base_subgroups
from cctbx.sgtbx import bravais_types,change_of_basis_op
from six.moves import zip
find_max_delta = sgtbx.lattice_symmetry_find_max_delta

#For LABELIT, the derive from the iotbx subgroup list such that the
#  input cell is NOT reduced to "minimum symmetry".  This allows alignment
#  of the basis with a previous data wedge.

class metric_subgroups(base_subgroups):
  def __init__(self,
        input_symmetry,
        max_delta,
        enforce_max_delta_for_generated_two_folds=True,
        bravais_types_only=True,force_minimum=False,
        best_monoclinic_beta=True,
        interest_focus="metric_symmetry"):
    adopt_init_args(self,locals())
    self.parse_reference()
    self.result_groups = []
    self.cb_op_best_cell_vector =[]
    self.minimum_symmetry = self.input_symmetry
    self.cb_op_inp_minimum = sgtbx.change_of_basis_op() #identity
    if force_minimum:
      self.change_input_to_minimum_cell()
    if interest_focus=="metric_symmetry":
      self.derive_result_group_list(group_of_interest=self.lattice_group_info())
    elif interest_focus=="input_symmetry":
      self.derive_result_group_list_original_input(group_of_interest=
        self.input_symmetry.space_group_info())

    #as a future reindexing reference, record the cb_op to best cell
    for i,j in zip(self.result_groups,self.cb_op_best_cell_vector):
      i['cb_op_best_cell']=j

  reasonable_cutoff = 10.0 # maximum direct-space mean-square deviation of two orientations
                           # measured in Angstrom^2; larger value indicates no match
                           # later, this could be redefined as a fraction of cell size

  def derive_result_group_list_original_input(self,group_of_interest):
    # more-or-less a capitulation to code duplication.  Don't want to copy code
    # from cctbx.sgtbx.lattice_symmetry but can't figure out a quick way around it.

    # Get list of sub-spacegroups
    subgrs = subgroups.subgroups(group_of_interest).groups_parent_setting()

    # Order sub-groups
    sort_values = flex.double()
    for group in subgrs:
      order_z = group.order_z()
      space_group_number = sgtbx.space_group_type(group, False).number()
      assert 1 <= space_group_number <= 230
      sort_values.append(order_z*1000+space_group_number)
    perm = flex.sort_permutation(sort_values, True)

    for i_subgr in perm:
      acentric_subgroup = subgrs[i_subgr]
      acentric_supergroup = metric_supergroup(acentric_subgroup)
      # Add centre of inversion to acentric lattice symmetry
      centric_group = sgtbx.space_group(acentric_subgroup)
      # Make symmetry object: unit-cell + space-group
      # The unit cell is potentially modified to be exactly compatible
      # with the space group symmetry.
      subsym = crystal.symmetry(
        unit_cell=self.minimum_symmetry.unit_cell(),
        space_group=centric_group,
        assert_is_compatible_unit_cell=False)
      supersym = crystal.symmetry(
        unit_cell=self.minimum_symmetry.unit_cell(),
        space_group=acentric_supergroup,
        assert_is_compatible_unit_cell=False)
      # Convert subgroup to reference setting
      cb_op_minimum_ref = subsym.space_group_info().type().cb_op()
      ref_subsym = subsym.change_basis(cb_op_minimum_ref)
      # Ignore unwanted groups
      if (self.bravais_types_only and
          not str(ref_subsym.space_group_info()) in bravais_types.centric):
        continue
      # Choose best setting for monoclinic and orthorhombic systems
      cb_op_best_cell = self.change_of_basis_op_to_best_cell(ref_subsym)
      best_subsym = ref_subsym.change_basis(cb_op_best_cell)
      # Total basis transformation
      cb_op_best_cell = change_of_basis_op(str(cb_op_best_cell),stop_chars='',r_den=144,t_den=144)
      cb_op_minimum_ref=change_of_basis_op(str(cb_op_minimum_ref),stop_chars='',r_den=144,t_den=144)
      self.cb_op_inp_minimum=change_of_basis_op(str(self.cb_op_inp_minimum),stop_chars='',r_den=144,t_den=144)
      cb_op_inp_best = cb_op_best_cell * (cb_op_minimum_ref * self.cb_op_inp_minimum)
      # Use identity change-of-basis operator if possible
      if (best_subsym.unit_cell().is_similar_to(self.input_symmetry.unit_cell())):
        cb_op_corr = cb_op_inp_best.inverse()
        try:
          best_subsym_corr = best_subsym.change_basis(cb_op_corr)
        except RuntimeError as e:
          if (str(e).find("Unsuitable value for rational rotation matrix.") < 0):
            raise
        else:
          if (best_subsym_corr.space_group() == best_subsym.space_group()):
            cb_op_inp_best = cb_op_corr * cb_op_inp_best
      """Note:  The following call does not work if n_ltr >1 for the space group"""
      if acentric_subgroup.n_ltr() == 1:
        m_a_d = find_max_delta(
                                  reduced_cell=self.minimum_symmetry.unit_cell(),
                                  space_group=acentric_subgroup)
      else:
        m_a_d = 0.0
      self.result_groups.append({'subsym':subsym,
                                 'supersym':supersym,
                                 'ref_subsym':ref_subsym,
                                 'best_subsym':best_subsym,
                                 'cb_op_inp_best':cb_op_inp_best,
                                 'max_angular_difference':m_a_d
                                })

  def change_of_basis_op_to_best_cell(self,ref_subsym):

    #plain algorithm for routine work.  Triclinic input symmetry is already
    #  in the minimum (reduced) form.  Determine the best cell for conventional
    #  monoclinic and orthorhombic crystal systems.
    subgroup_cb_op_best_cell = ref_subsym.change_of_basis_op_to_best_cell(
      best_monoclinic_beta=self.best_monoclinic_beta)

    #fancy algorithm for reindexing work.  Triclinic symmetry is aligned
    #  with a reference orientation, but is not necessarily in reduced
    #  form.  Guarantee that the best cell
    #  conforms to the one that has already been determined for a reference
    #  dataset collected on the same cystal.  This option is only
    #  executed if the reference information is present in a pickled file.
    if self.reference_subgroups != None: #subsequent code only for reindexing cases
      import inspect
        #admittedly this is a private interface.  Rely on having this method always
        #  called by derive_result_group_*.  Otherwise this breaks.
        #do not save any frame references, thus avoiding reference circles:
      current_supersym = inspect.currentframe().f_back.f_locals['supersym']
      current_bravais = str(bravais_types.bravais_lattice(sgtbx.space_group_info(
          group=current_supersym.space_group()).type().number()))
      current_cb_op_minimum_ref = inspect.currentframe().f_back.f_locals['cb_op_minimum_ref']
      guess1_cb_op_inp_best = subgroup_cb_op_best_cell * current_cb_op_minimum_ref *\
                              self.cb_op_inp_minimum
        #assume this frame is called in a very specific way such that four frames
        # back gives the lepage module, class character.  It's the only way to
        # pull out the current triclinic orientation.
      guess1_orient = self.current_orientation.change_basis(matrix.sqr(
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
        guess2_orient = self.current_orientation.change_basis(matrix.sqr(
          guess2_cb_op_inp_best.c().as_double_array()[0:9]).transpose().elems)
        dmsd = guess2_orient.direct_mean_square_difference(refitem['orient'])
        if dmsd < self.reasonable_cutoff:
          subgroup_cb_op_best_cell = refitem['cb_op_best_cell']
          break

    self.cb_op_best_cell_vector.append(subgroup_cb_op_best_cell)
    return subgroup_cb_op_best_cell

  def parse_reference(self):
    try:
      o_file = ('crystal_orientation')
      f = open(o_file,'r')
      reference_ori = pickle.load(f)
      self.reference_subgroups = pickle.load(f)
      f.close()

    except Exception:
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

  def sg_digest(self,add_inv=False):
    return ' | '.join([
      self['cb_op_inp_best'].as_xyz(),
      sgtbx.space_group_info(group=self['reduced_group']).type().lookup_symbol(),
      ])

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
