from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import adptbx
from libtbx import table_utils
import cctbx.sgtbx.lattice_symmetry
import cctbx.sgtbx.cosets
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
import mmtbx.scaling
from mmtbx.scaling import absolute_scaling
from mmtbx.scaling import matthews, twin_analyses
from mmtbx.scaling import basic_analyses, data_statistics
import libtbx.phil.command_line
from cStringIO import StringIO
from scitbx.python_utils import easy_pickle
import sys, os


class reindexing(object):
  __doc__=""" Reindexing matrices """
  def __init__(self,
               set_a,
               set_b,
               out=None,
               relative_length_tolerance=0.05,
               absolute_angle_tolerance=3.0,
               lattice_symmetry_max_delta=3.0):

    self.set_a = set_a
    self.set_b = set_b

    ##----------
    self.change_of_basis_op_to_minimum_cell_a=\
      set_a.change_of_basis_op_to_minimum_cell()
    self.change_of_basis_op_to_minimum_cell_b=\
      set_b.change_of_basis_op_to_minimum_cell()
    ##----------
    self.minimum_cell_symmetry_a = crystal.symmetry.change_basis(
      set_a,
      cb_op=self.change_of_basis_op_to_minimum_cell_a)
    self.minimum_cell_symmetry_b = crystal.symmetry.change_basis(
      set_b,
      cb_op=self.change_of_basis_op_to_minimum_cell_b)
    ##----------
    self.lattice_group_a = sgtbx.lattice_symmetry.group(
      self.minimum_cell_symmetry_a.unit_cell(),
      max_delta=lattice_symmetry_max_delta)
    self.lattice_group_a.expand_inv(sgtbx.tr_vec((0,0,0)))
    self.lattice_group_a.make_tidy()

    self.lattice_group_b = sgtbx.lattice_symmetry.group(
      self.minimum_cell_symmetry_b.unit_cell(),
      max_delta=lattice_symmetry_max_delta)
    self.lattice_group_b.expand_inv(sgtbx.tr_vec((0,0,0)))
    self.lattice_group_b.make_tidy()
    ##----------
    self.lattice_symmetry_a = crystal.symmetry(
      unit_cell=self.minimum_cell_symmetry_a.unit_cell(),
      space_group_info=sgtbx.space_group_info(group=self.lattice_group_a),
      assert_is_compatible_unit_cell=False)

    self.lattice_symmetry_b = crystal.symmetry(
      unit_cell=self.minimum_cell_symmetry_b.unit_cell(),
      space_group_info=sgtbx.space_group_info(group=self.lattice_group_b),
      assert_is_compatible_unit_cell=False)
    ##----------
    self.intensity_symmetry_a = \
       self.minimum_cell_symmetry_a.reflection_intensity_symmetry(
         anomalous_flag=set_a.anomalous_flag())

    self.intensity_symmetry_b = \
       self.minimum_cell_symmetry_b.reflection_intensity_symmetry(
         anomalous_flag=set_b.anomalous_flag())
    ##----------
    c_inv_rs = self.minimum_cell_symmetry_a.unit_cell().\
      similarity_transformations(
        other=self.minimum_cell_symmetry_b.unit_cell(),
        relative_length_tolerance=relative_length_tolerance,
        absolute_angle_tolerance=absolute_angle_tolerance)


    min_bases_msd = None
    similarity_cb_op = None

    for c_inv_r in c_inv_rs:
      c_inv = sgtbx.rt_mx(sgtbx.rot_mx(c_inv_r))
      cb_op = sgtbx.change_of_basis_op(c_inv).inverse()
      bases_msd = self.minimum_cell_symmetry_a.unit_cell() \
        .bases_mean_square_difference(
          other=cb_op.apply(self.minimum_cell_symmetry_b.unit_cell()))
      if (min_bases_msd is None
          or min_bases_msd > bases_msd):
        min_bases_msd = bases_msd
        similarity_cb_op = cb_op
    if (similarity_cb_op is None): return []

    common_lattice_group = sgtbx.space_group(self.lattice_group_a)
    for s in self.lattice_group_b.build_derived_acentric_group() \
               .change_basis(similarity_cb_op):
      try: common_lattice_group.expand_smx(s)
      except RuntimeError: return []
    common_lattice_group.make_tidy()
    result = []
    for s in sgtbx.cosets.double_unique(
               g=common_lattice_group,
               h1=self.intensity_symmetry_a.space_group()
                   .build_derived_acentric_group()
                   .make_tidy(),
               h2=self.intensity_symmetry_b.space_group()
                   .build_derived_acentric_group()
                   .change_basis(similarity_cb_op)
                   .make_tidy()):
      if (s.r().determinant() > 0):
        result.append(sgtbx.change_of_basis_op(s) * similarity_cb_op)
    self.matrices = result
    self.cc_values= []
    self.matches = []
    self.table=None
    self.analyse()

  def combined_cb_op(self, cb_op):
    s = self.change_of_basis_op_to_minimum_cell_a
    o = self.change_of_basis_op_to_minimum_cell_b
    return s.inverse() * cb_op.new_denominators(s) * o

  def analyse(self):
    ## As the coset decompision is carried out on the minimum cell
    ## The re-indexing laws need to be transform to the original
    ## spacegroup
    table_data=[]
    for cb_op in self.matrices:
      cb_op_comb = self.combined_cb_op(cb_op)
      ## FIX ASU MAPPING HERE

      tmp_set_b = self.set_b.change_basis(cb_op_comb).map_to_asu()
      tmp_set_a, tmp_set_b = self.set_a.map_to_asu().common_sets(
        tmp_set_b,
        assert_is_similar_symmetry=False)
      tmp_cc = tmp_set_a.correlation(
        tmp_set_b,
        assert_is_similar_symmetry=False)
      ## STore the cc values
      self.cc_values.append(  tmp_cc.coefficient()  )
      self.matches.append(
        float(tmp_set_a.indices().size())/float(self.set_a.indices().size()))




  def select_and_transform(self,
                           out=None,
                           matches_cut_off=0.75,
                           input_array=None):
    if out is None:
      out = sys.stdout
    ## hopsa
    max_cc=-1.0
    location = 0
    table_data=[]
    for ii in range(len(self.matrices)):
      table_data.append(
        [self.matrices[ii].as_hkl(),
         "%4.3f"%(self.cc_values[ii]),
         "%4.3f"%(self.matches[ii]),
         '   ']
        )

      if self.matches[ii]>=matches_cut_off:
        if max_cc<self.cc_values[ii]:
          max_cc = self.cc_values[ii]
          location = ii

    legend = ('Operator', 'Correlation', 'matches (%)', 'choice')
    table_data[location][3]=' <--- '
    self.table = table_utils.format([legend]+table_data,
                                       comments=None,
                                       has_header=True,
                                       separate_rows=False,
                                       prefix='| ',
                                       postfix=' |')

    print >> out, self.table
    
    if input_array is not None:
      return input_array.change_basis( self.matrices[location] ).map_to_asu()
