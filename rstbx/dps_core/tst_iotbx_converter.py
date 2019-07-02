from __future__ import absolute_import, division, print_function
from six.moves import cPickle as pickle
from six.moves import cStringIO as StringIO
from cctbx import crystal,sgtbx,uctbx
from cctbx.sgtbx import lattice_symmetry
from rstbx.dps_core.lepage import iotbx_converter

class subgroup_comparator:

  def __init__(self,symmetry):
    self.symmetry = symmetry

  def get_input_symmetry_subgroups(self):
    cb_op_input_to_primitive = self.symmetry.change_of_basis_op_to_primitive_setting()
    primitive_symmetry = self.symmetry.primitive_setting()
    cb_op_input_to_minimum = primitive_symmetry.change_of_basis_op_to_minimum_cell() * cb_op_input_to_primitive

    subgroup_list = iotbx_converter(
      self.symmetry.change_basis(cb_op_input_to_minimum).unit_cell(),max_delta=3.0,
      bravais_types_only=False,
      space_group_symbol= str(self.symmetry.change_basis(cb_op_input_to_minimum).space_group_info()),
      force_minimum=False,interest_focus="input_symmetry")

    return subgroup_list

  def get_metric_symmetry_subgroups(self):
    subgroup_list = lattice_symmetry.metric_subgroups(
      self.symmetry,3.0,bravais_types_only=False,
      best_monoclinic_beta=False).result_groups

    return subgroup_list

  def show_list(self,any_subgroup_list):
    for subgroup in any_subgroup_list:
      subgroup['best_subsym'].show_summary()
      print()

  def get_best_subsym_list(self, any_subgroup_list):
    return [subgroup['best_subsym'] for subgroup in any_subgroup_list]

  def get_list_as_string(self,any_subgroup_list):
    S = StringIO()
    pickle.dump( [subgroup['best_subsym'] for subgroup in any_subgroup_list], S)
    return S.getvalue()

  def get_string_as_list(self,any_subgroup_string):
    from io import BytesIO
    S = BytesIO(any_subgroup_string.encode("ascii"))
    subgroups = pickle.load( S )
    return subgroups

  def compare_lists(self,a,b):
    assert len(a) == len(b)
    all_comparisons=[]
    for item in a:
      all_comparisons.append(False)
      for ritem in b:
        if item.is_similar_symmetry(ritem):
          all_comparisons[-1]=True
    assert False not in all_comparisons



def get_one_example():
  return subgroup_comparator(symmetry =
  # symmetry from 3ged
   crystal.symmetry(
    unit_cell=uctbx.unit_cell((124.287,124.287,162.608,90.0,90.0,90.0)),
    space_group=sgtbx.space_group_info('I 41 2 2').group())
  )

expected_output = """Unit cell: (124.287, 124.287, 162.608, 90, 90, 90)
Space group: I 41 2 2 (No. 98)

Unit cell: (124.287, 124.287, 162.608, 90, 90, 90)
Space group: I 41 (No. 80)

Unit cell: (124.287, 124.287, 162.608, 90, 90, 90)
Space group: I 21 21 21 (No. 24)

Unit cell: (162.608, 175.768, 175.768, 90, 90, 90)
Space group: F 2 2 2 (No. 22)

Unit cell: (124.287, 162.608, 124.287, 90, 90, 90)
Space group: I 1 2 1 (No. 5)

Unit cell: (119.725, 175.768, 119.725, 90, 94.4545, 90)
Space group: I 1 2 1 (No. 5)

Unit cell: (124.287, 124.287, 162.608, 90, 90, 90)
Space group: I 1 2 1 (No. 5)

Unit cell: (124.287, 124.287, 162.608, 90, 90, 90)
Space group: I 1 2 1 (No. 5)

Unit cell: (119.725, 175.768, 119.725, 90, 94.4545, 90)
Space group: I 1 2 1 (No. 5)

Unit cell: (119.725, 119.725, 119.725, 94.4545, 117.462, 117.462)
Space group: P 1 (No. 1)

Unit cell: (124.287, 124.287, 162.608, 90, 90, 90)
Space group: I 4/m m m (No. 139)

Unit cell: (124.287, 124.287, 162.608, 90, 90, 90)
Space group: I 4/m (No. 87)

Unit cell: (124.287, 124.287, 162.608, 90, 90, 90)
Space group: I m m m (No. 71)

Unit cell: (162.608, 175.768, 175.768, 90, 90, 90)
Space group: F m m m (No. 69)

Unit cell: (175.768, 162.608, 124.287, 90, 135, 90)
Space group: C 1 2/m 1 (No. 12)

Unit cell: (162.608, 175.768, 119.725, 90, 132.773, 90)
Space group: C 1 2/m 1 (No. 12)

Unit cell: (204.667, 124.287, 124.287, 90, 127.392, 90)
Space group: C 1 2/m 1 (No. 12)

Unit cell: (204.667, 124.287, 124.287, 90, 127.392, 90)
Space group: C 1 2/m 1 (No. 12)

Unit cell: (162.608, 175.768, 119.725, 90, 132.773, 90)
Space group: C 1 2/m 1 (No. 12)

Unit cell: (119.725, 119.725, 119.725, 117.462, 117.462, 94.4545)
Space group: P -1 (No. 2)

"""

expected_subgroup_lists=[
"""(lp0
ccopy_reg
_reconstructor
p1
(ccctbx.crystal
symmetry
p2
c__builtin__
object
p3
Ntp4
Rp5
(dp6
S'_unit_cell'
p7
ccctbx_uctbx_ext
unit_cell
p8
((F124.28700000000005
F124.28700000000005
F162.608
F90.0
F90.0
F90.0
tp9
tp10
Rp11
sS'_space_group_info'
p12
g1
(ccctbx.sgtbx
space_group_info
p13
g3
Ntp14
Rp15
(ccctbx_sgtbx_ext
space_group
p16
(S' I 4bw 2bw'
p17
tp18
Rp19
tp20
bsbag1
(g2
g3
Ntp21
Rp22
(dp23
g7
g8
((F124.28700000000005
F124.28700000000005
F162.608
F90.0
F90.0
F90.0
tp24
tp25
Rp26
sg12
g1
(g13
g3
Ntp27
Rp28
(g16
(S' I 4bw'
p29
tp30
Rp31
tp32
bsbag1
(g2
g3
Ntp33
Rp34
(dp35
g7
g8
((F124.28700000000003
F124.28700000000005
F162.60800000000003
F90.0
F90.0
F90.0
tp36
tp37
Rp38
sg12
g1
(g13
g3
Ntp39
Rp40
(g16
(S' I 2b 2c'
p41
tp42
Rp43
tp44
bsbag1
(g2
g3
Ntp45
Rp46
(dp47
g7
g8
((F162.608
F175.76836102666491
F175.76836102666496
F90.0
F90.0
F90.0
tp48
tp49
Rp50
sg12
g1
(g13
g3
Ntp51
Rp52
(g16
(S' F 2 2'
p53
tp54
Rp55
tp56
bsbag1
(g2
g3
Ntp57
Rp58
(dp59
g7
g8
((F124.28700000000003
F162.608
F124.28700000000005
F90.0
F90.000000000000014
F90.0
tp60
tp61
Rp62
sg12
g1
(g13
g3
Ntp63
Rp64
(g16
(S' C 2y (x,y,-x+z)'
p65
tp66
Rp67
tp68
bsbag1
(g2
g3
Ntp69
Rp70
(dp71
g7
g8
((F119.72455721571909
F175.76836102666491
F119.72455721571913
F90.0
F94.454526861950484
F90.0
tp72
tp73
Rp74
sg12
g1
(g13
g3
Ntp75
Rp76
(g16
(S' C 2y (x,y,-x+z)'
p77
tp78
Rp79
tp80
bsbag1
(g2
g3
Ntp81
Rp82
(dp83
g7
g8
((F124.28700000000005
F124.28700000000003
F162.60800000000003
F90.0
F90.000000000000014
F90.0
tp84
tp85
Rp86
sg12
g1
(g13
g3
Ntp87
Rp88
(g16
(S' C 2y (x,y,-x+z)'
p89
tp90
Rp91
tp92
bsbag1
(g2
g3
Ntp93
Rp94
(dp95
g7
g8
((F124.28700000000011
F124.28700000000003
F162.60800000000003
F90.0
F90.000000000000057
F90.0
tp96
tp97
Rp98
sg12
g1
(g13
g3
Ntp99
Rp100
(g16
(S' C 2y (x,y,-x+z)'
p101
tp102
Rp103
tp104
bsbag1
(g2
g3
Ntp105
Rp106
(dp107
g7
g8
((F119.72455721571912
F175.76836102666496
F119.72455721571913
F90.0
F94.454526861950498
F90.0
tp108
tp109
Rp110
sg12
g1
(g13
g3
Ntp111
Rp112
(g16
(S' C 2y (x,y,-x+z)'
p113
tp114
Rp115
tp116
bsbag1
(g2
g3
Ntp117
Rp118
(dp119
g7
g8
((F119.72455721571913
F119.72455721571913
F119.72455721571913
F94.454526861950498
F117.4623774424071
F117.4623774424071
tp120
tp121
Rp122
sg12
g1
(g13
g3
Ntp123
Rp124
(g16
(S' P 1'
p125
tp126
Rp127
tp128
bsba.""",
"""(lp0
ccopy_reg
_reconstructor
p1
(ccctbx.crystal
symmetry
p2
c__builtin__
object
p3
Ntp4
Rp5
(dp6
S'_unit_cell'
p7
ccctbx_uctbx_ext
unit_cell
p8
((F124.28700000000005
F124.28700000000005
F162.608
F90.0
F90.0
F90.0
tp9
tp10
Rp11
sS'_space_group_info'
p12
g1
(ccctbx.sgtbx
space_group_info
p13
g3
Ntp14
Rp15
(ccctbx_sgtbx_ext
space_group
p16
(S'-I 4 2'
p17
tp18
Rp19
tp20
bsbag1
(g2
g3
Ntp21
Rp22
(dp23
g7
g8
((F124.28700000000005
F124.28700000000005
F162.608
F90.0
F90.0
F90.0
tp24
tp25
Rp26
sg12
g1
(g13
g3
Ntp27
Rp28
(g16
(S'-I 4'
p29
tp30
Rp31
tp32
bsbag1
(g2
g3
Ntp33
Rp34
(dp35
g7
g8
((F124.28700000000003
F124.28700000000005
F162.60800000000003
F90.0
F90.0
F90.0
tp36
tp37
Rp38
sg12
g1
(g13
g3
Ntp39
Rp40
(g16
(S'-I 2 2'
p41
tp42
Rp43
tp44
bsbag1
(g2
g3
Ntp45
Rp46
(dp47
g7
g8
((F162.608
F175.76836102666491
F175.76836102666496
F90.0
F90.0
F90.0
tp48
tp49
Rp50
sg12
g1
(g13
g3
Ntp51
Rp52
(g16
(S'-F 2 2'
p53
tp54
Rp55
tp56
bsbag1
(g2
g3
Ntp57
Rp58
(dp59
g7
g8
((F175.76836102666491
F162.608
F124.28700000000003
F90.0
F135.0
F90.0
tp60
tp61
Rp62
sg12
g1
(g13
g3
Ntp63
Rp64
(g16
(S'-C 2y'
p65
tp66
Rp67
tp68
bsbag1
(g2
g3
Ntp69
Rp70
(dp71
g7
g8
((F162.608
F175.76836102666491
F119.72455721571909
F90.0
F132.77273656902477
F90.0
tp72
tp73
Rp74
sg12
g1
(g13
g3
Ntp75
Rp76
(g16
(S'-C 2y'
p77
tp78
Rp79
tp80
bsbag1
(g2
g3
Ntp81
Rp82
(dp83
g7
g8
((F204.66709562848644
F124.28700000000003
F124.28700000000005
F90.0
F127.39194857784936
F90.0
tp84
tp85
Rp86
sg12
g1
(g13
g3
Ntp87
Rp88
(g16
(S'-C 2y'
p89
tp90
Rp91
tp92
bsbag1
(g2
g3
Ntp93
Rp94
(dp95
g7
g8
((F204.66709562848644
F124.28700000000003
F124.28700000000005
F90.0
F127.39194857784936
F90.0
tp96
tp97
Rp98
sg12
g1
(g13
g3
Ntp99
Rp100
(g16
(S'-C 2y'
p101
tp102
Rp103
tp104
bsbag1
(g2
g3
Ntp105
Rp106
(dp107
g7
g8
((F162.608
F175.76836102666496
F119.72455721571912
F90.0
F132.77273656902477
F90.0
tp108
tp109
Rp110
sg12
g1
(g13
g3
Ntp111
Rp112
(g16
(S'-C 2y'
p113
tp114
Rp115
tp116
bsbag1
(g2
g3
Ntp117
Rp118
(dp119
g7
g8
((F119.72455721571913
F119.72455721571913
F119.72455721571913
F117.4623774424071
F117.4623774424071
F94.454526861950498
tp120
tp121
Rp122
sg12
g1
(g13
g3
Ntp123
Rp124
(g16
(S'-P 1'
p125
tp126
Rp127
tp128
bsba."""]

if __name__=="__main__":
  EX = get_one_example()
  input_sym_subgroups = EX.get_input_symmetry_subgroups()
  metric_sym_subgroups = EX.get_metric_symmetry_subgroups()

  #EX.show_list(input_sym_subgroups)
  #EX.show_list(metric_sym_subgroups)

  #print  EX.get_list_as_string(input_sym_subgroups)
  #print  EX.get_list_as_string(metric_sym_subgroups)

  # list the subgroups of a particular space group & compare to expected reference
  input_sym = expected_subgroup_lists[0]
  input_sym_list = EX.get_string_as_list(input_sym)
  EX.compare_lists(input_sym_list, EX.get_best_subsym_list(input_sym_subgroups))

  # list the subgroups of a metric symmetry & compare to expected reference
  metric_sym = expected_subgroup_lists[1]
  metric_sym_list = EX.get_string_as_list(metric_sym)
  EX.compare_lists(metric_sym_list, EX.get_best_subsym_list(metric_sym_subgroups))

  print("OK")
