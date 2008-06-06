import copy
import scitbx.math
from scitbx import matrix
from scitbx.array_family import flex
from cctbx import sgtbx, crystal
from rstbx.symmetry.subgroup import metric_subgroups, MetricSubgroup, group_classification as group_classif

def group_classif_bravais():
  bravais_only_keys = [1,3,5,16,21,22,23,89,97,155,177,207,211,209]
  rd = {}
  for key in bravais_only_keys: rd[key]=group_classif[key]
  return rd

def echelon_constraints(group,reciprocal_space = 1):
  #    direct space : XT G X  = G
  #reciprocal space : X G* XT = G*
  #as tuple: g = g00,g01,g02,g11,g12,g22 = A F E B D C
  '''  G = g00 g01 g02 = A F E
           g01 g11 g12   F B D
           g02 g12 g22   E D C'''
  n0 = 6*len(group)
  n1 = 6
  m = flex.int(flex.grid(n0,n1))
  i = 0
  for x in group:
    if reciprocal_space==1:
      Rm = matrix.sqr(x.as_double_array()[0:9])
    else:
      Rm = matrix.sqr(x.as_double_array()[0:9]).transpose()
    R = Rm.elems
    Rt = Rm.transpose()
    e = (R[0]*R[0]-1, 2*R[0]*R[1],            2*R[0]*R[2],            R[1]*R[1],2*R[1]*R[2],           R[2]*R[2])
    for x in xrange(6): m[i] = int(e[x]); i+=1
    e = (R[0]*R[3],     R[1]*R[3]+R[0]*R[4]-1,  R[2]*R[3]+R[0]*R[5],  R[1]*R[4],  R[2]*R[4]+R[1]*R[5], R[2]*R[5])
    for x in xrange(6): m[i] = int(e[x]); i+=1
    e = (R[0]*R[6],     R[1]*R[6]+R[0]*R[7],    R[2]*R[6]+R[0]*R[8]-1,R[1]*R[7],  R[2]*R[7]+R[1]*R[8], R[2]*R[8])
    for x in xrange(6): m[i] = int(e[x]); i+=1
    e = (R[3]*R[3],   2*R[3]*R[4],            2*R[3]*R[5],            R[4]*R[4]-1,2*R[4]*R[5],         R[5]*R[5])
    for x in xrange(6): m[i] = int(e[x]); i+=1
    e = (R[3]*R[6],     R[4]*R[6]+R[3]*R[7],    R[5]*R[6]+R[3]*R[8],  R[4]*R[7],  R[5]*R[7]+R[4]*R[8]-1,R[5]*R[8])
    for x in xrange(6): m[i] = int(e[x]); i+=1
    e = (R[6]*R[6],   2*R[6]*R[7],            2*R[6]*R[8],            R[7]*R[7],2*R[7]*R[8],         R[8]*R[8]-1)
    for x in xrange(6): m[i] = int(e[x]); i+=1
  mnew = scitbx.math.row_echelon_form(m)
  i = 0
  #Rearrange row echelon changing coefficient order AFEBDC to ABCDEF
  n0 = mnew
  i=0
  C = flex.int(flex.grid(n0,n1))
  for x in xrange(mnew):
    C[i]=m[i]; C[i+1]=m[i+3]; C[i+2]=m[i+5]; C[i+3]=m[i+4]; C[i+4]=m[i+2]; C[i+5]=m[i+1]
    i+=6
  i=0
  return C

def equal(A,B,tolerance=0.99):
  return abs(A-B) < 1.- tolerance

def bestcmp(a,b):
  if equal(a['max_angular_difference'], b['max_angular_difference']):
    if a.number() > b.number(): return -1
    if a.number() == b.number(): return 0
    else: return 1
  if a['max_angular_difference'] > b['max_angular_difference']: return -1
  return 1

def iotbx_converter(unit_cell,max_delta,bravais_types_only=True,
    space_group_symbol="P 1",force_minimum=False):
    # with regard to "force_minimum": when autoindexing, the orientation
    # matrix may be derived from comparison to a previously indexed case;
    # the setting may be non-standard; therefore we do not want to
    # convert to the reduced cell when calculating metric subgroups.
  input_symmetry = crystal.symmetry(unit_cell=unit_cell,
    space_group_symbol=space_group_symbol)

  M = metric_subgroups(input_symmetry,max_delta,
                       enforce_max_delta_for_generated_two_folds=True,
                       bravais_types_only=bravais_types_only,
                       force_minimum=force_minimum)
  M.labelit_style = []
  for subgroup in M.result_groups:
    # required keys for subgroup:
    #   max_angular_difference
    #   subsym: the centrosymmetric group, referred to the input_cell basis
    #   cb_op_inp_best:  change of basis from the input cell to the best reference cell

    # derived keys added in other frames:
    #   orient: the orientation matrix, in the reference setting

    # methods:
    #   to_reference_setting_as_double_array_transpose (formerly 'matrix')
    #   number: the group number of subsym

    # other attributes:
    #   reduced_group: the acentric group, expressed in input_cell basis
    #   supersym: acentric metric supergroup, input_cell basis

    subgroup.update (  # Get 'bravais' and 'system' keys
      group_classif[sgtbx.space_group_info(
          group=subgroup['supersym'].space_group()).type().number()])
    subgroup['best_group']=subgroup['best_subsym'].space_group()
    #special procedure to get non-centrosymmetric group
    subgroup['reduced_group']=\
      subgroup['subsym'].space_group().build_derived_acentric_group()
    subgroup['constraints']=echelon_constraints(subgroup['reduced_group'])
    M.labelit_style.append(MetricSubgroup().import_iotbx_style(subgroup))
  M.labelit_style.sort(bestcmp)
  return M.labelit_style
