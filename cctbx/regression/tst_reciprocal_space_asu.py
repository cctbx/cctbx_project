from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import miller
from cctbx.development import debug_utils
import sys
from six.moves import range

def get_all_incides(space_group, abs_range):
  result = {}
  for h0 in range(-abs_range,abs_range+1):
    for h1 in range(-abs_range,abs_range+1):
      for h2 in range(-abs_range,abs_range+1):
        h = (h0, h1, h2)
        if (h == (0,0,0)): continue
        equiv = miller.sym_equiv_indices(space_group, h).indices()
        for h_e in equiv:
          h_s = h_e.h()
          assert h_s == h_e.hr()
          result[h_s] = 0
  return flex.miller_index(list(result.keys()))

cond_dict = {}

def filter_asu(space_group_info, indices, verbose):
  asu = space_group_info.reciprocal_space_asu()
  assert asu.is_reference()
  cond = asu.reference_as_string()
  if (verbose): print(cond)
  cond_dict[cond] = 0
  is_centric = space_group_info.group().is_centric
  result = flex.miller_index()
  for hkl in indices:
    h,k,l = hkl
    vfy = eval(cond)
    assert vfy == asu.is_inside(hkl)
    if (asu.is_inside(hkl)):
      result.append(hkl)
      if (not is_centric(hkl)):
        result.append([-e for e in hkl])
  return result

def expand_indices(space_group, asu_indices):
  result = {}
  for hkl in asu_indices:
    equiv = miller.sym_equiv_indices(space_group, hkl).indices()
    for h_e in equiv:
      h_s = h_e.h()
      assert h_s == h_e.hr()
      assert h_s not in result
      result[h_s] = 0
  return flex.miller_index(list(result.keys()))

def exercise(space_group_info, verbose=0):
  if (not space_group_info.reciprocal_space_asu().is_reference()): return
  all_indices = get_all_incides(space_group_info.group(), abs_range=4)
  if (verbose): print("all_indices.size():", all_indices.size())
  asu_indices = filter_asu(space_group_info, all_indices, verbose)
  if (verbose): print("asu_indices.size():", asu_indices.size())
  exp_indices = expand_indices(space_group_info.group(), asu_indices)
  if (verbose): print("exp_indices.size():", exp_indices.size())
  assert all_indices.size() == exp_indices.size()

def run_call_back(flags, space_group_info):
  exercise(space_group_info, verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  if (0):
    print("len(cond_dict):", len(cond_dict))
    for k in cond_dict.keys():
      print(k)

if (__name__ == "__main__"):
  run()
