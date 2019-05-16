from __future__ import absolute_import, division, print_function
from six.moves import range

def intify(a):
    return tuple([int(round(val)) for val in a])

def reference_map(sg, mi):
  from cctbx import sgtbx
  asu = sgtbx.reciprocal_space_asu(sg.type())
  isym_ = []
  mi_ = []

  for hkl in mi:
    found = False
    for i_inv in range(sg.f_inv()):
      for i_smx in range(sg.n_smx()):
        rt_mx = sg(0, i_inv, i_smx)
        hkl_ = intify(hkl * rt_mx.r())
        if asu.is_inside(hkl_):
          mi_.append(hkl_)
          if i_inv:
            isym_.append(- i_smx)
          else:
            isym_.append(i_smx)
          found = True
          break

    if found:
      continue
    else:
      assert(not sg.is_centric())

    for i_inv in range(sg.f_inv()):
      for i_smx in range(sg.n_smx()):
        rt_mx = sg(0, i_inv, i_smx)
        _hkl = [-h for h in hkl]
        mhkl_ = intify(_hkl * rt_mx.r())
        if asu.is_inside(mhkl_):
          mi_.append(mhkl_)
          isym_.append(- i_smx)
          found = True
          break

  return mi_, isym_

def tst_map_to_asu_isym(anomalous_flag):
  from cctbx import sgtbx
  from cctbx.miller import map_to_asu_isym
  from cctbx.array_family import flex

  mi = flex.miller_index()
  i = flex.int()

  import random

  nhkl = 1000

  for j in range(nhkl):
    hkl = [random.randint(-10, 10) for j in range(3)]
    mi.append(hkl)
    i.append(0)

  spacegroup = sgtbx.space_group_symbols(195).hall()
  sg = sgtbx.space_group(spacegroup)
  mi_, isym_ = reference_map(sg, mi)
  map_to_asu_isym(sg.type(), anomalous_flag, mi, i)

  for j in range(nhkl):
    assert(i[j] == isym_[j])

if __name__ == '__main__':
  tst_map_to_asu_isym(True)
  tst_map_to_asu_isym(False)
  print('OK')
