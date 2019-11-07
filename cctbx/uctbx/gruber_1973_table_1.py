from __future__ import absolute_import, division, print_function
from scitbx import matrix

raw_mk2 = """
N
100/110/00m
N
110/010/00m
N
100/1m0/001
110/100/00m
N
100/0m0/101
N
100/010/10m
N
100/0m0/101
N
100/0m0/101
N
100/0m0/101
N
100/1m0/001
N
100/110/00m
N
100/0m0/0mm
100/0m0/101
N
100/110/00m
N
100/0mm/00m
100/110/00m
N
100/0m0/0mm
N
100/0mm/00m
N
100/010/mmm
N
100/mmm/001
N
100/010/01m
100/0mm/0m0
N
100/mmm/010
100/010/mmm
N
100/010/01m
100/0mm/0m0
100/110/111
100/110/00m
N
100/010/01m
N
100/010/mmm
N
100/010/01m
100/110/111
100/110/00m
N
100/110/00m
N
100/0m0/0mm
N
100/0m0/101
N
100/0m0/101
N
100/110/00m
"""

raw_type_conditions = """
01:a=b<=c|a<b=c:q,q,a
02:a=b<=c:q/2,q,a
03:a=b<=c:q-p/2,q,a
04:a<=b<c:q,a,q
05:q<a=b<c|a<b<c:q-p/2,a,q
06:a=b<c:q,a,q/2
07:a=b<c:q,a,p/2
08:a=b<c:q,a,q-p/2
09:a<b=c|q<a<b<c:q-p/2,q,a
10:a<b<=c:b+m*(p-q),p,a
11:q=a<b=c|a<b<c:b,a,q
12:a<b=c:n*p+q,q,a
13:a<b=c:b,a/2,a
14:q<a<b=c:b,q,q
15:q<a<b=c:b,q/2,q
16:q<a<b=c:-b+a-2*q/3,-a+q/3,-a+q/3
17:q<a<b=c:-b+a/2-q/6,-a/2-q/6,-a+q/3
18:q<a<b=c:b,q-p/2,q
19:q<a<b=c:-b+a-p/6-q/2,-a-p/6+q/2,-a+p/3
20:q<a<b=c:b,a-q/2,a
21:q<a<b<c:b,q-p/2,q
22:q<a<b<c:-b+a+p/2-q,-a+p/2,-a-p+q
23:q<a<b<c:b,a-q/2,a
24:a<b<c:b,a/2,a
25:a<b<c:b,a+p-q,p
26:a<b<c:n*p+q,a,q
27:a<b<c:b+m*(p-q),a,p
28:a<b<c:b+p-m*q,a+p-q,a
"""

def get_mk2_sets():
  sets = {}
  k = 0
  set = []
  for line in raw_mk2.split():
    if (line == "N"):
      if (k): sets[k] = set
      k += 1
      set = []
      continue
    assert line[3] == "/"
    assert line[7] == "/"
    m = []
    for e in line.replace("/",""):
      if   (e == "0"): m.append(0)
      elif (e == "1"): m.append(1)
      elif (e == "m"): m.append(-1)
      else:
        raise RuntimeError("Corrupt internal table.")
    m = matrix.sqr(m)
    assert abs(m.determinant()) == 1, (line, m.elems, m.determinant())
    set.append(m)
  assert k == 28
  sets[k] = set
  return sets

class type_conditions(object):

  def __init__(self, ck_types, defks):
    self.ck_types = ck_types
    self.defks = defks

def get_type_conditions():
  result = {}
  k = 0
  for line in raw_type_conditions.split():
    k += 1
    k_ctrl,cks,defks = line.split(":")
    assert int(k_ctrl) == k
    result[k] = type_conditions(cks.split("|"), defks)
  assert k == 28
  return result
