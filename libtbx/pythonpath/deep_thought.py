from __future__ import division
from __future__ import print_function
from builtins import range
s = set()
for i in range(1,11):
  for j in range(1,11):
    s.add(i*j)
print(len(s))
