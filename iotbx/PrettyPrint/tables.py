from __future__ import division
from __future__ import print_function
def binTable(binners, labels):
  pass
  size = binners[0].size()
  print(size,"bins:")
  length = len(labels)

  print("             ", end=' ')
  for x in xrange(length):
   if len(labels)>=x+1 and labels[x] is not None:
     print("%12s"%labels[x], end=' ')
   else:
     print("%12s"%None, end=' ')
  print()

  binlimits=binners[0].getbinlimits()
  for x in xrange(size):
    if x == 0:
      print("       %6.2f"%binlimits[0], end=' ')
    else:
      print("%6.2f %6.2f"%(binlimits[x-1],binlimits[x]), end=' ')

    for list in binners:
      value = list[x]

      if value is not None:
        if -2.0<value<2.0:
          print("%12.5f"%value, end=' ')
        else:
          print("%12.2f"%value, end=' ')
      else:
        print("            ", end=' ')
    print()
