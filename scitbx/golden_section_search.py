from __future__ import absolute_import, division, print_function
import sys
from six.moves import range

def print_progress_dots(out,ii,n_monitor=10,width=80,offset=5, is_final=False):
  if ii-1 == 0:
    print(" "*offset+"Starting search", file=out)

  if (ii-1)%width==0:
    print("\n"+offset*" ", end=' ', file=out)
  if (ii-1)%n_monitor==0:
    print(".", end=' ', file=out)
    out.flush()
  if is_final:
    print(file=out)
    print(file=out)
    print("    Done after %s iterations"%ii, file=out)
    print(file=out)


def gss(f,a,b,eps=1e-7,N=1000,out=None,monitor_progress=False,n_monitor=1,offset=5):
  if out is None:
    out = sys.stdout

  c = 0.61803398875 #golden ratio
  x1 = c*a + (1-c)*b
  fx1 = f(x1)
  x2 = (1-c)*a + c*b
  fx2 = f(x2)
  for i in range(1,N-2):
     if monitor_progress:
       print_progress_dots(out,i,n_monitor,width=40,offset=offset)

     if fx1 < fx2:
        b = x2
        x2 = x1
        fx2 = fx1
        x1 = c*a + (1-c)*b
        fx1 = f(x1)
     else:
        a = x1
        x1 = x2
        fx1 = fx2
        x2 = (1-c)*a + c*b
        fx2 = f(x2)
     if (abs(b-a) < eps):
       if fx1 < fx2:
         if monitor_progress: print_progress_dots(out,i,n_monitor,width=40,offset=offset, is_final=True)
         return x1
       else:
        if monitor_progress: print_progress_dots(out,i,n_monitor,width=40,offset=offset, is_final=True)
        return x2
  if fx1 < fx2:
    if monitor_progress: print_progress_dots(out,i,n_monitor,width=40,offset=offset, is_final=True)
    return x1
  else:
    if monitor_progress: print_progress_dots(out,i,n_monitor,width=40,offset=offset, is_final=True)
    return x2
