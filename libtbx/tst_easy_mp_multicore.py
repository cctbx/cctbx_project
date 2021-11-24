from __future__ import absolute_import, division, print_function

def factorial(n):
  f = 1
  while n>1:
    f*=n
    n-=1
  return f

def main(nproc=2):
  from libtbx import easy_mp
  argss = []
  for i in range(7):
    argss.append([i+1])
  for args, res, err_str in easy_mp.multi_core_run(factorial,
                                                   tuple(argss),
                                                   nproc,
                                                   ):
    if res:
      print (args, res)
    assert not err_str, 'Error: %s' % err_str

if __name__ == '__main__':
  main()
