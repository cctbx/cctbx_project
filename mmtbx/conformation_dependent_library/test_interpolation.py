from __future__ import division
from __future__ import print_function
from builtins import range
import sys


def run():
  outl = ""
  for i in range(90,120):
    for j in range(-180,181):
      grid = get_grid_values("Gly_nonxpro", float(i),float(j))
      print_grid(grid, i, j)
      index = get_index(float(i), float(j))
      print(index, end=' ')
      r = interpolate_2d(grid, index)
      print(r)
      outl += " %f %f\n" % (j, r)
    print(outl)
    outl += "\n"
    for j in range(-180, 181, 10):
      outl += " %f %f\n" % (j,
                          get_db_result(cdl_database["Gly_nonxpro"], (90, j), 2))
    f=file("junk.dat", "wb")
    f.write(outl)
    f.close()
    assert 0

def exercise () :
  values = [
    [0.1, 1.1, 2.1, 1.2],
    [0.2, 2.2, 2.3, 1.5],
    [1.4, 1.5, 3.1, 1.6],
    [1.7, 2, 3, 0] ]
  print(values)
  answers = [2.2,1.5,2.3,3.1]
  i=0
  for x in range(2):
    for y in range(2):
      #x, y = (0, 0)
      r = interpolate_2d(values, (x, y))
      print(x,y,r,answers[i])
      assert abs(r-answers[i])<0.01
      i+=1

if __name__=="__main__":
  exercise()
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
