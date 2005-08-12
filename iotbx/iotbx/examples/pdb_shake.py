from iotbx import pdb
import iotbx.pdb.interpretation
from cctbx.array_family import flex
import sys

def run():
  stage_1 = pdb.interpretation.stage_1(file_name=sys.argv[1])

  # add random numbers [-0.5,0.5) to coordinates
  new_sites_cart = stage_1.get_sites_cart().as_double()
  new_sites_cart += flex.random_double(size=new_sites_cart.size())-0.5
  new_sites_cart = flex.vec3_double(new_sites_cart)

  # reset B-factors (min=0.01, max=0.21)
  new_u_iso = flex.random_double(size=new_sites_cart.size(), factor=0.2) + 0.01

  stage_1.write_modified(
    out=sys.stdout,
    new_sites_cart=new_sites_cart,
    new_u_iso=new_u_iso)

if (__name__ == "__main__"):
  run()
