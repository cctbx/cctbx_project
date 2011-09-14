"""
The rhombohedral space groups have two commonly used settings, using
two different basis systems: "hexagonal basis", and "rhombohedral
basis". This confused me (rwgk) for the longest time.

The best way to think about this is to compare to, e.g. the
tetragonal case, where you have two "Bravais types", "tetragonal P"
and "tetragonal I". Now, you can transform the tetragonal I groups to
a primitive setting. It is still the "tetragonal I Bravais type", just
in a primitive setting, which leads to "weird" unit cell parameters
that nobody wants to use.

In the trigonal system you'll find "trigonal P" and "trigonal R".
As in the case of the tetragonal I, you can convert the trigonal R
groups to a primitive setting. However, in contrast to the tetragonal
case, you get "nice" unit cell parameters.

See also: http://en.wikipedia.org/wiki/Bravais_lattice

Note that "trigonal P" and "hexagonal P" are the same Bravais type.
"""

from cctbx import crystal
import sys

def run(args):
  assert len(args) == 0
  #
  print 'This is the "tetragonal I" setting humans like:'
  tetragonal_i = crystal.symmetry(
    unit_cell=(10,10,13,90,90,90),
    space_group_symbol="I4")
  tetragonal_i.show_summary()
  print
  print 'Exact same symmetry, but using a basis system that is not very'
  print 'accessible to humans:'
  cb_op = tetragonal_i.change_of_basis_op_to_primitive_setting()
  print 'Change of basis:', cb_op
  tetragonal_i.change_basis(cb_op=cb_op).show_summary()
  print
  #
  print 'This is the "trigonal R" setting most humans like best:'
  trigonal_r = crystal.symmetry(
    unit_cell=(10,10,13,90,90,120),
    space_group_symbol="R3")
  trigonal_r.show_summary()
  print
  print 'Same symmetry, sometimes preferred:'
  cb_op = trigonal_r.change_of_basis_op_to_primitive_setting()
  print 'Change of basis:', cb_op
  trigonal_r.change_basis(cb_op).show_summary()
  print

if (__name__ == "__main__"):
  run(sys.argv[1:])
