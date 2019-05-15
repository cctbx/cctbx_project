from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import libtbx.load_env
import sys
from six.moves import range

def extend_sys_path():
  sys.path.insert(0,
    libtbx.env.under_build("scitbx/array_family/boost_python"))

def exercise_std_vector_conversions(verbose=0):
  if (verbose): print('Checking std::vector conversions')
  assert rt.std_vector((8,9,10)) == 27
  assert rt.std_vector([7,9,10]) == 26
  assert rt.std_vector(iter(range(3))) == 3
  assert rt.std_vector(iter(range(10,13))) == 33
  assert rt.std_vector((8,9,10,11)) == 38
  try: rt.std_vector((0,1,"c"))
  except TypeError: pass
  else: raise RuntimeError("TypeError expected.")
  if (verbose): print('OK')

def exercise_std_list_conversions(verbose=0):
  if (verbose): print('Checking std::list conversions')
  assert rt.std_list((8,9,10)) == 27
  assert rt.std_list([7,9,10]) == 26
  assert rt.std_list(iter(range(3))) == 3
  assert rt.std_list(iter(range(10,13))) == 33
  assert rt.std_list((8,9,10,11)) == 38
  try: rt.std_list((0,1,"c"))
  except TypeError: pass
  else: raise RuntimeError("TypeError expected.")
  if (verbose): print('OK')

def exercise_boost_array_conversions(verbose=0):
  if (verbose): print('Checking boost:array conversions')
  assert rt.boost_array((8,9,10)) == 27
  assert rt.boost_array([7,9,10]) == 26
  assert rt.boost_array(range(3)) == 3
  assert rt.boost_array(range(10,13)) == 33
  assert rt.boost_array((8,9,10,11)) == 38
  try: rt.boost_array((8,9,10,11,12))
  except TypeError: pass
  else: raise RuntimeError("TypeError expected.")
  try: rt.boost_array((0, 1, "c"))
  except TypeError: pass
  else: raise RuntimeError("TypeError expected.")
  if (verbose): print('OK')

def exercise_small_conversions(verbose=0):
  if (verbose): print('Checking af::small conversions')
  assert rt.small((8,9,10)) == 27
  assert rt.small([7,9,10]) == 26
  assert rt.small(iter(range(3))) == 3
  assert rt.small(iter(range(10,14))) == 46
  assert rt.small((8,9,10,11)) == 38
  try: rt.small((0,1,"c"))
  except TypeError: pass
  else: raise RuntimeError("TypeError expected.")
  if (verbose): print('OK')

def exercise_shared_flex_conversions(verbose=0):
  if (verbose): print('Checking shared<->flex conversions')
  assert tuple(rt.make_shared()) == (3,1,2)
  s = rt.make_shared()
  assert rt.use_shared(s) == 6
  # rt.modify_shared(s) # XXX needs lvalue converter
  # assert tuple(s) == (6,2,4)
  assert rt.use_shared((8,9,10)) == 27
  if (verbose): print('OK')

def exercise_ref_flex_conversions(verbose=0):
  if (verbose): print('Checking flex->ref conversions')
  s = rt.make_shared()
  assert rt.use_const_ref(s) == 6
  rt.modify_ref(s)
  assert tuple(s) == (6,2,4)
  assert rt.use_const_ref(None) == 0
  rt.modify_ref(None)
  if (verbose): print('OK')

def exercise_ref_flex_grid_flex_conversions(verbose=0):
  if (verbose): print('Checking flex->ref_flex_grid conversions')
  a = flex.double(flex.grid((2, 3)))
  assert a.accessor() == rt.use_const_ref_flex_grid(a)
  assert a.all() == rt.use_const_ref_flex_grid(a).all()
  if (verbose): print('OK')

def exercise_c_grid_conversions(verbose=0):
  if (verbose): print('Checking flex->ref_c_grid conversions')
  a = flex.double(flex.grid((2, 3)))
  assert a.accessor() == rt.use_const_ref_c_grid_2(a)
  assert a.all() == rt.use_const_ref_c_grid_2(a).all()
  assert a.accessor().all() == rt.use_const_ref_c_grid_padded_2(a).all()
  assert a.accessor().all() == rt.use_const_ref_c_grid_padded_2(a).focus()
  a = flex.double(flex.grid((2, 3, 4)))
  assert a.accessor() == rt.use_const_ref_c_grid_3(a)
  assert a.all() == rt.use_const_ref_c_grid_3(a).all()
  assert a.accessor().all() == rt.use_const_ref_c_grid_padded_3(a).all()
  assert a.accessor().all() == rt.use_const_ref_c_grid_padded_3(a).focus()
  a = flex.double(flex.grid((2, 3, 5)).set_focus((2, 3, 4)))
  assert a.accessor().all() == rt.use_const_ref_c_grid_padded_3(a).all()
  assert a.accessor().focus() == rt.use_const_ref_c_grid_padded_3(a).focus()
  assert a.is_0_based()
  assert a.is_padded()
  b = flex.double(flex.grid((1,2,3), (4,6,8)))
  assert not b.is_0_based()
  assert not b.is_padded()
  for c in [a,b]:
    try: rt.use_const_ref_c_grid_3(c)
    except KeyboardInterrupt: raise
    except Exception as e:
      assert str(e).startswith("Python argument types in")
    else: raise RuntimeError("Boost.Python.ArgumentError expected.")
  if (verbose): print('OK')

def exercise_to_tuple(verbose=0):
  if (verbose): print('Checking to_tuple conversions')
  assert rt.make_boost_int_2(3, 5) == (3, 5)
  assert rt.make_boost_int_2(3) == (3, 2)
  assert rt.make_boost_int_2() == (7, 2)
  if (verbose): print('OK')

def run(args):
  iterations = 100
  if (len(args) > 0):
    iterations = int(args[0])
  verbose = 0
  if (iterations and iterations < 4): verbose = 1
  i = 0
  while (iterations == 0 or i < iterations):
    exercise_std_vector_conversions(verbose)
    exercise_std_list_conversions(verbose)
    exercise_boost_array_conversions(verbose)
    exercise_small_conversions(verbose)
    exercise_shared_flex_conversions(verbose)
    exercise_ref_flex_conversions(verbose)
    exercise_ref_flex_grid_flex_conversions(verbose)
    exercise_c_grid_conversions(verbose)
    exercise_to_tuple(verbose)
    i += 1
  if (not verbose): print('OK')

if (__name__ == "__main__"):
  extend_sys_path()
  import regression_test_ext as rt
  import sys
  run(sys.argv[1:])
