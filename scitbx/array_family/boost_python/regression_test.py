def extend_sys_path():
  import sys, os.path
  libtbx_build = os.environ["LIBTBX_BUILD"]
  sys.path.insert(0, os.path.normpath(os.path.join(libtbx_build,
    "scitbx/array_family/boost_python")))

def exercise_std_vector_conversions(verbose=0):
  if (verbose): print 'Checking std::vector conversions'
  assert rt.std_vector((8,9,10)) == 27
  assert rt.std_vector([7,9,10]) == 26
  assert rt.std_vector(iter(xrange(3))) == 3
  assert rt.std_vector(iter(xrange(10,13))) == 33
  assert rt.std_vector((8,9,10,11)) == 38
  try: rt.std_vector((0,1,"c"))
  except TypeError: pass
  else: raise RuntimeError, "TypeError expected."
  if (verbose): print 'OK'

def exercise_std_list_conversions(verbose=0):
  if (verbose): print 'Checking std::list conversions'
  assert rt.std_list((8,9,10)) == 27
  assert rt.std_list([7,9,10]) == 26
  assert rt.std_list(iter(xrange(3))) == 3
  assert rt.std_list(iter(xrange(10,13))) == 33
  assert rt.std_list((8,9,10,11)) == 38
  try: rt.std_list((0,1,"c"))
  except TypeError: pass
  else: raise RuntimeError, "TypeError expected."
  if (verbose): print 'OK'

def exercise_boost_array_conversions(verbose=0):
  if (verbose): print 'Checking boost:array conversions'
  assert rt.boost_array((8,9,10)) == 27
  assert rt.boost_array([7,9,10]) == 26
  assert rt.boost_array(xrange(3)) == 3
  assert rt.boost_array(xrange(10,13)) == 33
  assert rt.boost_array((8,9,10,11)) == 38
  try: rt.boost_array((8,9,10,11,12))
  except TypeError: pass
  else: raise RuntimeError, "TypeError expected."
  try: rt.boost_array((0, 1, "c"))
  except TypeError: pass
  else: raise RuntimeError, "TypeError expected."
  if (verbose): print 'OK'

def exercise_small_conversions(verbose=0):
  if (verbose): print 'Checking af::small conversions'
  assert rt.small((8,9,10)) == 27
  assert rt.small([7,9,10]) == 26
  assert rt.small(iter(xrange(3))) == 3
  assert rt.small(iter(xrange(10,14))) == 46
  assert rt.small((8,9,10,11)) == 38
  try: rt.small((0,1,"c"))
  except TypeError: pass
  else: raise RuntimeError, "TypeError expected."
  if (verbose): print 'OK'

def exercise_shared_flex_conversions(verbose=0):
  if (verbose): print 'Checking shared<->flex conversions'
  assert tuple(rt.return_shared()) == (3,1,2)
  s = rt.return_shared()
  assert rt.use_shared(s) == 6
  # rt.modify_shared(s) # XXX needs lvalue converter
  # assert tuple(s) == (6,2,4)
  assert rt.use_shared((8,9,10)) == 27
  if (verbose): print 'OK'

def exercise_ref_flex_conversions(verbose=0):
  if (verbose): print 'Checking flex->ref conversions'
  s = rt.return_shared()
  assert rt.use_const_ref(s) == 6
  rt.modify_ref(s)
  assert tuple(s) == (6,2,4)
  if (verbose): print 'OK'

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
    i += 1
  if (not verbose): print 'OK'

if (__name__ == "__main__"):
  extend_sys_path()
  import regression_test_ext as rt
  import sys
  run(sys.argv[1:])
