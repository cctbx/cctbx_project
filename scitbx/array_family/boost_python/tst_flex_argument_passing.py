from scitbx.array_family import flex

def exercise():
  a = flex.double()
  test = flex.flex_argument_passing()
  try:
    a = flex.double()
    test.shared_as_reference_fails(a)
  except Exception, ex:
    assert str(type(ex)) == "<class 'Boost.Python.ArgumentError'>"
  a = flex.double()
  test.shared_as_value_fails(a)
  assert len(a) == 0
  a = flex.double()
  test.versa_flex_grid_as_value_fails(a)
  assert len(a) == 0

  a = flex.double()
  test.versa_flex_grid_as_reference_succeeds(a)
  assert list(a) == [1.5, 2.5, 3.5]

  a = flex.double()
  test.easy_versa_flex_grid_as_reference(a)
  assert list(a) == [1.5, 0.5, 2.5, 3.5, 4.5]

def run():
  exercise()
  print 'OK'

if __name__ == '__main__':
  run()
