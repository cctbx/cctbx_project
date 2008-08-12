from libtbx import itertbx

def exercise_count():
  it = itertbx.count()
  l = []
  for i in xrange(4): l.append(it.next())
  assert l == list(xrange(4))

  it = itertbx.count(4)
  l = []
  for i in xrange(8): l.append(it.next())
  assert l == list(xrange(4,4+8))

def exercise_islice():
  l = [1,6,8,9,3,7,2,0]
  it = itertbx.islice(l, 3)
  assert list(it) == l[0:3]
  it = itertbx.islice(l, 2, 5)
  assert list(it) == l[2:5]
  it = itertbx.islice(l, 1, 7, 2)
  assert list(it) == [ l[1], l[3], l[5] ]

def exercise_step():
  it = itertbx.step(increment=2)
  assert list(itertbx.islice(it,5)) == [0,2,4,6,8]
  it = itertbx.step(2,3)
  assert list(itertbx.islice(it,4)) == [2,5,8,11]

def exercise_izip():
  l1, l2, l3 = [1,2], [3,4,5], [6,7,8,9]
  i = itertbx.izip(l1,l2,l3)
  assert i.next() == (1,3,6)
  assert i.next() == (2,4,7)
  try: i.next()
  except StopIteration: pass
  else: raise AssertionError

def run():
  exercise_count()
  exercise_islice()
  exercise_step()
  exercise_izip()
  print 'OK'

if __name__ == '__main__':
  run()
