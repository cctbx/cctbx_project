import boost.tuple

def exercise():
  assert boost.tuple.exercise(1) == (2, 0.5)
  assert boost.tuple.exercise(2) == (4, 1)

def run():
  exercise()
  print 'OK'

if __name__ == '__main__':
  run()
