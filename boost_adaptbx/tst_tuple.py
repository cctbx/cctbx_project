import re
import boost.tuple

def exercise():
  doc = boost.tuple.exercise.__doc__
  doc = re.sub("\s+", "", doc)
  assert re.match("exercise\(\(int\)arg1\)->object:"
                  "C\+\+signature:"
                  "boost::tuples::tuple<int,double,"
                  "(boost::tuples::null_type,)+ boost::tuples::null_type >"
                  "exercise\(int\)",
                  doc,
                  re.X)
  assert boost.tuple.exercise(1) == (2, 0.5)
  assert boost.tuple.exercise(2) == (4, 1)

def run():
  exercise()
  print 'OK'

if __name__ == '__main__':
  run()
