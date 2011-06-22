path_random_picks_pickle = (
  "/net/anaconda/scratch1/rwgk/pdbtools/minidb/random_picks.pickle")

import sys, os, time, cPickle

def get_primes(n):
  n = abs(n)
  primes = {}
  if (n != 0):
    while (n % 2 == 0):
      try: primes[2] += 1
      except Exception: primes[2] = 1
      n /= 2
    d = 3
    while (n != 1):
      while (n % d == 0):
        try: primes[d] += 1
        except Exception: primes[d] = 1
        n /= d
      d += 2
  return primes

def are_good_primes(primes):
  for p in primes.keys():
    if (not p in (2,3,5)): return 0
  return 1

def next_best(n):
  while 1:
    primes = get_primes(n)
    if (are_good_primes(primes)): return n
    n += 1

def get_grid(ucparams, resolution):
  return [next_best(int(round(u * 3 / resolution))) for u in ucparams[:3]]

def memory(grid, sizeof_FloatType = 8):
  Nz_complex = grid[2]/2+1
  return grid[0] * grid[1] * 2 * Nz_complex * sizeof_FloatType

def run_one(memlimit, package, iter, type_and_dir):
  f = open(path_random_picks_pickle, "r")
  Records = cPickle.load(f)
  f.close()
  for idCode, data in Records.items():
    depDate, ucparams, resolution = data
    grid = get_grid(ucparams, resolution)
    memuse = memory(grid)
    if (memlimit == 0 or memuse <= memlimit):
      t0 = time.time()
      os.system("tst3d %s %s %d %d %d %d" % (
       (package, type_and_dir) + tuple(grid) + (iter,)))
      t = time.time() - t0
      print idCode, ucparams, resolution, grid, memuse, iter, t
      sys.stdout.flush()

def run(type_and_dir, memlimit = 0):
  for package in ("fftw", "fftpack"):
    for iter in (0, 1):
      t0 = time.time()
      run_one(memlimit, package, iter, type_and_dir)
      print "Time %s %d:" % (package, iter), time.time() - t0
      sys.stdout.flush()

if (__name__ == "__main__"):
  type_and_dir = "rb"
  memlimit = 0
  if (len(sys.argv) > 1):
    type_and_dir = sys.argv[1]
  if (len(sys.argv) > 2):
    memlimit = int(sys.argv[2]) * (1024 * 1024)
  assert type_and_dir in ("cf", "cb", "rf", "rb")
  run(type_and_dir, memlimit)
