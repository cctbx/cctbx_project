""" Generator for mmtbx.refinement.llgi_phi_table (developer script; needs
scipy). Information fraction phi(sigmaA, sigma) = (E[<E>^2] - sigmaA^2)/(1 - sigmaA^2),
expectation over Ec (Wilson), true E given Ec (sigmaA model), I = |E|^2 + sigma z,
with <E> the exact posterior mean along the model phase. """
import numpy as np
from scipy.special import roots_legendre, roots_hermitenorm, ive, i0e
from cctbx.array_family import flex
from cctbx.xray import ext

def gl(a, b, n):
    x, w = roots_legendre(n); return a + (b - a)*(x + 1)/2, w*(b - a)/2

def composite(breaks, n):
    xs, ws = [], []
    for a, b in zip(breaks[:-1], breaks[1:]):
        x, w = gl(a, b, n); xs.append(x); ws.append(w)
    return np.concatenate(xs), np.concatenate(ws)

def ratio(x, centric):
    return np.tanh(x) if centric else ive(1, x)/ive(0, x)

def phi(sA, sig, centric, nr=10, nt=40, nz=20):
    if sA <= 0: return 0.0
    v = 1 - sA*sA
    b = 0.5 if centric else 1.0
    brk = [0, 1e-3, 1e-2, 3e-2, 0.1, 0.3, 0.6, 1, 1.5, 2, 3, 4.5, 6.5]
    r, wr = composite(brk, nr)
    if centric:
        wr = wr*2*np.exp(-r*r/2)/np.sqrt(2*np.pi)   # c >= 0, both signs
    else:
        wr = wr*2*r*np.exp(-r*r)
    z, wz = roots_hermitenorm(nz); wz = wz/np.sqrt(2*np.pi)
    R_list, T_list, W_list = [], [], []
    for rk, wk in zip(r, wr):
        mean = sA*rk
        if centric:
            t, wt = gl(mean - 9*np.sqrt(v), mean + 9*np.sqrt(v), nt)
            pt = np.exp(-(t - mean)**2/(2*v))/np.sqrt(2*np.pi*v)
        else:
            t, wt = gl(max(0.0, mean - 9*np.sqrt(v)), mean + 9*np.sqrt(v), nt)
            a = 2*sA*rk*t/v
            pt = 2*t/v*np.exp(-(t - mean)**2/v)*i0e(a)
        R_list.append(np.full(nt, rk)); T_list.append(t); W_list.append(wk*wt*pt)
    rr = np.concatenate(R_list); tt = np.concatenate(T_list); ww = np.concatenate(W_list)
    norm = ww.sum()
    if sig == 0:
        kap = 2*b*sA*rr/v
        e = np.abs(tt)*ratio(kap*np.abs(tt), centric)
        return (np.sum(ww*e*e)/norm - sA*sA)/v, norm
    I = (tt[:, None]**2 + sig*z[None, :]).ravel()
    W = (ww[:, None]*wz[None, :]).ravel()
    Ec = np.repeat(rr, nz)
    n = I.size
    res = ext.llgi_exact_evaluate(flex.double(I), flex.double(n, sig), flex.double(Ec),
      flex.double(n, sA), flex.bool(n, centric), flex.double(n, 0.0))
    e = res.e_expected.as_numpy_array()
    return (np.sum(W*e*e)/norm - sA*sA)/v, norm

SIGMAA_GRID = [0.01] + [round(0.05*i, 2) for i in range(1, 20)] + [0.97, 0.98, 0.99, 0.995]
SIGMA_GRID = [0.0, 0.01, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7,
  1.0, 1.4, 2.0, 2.8, 4.0, 5.6, 8.0, 11.0, 16.0, 23.0, 32.0]

def _job(args):
  c, i, j = args
  return (c, i, j, phi(SIGMAA_GRID[i], SIGMA_GRID[j], bool(c))[0])

if __name__ == "__main__":
  # Writes phi_table.json (about 1 hour of CPU; parallelised over 10
  # processes). mfix.py/phi_ref.py reference values are reproduced to
  # <1e-5 (centric) and <3e-4 relative (acentric; phi_ref's own acentric
  # grid is the less accurate of the two -- see the large-sigma limit).
  import json
  from multiprocessing import Pool
  jobs = [(c, i, j) for c in (0, 1) for i in range(len(SIGMAA_GRID))
          for j in range(len(SIGMA_GRID))]
  with Pool(10) as p: res = p.map(_job, jobs, chunksize=4)
  tab = np.zeros((2, len(SIGMAA_GRID), len(SIGMA_GRID)))
  for c, i, j, v in res: tab[c, i, j] = v
  json.dump({"sigmaa": SIGMAA_GRID, "sigma": SIGMA_GRID,
    "phi_acentric": tab[0].tolist(), "phi_centric": tab[1].tolist()},
    open("phi_table.json", "w"))
