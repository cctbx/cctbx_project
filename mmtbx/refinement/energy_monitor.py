from __future__ import absolute_import, division, print_function

to_kcal_mol = {'ev':23.0609,
  }

def _print_energy_in_kcal(e, units):
    if units.lower() in to_kcal_mol:
      return '%15.3f %s' % (e*to_kcal_mol[units.lower()], 'kcal/mol')
    else:
      return '%15.3f %s' % (e, units)

def print_energy_in_kcal(ga):
  s=[]
  if ga is None: return s
  for d, e, l, b in ga.energies:
    units=ga.units.lower()
    if d in ['opt', 'bound']: atoms=b
    elif d in ['energy', 'strain']: atoms=l
    s.append('%-12s %s (atoms %4d)  ' % (d,
                                          _print_energy_in_kcal(e, units), atoms))
  return s

class energies(list):
  def __init__(self):
    pass

  def as_string(self, verbose=False):
    # from libtbx import easy_pickle
    # easy_pickle.dump('ga.pickle', self)
    s='QM energies\n'
    for i, gas in enumerate(self):
      t=''
      for j, ga in enumerate(gas):
        rc = print_energy_in_kcal(ga)
        if rc:
          for line in rc:
            t += '    %s\n' % line
      if verbose: print('macro_cycle %d %s' % (i+1,t))
      if i:
        def _add_dE(e1, e2, units):
          s=''
          if e1 and e2:
            if e1[0]==e2[0]:
              if e1[0] in ['opt']:
                b1=e1[3]==e2[3]
              if e1[0] in ['strain', 'energy']:
                b1=e1[2]==e2[2]
              if b1:
                de = e2[1]-e1[1]
                s+='    %-12s %s\n' % ('%s dE' % e2[0],
                                      _print_energy_in_kcal(e2[1]-e1[1],units))
          return s
        e1=e2=None
        for k in range(2):
          if gas[k] and first[k]:
            e2=gas[k].energies
            e1=first[k].energies
            for f1, f2 in zip(e1,e2):
              t+= _add_dE(f1,f2,gas[k].units)
              if verbose: print(i+1,f1,f2,s)
      else:
        first=gas
      if t:
        s+='  Macro cycle %d\n' % (i+1)
        s+=t
    return s

if __name__ == '__main__':
  from libtbx import easy_pickle
  e=easy_pickle.load('ga.pickle')
  rc=e.as_string()
  print(rc)
