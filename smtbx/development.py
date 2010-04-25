from cctbx.development import random_structure
from cctbx import sgtbx

class random_xray_structure(random_structure.xray_structure):

  def __init__(self, space_group_info, u_iso_xor_u_aniso=True, **kwds):
    super(random_xray_structure, self).__init__(space_group_info, **kwds)
    if u_iso_xor_u_aniso: return
    if kwds['use_u_iso'] and kwds['use_u_aniso']:
      for sc in self.scatterers():
        sc.flags.set_use_u_iso(True).set_use_u_aniso(True)


class test_case(object):

  class __metaclass__(type):
    def __init__(cls, classname, bases, classdict):
      exercises = []
      for base in bases:
        try:
          exercises.extend(base.exercises)
        except AttributeError:
          pass
      for name, attr in classdict.items():
        if callable(attr) and name.startswith('exercise'):
          exercises.append(attr)
      dsu = [ (ex.__name__, ex) for ex in exercises ]
      dsu.sort()
      cls.exercises = [ ex for foo, ex in dsu ]

  def run(cls, verbose=False, *args, **kwds):
    if verbose: print cls.__name__
    for exercise in cls.exercises:
      if verbose: print "\t%s ... " % exercise.__name__,
      o = cls(*args, **kwds)
      exercise(o)
      if verbose: print "OK"
  run = classmethod(run)
