import random

def shaked_structure(xs0, thermal_shift, site_shift):
  xs = xs0.deep_copy_scatterers()
  for a in xs.scatterers():
    if(a.flags.grad_site()):
      a.site = [ x + site_shift for x in a.site ]
    elif(a.flags.use_u_iso() and a.flags.grad_u_iso()):
      a.u_iso += (thermal_shift * random.random())
    elif(a.flags.use_u_aniso() and a.flags.grad_u_aniso()):
      a.u_star = [ u + thermal_shift * random.random() for u in a.u_star ]
  return xs
