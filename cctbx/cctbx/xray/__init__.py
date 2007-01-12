from cctbx.xray.ext import *
from cctbx.xray.observation_types import *
from cctbx.xray.scatterer import *
from cctbx.xray.structure import *
from cctbx.xray import structure_factors
from cctbx.xray.target_functors import *
from cctbx.xray.minimization import *
from cctbx.xray import ext

def set_selected_scatterer_grad_flags(scatterers,
                                      site      = None,
                                      u_iso     = None,
                                      u_aniso   = None,
                                      occupancy = None,
                                      fp        = None,
                                      fdp       = None):
  ext.set_selected_scatterer_grad_flags(scatterers = scatterers,
                                        site       = site,
                                        u_iso      = u_iso,
                                        u_aniso    = u_aniso,
                                        occupancy  = occupancy,
                                        fp         = fp,
                                        fdp        = fdp)
