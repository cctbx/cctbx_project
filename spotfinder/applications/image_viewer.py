from __future__ import absolute_import, division, print_function
from spotfinder.applications.wrappers import DistlOrganizer

class Empty: pass

"Later go back and refactor this module and signal_strength to avoid code duplication."
class run_signal_strength_class(DistlOrganizer):

  def __init__(self,params):
    E = Empty()
    E.argv=['Empty']
    E.argv.append(params.distl.image)

    self.verbose = params.distl.verbose
    if params.distl.res.inner!=None:
      params.distl_lowres_limit = params.distl.res.inner
    if params.distl.res.outer!=None:
      params.force_method2_resolution_limit = params.distl.res.outer
      params.distl_highres_limit = params.distl.res.outer

    params.distl_force_binning = False
    params.distl_permit_binning = False
    params.wedgelimit = len(E.argv)
    params.spotfinder_header_tests = False
    DistlOrganizer.__init__(self,verbose = True, argument_module=E,
                         phil_params=params)
    self.S = None # need to initialize determined by class SpotFrame

  def view(self):
    from rstbx.viewer.spotfinder_wrap import spot_wrapper
    spot_wrapper(working_phil=self.phil_params).display(path = self.phil_params.distl.image,
                                     organizer = self)
