from labelit.command_line.stats_index import AutoIndexOrganizer
from labelit.command_line.stats_index import best_character_to_IndexPrinter
from labelit.command_line.default_param import establish_dict_for_refinement
from labelit.dptbx.autoindex import index_and_refine

class new_horizons_state:
  def __init__(self,horizons_phil,args):
    self.horizons_phil = horizons_phil
    self.organizer = AutoIndexOrganizer(
      verbose = self.horizons_phil.distl.bins.verbose,
      argument_module = args,
      horizons_phil = horizons_phil,
      delegate = self.index_and_integrate)
    self.organizer.process()

  def index_and_integrate(self,frames,files,spotfinder_results):
    self.pd = establish_dict_for_refinement(frames,spotfinder_results)
    #------------------------------------------------------------
    ai,P = index_and_refine(pd = self.pd,
                            rawframes = files,
                            spotfinder_results = spotfinder_results,
                            verbose = False,
                            horizon_phil = self.horizons_phil)
    #------------------------------------------------------------
    if self.horizons_phil.compatibility_allow==False:
      M = best_character_to_IndexPrinter(ai,P,self.pd,True,self.horizons_phil)
    else:
      from labelit.diffraction.compatibility import best_compatibility_to_IndexPrinter
      M = best_compatibility_to_IndexPrinter(ai,P,self.pd,files,
          spotfinder_results)
    #------------------------------------------------------------
    if self.horizons_phil.__dict__.has_key("writer"):
      self.horizons_phil.writer.make_image_plots_detail(
        ai=ai,pd=self.pd,inframes=files,spotfinder_results=spotfinder_results)

    if not self.horizons_phil.index_only:
     if 0:
      from labelit.dps import IntegrateCharacters
      IC = IntegrateCharacters(M,self.pd)
      IC.write_mosflm_matrices()
      IC.find_best()
      IC.show()
     if 1:
      from rstbx.new_horizons.oscillation_shots import IntegrateCharacters
      IC = IntegrateCharacters(M,self.pd,self.horizons_phil,files,
        spotfinder_results)
      IC.find_best()
      IC.show()
    return self.pd

def pre_indexing_validation(horizons_phil):
  if horizons_phil.indexing.data!=[]:
    horizons_phil.wedgelimit=len(horizons_phil.indexing.data)
    # for abutting or near-abutting images it is still necessary to set
    # codecamp.maxcell

def pack_names(horizons_phil):
  class Empty:pass
  E = Empty()
  E.argv=['Empty']
  for x in horizons_phil.indexing.data:
    E.argv.append(x)
  return E

def run_index(horizons_phil):
  pre_indexing_validation(horizons_phil)
  imagefile_arguments = pack_names(horizons_phil)
  S = new_horizons_state(horizons_phil,imagefile_arguments)
