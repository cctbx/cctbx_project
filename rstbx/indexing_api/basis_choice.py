from __future__ import absolute_import, division, print_function
import exceptions
from libtbx.utils import Sorry

from cctbx.uctbx.reduction_base import iteration_limit_exceeded as KGerror
from rstbx.dps_core.cell_assessment import SmallUnitCellVolume
from rstbx.dps_core.basis_choice import HandleCombo as HandleComboBase
from rstbx.indexing_api.tools import AbsenceHandler

class FewSpots(exceptions.Exception): pass

unphysical_cell_cutoff_macromolecular_regime = 100. # Angstrom^3

class HandleCombo(HandleComboBase):
  def __init__(self,ai,combo):
    HandleComboBase.__init__(self,ai,combo,unphysical_cell_cutoff_macromolecular_regime)

  def handle_absences(self):
      Abs = AbsenceHandler()
      while Abs.absence_detected(self.ai.hklobserved()):
        newmat = Abs.correct(self.ai.getOrientation())
        self.ai.setOrientation(newmat)
        self.ai.niggli()

def select_best_combo_of(ai,better_than=0.36,candidates=20,basis=15,spot_sep=0.4,opt_inputs=None):
  """Take the first few candidates of ai.combos().  Search for all combos
     with hkl obs-calc better than a certain fraction limit, then
     handle absences, and return the best one"""
  best_combo = None
  #Establish a minimum acceptable likelihood requirement for the best_combo.
  # This parameter used to be 0.0; but was adjusted upward to 0.30 when preparing
  # figure 4, to force a wider search for good indexing solutions.
  # Explanation:
  #    The ai.model_likelihood function computes what percentage of observed
  #    spots have corresponding predicted spots within 0.5 x random rmsd.
  #    This will tend to be lower than 100% because:
  #      a) the default model mosaicity (0.1deg) is too low
  # If indexing fails on a good image either the best_likelihood must
  #  be lowered or the mosaicity must be increased.  Not yet decided how
  #  to do this for automated operations.--10/09/2003

  base_likelihood = 0.30
  best_likelihood = base_likelihood

  C = ai.combos(basis)
  maxtry = min(candidates,len(C))
  try_counter = 0
  solutions = SolutionTracker()

  for combo in C:
    #print "COMBO: (%d,%d,%d)"%(combo[0],combo[1],combo[2])
    try:
      HC = HandleCombo(ai,combo)

      dev = ai.rmsdev()
      if dev < better_than:
        HC.handle_absences()

        model_likelihood = ai.model_likelihood(spot_sep)
        #printcombo(ai,combo,model_likelihood)

        #  XXXXXX come back to this later
        #while ((best_likelihood<=base_likelihood) and
        #       ai.getMosaicity()<1.5 and
        #       model_likelihood <= best_likelihood):
        #         ai.setMosaicity(ai.getMosaicity() + 0.1)
        #         model_likelihood = ai.model_likelihood(spot_sep)
        #  XXXXXX

        if model_likelihood > best_likelihood:
          best_likelihood = model_likelihood
        this_solution = {'combo':combo,'model_likelihood':model_likelihood,
                          'volume':ai.getOrientation().unit_cell().volume(),"rmsdev":dev,
                          'lattice_likelihood':0.}

        solutions.append(this_solution)
      try_counter+=1
    except (FewSpots) as f:
      #print "COMBO: (%d,%d,%d) rejected on too few spots"%(combo[0],combo[1],combo[2])
      #printcombo(ai,combo)
      continue
    except (SmallUnitCellVolume) as f:
      #print "COMBO: (%d,%d,%d) rejected on small cell volume"%(combo[0],combo[1],combo[2])
      continue
    except (KGerror) as f:
      #print "COMBO: (%d,%d,%d) rejected on Krivy-Gruber iterations"%(combo[0],combo[1],combo[2])
      #printcombo(ai,combo)
      continue
    except (RuntimeError) as f:
      if str(f).find("Iteration limit exceeded")>0: continue
      print("Report this problem to DIALS developers:")
      print("COMBO: (%d,%d,%d) rejected on C++ runtime error"%(combo[0],combo[1],combo[2]))
      #printcombo(ai,combo)
      continue
    except Exception:
      raise
    if solutions.halts():
      return solutions
    if try_counter == maxtry: break
  return solutions

class SolutionTracker:
  def __init__(self):
    self.all_solutions=[]
    self.volume_filtered=[]
  def append(self,item):
    item['serial']=len(self.all_solutions)
    self.all_solutions.append(item)
    self.update_analysis()
  def update_analysis(self):
    self.best_likelihood = max([i['model_likelihood'] for i in self.all_solutions])
    self.close_solutions = [
      i for i in self.all_solutions
        if i['model_likelihood']>=0.9*self.best_likelihood ]
    '''added the following heuristic to accomodate poor diffraction such as
    19198 & 19199:  combos tend to form a big cluster in terms of
    model_likelihood.  In this top cluster, choose a good one that has the
    least primitive unit cell volume.  The other primitiveness test
    (in the paper) doesn't always work if there are too many close spots'''
    self.min_volume = min([i['volume'] for i in self.close_solutions])
      # new method, Jan 2005.  Rely on likelihood, except filter out solutions
      #  where the volume is significantly higher than the minimum:
    self.volume_filtered = [i for i in self.close_solutions if i['volume']<1.25*self.min_volume]
    self.best_volume_filtered_likelihood = max([i['model_likelihood'] for i in self.volume_filtered])
  def best_combo(self):
    if len(self.all_solutions)==0: return None
    combo = [i for i in self.volume_filtered if i['model_likelihood']==
            self.best_volume_filtered_likelihood][0]
    #print "In this round the best combo was",combo,"with likelihood",self.best_volume_filtered_likelihood
    return combo
  def best_lattice(self):
    if len(self.volume_filtered)==0: return None
    max_lattice = max([ele['lattice_likelihood'] for ele in self.volume_filtered])
    combo = [i for i in self.volume_filtered if i['lattice_likelihood']==
            max_lattice][0]
    return combo
  def halts(self):
    return len(self.volume_filtered)>=5

def increase_mosaicity(pd,ai,verbose=1):
  '''the parameter dictionary must be revised if the combo search determined
      that the minimum necessary mosaicity needed to be raised'''
  if 'mosaicity' in pd:
    old_mosaicity = float(pd['mosaicity'])
    new_mosaicity = ai.getMosaicity()
    if new_mosaicity>old_mosaicity:
      if verbose:
        print('input mosaicity %4.1f; new value %4.1f'%(old_mosaicity,new_mosaicity))
      pd['mosaicity']='%f'%new_mosaicity


### Limitations of this quick approach
# 1) local optimization of direct beam; no search of nearby local minima
# 2) no optimization of mosaicity
# 3) no rejection of second lattice
# 4) no check against the raw data for likelihood
# 5)  doesn't allow opt_inputs
# 6) target cell not allowed at present (FIXED)
# 7) no quick refinement of Direction vectors; requires laborious reindexing at present
# 8) no analytical approximation of fringe functions; could use LBFGS if this could be achieved
# 9) tied to refinement of beam vector instead of d0 origin vector (FIXED)
# 10) current supports single panel only (FIXED)

class SelectBasisMetaprocedure:
  def __init__(self,input_index_engine,input_dictionary,horizon_phil,
    opt_rawframes=None,opt_target=False,reduce_target=True):

    from libtbx import adopt_init_args
    adopt_init_args(self,locals())

    if self.horizon_phil.target_cell!=None:

      # change target to primitive centering type
      from cctbx import crystal
      input_symmetry = crystal.symmetry(
        unit_cell=self.horizon_phil.target_cell,
        space_group_symbol="Hall: %s 1" % self.horizon_phil.target_cell_centring_type)
      from cctbx.sgtbx.lattice_symmetry import metric_subgroups
      groups = metric_subgroups(input_symmetry, 0.0,
        enforce_max_delta_for_generated_two_folds=True)
      #groups.show()
      primitive_target_cell = groups.result_groups[-1]["best_subsym"].unit_cell()

      from rstbx.indexing_api.force_cell import force_cell
      best = force_cell(self.input_index_engine,primitive_target_cell)
      try:
        #print "Best score %.1f, triangle %12s"%(best["score"],str(best["triangle"])),best["orientation"].unit_cell()
        self.input_index_engine.setOrientation(best["orientation"])
      except Exception:
        raise Sorry("""Cannot index with the target_cell.  It is possible the target cell is wrong; try indexing
       without one.  It may be necessary to change the beam position, distance, or two theta angle on the
       command line.  See the http://cci.lbl.gov/labelit usage primer.""")
      # originally implemented with default conversion to niggli cell (reduce_target=True).
      # Added this as a configurable option in the context of indexing for sparse
      # nanocrystal stills, since we want to restrain to the originally-input target
      # setting, not necessarily the reduced cell:
      if reduce_target: self.input_index_engine.niggli()
      return

    #initial search
    all_sol = select_best_combo_of(input_index_engine,
                                      better_than=0.36,
                                      candidates=25,
                                      opt_inputs=(self.input_dictionary,opt_rawframes))

    best_combo = all_sol.best_combo()
    if best_combo!=None:
      self.evaluate_combo(best_combo)
      # XXX revisit the question of codecamp maxcell and whether this test should be functional
      #if self.horizon_phil.codecamp.maxcell != None: return
      if opt_rawframes == None or best_combo['lattice_likelihood'] > 3.0:
        return # quick abort to test out indexing 10/16/13

    raise Exception("""No autoindexing solution.
           Possible: incorrect beam center; multiple lattices; too few spots.""")

  def evaluate_combo(self,best_combo):
    increase_mosaicity(self.input_dictionary,self.input_index_engine,verbose=0)
    #print "best combo",best_combo
    HC = HandleCombo(self.input_index_engine,best_combo['combo'])
    HC.handle_absences()
