import exceptions,math
from scitbx import matrix
from cctbx.uctbx.reduction_base import iteration_limit_exceeded as KGerror
from rstbx.dps_core.cell_assessment import unit_cell_too_small,SmallUnitCellVolume

class AbsenceHandler:
  def __init__(self):
    self.recursion_limit=8

  def absence_detected(self,hkllist):
    self.hkl = hkllist
    self.N   = self.hkl.size()
    self.flag = None

    from cctbx.sgtbx.sub_lattice_tools import generate_matrix
    allGenerators=[]
    # include identity for on_means calculation
    for mod in [6,5,4,3,2]:
        for matS in generate_matrix(mod):
          allGenerators.append(matS)
    self.allGenerators = allGenerators

    for idx,matS in enumerate(allGenerators):
        invS = matS.inverse()

        idx_possible = True
        for miller in [matrix.row(i) for i in hkllist]:
          transformed_miller = miller*invS.transpose()
          #print transformed_miller
          for element in (transformed_miller).elems:
            if element.denominator() > 1:
              idx_possible = False

        #print "transformation",invS.transpose().elems,{True:"is",False:"not"}[idx_possible],"possible"
        if idx_possible:
          print "There are systematic absences. Applying transformation",invS.transpose().elems
          self.cb_op = invS.transpose().inverse()  # not sure if transpose.inverse or just inverse
          self.flag = True
          return 1

    return 0

  def correct(self,orientation):
    print "before", orientation.unit_cell(),orientation.unit_cell().volume()
    print [float(i) for i in self.cb_op.elems]
    orientation.change_basis([float(i) for i in self.cb_op.elems])
    print "after", orientation.unit_cell(),orientation.unit_cell().volume()
    unit_cell_too_small(orientation.unit_cell(),cutoff=25.)
    return orientation

class FewSpots(exceptions.Exception): pass

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
    self.min_volume = min([i['volume'] for i in self.close_solutions])
      # new method, Jan 2005.  Rely on likelihood, except filter out solutions
      #  where the volume is significantly higher than the minimum:
    self.volume_filtered = [i for i in self.close_solutions if i['volume']<1.25*self.min_volume]
    self.best_volume_filtered_likelihood = max([i['model_likelihood'] for i in self.volume_filtered])
  def best_combo(self):
    print "There are %d combos"%(len(self.all_solutions))
    if len(self.all_solutions)==0: return None
    combo = [i for i in self.volume_filtered if i['model_likelihood']==
            self.best_volume_filtered_likelihood][0]
    #print "In this round the best combo was",combo,"with likelihood",self.best_volume_filtered_likelihood
    return combo
  def halts(self):
    return len(self.volume_filtered)>=20

class HandleCombo:
  def __init__(self,ai,combo):
    solns=[ai[combo[i]] for i in [0,1,2]]
    ai.setA(solns)
    unit_cell_too_small(ai.getOrientation().unit_cell(),cutoff=25.)  #backported from Fig 2
    ai.niggli() # reduce cell
    self.ai = ai
    self.combo = combo

  def handle_absences(self):
      Abs = AbsenceHandler()
      if Abs.absence_detected(self.ai.hklobserved()):
        newmat = Abs.correct(self.ai.getOrientation())
        self.ai.setOrientation(newmat)
        self.ai.niggli(cutoff=25.)

def select_best_combo_of(ai,better_than=0.15,candidates=20,basis=15):
  """Take the first few candidates of ai.combos().  Search for all combos
     with hkl obs-calc better than a certain fraction limit, then
     handle absences, and return the best one"""
  best_combo = None

  base_likelihood = 0.30
  best_likelihood = base_likelihood
  C = ai.combos(basis)
  maxtry = min(candidates,len(C))
  try_counter = 0
  solutions = SolutionTracker()
  for x in xrange(ai.n_candidates()):
    D = ai[x]
    print "BC%d"%x,"%.4f %8.2f %8.2f kmax=%2d kval=%5.1f kval2=%5.1f kval3=%5.1f"%(D.real,
    180*D.psi/math.pi, 180.*D.phi/math.pi, D.kmax, D.kval,D.kval2,D.kval3)

  for combo in C:
    #print "COMBO: (%d,%d,%d)"%(combo[0],combo[1],combo[2])
    try:
      HC = HandleCombo(ai,combo)
      #HC.handle_absences()  #might need to add this in later
      if ai.rmsdev() < better_than:
        model_likelihood = 1. - ai.rmsdev() # provisional expression for likelihood
        if model_likelihood > best_likelihood:
          best_likelihood = model_likelihood
        this_solution = {'combo':combo,'model_likelihood':model_likelihood,
                          'volume':ai.getOrientation().unit_cell().volume()}
        solutions.append(this_solution)
      try_counter+=1
    except (FewSpots),f:
      #print "COMBO: (%d,%d,%d) rejected on too few spots"%(combo[0],combo[1],combo[2])
      #printcombo(ai,combo)
      continue
    except (SmallUnitCellVolume),f:
      #print "Small COMBO: (%d,%d,%d) rejected on small cell volume"%(combo[0],combo[1],combo[2])
      continue
    except (KGerror),f:
      #print "KG COMBO: (%d,%d,%d) rejected on Krivy-Gruber iterations"%(combo[0],combo[1],combo[2])
      #printcombo(ai,combo)
      continue
    except (RuntimeError),f:
      if str(f).find("Iteration limit exceeded")>0: continue
      if str(f).find("Corrupt metrical matrix")>=0: continue # colinear or coplanar
      print "Report this problem to LABELIT developers:"
      print "COMBO: (%d,%d,%d) rejected on C++ runtime error"%(combo[0],combo[1],combo[2])
      #printcombo(ai,combo)
      continue
    except:
      raise
    if solutions.halts():
      return solutions
    if try_counter == maxtry: break
  return solutions

class SelectBasisMetaprocedure:
  def __init__(self,input_index_engine):
    self.input_index_engine = input_index_engine

    #initial search
    all_sol = select_best_combo_of(input_index_engine,
                                   better_than=0.36,
                                   candidates=25)
    best_combo = all_sol.best_combo()
    print "Best combo",best_combo
    if best_combo!=None:
      self.evaluate_combo(best_combo)
      return

    raise exceptions.Exception("AutoIndexing Failed to Select Basis")

  def evaluate_combo(self,best_combo,verbose=False):
    #print "best combo",best_combo
    HC = HandleCombo(self.input_index_engine,best_combo['combo'])
    HC.handle_absences()

  def show_rms(self):
    print "+++++++++++++++++++++++++++++++"
    print self.input_index_engine.getOrientation().unit_cell()
    print "RMSDEV:",self.input_index_engine.rmsdev()
    print "-------------------------------"

    for hkl,obs in zip(self.input_index_engine.hklobserved(),self.input_index_engine.observed()):
      displace = matrix.col(hkl) - matrix.col(obs)
      diff = math.sqrt(displace.dot(displace))
      print hkl,diff
