import math
from scitbx.array_family import flex
from cctbx import miller
from cctbx import maptbx
from libtbx.utils import Sorry

def Hn(m):
  m_ = m
  sc = math.log(m_.size())
  s = m_>0
  m_ = m_.select(s.iselection())
  m_ = m_/flex.sum(m_)
  return -flex.sum(m_*flex.log(m_))/sc

def Hw(m):
  s = m>0
  m_ = m
  m_ = m_.select(s.iselection())
  return -flex.sum(m_*flex.log(m_))

def Q_X(x,y):
  x = abs(x).data()
  y = abs(y).data()
  delta = flex.abs(x-y)
  return flex.sum(delta*delta)

def scale(x, y):
  assert type(x) == type(y)
  if(type(x) == miller.array):
    x = x.data()
    y = y.data()
  x = flex.abs(x)
  y = flex.abs(y)
  d = flex.sum(y*y)
  if d == 0: return 1
  else:
    return flex.sum(x*y)/d

def r_factor(x,y, use_scale):
  try:
    x = flex.abs(x.data())
    y = flex.abs(y.data())
  except: pass
  sc=1
  if(use_scale): sc = scale(x,y)
  return flex.sum(flex.abs(x-sc*y))/flex.sum(x)

def show_map_stat(m, prefix):
  print "%s (min/max/mean/sum, Hw, Hn): %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f"%(
    prefix, flex.min(m), flex.max(m), flex.mean(m), flex.sum(m), Hw(m), Hn(m))

class run(object) :
  def __init__ (self,
                f,
                f_000             = None,
                lam               = None,
                start_map         = "lde",
                resolution_factor = 0.25,
                mean_density      = 0.375,
                max_iterations    = 2000,
                beta              = 0.9,
                verbose           = False):
    self.start_map      = start_map
    assert start_map in ["flat", "lde", "min_shifted"]
    self.lam            = lam
    self.max_iterations = max_iterations
    self.f_000          = f_000
    self.verbose        = verbose
    self.beta           = beta
    self.f              = f
    # current monitor and optimized functional values
    self.cntr         = None
    self.r            = None
    self.h_n          = None
    self.h_w          = None
    self.Z            = None
    self.scale_kc     = None
    self.a_gd         = None
    self.q_x          = None
    self.q_tot        = None
    self.tp           = None
    self.rho_tilda    = None
    self.f_mem        = None
    self.header_shown = False
    #
    self.crystal_gridding = self.f.crystal_gridding(
      d_min                   = self.f.d_min(),
      resolution_factor       = resolution_factor,
      grid_step               = None,
      symmetry_flags          = None,
      mandatory_factors       = None,
      max_prime               = 5,
      assert_shannon_sampling = True)
    self.n_real = self.crystal_gridding.n_real()
    max_index = [int((i-1)/2.) for i in self.n_real]
    self.N = self.n_real[0]*self.n_real[1]*self.n_real[2]
    self.full_set = self.f.complete_set(max_index=max_index)
    if(verbose):
      print "Resolution factor: %-6.4f"%resolution_factor
      print "  N, n1,n2,n3:",self.N,self.n_real[0],self.n_real[1],self.n_real[2]
      print "Box: "
      print "  resolution: %6.4f"%self.full_set.d_min()
      print "  max. index |h|,|k|,|l|<nreal/2:", max_index
      print "  n.refl.:", self.full_set.indices().size()
    # STEP 1
    if(self.f_000 is None): self.f_000=mean_density*self.f.unit_cell().volume()
    Cobs = 1
    Ca = 0.37
    self.Agd = Ca/self.N
    if(verbose):
      print "unit cell:", self.f.unit_cell().parameters()
      print "Cobs, Ca, Agd, f_000:", Cobs, Ca, self.Agd, "%6.3f"%self.f_000
      print "Cobs/(N*f_000): ",Cobs/(self.N*self.f_000)
      print "memory factor (beta):", self.beta
    self.f = self.f.customized_copy(data=self.f.data()*Cobs/(self.N*self.f_000))
    fft_map = miller.fft_map(
      crystal_gridding     = self.crystal_gridding,
      fourier_coefficients = self.f)
    self.rho_obs = fft_map.real_map_unpadded()
    if(verbose):
      show_map_stat(m = self.rho_obs, prefix = "rho_obs (formula #13)")
    # STEP 2
    self.rho = self.normalize_start_map()
    if(verbose):
      show_map_stat(m = self.rho, prefix = "rho_0 (initial approximation)")
    # STEP 3
    self.iterations()

  def normalize_start_map(self):
    rho = self.rho_obs.deep_copy()
    if(self.start_map == "flat"):
      rho = flex.double(flex.grid(self.n_real), 1./self.N)
    elif(self.start_map == "lde"):
      eps = flex.max(rho)/100.
      selection_nonpositive = rho <= eps
      rho = rho.set_selected(selection_nonpositive, eps)
    elif(self.start_map == "min_shifted"):
      rho = rho - flex.min(rho)
    else: raise Sorry("Invalid initial map modification choice.")
    return rho / flex.sum(rho)

  def write_mtz_file(self, file_name = "me.mtz", column_root_label = "MEM"):
    f_mem = self.full_set.structure_factors_from_map(
      map            = self.rho,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)#.resolution_filter(d_min=1.0)
    mtz_dataset = f_mem.as_mtz_dataset(column_root_label=column_root_label)
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name = file_name)

  def update_metrics(self):
    self.r        = r_factor(self.f, self.f_mem, use_scale=False)
    self.h_n      = Hn(self.rho.deep_copy())
    self.h_w      = Hw(self.rho.deep_copy())
    self.scale_kc = scale(self.f, self.f_mem)
    self.q_x      = Q_X(self.f_mem, self.f)
    self.q_tot    = Hw(self.rho)-self.lam/2*self.q_x*self.N
    self.tp       = flex.sum(self.rho)

  def show(self, verbose=True):
    if(not self.header_shown):
      print "lam:", self.lam
      print "  nit    TP       Hw        Q_X              Qtot         Agd    scFc       Hn     R   *   tpgd      hngd     romingd   romaxgd"
      self.header_shown = True
    self.update_metrics()
    fs = " ".join(["%5d"%self.cntr, "%8.6f"%self.tp, "%9.6f"%self.h_w,\
         "%12.8e"%self.q_x, "%11.6f"%self.q_tot, "%10.6f"%self.Agd,\
         "%6.4f"%self.scale_kc, "%9.6f"%self.h_n, "%6.4f"%self.r,\
         "*", "%8.6f"%self.Z, "%9.6f"%Hn(self.rho_tilda),\
         "%9.6f %9.6f"%(flex.min(self.rho_tilda), flex.max(self.rho_tilda))])
    if(verbose): print fs
    return fs

  def iterations(self):
    self.cntr = 0
    while self.cntr < self.max_iterations:
      self.f_mem = self.f.structure_factors_from_map(
        map            = self.rho,
        use_scale      = False,
        anomalous_flag = False,
        use_sg         = False)
      self.f_mem = self.f_mem.customized_copy(data = self.f_mem.data()/self.N)
      fft_map = miller.fft_map(
        crystal_gridding     = self.crystal_gridding,
        fourier_coefficients = self.f_mem)
      rho_mod = fft_map.real_map_unpadded()
      delta = rho_mod - self.rho_obs #/scale(rho_mod, self.rho_obs)
      self.rho_tilda = self.rho.deep_copy()
      maptbx.compute_mem_iteration(
        rho   = self.rho_tilda,
        delta = delta,
        lam   = self.lam,
        n     = 1,#self.N, # coupled with lam: see iteration formula
        a_gd   = self.Agd)
      self.Z = flex.sum(self.rho_tilda)
      self.rho = (1-self.beta)*self.rho + self.beta*self.rho_tilda
      if(self.verbose and self.cntr==0 or self.cntr%10==0):
        self.show()
      if(self.cntr%25==0): self.Agd = self.Agd/self.Z
      self.cntr += 1
    self.update_metrics()
