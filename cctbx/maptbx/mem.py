from __future__ import absolute_import, division, print_function
import math
from scitbx.array_family import flex
from cctbx import miller
from cctbx import maptbx
from libtbx.utils import Sorry
import sys
from six.moves import range

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
  except Exception: pass
  sc=1
  if(use_scale): sc = scale(x,y)
  return flex.sum(flex.abs(x-sc*y))/flex.sum(x)

def show_map_stat(m, prefix, out=sys.stdout):
  print("%s (min/max/mean/sum, Hw, Hn): %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f"%(
    prefix, flex.min(m), flex.max(m), flex.mean(m), flex.sum(m), Hw(m), Hn(m)), file=out)

class run(object):
  def __init__(self,
                f,
                f_000             = None,
                lam               = None,
                start_map         = "lde",
                resolution_factor = 0.25,
                mean_density      = 0.375,
                max_iterations    = 2000,
                beta              = 0.9,
                use_modification  = True,
                xray_structure    = None,
                verbose           = False,
                lambda_increment_factor = None,
                convergence_at_r_factor = 0,
                convergence_r_threshold = 0.1,
                detect_convergence = True,
                crystal_gridding  = None,
                use_scale         = True,
                log               = None):
    if (log is None) : log = sys.stdout
    self.log              = log
    self.start_map        = start_map
    assert start_map in ["flat", "lde", "min_shifted"]
    self.lam              = lam
    self.max_iterations   = max_iterations
    self.f_000            = f_000
    self.verbose          = verbose
    self.beta             = beta
    self.f                = f
    self.use_modification = use_modification
    self.meio_obj         = None
    self.xray_structure   = xray_structure
    self.lambda_increment_factor = lambda_increment_factor
    self.convergence_at_r_factor = convergence_at_r_factor
    self.detect_convergence      = detect_convergence
    self.crystal_gridding        = crystal_gridding
    self.use_scale               = use_scale
    self.convergence_r_threshold = convergence_r_threshold
    #
    if(self.f.anomalous_flag()):
       merged = self.f.as_non_anomalous_array().merge_equivalents()
       self.f = merged.array().set_observation_type( self.f )
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
    self.f_mem        = None
    self.header_shown = False
    self.cc           = None
    self.r_factors    = flex.double()
    self.cc_to_answer = flex.double()
    #
    if(self.crystal_gridding is None):
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
    self.f_calc = None
    if(self.xray_structure is not None):
      self.f_calc = self.full_set.structure_factors_from_scatterers(
        xray_structure = self.xray_structure).f_calc()
    if(verbose):
      print("Resolution factor: %-6.4f"%resolution_factor, file=self.log)
      print("  N, n1,n2,n3:",self.N,self.n_real[0],self.n_real[1],self.n_real[2], file=self.log)
      print("Box: ", file=self.log)
      print("  resolution: %6.4f"%self.full_set.d_min(), file=self.log)
      print("  max. index |h|,|k|,|l|<nreal/2:", max_index, file=self.log)
      print("  n.refl.:", self.full_set.indices().size(), file=self.log)
    # STEP 1
    if(self.f_000 is None): self.f_000=mean_density*self.f.unit_cell().volume()
    Cobs = 1
    Ca = 0.37
    self.Agd = Ca/self.N
    if(verbose):
      print("Cobs, Ca, Agd, f_000:", Cobs, Ca, self.Agd, "%6.3f"%self.f_000, file=self.log)
      print("Cobs/(N*f_000): ",Cobs/(self.N*self.f_000), file=self.log)
      print("memory factor (beta):", self.beta, file=self.log)
    self.f = self.f.customized_copy(data=self.f.data()*Cobs/(self.N*self.f_000))
    fft_map = miller.fft_map(
      crystal_gridding     = self.crystal_gridding,
      fourier_coefficients = self.f)
    self.rho_obs = fft_map.real_map_unpadded()
    if(verbose):
      show_map_stat(m = self.rho_obs, prefix = "rho_obs (formula #13)",
        out=self.log)
    # STEP 2
    self.rho = self.normalize_start_map()
    if(verbose):
      show_map_stat(m = self.rho, prefix = "rho_0 (initial approximation)",
        out=self.log)
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

  def map_coefficients(self, d_min=None):
    o=maptbx.non_linear_map_modification_to_match_average_cumulative_histogram(
      map_1 = self.rho, map_2 = self.rho_obs)
    # XXX p1 must be equal to p2 very accurately; sometimes not the case
    # currently
    for x in [0.5,1,1.5]:
      p1 = (o.map_1()>x).count(True)*1./o.map_1().size()
      p2 = (o.map_2()>x).count(True)*1./o.map_2().size()
      print("Cumulative histograms match for two maps: %9.6f %9.6f" %(
        p1,p2), file=self.log)
    f1 = self.full_set.structure_factors_from_map(
      map            = o.map_1(),
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
    f2 = self.f.structure_factors_from_map(
      map            = o.map_2(),
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
    if(d_min is not None):
      f1 = f1.resolution_filter(d_min=d_min)
      f2 = f2.resolution_filter(d_min=d_min)
    return f1, f2

  def write_mtz_file(self, file_name = "me.mtz", column_root_label = "MEM",
      d_min=None):
    f1, f2 = self.map_coefficients(d_min=d_min)
    mtz_dataset = f1.as_mtz_dataset(column_root_label=column_root_label)
    mtz_dataset.add_miller_array(miller_array=f2, column_root_label="ORIG")
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name = file_name)

  def update_metrics(self):
    self.r        = r_factor(self.f, self.f_mem, use_scale=False)
    self.h_n      = self.meio_obj.hn()
    self.h_w      = self.meio_obj.hw()
    self.scale_kc = scale(self.f, self.f_mem)
    self.q_x      = Q_X(self.f_mem, self.f)
    self.q_tot    = self.h_w-self.lam/2*self.q_x*self.N
    self.tp       = self.meio_obj.tp()

  def show(self, verbose=True):
    if(self.verbose and not self.header_shown):
      print("lam:", self.lam, file=self.log)
      print("  nit    TP       Hw        Q_X              Qtot         Agd    scFc       Hn     R    * tpgd       lam        cc", file=self.log)
      self.header_shown = True
    self.update_metrics()
    fs = " ".join(["%5d"%self.cntr, "%8.6f"%self.tp, "%9.6f"%self.h_w,\
         "%12.8e"%self.q_x, "%11.6f"%self.q_tot, "%10.6f"%self.Agd,\
         "%7.4f"%self.scale_kc, "%9.6f"%self.h_n, "%6.4f"%self.r,\
         "*", "%8.6f"%self.Z, "%12.6f"%self.lam])
    if(self.cc is not None): cc = "%7.5f"%self.cc
    else: cc = self.cc
    print(fs, cc, file=self.log)
    return fs

  def is_converged(self, rho_trial):
    result = False
    r = r_factor(self.f, self.f_mem, use_scale=False)
    if(r < self.convergence_r_threshold):
      if(self.xray_structure is None):
        self.r_factors.append(r)
        size = self.r_factors.size()
        if(size>=3):
          tmp = flex.mean(self.r_factors[size-3:])
          if(tmp <= r or r <= self.convergence_at_r_factor): result = True
      else:
        f_mem = self.full_set.structure_factors_from_map(
          map            = rho_trial,
          use_scale      = False,
          anomalous_flag = False,
          use_sg         = False)
        self.cc = f_mem.map_correlation(other = self.f_calc)
        self.cc_to_answer.append(self.cc)
        def max_change_so_far(x):
          result = flex.double()
          if(self.cc_to_answer.size()):
            for i in range(self.cc_to_answer.size()):
              if(i>0):
                result.append(self.cc_to_answer[i]-self.cc_to_answer[i-1])
          return flex.max(result)
        size = self.cc_to_answer.size()
        if(size>=3):
          mcsf = max_change_so_far(x = self.cc_to_answer)
          tmp = flex.mean(self.cc_to_answer[size-3:])
          if(tmp >= self.cc-1.e-6 or mcsf/5 >= self.cc-tmp):
            result = True
    else: self.cc=None
    return result

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
      rho_trial = self.rho.deep_copy()
      self.meio_obj = maptbx.mem_iteration(
        rho_mod,
        self.rho_obs,
        rho_trial,
        self.lam*self.N,
        self.n_real,
        self.Agd,
        self.beta,
        self.use_scale)
      if(self.detect_convergence and self.is_converged(rho_trial=rho_trial)):
        break
      else: self.rho = rho_trial
      self.Z = self.meio_obj.z()
      if(self.verbose): self.show()
      if(self.cntr%25==0): self.Agd = self.Agd/self.Z
      self.cntr += 1
      if(self.lambda_increment_factor is not None):
        self.lam *= self.lambda_increment_factor
    self.update_metrics()
