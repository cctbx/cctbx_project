from __future__ import absolute_import, division, print_function
import sys, math
from cctbx import xray
from cctbx import adptbx
from iotbx import pdb
from libtbx import adopt_init_args
from cctbx.array_family import flex
import cctbx.eltbx.xray_scattering
from cctbx import eltbx
import iotbx.phil
from cctbx import maptbx
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry
from six.moves import zip
from six.moves import range

ias_master_params = iotbx.phil.parse("""\
  b_iso_max = 100.0
    .type = float
  occupancy_min = -1.0
    .type = float
  occupancy_max = 1.5
    .type = float
  ias_b_iso_max = 100.0
    .type = float
  ias_b_iso_min = 0.0
    .type = float
  ias_occupancy_min = 0.01
    .type = float
  ias_occupancy_max = 3.0
    .type = float
  initial_ias_occupancy = 1.0
    .type = float
  build_ias_types = L R B BH
    .type = strings
    .optional = False
  ring_atoms = None
    .type = strings
    .multiple = True
  lone_pair
    .multiple = True
  {
    atom_x = CA
      .type = str
    atom_xo = C
      .type = str
    atom_o = O
      .type = str
  }
  use_map = True
    .type = bool
  build_only = False
    .type = bool
  file_prefix = None
    .type = str
  peak_search_map {
     map_type = *Fobs-Fmodel mFobs-DFmodel
       .type=choice(multi=False)
     grid_step = 0.25
       .type = float
     scaling = *volume sigma
       .type=choice(multi=False)
  }
""")

ias_scattering_dict = {
    "IS6" : eltbx.xray_scattering.gaussian([0.078250],[1.018360],0),
    "IS8" : eltbx.xray_scattering.gaussian([0.113700],[1.007725],0),
    "IS7" : eltbx.xray_scattering.gaussian([0.066930],[0.927570],0),
    "IS4" : eltbx.xray_scattering.gaussian([0.141092],[1.292696],0),
    "IS1" : eltbx.xray_scattering.gaussian([0.224536],[1.314207],0),
    "IS2" : eltbx.xray_scattering.gaussian([0.156768],[1.105028],0),
    "IS3" : eltbx.xray_scattering.gaussian([0.094033],[0.922350],0),
    "IS5" : eltbx.xray_scattering.gaussian([0.198387],[1.375240],0),
    "IS9" : eltbx.xray_scattering.gaussian([-0.29830],[4.232900],0),
    "IS10": eltbx.xray_scattering.gaussian([0.0465],[0.3407],0) }

# IAS types: B = bond, R - ring center, L - lone pair
BH_type = ["IS1","IS2","IS3"]
B_type  = ["IS4","IS5","IS6","IS7","IS8"]
R_type  = ["IS9"]
L_type  = ["IS10"]

class ias_counters(object):
  def __init__(self, iass):
    self.n_ias         = len(iass)
    self.n_ias_b       = 0
    self.n_ias_bh      = 0
    self.n_ias_r       = 0
    self.n_ias_l       = 0
    self.n_bonds_non_h = 0
    for ias in iass:
      if(ias.type == "B"):  self.n_ias_b  += 1
      if(ias.type == "BH"): self.n_ias_bh += 1
      if(ias.type == "L"):  self.n_ias_l  += 1
      if(ias.type == "R"):  self.n_ias_r  += 1

  def show(self, out = None):
    if(out is None): out = sys.stdout
    print("Number of IAS built:", file=out)
    print("   total:             ", self.n_ias, file=out)
    print("   bond (X-X):        ", self.n_ias_b, file=out)
    print("   bond (X-H):        ", self.n_ias_bh, file=out)
    print("   ring centers:      ", self.n_ias_r, file=out)
    print("   lone pairs:        ", self.n_ias_l, file=out)
    assert self.n_ias == self.n_ias_b+self.n_ias_r+self.n_ias_l+self.n_ias_bh

class atom(object):
  def __init__(self, site_cart = None,
                     q         = None,
                     b_iso     = None,
                     name      = None,
                     i_seq     = None,
                     element   = None):
    adopt_init_args(self, locals())

class bond_density(object):
  def __init__(self, data                    = None,
                     dist                    = None,
                     ideal_data              = None,
                     one_dim_point           = None,
                     one_dim_point_predicted = None,
                     peak_data               = None,
                     peak_dist               = None):
    self.data = data
    self.dist = dist
    self.ideal_data = ideal_data
    self.one_dim_point = one_dim_point
    self.one_dim_point_predicted = one_dim_point_predicted
    self.peak_data = peak_data
    self.peak_dist = peak_dist

class ias(object):
  def __init__(self, atom_1              = None,
                     atom_2              = None,
                     atom_3              = None,
                     site_cart_predicted = None,
                     peak_value          = None,
                     peak_position_cart  = None,
                     q                   = None,
                     b_iso               = None,
                     name                = None,
                     type                = None,
                     bond_density        = None,
                     status              = None):
    adopt_init_args(self, locals())

  def deep_copy(self):
    if(self.bond_density is not None):
       bond_density = bond_density(
           data                    = self.bond_density.data,
           dist                    = self.bond_density.dist,
           ideal_data              = self.bond_density.ideal_data,
           one_dim_point           = self.bond_density.one_dim_point,
           one_dim_point_predicted = self.bond_density.one_dim_point_predicted,
           peak_data               = self.bond_density.peak_data,
           peak_dist               = self.bond_density.peak_dist)
    else: bond_density = None
    return ias(atom_1              = self.atom_1,
               atom_2              = self.atom_2,
               atom_3              = self.atom_3,
               site_cart_predicted = self.site_cart_predicted,
               peak_value          = self.peak_value,
               peak_position_cart  = self.peak_position_cart,
               q                   = self.q,
               b_iso               = self.b_iso,
               name                = self.name,
               type                = self.type,
               bond_density        = bond_density,
               status              = self.status)

class extract_ias(object):
  def __init__(self,
               xray_structure,
               bond_proxies_simple,
               pdb_atoms,
               planarity_proxies,
               params,
               log = None):
    if(log is None): self.log = sys.stdout
    self.params = params
    assert xray_structure.scatterers().size() == pdb_atoms.size()
    self.pdb_atoms = pdb_atoms
    self.sites_cart = xray_structure.sites_cart()
    self.u_iso = xray_structure.extract_u_iso_or_u_equiv()
    self.b_iso = self.u_iso*math.pi**2*8.
    self.occupancies = xray_structure.scatterers().extract_occupancies()
    self.iass = []
    for proxy in bond_proxies_simple:
      i, j = proxy.i_seqs
      atom_i = self._get_atom(i)
      atom_j = self._get_atom(j)
      if(atom_i.element in ["H", "D"] or atom_j.element in ["H", "D"]):
         type = "BH"
      else: type = "B"
      self.iass.append( ias(atom_1 = atom_i,
                            atom_2 = atom_j,
                            type   = type) )
    for i_proxy, proxy in enumerate(planarity_proxies):
      if(len(proxy.i_seqs) >= 6):
         i, j = self._is_phe_tyr_ring(i_seqs = proxy.i_seqs)
         if([i,j].count(None)==0):
            self.iass.append( ias(atom_1 = self._get_atom(i),
                                  atom_2 = self._get_atom(j),
                                  type   = "R") )
         elif(len(self.params.ring_atoms)>0):
           for ring_atoms in self.params.ring_atoms:
             i, j = self._is_any_ring(i_seqs=proxy.i_seqs, ring_atoms=ring_atoms)
             if([i,j].count(None)==0):
               self.iass.append( ias(atom_1 = self._get_atom(i),
                                     atom_2 = self._get_atom(j),
                                     type   = "R") )
      if(len(proxy.i_seqs) >= 4):
         i, j, k = self._is_peptide_plane(i_seqs = proxy.i_seqs,
           params = params.lone_pair)
         if([i,j,k].count(None)==0):
            self.iass.append( ias(atom_1 = self._get_atom(i),
                                  atom_2 = self._get_atom(j),
                                  atom_3 = self._get_atom(k),
                                  type   = "L") )

  def _get_atom(self, i_seq):
    pdb_atom = self.pdb_atoms[i_seq]
    return atom(site_cart = self.sites_cart[i_seq],
                q         = self.occupancies[i_seq],
                b_iso     = self.b_iso[i_seq],
                name      = pdb_atom.name.strip(),
                i_seq     = i_seq,
                element   = pdb_atom.element.strip())

  def _is_any_ring(self, i_seqs, ring_atoms):
    counter = 0
    i_best,j_best = None,None
    for i_seq in i_seqs:
      atom_i_name = self.pdb_atoms[i_seq].name.strip()
      if(atom_i_name in ring_atoms):
        counter += 1
    if(counter != len(ring_atoms)): return i_best,j_best
    dist = 0.0
    for i in i_seqs:
      ri = self.pdb_atoms[i].xyz
      atom_i_name = self.pdb_atoms[i].name.strip()
      for j in i_seqs:
        rj = self.pdb_atoms[j].xyz
        atom_j_name = self.pdb_atoms[j].name.strip()
        d_ = math.sqrt((ri[0]-rj[0])**2+(ri[1]-rj[1])**2+(ri[2]-rj[2])**2)
        if(d_ > dist and atom_i_name in ring_atoms and atom_j_name in ring_atoms):
          dist = d_
          i_best,j_best = i,j
    #if 0: print self.pdb_atoms[i_best].name, self.pdb_atoms[j_best].name
    return i_best,j_best

  def _is_phe_tyr_ring(self, i_seqs):
    ring_atoms = ["CZ","CE1","CE2","CD1","CD2","CG"]
    counter = 0
    cz_i_seq, cg_i_seq = None, None
    for i_seq in i_seqs:
        atom_i_name = self.pdb_atoms[i_seq].name.strip()
        if(atom_i_name in ring_atoms):
           counter += 1
        if(atom_i_name == "CZ"): cz_i_seq = i_seq
        if(atom_i_name == "CG"): cg_i_seq = i_seq
    if(counter == 6):
       return cg_i_seq, cz_i_seq
    else:
       return None, None

  def _is_peptide_plane(self, i_seqs, params):
    counter = 0
    ca_i_seq, c_i_seq, o_i_seq = None, None, None
    i_pg = None
    for i_seq in i_seqs:
      atom_i_name = self.pdb_atoms[i_seq].name.strip()
      for i_pg_, pg in enumerate(params):
        if(atom_i_name in [pg.atom_x,pg.atom_xo,pg.atom_o]):
          counter += 1
          #if(i_pg is None): i_pg = i_pg_
          #else: assert i_pg == i_pg_, [i_pg,i_pg_]
        if(atom_i_name == pg.atom_x ): ca_i_seq = i_seq
        if(atom_i_name == pg.atom_xo): c_i_seq  = i_seq
        if(atom_i_name == pg.atom_o ): o_i_seq  = i_seq
    if(counter >= 3):
       return ca_i_seq, c_i_seq, o_i_seq
    else:
       return None, None, None

def ias_site_position(site_i, site_j, alp):
  alp1 = alp+1.0
  return ((site_i[0] + site_j[0]*alp)/alp1,
          (site_i[1] + site_j[1]*alp)/alp1,
          (site_i[2] + site_j[2]*alp)/alp1)

def ias_position_at_lone_pairs(ias, params):
  site_c, site_ca, site_o = None, None, None
  ias_sites, label = None, None
  for a in [ias.atom_1, ias.atom_2, ias.atom_3]:
    for pg in params:
      if(a.name == pg.atom_x ): site_ca = a.site_cart
      if(a.name == pg.atom_xo): site_c  = a.site_cart
      if(a.name == pg.atom_o ): site_o  = a.site_cart
  assert [site_c,site_ca,site_o].count(None) == 0
  dist_co = math.sqrt( (site_c[0]-site_o[0])**2 +
                       (site_c[1]-site_o[1])**2 +
                       (site_c[2]-site_o[2])**2 )
  ias_sites = add_lone_pairs_for_peptyde_o(site_c  = site_c,
                                           site_o  = site_o,
                                           site_ca = site_ca,
                                           dist_co = dist_co)
  return ias_sites

def set_ias_name_and_predicted_position(iass, params):
  phe_ring = ["CZ","CE1","CE2","CD1","CD2","CG"]
  if(params.ring_atoms is not None):
    for ra in params.ring_atoms:
      if(ra is not None):
        phe_ring += ra
  elbow    = ["CA","CB","CG"]
  main_cn  = ["C","N"]
  main_can = ["CA","N","CD"]
  main_cod = ["C","O","OXT"]
  main_cos = ["CZ","OH", "CB","OG1"]
  main_cac = ["CA","C"]
  any_ch = ["H","C"]
  any_nh = ["H","N"]
  any_oh = ["H","O"]
  new_iass = []
  for ias_ in iass:
    ias_site, ias_sites, label = None, None, None
    name_i, name_j = ias_.atom_1.name, ias_.atom_2.name
    site_i, site_j = ias_.atom_1.site_cart, ias_.atom_2.site_cart
    if(ias_.type == "B" or ias_.type == "BH"):
       if(name_i in phe_ring and name_j in phe_ring):
          label = "IS5"
          ias_site = ias_site_position(site_i, site_j, 1.0)
       elif(name_i in elbow and name_j in elbow):
          label = "IS4"
          ias_site = ias_site_position(site_i, site_j, 1.0)
       elif(name_i in main_cn and name_j in main_cn):
          label = "IS6"
          ias_site = ias_site_position(site_i, site_j, 1.0)
       elif(name_i in main_can and name_j in main_can):
          label = "IS6"
          alp = 0.773960 / 0.660850
          if(name_j[0] == "N"):
             ias_site = ias_site_position(site_i, site_j, alp)
          else:
             assert name_j[0] == "C"
             ias_site = ias_site_position(site_j, site_i, alp)
       elif(name_i in main_cod and name_j in main_cod):
          label = "IS8"
          alp = 0.703388 / 0.520963
          if(name_j[0] == "O"):
             ias_site = ias_site_position(site_j, site_i, alp)
          else:
             assert name_j[0] == "C"
             ias_site = ias_site_position(site_i, site_j, alp)
       elif(name_i in main_cos and name_j in main_cos):
          label = "IS7"
          alp = 0.657100 / 0.681920
          if(name_j[0] == "O"):
             ias_site = ias_site_position(site_j, site_i, alp)
          else:
             assert name_j[0] == "C"
             ias_site = ias_site_position(site_i, site_j, alp)
       elif(name_i in main_cac and name_j in main_cac):
          label = "IS4"
          ias_site = ias_site_position(site_i, site_j, 1.0)
       elif(name_i[0] in any_ch and name_j[0] in any_ch and
          [name_i[0],name_j[0]].count("H")==1 ):
          label = "IS1"
          alp = 0.856968 / 0.244483
          if(name_j[0] == "H"):
             ias_site = ias_site_position(site_i, site_j, alp)
          else:
             ias_site = ias_site_position(site_j, site_i, alp)
       elif(name_i[0] in any_nh and name_j[0] in any_nh):
          label = "IS2"
          alp = 0.760400 / 0.267816
          if(name_j[0] == "H"):
             ias_site = ias_site_position(site_i, site_j, alp)
          else:
             ias_site = ias_site_position(site_j, site_i, alp)
       elif(name_i[0] in any_oh and name_j[0] in any_oh):
          label = "IS3"
          alp = 0.716517 / 0.288733
          if(name_j[0] == "H"):
             ias_site = ias_site_position(site_i, site_j, alp)
          else:
             ias_site = ias_site_position(site_j, site_i, alp)
       else:
          label = "IS4"
          ias_site = ias_site_position(site_i, site_j, 1.0)
    if(ias_.type == "L"):
       label = "IS10"
       ias_sites = ias_position_at_lone_pairs(ias_, params.lone_pair)
    if(ias_.type == "R"):
       label = "IS9"
       ias_site = ias_site_position(site_i, site_j, 1.0)
    if([label, ias_site].count(None)==0):
      assert ias_sites is None
      ias_.site_cart_predicted = ias_site
      ias_.name = label
      new_iass.append(ias_)
    if([label, ias_sites].count(None)==0):
      assert ias_site is None
      ias_.site_cart_predicted = ias_sites[0]
      ias_.name = label
      new_iass.append(ias_)
      ias_dc = ias_.deep_copy()
      ias_dc.site_cart_predicted = ias_sites[1]
      ias_dc.name = label
      new_iass.append(ias_dc)
  return new_iass

def set_peaks(iass, fmodel, grid_step, map_type, scaling):
  assert grid_step > 0
  fft_map = fmodel.electron_density_map().fft_map(
    map_type = map_type,
    resolution_factor = 1./(int(fmodel.f_obs().d_min()/grid_step)+1))
  if(scaling == "volume"):
     fft_map.apply_volume_scaling()
  if(scaling == "sigma"):
     fft_map.apply_sigma_scaling()
  fft_map_data = fft_map.real_map_unpadded()
  unit_cell = fmodel.f_obs().unit_cell()
  for ias in iass:
    if(ias.type in ["B", "BH"] ):
       bp = find_peak_at_bond(map_data    = fft_map_data,
                              unit_cell   = unit_cell,
                              label       = ias.name,
                              site_cart_1 = ias.atom_1.site_cart,
                              site_cart_2 = ias.atom_2.site_cart)
       ias.peak_value         = bp.peak_value
       ias.peak_position_cart = bp.peak_site_cart
       ias.q                  = bp.q_estimated
       ias.b_iso              = bp.b_estimated
       site_frac_1 = unit_cell.fractionalize(ias.atom_1.site_cart)
       site_frac_predicted = unit_cell.fractionalize(ias.site_cart_predicted)
       one_dim_point_predicted = \
                           unit_cell.distance(site_frac_1, site_frac_predicted)
       ias.bond_density = bond_density(
                             data                    = bp.data,
                             dist                    = bp.dist,
                             ideal_data              = None,
                             one_dim_point           = bp.one_dim_point,
                             one_dim_point_predicted = one_dim_point_predicted,
                             peak_data               = bp.peak_data,
                             peak_dist               = bp.peak_dist)

       ias.status = bp.status

class find_peak_at_bond(object):
  def __init__(self,map_data,unit_cell,label,site_cart_1,site_cart_2,step=0.005):
    x1,y1,z1 = site_cart_1
    x2,y2,z2 = site_cart_2
    self.one_dim_point  = None
    self.peak_value     = None
    self.peak_site_cart = None
    self.status = None
    self.bond_length = math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    alp = 0
    self.data = flex.double()
    self.dist = flex.double()
    self.peak_sites = flex.vec3_double()
    i_seq = 0
    while alp <= 1.0+1.e-6:
      xp = x1+alp*(x2-x1)
      yp = y1+alp*(y2-y1)
      zp = z1+alp*(z2-z1)
      site_frac = unit_cell.fractionalize((xp,yp,zp))
      ed_ = maptbx.eight_point_interpolation(map_data, site_frac)
      self.dist.append( math.sqrt((x1-xp)**2+(y1-yp)**2+(z1-zp)**2) )
      self.data.append(ed_)
      self.peak_sites.append(unit_cell.orthogonalize(site_frac))
      alp += step
      i_seq += 1
    i_seq_left, i_seq_right, max_peak_i_seq = self.find_peak()
    self.b_estimated, self.q_estimated = None, None
    self.a, self.b = None, None
    self.peak_data, self.peak_dist = None, None
    if([i_seq_left, i_seq_right].count(None) == 0):
      self.one_dim_point  = self.dist[max_peak_i_seq]
      self.peak_value     = self.data[max_peak_i_seq]
      self.peak_site_cart = self.peak_sites[max_peak_i_seq]
      self.peak_data = self.data[i_seq_left:i_seq_right+1]
      self.peak_dist = self.dist[i_seq_left:i_seq_right+1]
      assert (self.peak_data < 0.0).count(True) == 0
      origin = self.dist[max_peak_i_seq]

      dist = (self.peak_dist - origin).deep_copy()
      sel = self.peak_data > 0.0
      data = self.peak_data.select(sel)
      dist = dist.select(sel)
      if(data.size() > 0):
         approx_obj = maptbx.one_gaussian_peak_approximation(
                                                data_at_grid_points    = data,
                                                distances              = dist,
                                                use_weights            = False,
                                                optimize_cutoff_radius = False)
         a_real = approx_obj.a_real_space()
         b_real = approx_obj.b_real_space()
         gof = approx_obj.gof()
         self.a = ias_scattering_dict[label].array_of_a()[0]
         self.b = ias_scattering_dict[label].array_of_b()[0]
         self.b_estimated = approx_obj.b_reciprocal_space()-self.b
         self.q_estimated = approx_obj.a_reciprocal_space()/self.a
         #print "%.2f %.2f"%(self.q_estimated, self.b_estimated)
         if(self.b_estimated <= 0.0):
            self.b_estimated = self.b
         if(self.q_estimated <= 0.0):
            self.q_estimated = self.a
    self.set_status()

  def find_peak(self, peak_width = 0.3, offset = 0.1):
    last_index = len(self.dist)-1
    assert approx_equal(self.bond_length, self.dist[last_index])
    peak_width_ = 0.0
    for i_seq, d in enumerate(self.dist):
      peak_width_ += d
      if(peak_width_ >= peak_width): break
    n_points = max(0, i_seq-1)//2+1
    tol = max(0, i_seq-1)
    i_seq_left, i_seq_right, max_peak_i_seq = None, None, None
    peak_i_seqs = []
    for i_seq, data in enumerate(self.data):
      if(i_seq >= n_points and i_seq <= last_index-n_points):
        if(self.dist[i_seq] > offset and
                                   self.bond_length-self.dist[i_seq] > offset):
          data_at_i_seq = self.data[i_seq]
          dist_at_i_seq = self.dist[i_seq]
          max_at_i_seq = True
          for i in range(-n_points, n_points+1):
            if(i != 0):
              j = i_seq+i
              if(j < 0): j=0
              if(j > last_index): j=last_index
              if(self.data[j] > data_at_i_seq):
                max_at_i_seq = False
                break
          if(max_at_i_seq): peak_i_seqs.append(i_seq)
    if(len(peak_i_seqs) > 0):
      max_peak_i_seq = peak_i_seqs[0]
      max_peak = self.data[max_peak_i_seq]
      peaks = []
      for i_seq in peak_i_seqs:
        i_seq_left = i_seq
        max_peak_left = self.data[i_seq]
        counter = 0
        while i_seq_left >= 0:
          if(self.data[i_seq_left] <= max_peak_left):
             max_peak_left = self.data[i_seq_left]
          else: counter += 1
          if(self.data[i_seq_left] <= 0.0 or counter == tol):
            i_seq_left += 1
            i_seq_left += counter
            break
          i_seq_left -= 1
        #assert i_seq_left < i_seq
        i_seq_right = i_seq
        max_peak_right =  self.data[i_seq]
        counter = 0
        while i_seq_right <= len(self.data)-1:
          if(self.data[i_seq_right] <= max_peak_right):
             max_peak_right = self.data[i_seq_right]
          else: counter += 1
          if(self.data[i_seq_right] <= 0.0 or counter == tol):
            i_seq_right -= 1
            i_seq_right -= counter
            break
          i_seq_right += 1
        #assert i_seq_right > i_seq
        if(i_seq_left < i_seq and i_seq_right > i_seq):
          peaks.append((i_seq, i_seq_left, i_seq_right))
      d=-1.e+6
      peak_sel = None
      for peak in peaks:
        d_ = abs(peak[0]-peak[1])+abs(peak[1]-peak[0])
        if(d_ > d):
          peak_sel = peak
          d=d_
      if(peak_sel is not None):
        i_seq_left, i_seq_right, max_peak_i_seq = \
                                          peak_sel[1], peak_sel[2], peak_sel[0]
        assert i_seq_right > max_peak_i_seq
        assert i_seq_left < max_peak_i_seq
        if(i_seq_left < 0): i_seq_left=0
        if(i_seq_right > len(self.data)-1): i_seq_right=len(self.data)-1
        if(i_seq_left >= i_seq_right):
          i_seq_left, i_seq_right, max_peak_i_seq = None, None, None
      else:
        i_seq_left, i_seq_right, max_peak_i_seq = None, None, None
    return i_seq_left, i_seq_right, max_peak_i_seq

  def set_status(self):
    self.status = True
    if(self.one_dim_point is not None):
      d_l = self.one_dim_point
      d_r = self.bond_length - self.one_dim_point
      if(d_l < 0.1 or d_r < 0.1):
         self.status = False
      if([self.a, self.b].count(None) == 0 and (self.a <= 0.0 or self.b <= 0.0)):
         self.status = False
      if([self.q_estimated, self.b_estimated].count(None) == 0 and
         (self.q_estimated <= 0.0 or self.b_estimated <= 0.0)):
         self.status = False
      if((self.data > 0).count(True) == 0):
         self.status = False
    else:
      self.status = False

def set_status(iass, params):
  if (params.build_ias_types is None):
    raise Sorry("build_ias_types must be specified.")
  for ias in iass:
    if(ias.status is not False and ias.type in ["B", "BH"] and
       ias.type in params.build_ias_types):
       is_ok = ias.atom_1.b_iso <= params.b_iso_max and \
               ias.atom_2.b_iso <= params.b_iso_max and \
               ias.atom_1.q <= params.occupancy_max and \
               ias.atom_2.q >= params.occupancy_min
       if(ias.b_iso is not None):
          is_ok = is_ok and ias.b_iso <= params.ias_b_iso_max and \
                            ias.b_iso >  params.ias_b_iso_min
       if(ias.q is not None):
          is_ok = is_ok and ias.q <= params.ias_occupancy_max and \
                            ias.q >  params.ias_occupancy_min
       if(not is_ok): ias.status = False
       else: ias.status = True
       if(ias.status and ias.peak_position_cart is not None):
         check_at_bond_vector(ias)
    elif(ias.status is not False and ias.type not in ["B", "BH"] and
       ias.type in params.build_ias_types):
       is_ok = ias.atom_1.b_iso <= params.b_iso_max and \
               ias.atom_2.b_iso <= params.b_iso_max and \
               ias.atom_1.q <= params.occupancy_max and \
               ias.atom_2.q >= params.occupancy_min
       assert ias.b_iso is None
       assert ias.q is None
       if(not is_ok): ias.status = False
       else: ias.status = True
    else:
       ias.status = False
    assert ias.status is not None

def check_at_bond_vector(ias):
  a1 = ias.atom_1.site_cart
  a2 = ias.atom_2.site_cart
  e  = ias.peak_position_cart
  d = math.sqrt((a1[0]-a2[0])**2 + (a1[1]-a2[1])**2 + (a1[2]-a2[2])**2)
  d1 = math.sqrt((e[0]-a2[0])**2 + (e[1]-a2[1])**2 + (e[2]-a2[2])**2)
  d2 = math.sqrt((a1[0]-e[0])**2 + (a1[1]-e[1])**2 + (a1[2]-e[2])**2)
  assert approx_equal(d, d1+d2, 1.e-4)

class manager(object):
  def __init__(self, geometry,
                     pdb_atoms,
                     xray_structure,
                     params    = None,
                     fmodel    = None,
                     file_name = None,
                     log       = None):
     adopt_init_args(self, locals())
     self.ias_selection = None
     if(log is None): self.log = sys.stdout
     if(self.params is None): self.params = ias_master_params.extract()
     self.geometry.pair_proxies(xray_structure.sites_cart())
     bond_proxies_simple, asu = self.geometry.get_covalent_bond_proxies()
     print("Total number of covalent bonds = ", len(bond_proxies_simple), file=self.log)
     iass = extract_ias(xray_structure       = self.xray_structure,
                        bond_proxies_simple  = bond_proxies_simple,
                        pdb_atoms            = self.pdb_atoms,
                        planarity_proxies    = self.geometry.planarity_proxies,
                        params               = self.params,
                        log                  = self.log).iass
     iass = set_ias_name_and_predicted_position(iass = iass, params = self.params)
     if(self.fmodel is not None):
       set_peaks(iass      = iass,
                 fmodel    = fmodel,
                 grid_step = self.params.peak_search_map.grid_step,
                 map_type  = self.params.peak_search_map.map_type,
                 scaling   = self.params.peak_search_map.scaling)
     if(file_name is not None):
       self.all_bonds(iass = iass, file_name = file_name)
     print("IAS considered: ", file=self.log)
     ias_counters(iass).show(out = self.log)
     set_status(iass, self.params)
     print("IAS selected: ", file=self.log)
     ias_counters(iass).show(out = self.log)
     self.ias_xray_structure = self.iass_as_xray_structure(iass)
     print("IAS scattering dictionary:", file=log)
     self.ias_xray_structure.scattering_type_registry().show(out = self.log)
     if(1):
        self.write_pdb_file(out=self.log)

  def set_ias_selection(self, ias_selection):
    self.ias_selection = ias_selection

  def setup_ias_selection(self):
    ias_size = self.ias_xray_structure.scatterers().size()
    tail = flex.bool(ias_size, True)
    self.ias_selection = flex.bool(self.xray_structure.scatterers().size(),False)
    self.ias_selection.extend(tail)

  def get_ias_selection(self):
    if self.ias_selection is None:
      self.setup_ias_selection()
    return self.ias_selection

  def iass_as_xray_structure(self, iass):
     ias_xray_structure = xray.structure(
                     crystal_symmetry = self.xray_structure.crystal_symmetry())
     for ias in iass:
       assert ias.status is not None
       if(ias.status):
          if(self.params.use_map):
             site_cart = ias.peak_position_cart
             b_iso     = ias.b_iso
             q         = ias.q
          else:
             site_cart = ias.site_cart_predicted
             b_iso     = (ias.atom_1.b_iso + ias.atom_2.b_iso)*0.5
             q         = self.params.initial_ias_occupancy
          if(b_iso is None):
             b_iso = (ias.atom_1.b_iso + ias.atom_2.b_iso)*0.5
          if(q is None):
             q = self.params.initial_ias_occupancy
          if(site_cart is None):
             site_cart = ias.site_cart_predicted
          assert [b_iso, q].count(None) == 0
          if(site_cart is not None):
            ias_scatterer = xray.scatterer(
              label           = ias.name,
              scattering_type = ias.name,
              site = self.xray_structure.unit_cell().fractionalize(site_cart),
              u               = adptbx.b_as_u(b_iso),
              occupancy       = q)
            ias_xray_structure.add_scatterer(ias_scatterer)
     ias_xray_structure.scattering_type_registry(
                                             custom_dict = ias_scattering_dict)
     return ias_xray_structure

  def write_pdb_file(self, out = None):
    if (out is None): out = sys.stdout
    sites_cart = self.ias_xray_structure.sites_cart()
    for i_seq, sc in enumerate(self.ias_xray_structure.scatterers()):
      a = pdb.hierarchy.atom_with_labels()
      a.hetero = True
      a.serial = i_seq+1
      a.name = sc.label[:4] # XXX
      a.resname = "IAS"
      a.resseq = i_seq+1
      a.xyz = sites_cart[i_seq]
      a.occ = sc.occupancy
      a.b = adptbx.u_as_b(sc.u_iso)
      a.element = sc.label[:2] # XXX
      print(a.format_atom_record_group(), file=out)

  def all_bonds(self, iass, file_name):
    sorted = []
    cc = []
    co = []
    cn = []
    xh = []
    other = []
    for ias in iass:
      e1 = ias.atom_1.element
      e2 = ias.atom_2.element
      if(e1 == "C" and e2 == "C"): cc.append(ias)
      elif(e1 in ["C","N"] and e2 in ["C","N"]): cn.append(ias)
      elif(e1 in ["C","O"] and e2 in ["C","O"]): co.append(ias)
      elif(e1 in ["H"] or e2 in ["H"]): xh.append(ias)
      else: other.append(ias)
    sorted = cc+co+cn+xh
    fout = open(file_name, "w")
    fmt1 = "*** %s_%s(%s)-%s_%s(%s)=%s"
    fmt2 = "%10.4f %10.4f"
    for ias in sorted:
      if(ias.bond_density is not None):
         name1  = ias.atom_1.name.strip()
         i_seq1 = str(ias.atom_1.i_seq).strip()
         b_iso1 = str("%.1f"%ias.atom_1.b_iso).strip()
         name2  = ias.atom_2.name.strip()
         i_seq2 = str(ias.atom_2.i_seq).strip()
         b_iso2 = str("%.1f"%ias.atom_2.b_iso).strip()
         s1 = ias.atom_1.site_cart
         s2 = ias.atom_2.site_cart
         dist = str("%.2f"%math.sqrt(
                   (s1[0]-s2[0])**2+(s1[1]-s2[1])**2+(s1[2]-s2[2])**2)).strip()
         print(fmt1%(name1,i_seq1,b_iso1,name2,i_seq2,b_iso2, dist), file=fout)
         for d, e in zip(ias.bond_density.dist, ias.bond_density.data):
           print(fmt2%(d, e), file=fout)
         if(ias.bond_density.peak_dist is not None):
           for d, e in zip(ias.bond_density.peak_dist, ias.bond_density.peak_data):
             print("peak= %10.4f %10.4f"%(d, e), file=fout)
         print("IAS: ", ias.bond_density.one_dim_point, ias.peak_value, file=fout)
         print("status= ", ias.status, file=fout)

def add_lone_pairs_for_peptyde_o(site_c, site_o, site_ca, dist_co,
                                 R = 0.35,
                                 ALPHA = 180.-113.364532,
                                 small = 1.e-6):
  X1,Y1,Z1 = site_c
  X2,Y2,Z2 = site_o
  XC,YC,ZC = site_ca
  R = R * R
  DE = dist_co
  Q1 = DE #Q1=1.24D0
  Q2 = math.cos(ALPHA*math.pi/180)
  A1 = X2-X1
  A2 = XC-X1
  A3 = Y2-Y1
  A4 = YC-Y1
  A5 = Z2-Z1
  A6 = ZC-Z1
  A  =   A3*A6-A4*A5
  B  = -(A1*A6-A2*A5)
  C  =   A1*A4-A2*A3
  D  =  -X1*A-Y1*B-Z1*C
  F1 = -A1*X2-A3*Y2-A5*Z2
  F2 = F1-Q1*Q2*math.sqrt(R)
  F  = F2 # F1
  #assert A != 0.0 and A1 != 0.0
  if(A != 0.0 and A1 != 0.0):
    tmp = (A3/A1-B/A)
    assert tmp != 0.0
    B1=(C/A-A5/A1)/tmp
    B2=(D/A-F/A1)/tmp
    B3=(C+B*B1)/A
    B4=(D+B*B2)/A
    D1=1+B1**2+B3**2
    D2=2.0*B3*(B4+X2)+2.0*B1*(B2-Y2)-2.0*Z2
    D3=(B4+X2)**2+(B2-Y2)**2+Z2**2-R
    P=D2**2-4.0*D1*D3
    assert P >= 0.0
    assert D1 != 0.0
    ZI1=(-D2+math.sqrt(P))/(2.0*D1)
    YI1=B2+B1*ZI1
    XI1=-B4-B3*ZI1
    ZI2=(-D2-math.sqrt(P))/(2.0*D1)
    YI2=B2+B1*ZI2
    XI2=-B4-B3*ZI2
    assert abs(A*XI2+B*YI2+C*ZI2+D) < small
    assert abs(A*XI1+B*YI1+C*ZI1+D) < small
    assert abs(A*X1+B*Y1+C*Z1+D)    < small
    assert abs(A*X2+B*Y2+C*Z2+D)    < small
    assert abs(A*XC+B*YC+C*ZC+D)    < small
    return (XI1,YI1,ZI1), (XI2,YI2,ZI2)
  else:
    return None, None
