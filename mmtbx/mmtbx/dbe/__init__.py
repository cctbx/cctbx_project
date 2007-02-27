import os, sys, math
from cctbx import xray
from cctbx import adptbx
from iotbx import pdb
from libtbx import adopt_init_args
from cctbx.array_family import flex
import cctbx.eltbx.xray_scattering
from cctbx import eltbx
import iotbx.phil
from cctbx import maptbx

dbe_master_params = iotbx.phil.parse("""\
  b_iso_max = 10.0
    .type = float
  occupancy_min = 0.7
    .type = float
  occupancy_max = 1.5
    .type = float
  initial_dbe_occupancy = 1.0
    .type = float
  build_dbe_types = L R B BH
    .type = strings
  number_of_macro_cycles = 10
    .type = int
  restraints {
    selection = D1 D2 D3 D4 D5 D6 D7 D8 D9 D10
      .type = strings
    weight_scale = 10.0
      .type = float
  }
""")

dbe_scattering_dict = {
    "D6" : eltbx.xray_scattering.gaussian([0.078250],[1.018360],0),
    "D8" : eltbx.xray_scattering.gaussian([0.113700],[1.007725],0),
    "D7" : eltbx.xray_scattering.gaussian([0.066930],[0.927570],0),
    "D4" : eltbx.xray_scattering.gaussian([0.141092],[1.292696],0),
    "D1" : eltbx.xray_scattering.gaussian([0.224536],[1.314207],0),
    "D2" : eltbx.xray_scattering.gaussian([0.156768],[1.105028],0),
    "D3" : eltbx.xray_scattering.gaussian([0.094033],[0.922350],0),
    "D5" : eltbx.xray_scattering.gaussian([0.198387],[1.375240],0),
    "D9" : eltbx.xray_scattering.gaussian([-0.29830],[4.232900],0),
    "D10": eltbx.xray_scattering.gaussian([0.0465],[0.3407],0) }

# DBE types: B = bond, R - ring center, L - lone pair
BH_type = ["D1","D2","D3"]
B_type  = ["D4","D5","D6","D7","D8"]
R_type  = ["D9"]
L_type  = ["D10"]

class dbe_counters(object):
  def __init__(self,
               n_dbe         = 0,
               n_bonds       = 0,
               n_dbe_b       = 0,
               n_dbe_r       = 0,
               n_dbe_l       = 0,
               n_bonds_non_h = 0):
     adopt_init_args(self, locals())

  def show(self, out = None):
    if(out is None): out = sys.stdout
    print >> out, "Total number of"
    print >> out, "   DBE built:            ", self.n_dbe
    print >> out, "      bond:              ", self.n_dbe_b
    print >> out, "      ring centers:      ", self.n_dbe_r
    print >> out, "      lone pairs:        ", self.n_dbe_l
    print >> out, "   covalent bonds:       ", self.n_bonds
    print >> out, "   non-H covalent bonds: ", self.n_bonds_non_h
    assert self.n_dbe == self.n_dbe_b + self.n_dbe_r + self.n_dbe_l

class manager(object):
  def __init__(self, geometry,
                     atom_attributes_list,
                     xray_structure,
                     params = dbe_master_params.extract(),
                     fmodel = None,
                     log    = None):
     adopt_init_args(self, locals())
     if(log is None): self.log = sys.stdout
     if(self.params is None): self.params = dbe_master_params.extract()
     self.dbe_counters = dbe_counters()
     self.sites_cart = self.xray_structure.sites_cart()
     self.u_iso = self.xray_structure.extract_u_iso_or_u_equiv()
     self.b_iso = self.u_iso*math.pi**2*8.
     self.q = self.xray_structure.scatterers().extract_occupancies()
     self.dbe_xray_structure = xray.structure(
                     crystal_symmetry = self.xray_structure.crystal_symmetry())
     bond_proxies_simple = self.geometry.pair_proxies().bond_proxies.simple
     self.dbe_counters.n_bonds = len(bond_proxies_simple)
     self.atom_indices = []
     self.target = None
     self.gradients = None
     ###> Locate DBE type "B" => loop over bond proxies:
     need_B = "B"  in self.params.build_dbe_types
     need_BH= "BH" in self.params.build_dbe_types
     if(need_B or need_BH):
        for proxy in bond_proxies_simple:
            i_seqs = proxy.i_seqs
            i,j = proxy.i_seqs
            atom_i = self.atom_attributes_list[i]
            atom_j = self.atom_attributes_list[j]
            if(atom_i.element.strip() not in ["H","D"] and
               atom_j.element.strip() not in ["H","D"]):
               self.dbe_counters.n_bonds_non_h += 1
            if(self._check_b_and_q(i) and self._check_b_and_q(j)):
               dbe_site, dbe_label = self._dbe_position_B(
                                      name_i = atom_i.name.strip(),
                                      name_j = atom_j.name.strip(),
                                      site_i = flex.double(atom_i.coordinates),
                                      site_j = flex.double(atom_j.coordinates))
               if(dbe_label is not None):
                  self.dbe_counters.n_dbe_b += 1
                  self._add_dbe_scatterer(anchor_i  = i,
                                          anchor_j  = j,
                                          dbe_site  = dbe_site,
                                          dbe_label = dbe_label)
     ###> Locate DBE type "R" and "L" => loop over planarity proxies:
     need_L = "L" in self.params.build_dbe_types
     need_R = "R" in self.params.build_dbe_types
     if(need_L or need_R):
        for proxy in self.geometry.planarity_proxies:
            if(need_R and self._is_phe_tyr_ring(i_seqs = proxy.i_seqs)):
               dbe_site, dbe_label, i, j = self._dbe_position_at_ring_center(
                                                         i_seqs = proxy.i_seqs)
               if(dbe_label is not None):
                  self.dbe_counters.n_dbe_r += 1
                  self._add_dbe_scatterer(anchor_i  = i,
                                          anchor_j  = j,
                                          dbe_site  = dbe_site,
                                          dbe_label = dbe_label)
            if(need_L and self._is_peptide_plane(i_seqs = proxy.i_seqs)):
               results = self._dbe_position_at_lone_pairs(
                                                         i_seqs = proxy.i_seqs)
               for result in results:
                   dbe_site, dbe_label, i = result
                   if(dbe_label is not None):
                      self.dbe_counters.n_dbe_l += 1
                      self._add_dbe_scatterer(anchor_i  = i,
                                              anchor_j  = i,
                                              dbe_site  = dbe_site,
                                              dbe_label = dbe_label)
     self.restraints_selection = flex.std_string(
                           self.dbe_xray_structure.scatterers().size(), "None")
     if(self.params.restraints.selection is not None):
        for i_seq, sc in enumerate(self.dbe_xray_structure.scatterers()):
            if(sc.label in self.params.restraints.selection):
               if(sc.label in B_type):  self.restraints_selection[i_seq] = "B"
               if(sc.label in BH_type): self.restraints_selection[i_seq] = "BH"
               if(sc.label in R_type):  self.restraints_selection[i_seq] = "R"
               if(sc.label in L_type):  self.restraints_selection[i_seq] = "L"
        self.need_restraints = self.restraints_selection.count("None") != \
                               self.restraints_selection.size()
     else: self.need_restraints = False
     self.dbe_counters.show(out = log)
     print >> log
     print >> log, "DBE scattering dictionaries:"
     self.dbe_xray_structure.scattering_type_registry(
                                             custom_dict = dbe_scattering_dict)
     self.dbe_xray_structure.scattering_type_registry().show()
     if(1):
        self.write_pdb_file()

  def _check_b_and_q(self, i_seq):
    return self.b_iso[i_seq] <= self.params.b_iso_max and     \
           self.q[i_seq]     <= self.params.occupancy_max and \
           self.q[i_seq]     >= self.params.occupancy_min

  def _add_dbe_scatterer(self, anchor_i, anchor_j, dbe_site, dbe_label):
    self.atom_indices.append([anchor_i, anchor_j, self.dbe_counters.n_dbe])
    self.dbe_counters.n_dbe += 1
    dbe_scatterer = xray.scatterer(
     label           = dbe_label,
     scattering_type = dbe_label,
     site            = self.xray_structure.unit_cell().fractionalize(dbe_site),
     u               = (self.u_iso[anchor_i] + self.u_iso[anchor_j])*0.5,
     occupancy       = self.params.initial_dbe_occupancy)
    self.dbe_xray_structure.add_scatterer(dbe_scatterer)

  def _is_phe_tyr_ring(self, i_seqs):
    ring_atoms = ["CZ","CE1","CE2","CD1","CD2","CG"]
    counter = 0
    for i_seq in i_seqs:
        atom_i = self.atom_attributes_list[i_seq]
        if(atom_i.name.strip() in ring_atoms and self._check_b_and_q(i_seq)):
           counter += 1
    if(counter == 6): return True
    else:             return False

  def _is_peptide_plane(self, i_seqs):
    peptide_atoms = ["O","C","CA","N"]
    counter = 0
    for i_seq in i_seqs:
        atom_i = self.atom_attributes_list[i_seq]
        if(atom_i.name.strip() in peptide_atoms and self._check_b_and_q(i_seq)):
           counter += 1
    if(counter == 4): return True
    else:             return False

  def _dbe_position_B(self, name_i, name_j, site_i, site_j):
    #
    if(self.fmodel is not None):
       fft_map = self.fmodel.electron_density_map(
                              map_type          = "k*Fobs-n*Fmodel",
                              k                 = 1,
                              n                 = 1,
                              resolution_factor = 1/5.)
       fft_map.apply_volume_scaling()
       fft_map_data = fft_map.real_map_unpadded()
    #
    phe_ring = ["CZ","CE1","CE2","CD1","CD2","CG"]
    elbow    = ["CA","CB","CG"]
    main_cn  = ["C","N"]
    main_can = ["CA","N"]
    main_cod = ["C","O","OXT"]
    main_cos = ["CZ","OH", "CB","OG1"]
    main_cac = ["CA","C"]
    any_ch = ["H","C"]
    any_nh = ["H","N"]
    any_oh = ["H","O"]
    dbe_site, label = None, None
    need_B = "B"  in self.params.build_dbe_types
    need_BH= "BH" in self.params.build_dbe_types
    if(name_i in phe_ring and name_j in phe_ring and need_B):
       label = "D5"
       dbe_site = self.dbe_site_position(site_i, site_j, 1.0)
    elif(name_i in elbow and name_j in elbow and need_B):
       label = "D4"
       dbe_site = self.dbe_site_position(site_i, site_j, 1.0)
    elif(name_i in main_cn and name_j in main_cn and need_B):
       label = "D6"
       dbe_site = self.dbe_site_position(site_i, site_j, 1.0)
    elif(name_i in main_can and name_j in main_can and need_B):
       label = "D6"
       alp = 0.773960 / 0.660850
       if(name_j[0] == "N"):
          dbe_site = self.dbe_site_position(site_i, site_j, alp)
       else:
          assert name_j[0] == "C"
          dbe_site = self.dbe_site_position(site_j, site_i, alp)
    elif(name_i in main_cod and name_j in main_cod and need_B):
       label = "D8"
       alp = 0.703388 / 0.520963
       if(name_j[0] == "O"):
          dbe_site = self.dbe_site_position(site_j, site_i, alp)
       else:
          assert name_j[0] == "C"
          dbe_site = self.dbe_site_position(site_i, site_j, alp)
    elif(name_i in main_cos and name_j in main_cos and need_B):
       label = "D7"
       alp = 0.657100 / 0.681920
       if(name_j[0] == "O"):
          dbe_site = self.dbe_site_position(site_j, site_i, alp)
       else:
          assert name_j[0] == "C"
          dbe_site = self.dbe_site_position(site_i, site_j, alp)
    elif(name_i in main_cac and name_j in main_cac and need_B):
       label = "D4"
       dbe_site = self.dbe_site_position(site_i, site_j, 1.0)
    elif(name_i[0] in any_ch and name_j[0] in any_ch and need_BH):
       label = "D1"
       alp = 0.856968 / 0.244483
       if(name_j[0] == "H"):
          dbe_site = self.dbe_site_position(site_i, site_j, alp)
       else:
          dbe_site = self.dbe_site_position(site_j, site_i, alp)
    elif(name_i[0] in any_nh and name_j[0] in any_nh and need_BH):
       label = "D2"
       alp = 0.760400 / 0.267816
       if(name_j[0] == "H"):
          dbe_site = self.dbe_site_position(site_i, site_j, alp)
       else:
          dbe_site = self.dbe_site_position(site_j, site_i, alp)
    elif(name_i[0] in any_oh and name_j[0] in any_oh and need_BH):
       label = "D3"
       alp = 0.716517 / 0.288733
       if(name_j[0] == "H"):
          dbe_site = self.dbe_site_position(site_i, site_j, alp)
       else:
          dbe_site = self.dbe_site_position(site_j, site_i, alp)
    if(label is not None and self.fmodel is not None):
       site_cart, gof = self.find_dbe(fft_map_data = fft_map_data,
                                      site_cart_1  = site_i,
                                      site_cart_2  = site_j,
                                      step         = 0.01,
                                      label        = label)
       dbe_site = site_cart
    return dbe_site, label

  def _dbe_position_at_ring_center(self, i_seqs):
    ring_atoms = ["CZ","CE1","CE2","CD1","CD2","CG"]
    dbe_site, label, i, j = None,None,None,None
    unit_cell = self.xray_structure.unit_cell()
    dbe_site, label = None, None
    for i_seq in i_seqs:
      for j_seq in i_seqs:
        if(i_seq != j_seq):
           atom_i = self.atom_attributes_list[i_seq]
           atom_j = self.atom_attributes_list[j_seq]
           if(atom_i.name.strip() in ring_atoms and
              atom_j.name.strip() in ring_atoms):
              site_i = atom_i.coordinates
              site_j = atom_j.coordinates
              site_i_frac = unit_cell.fractionalize(atom_i.coordinates)
              site_j_frac = unit_cell.fractionalize(atom_j.coordinates)
              dist = unit_cell.distance(site_i_frac, site_j_frac)
              if(abs(dist - 2.8) < 0.2):
                 i, j = i_seq, j_seq
                 label = "D9"
                 dbe_site = self.dbe_site_position(site_i, site_j, 1.0)
                 break
    return dbe_site, label, i, j

  def _dbe_position_at_lone_pairs(self, i_seqs):
    peptide_plane_atoms = ["O","C","CA","N"]
    unit_cell = self.xray_structure.unit_cell()
    site_c,site_ca,site_o = None,None,None
    dbe_sites, label = None, None
    for i_seq in i_seqs:
        atom_i = self.atom_attributes_list[i_seq]
        if(atom_i.name.strip() == "CA" and site_ca is None):
           site_ca = atom_i.coordinates
        if(atom_i.name.strip() == "C" and site_c is None):
           site_c = atom_i.coordinates
        if(atom_i.name.strip() == "O" and site_o is None):
           site_o = atom_i.coordinates
           i_seq_o = i_seq
    if([site_c,site_ca,site_o].count(None) == 0):
       site_c_frac = unit_cell.fractionalize(site_c)
       site_o_frac = unit_cell.fractionalize(site_o)
       dist = unit_cell.distance(site_c_frac, site_o_frac)
       if(dist < 1.3 and dist > 1.1):
          label = "D10"
          dbe_sites = add_lone_pairs_for_peptyde_o(site_c  = site_c,
                                                   site_o  = site_o,
                                                   site_ca = site_ca,
                                                   dist_co = dist)
    if([dbe_sites, label].count(None) == 0):
       return (dbe_sites[0], label, i_seq_o), (dbe_sites[1], label, i_seq_o)
    else:
       return (None,None,None), (None,None,None)

  def dbe_site_position(self, site_i, site_j, alp):
    alp1 = alp+1.0
    return ((site_i[0] + site_j[0]*alp)/alp1,
            (site_i[1] + site_j[1]*alp)/alp1,
            (site_i[2] + site_j[2]*alp)/alp1)

  def write_pdb_file(self, out = None):
    if(out is None):
       out = sys.stdout
    sites_cart = self.dbe_xray_structure.sites_cart()
    for i_seq, sc in enumerate(self.dbe_xray_structure.scatterers()):
        print >> out, pdb.format_atom_record(
                                    record_name = "HETATM",
                                    serial      = i_seq+1,
                                    name        = sc.label,
                                    resName     = "DBE",
                                    resSeq      = i_seq+1,
                                    site        = sites_cart[i_seq],
                                    occupancy   = sc.occupancy,
                                    tempFactor  = adptbx.u_as_b(sc.u_iso),
                                    element     = sc.label)

  def target_and_gradients(self, sites_cart, dbe_selection):
    dbe_scatterers = self.dbe_xray_structure.scatterers()
    mac = sites_cart.select(~dbe_selection)
    dbe = sites_cart.select(dbe_selection)
    mac_offset = mac.size()
    self.target = 0.0
    self.gradients = flex.vec3_double(sites_cart.size())
    self.number_of_restraints = 0
    if(self.need_restraints):
       for index in self.atom_indices:
         if(self.restraints_selection[index[2]] == "B"):
            sa = mac[index[0]]
            sb = mac[index[1]]
            sd = dbe[index[2]]
            d0 = math.sqrt((sa[0]-sb[0])**2+(sa[1]-sb[1])**2+(sa[2]-sb[2])**2)
            d1 = math.sqrt((sa[0]-sd[0])**2+(sa[1]-sd[1])**2+(sa[2]-sd[2])**2)
            d2 = math.sqrt((sb[0]-sd[0])**2+(sb[1]-sd[1])**2+(sb[2]-sd[2])**2)
            delta = d1 + d2 - d0
            self.target += delta**2
            g1 = -2.*delta*( (sa[0]-sd[0])/d1 + (sb[0]-sd[0])/d2 )
            g2 = -2.*delta*( (sa[1]-sd[1])/d1 + (sb[1]-sd[1])/d2 )
            g3 = -2.*delta*( (sa[2]-sd[2])/d1 + (sb[2]-sd[2])/d2 )
            self.gradients[mac_offset+index[2]] = [g1,g2,g3]
            self.number_of_restraints += 1
         if(self.restraints_selection[index[2]] == "BH"):
            sa = mac[index[0]]
            sb = mac[index[1]]
            sd = dbe[index[2]]
            d0 = math.sqrt((sa[0]-sb[0])**2+(sa[1]-sb[1])**2+(sa[2]-sb[2])**2)
            d1 = math.sqrt((sa[0]-sd[0])**2+(sa[1]-sd[1])**2+(sa[2]-sd[2])**2)
            d2 = math.sqrt((sb[0]-sd[0])**2+(sb[1]-sd[1])**2+(sb[2]-sd[2])**2)
            delta = d1 + d2 - d0
            self.target += delta**2
            g1 = -2.*delta*( (sa[0]-sd[0])/d1 + (sb[0]-sd[0])/d2 )
            g2 = -2.*delta*( (sa[1]-sd[1])/d1 + (sb[1]-sd[1])/d2 )
            g3 = -2.*delta*( (sa[2]-sd[2])/d1 + (sb[2]-sd[2])/d2 )
            self.gradients[mac_offset+index[2]] = [g1,g2,g3]
            self.number_of_restraints += 1
         if(self.restraints_selection[index[2]] == "R"):
            anchor_site = (flex.double(mac[index[0]]) + flex.double(mac[index[1]]))/2.
            sd = dbe[index[2]]
            delta1 = (anchor_site[0]-sd[0])
            delta2 = (anchor_site[1]-sd[1])
            delta3 = (anchor_site[2]-sd[2])
            self.target += (delta1**2 + delta2**2 + delta3**2)
            g1 = -2.*delta1
            g2 = -2.*delta2
            g3 = -2.*delta3
            self.gradients[mac_offset+index[2]] = [g1,g2,g3]
            self.number_of_restraints += 1
         if(self.restraints_selection[index[2]] == "L"):
            assert index[0] == index[1]
            sa = mac[index[0]]
            sd = dbe[index[2]]
            d0 = math.sqrt((sa[0]-sd[0])**2+(sa[1]-sd[1])**2+(sa[2]-sd[2])**2)
            delta = 0.35 - d0
            self.target += delta**2
            g1 =  2.*delta*(sa[0]-sd[0])/d0
            g2 =  2.*delta*(sa[1]-sd[1])/d0
            g3 =  2.*delta*(sa[2]-sd[2])/d0
            self.gradients[mac_offset+index[2]] = [g1,g2,g3]
            self.number_of_restraints += 1
       self.target *= (1./self.number_of_restraints)
       self.gradients *= (1./self.number_of_restraints)

  def find_dbe(self, fft_map_data, site_cart_1, site_cart_2, label, step = 0.01):
    x1,y1,z1 = site_cart_1
    x2,y2,z2 = site_cart_2
    ed = -1.e+6
    result = None
    d = math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    if(label in B_type):
       d_lim = 0.4
    if(label in BH_type):
       d_lim = 0.2
    alp = d_lim / d
    alp_end = (d - d_lim) / d
    data = flex.double()
    dist = flex.double()
    i_seq = 0
    i_max = 0
    gof = 1.e+6
    while alp < alp_end:
      xp = x1+alp*(x2-x1)
      yp = y1+alp*(y2-y1)
      zp = z1+alp*(z2-z1)
      site_frac = self.xray_structure.unit_cell().fractionalize((xp,yp,zp))
      ed_ = maptbx.eight_point_interpolation(fft_map_data, site_frac)
      dist.append( math.sqrt((x1-xp)**2+(y1-yp)**2+(z1-zp)**2) )
      data.append(ed_)
      if(ed_ > ed):
         ed = ed_
         result = self.xray_structure.unit_cell().orthogonalize(site_frac)
         i_max = i_seq
      #print "%10.4f %10.4f" % (alp, ed_)
      alp += step
      i_seq += 1
    origin = dist[i_max]
    dist = dist - origin
    sel = data > 0.0
    data = data.select(sel)
    dist = dist.select(sel)
    if(data.size() > 20):
       approx_obj = maptbx.one_gaussian_peak_approximation(
                                              data_at_grid_points    = data,
                                              distances              = dist,
                                              use_weights            = False,
                                              optimize_cutoff_radius = False)
       a_real = approx_obj.a_real_space()
       b_real = approx_obj.b_real_space()
       gof = approx_obj.gof()
       print "%8.3f %8.3f %6.2f" % (a_real, b_real, gof)
    else:
       print "No peak"
    return result, gof

  def refresh_dbe_positions(self, fmodel, dbe_selection):
    assert fmodel.xray_structure is self.xray_structure
    fmodel_dc = fmodel.deep_copy()
    sites_cart = fmodel_dc.xray_structure.sites_cart()
    unit_cell = self.xray_structure.unit_cell()
    mac_offset = sites_cart.select(~dbe_selection).size()
    scat_types = self.xray_structure.scatterers().extract_scattering_types()
    n_shifted = 0
    fmodel_dc.update_xray_structure(
              xray_structure = fmodel_dc.xray_structure.select(~dbe_selection),
              update_f_calc  = True)
    fft_map = fmodel_dc.electron_density_map(
                                         map_type          = "k*Fobs-n*Fmodel",
                                         k                 = 1,
                                         n                 = 1,
                                         resolution_factor = 1/5.)
    fft_map.apply_volume_scaling()
    fft_map_data = fft_map.real_map_unpadded()
    for index in self.atom_indices:
        dbe_scat_type = scat_types[mac_offset+index[2]]
        if(dbe_scat_type in B_type):
           d_lim = 0.4
        if(dbe_scat_type in BH_type):
           d_lim = 0.2
        if(dbe_scat_type in B_type+BH_type):
           i_atom_i, i_atom_j, i_dbe = index[0], index[1], index[2]
           dbe_site = sites_cart[mac_offset+i_dbe]
           dbe_site_frac = unit_cell.fractionalize(dbe_site)
           site_i_frac = unit_cell.fractionalize(sites_cart[i_atom_i])
           site_j_frac = unit_cell.fractionalize(sites_cart[i_atom_j])
           d1s = unit_cell.distance(site_i_frac, dbe_site_frac)
           d2s = unit_cell.distance(site_j_frac, dbe_site_frac)
           d12 = unit_cell.distance(site_i_frac, site_j_frac)
           new_dbe_site,gof = self.find_dbe(
                                        fft_map_data = fft_map_data,
                                        site_cart_1  = sites_cart[i_atom_i],
                                        site_cart_2  = sites_cart[i_atom_j],
                                        step         = 0.01,
                                        label        = dbe_scat_type)
           new_dbe_site_frac = unit_cell.fractionalize(new_dbe_site)
           d1 = unit_cell.distance(site_i_frac, new_dbe_site_frac)
           d2 = unit_cell.distance(site_j_frac, new_dbe_site_frac)
           shift = unit_cell.distance(new_dbe_site_frac, dbe_site_frac)
           if((d1 > d_lim and d2 > d_lim and gof < 25.0) and (d1s+d2s-d12 > 0.1 or d1s < 0.2 or d2s < 0.2)):
              n_shifted += 1
              fmodel.xray_structure.scatterers()[mac_offset+i_dbe].site = \
                    self.xray_structure.unit_cell().fractionalize(new_dbe_site)
           else:
              if(d1s+d2s-d12 > 0.2 or d1s < 0.15 or d2s < 0.15):
                 n_shifted += 1
                 new_dbe_site = self.dbe_site_position(sites_cart[i_atom_i],
                                                       sites_cart[i_atom_j],1.)
                 fmodel.xray_structure.scatterers()[mac_offset+i_dbe].site = \
                                          unit_cell.fractionalize(new_dbe_site)
    print >> self.log, "n_shifted = ", n_shifted
    fmodel.update_xray_structure(update_f_calc = True)


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
  assert A != 0.0 and A1 != 0.0
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
