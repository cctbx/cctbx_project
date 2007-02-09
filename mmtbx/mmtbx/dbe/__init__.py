import os, sys, math
from cctbx import xray
from cctbx import adptbx
from iotbx import pdb
from libtbx import adopt_init_args
from cctbx.array_family import flex
import cctbx.eltbx.xray_scattering
from cctbx import eltbx
import iotbx.phil

dbe_master_params = iotbx.phil.parse("""\
  b_iso_max = 10.0
    .type = float
  occupancy_min = 0.8
    .type = float
  occupancy_max = 1.5
    .type = float
  initial_dbe_occupancy = 0.5
    .type = float
  restraints {
    selection = D1 D2 D3 D4 D5 D6 D7 D8
      .type = strings
    weight_scale = 10.0
      .type = float
  }
""")

all = {
    "D6": eltbx.xray_scattering.gaussian([0.078250],[1.018360],0),
    "D8": eltbx.xray_scattering.gaussian([0.113700],[1.007725],0),
    "D7": eltbx.xray_scattering.gaussian([0.066930],[0.927570],0),
    "D4": eltbx.xray_scattering.gaussian([0.141092],[1.292696],0),
    "D1": eltbx.xray_scattering.gaussian([0.224536],[1.314207],0),
    "D2": eltbx.xray_scattering.gaussian([0.156768],[1.105028],0),
    "D3": eltbx.xray_scattering.gaussian([0.094033],[0.922350],0),
    "D5": eltbx.xray_scattering.gaussian([0.198387],[1.375240],0) }

class manager(object):
  def __init__(self, geometry,
                     atom_attributes_list,
                     xray_structure,
                     params = dbe_master_params.extract(),
                     log    = None):
     adopt_init_args(self, locals())
     if(log is None): self.log = sys.stdout
     self.sites_cart = self.xray_structure.sites_cart()
     self.u_iso = self.xray_structure.extract_u_iso_or_u_equiv()
     self.b_iso = self.u_iso*math.pi**2*8.
     self.q = self.xray_structure.scatterers().extract_occupancies()
     self.dbe_xray_structure = xray.structure(
                     crystal_symmetry = self.xray_structure.crystal_symmetry())
     bond_proxies_simple = self.geometry.pair_proxies().bond_proxies.simple
     n_bonds = len(bond_proxies_simple)
     n_dbe = 0
     n_non_h = 0
     self.atom_indices = []
     self.target = None
     self.gradients = None
     for proxy in bond_proxies_simple:
         i_seqs = proxy.i_seqs
         i,j = proxy.i_seqs
         atom_i = self.atom_attributes_list[i]
         atom_j = self.atom_attributes_list[j]
         if(atom_i.element.strip() not in ["H","D"] and
            atom_j.element.strip() not in ["H","D"]):
            n_non_h += 1
         if(self.b_iso[i] <= self.params.b_iso_max and
            self.b_iso[j] <= self.params.b_iso_max and
            self.q[i] <= self.params.occupancy_max and
            self.q[j] <= self.params.occupancy_max and
            self.q[i] >= self.params.occupancy_min and
            self.q[j] >= self.params.occupancy_min):
            dbe_site, label = self.dbe_position(
                                      name_i = atom_i.name.strip(),
                                      name_j = atom_j.name.strip(),
                                      site_i = flex.double(atom_i.coordinates),
                                      site_j = flex.double(atom_j.coordinates))
            if(label is not None):
               self.atom_indices.append([i,j,n_dbe])
               n_dbe += 1
               dbe_site = self.xray_structure.unit_cell().fractionalize(
                                         (dbe_site[0],dbe_site[1],dbe_site[2]))
               dbe_u_iso = (self.u_iso[i] + self.u_iso[j])*0.5
               dbe_scatterer = xray.scatterer(
                           label           = label,
                           scattering_type = label,
                           site            = dbe_site,
                           u               = dbe_u_iso,
                           occupancy       = self.params.initial_dbe_occupancy)
               self.dbe_xray_structure.add_scatterer(dbe_scatterer)
     self.restraints_selection = flex.bool(
                            self.dbe_xray_structure.scatterers().size(), False)
     for i_seq, sc in enumerate(self.dbe_xray_structure.scatterers()):
         if(sc.label in self.params.restraints.selection):
            self.restraints_selection[i_seq] = True
     print self.restraints_selection.count(True),self.restraints_selection.count(False)
     print >> log, "Total number of"
     print >> log, "   DBE built:            ", n_dbe
     print >> log, "   covalent bonds:       ", n_bonds
     print >> log, "   non-H covalent bonds: ", n_non_h
     print >> log
     print >> log, "DBE scattering dictionaries:"
     self.dbe_xray_structure.scattering_type_registry(custom_dict = all)
     self.dbe_xray_structure.scattering_type_registry().show()
     if(1):
        self.write_pdb_file()

  def dbe_position(self, name_i, name_j, site_i, site_j):
    phe_ring = ["CZ","CE1","CE2","CD1","CD2","CG"]
    elbow    = ["CA","CB","CG"]
    main_cn  = ["C","N"]
    main_can = ["CA","N"]
    main_cod = ["C","O","OXT"]
    main_cos = ["CZ","OH"]
    main_cac = ["CA","C"]
    any_ch = ["H","C"]
    any_nh = ["H","N"]
    any_oh = ["H","O"]
    label = None
    dbe_site = None
    if(name_i in phe_ring and name_j in phe_ring):
       label = "D5"
       dbe_site = self.dbe_site_position(site_i, site_j, 1.0)
    elif(name_i in elbow and name_j in elbow):
       label = "D4"
       dbe_site = self.dbe_site_position(site_i, site_j, 1.0)
    elif(name_i in main_cn and name_j in main_cn):
       label = "D6"
       dbe_site = self.dbe_site_position(site_i, site_j, 1.0)
    elif(name_i in main_can and name_j in main_can):
       label = "D6"
       alp = 0.773960 / 0.660850
       if(name_j[0] == "N"):
          dbe_site = self.dbe_site_position(site_i, site_j, alp)
       else:
          assert name_j[0] == "C"
          dbe_site = self.dbe_site_position(site_j, site_i, alp)
    elif(name_i in main_cod and name_j in main_cod):
       label = "D8"
       alp = 0.703388 / 0.520963
       if(name_j[0] == "O"):
          dbe_site = self.dbe_site_position(site_j, site_i, alp)
       else:
          assert name_j[0] == "C"
          dbe_site = self.dbe_site_position(site_i, site_j, alp)
    elif(name_i in main_cos and name_j in main_cos):
       label = "D7"
       alp = 0.657100 / 0.681920
       if(name_j[0] == "O"):
          dbe_site = self.dbe_site_position(site_j, site_i, alp)
       else:
          assert name_j[0] == "C"
          dbe_site = self.dbe_site_position(site_i, site_j, alp)
    elif(name_i in main_cac and name_j in main_cac):
       label = "D4"
       dbe_site = self.dbe_site_position(site_i, site_j, 1.0)
    elif(name_i[0] in any_ch and name_j[0] in any_ch):
       label = "D1"
       alp = 0.856968 / 0.244483
       if(name_j[0] == "H"):
          dbe_site = self.dbe_site_position(site_i, site_j, alp)
       else:
          dbe_site = self.dbe_site_position(site_j, site_i, alp)
    elif(name_i[0] in any_nh and name_j[0] in any_nh):
       label = "D2"
       alp = 0.760400 / 0.267816
       if(name_j[0] == "H"):
          dbe_site = self.dbe_site_position(site_i, site_j, alp)
       else:
          dbe_site = self.dbe_site_position(site_j, site_i, alp)
    elif(name_i[0] in any_oh and name_j[0] in any_oh):
       label = "D3"
       alp = 0.716517 / 0.288733
       if(name_j[0] == "H"):
          dbe_site = self.dbe_site_position(site_i, site_j, alp)
       else:
          dbe_site = self.dbe_site_position(site_j, site_i, alp)
    return dbe_site, label

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
    for index in self.atom_indices:
      if(self.restraints_selection[index[2]]):
         sa = mac[index[0]]
         sb = mac[index[1]]
         sd = mac[index[2]]
         d0 = math.sqrt((sa[0]-sb[0])**2 + (sa[1]-sb[1])**2 + (sa[2]-sb[2])**2)
         dad= math.sqrt((sa[0]-sd[0])**2 + (sa[1]-sd[1])**2 + (sa[2]-sd[2])**2)
         ddb= math.sqrt((sb[0]-sd[0])**2 + (sb[1]-sd[1])**2 + (sb[2]-sd[2])**2)
         delta = dad + ddb - d0
         self.target = delta**2
         g1 = -2.*(sa[0]+sb[0]-2.*sd[0])
         g2 = -2.*(sa[1]+sb[1]-2.*sd[1])
         g3 = -2.*(sa[2]+sb[2]-2.*sd[2])
         self.gradients[mac_offset+index[2]] = [g1,g2,g3]
         self.number_of_restraints += 1
    self.target *= (1./self.number_of_restraints)
    self.gradients *= (1./self.number_of_restraints)
