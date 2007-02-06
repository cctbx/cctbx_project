import os, sys, math
from cctbx import xray
from cctbx import adptbx
from iotbx import pdb
from libtbx import adopt_init_args
from cctbx.array_family import flex
import cctbx.eltbx.xray_scattering
from cctbx import eltbx

# experimental
all = {
    # peptide
    "MCAN": eltbx.xray_scattering.gaussian([0.08],[1.22],0),
    "MCN":  eltbx.xray_scattering.gaussian([0.15],[1.40],0),
    "MCO":  eltbx.xray_scattering.gaussian([0.14],[1.27],0),
    "MCAC": eltbx.xray_scattering.gaussian([0.20],[1.75],0),
    "HHH":  eltbx.xray_scattering.gaussian([0.20],[1.25],0),
    # ring
    "RL": eltbx.xray_scattering.gaussian([0.17],[1.64],0),
    "RIN": eltbx.xray_scattering.gaussian([0.26],[1.90],0) }

class manager(object):
  def __init__(self, geometry,
                     atom_attributes_list,
                     xray_structure,
                     dbe_parameters = None,
                     b_iso_max = 10.0):
     adopt_init_args(self, locals())
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
         if(self.b_iso[i] <= b_iso_max and self.b_iso[j] <= b_iso_max):
            dbe_site, label = self.dbe_position(
                                      name_i = atom_i.name.strip(),
                                      name_j = atom_j.name.strip(),
                                      site_i = flex.double(atom_i.coordinates),
                                      site_j = flex.double(atom_j.coordinates))
            if(label is not None):
               self.atom_indices.append([i,j,n_dbe])
               n_dbe += 1
               dbe_site = (dbe_site[0],dbe_site[1],dbe_site[2])
               dbe_site = self.xray_structure.unit_cell().fractionalize(dbe_site)
               dbe_u_iso = (self.u_iso[i] + self.u_iso[j])*0.5
               dbe_scatterer = xray.scatterer(label     = label,
                                              scattering_type = label,
                                              site      = dbe_site,
                                              u         = dbe_u_iso,
                                              occupancy = 0.5)
               self.dbe_xray_structure.add_scatterer(dbe_scatterer)
         if(1):
            if(label is not None):
               print atom_i.name.strip(), atom_j.name.strip(), label
            else:
               print atom_i.name.strip(), atom_j.name.strip(), None
     print "Number of DBE built:         ", n_dbe
     print "Total number of bonds:       ", n_bonds
     print "Total number of non-H bonds: ", n_non_h

     self.dbe_xray_structure.scattering_type_registry(custom_dict = all)

     reg = self.dbe_xray_structure.scattering_type_registry().as_type_gaussian_dict()
     self.dbe_xray_structure.scattering_type_registry().show()
     print

     reg = self.xray_structure.scattering_type_registry().as_type_gaussian_dict()
     self.xray_structure.scattering_type_registry().show()
     print

  def dbe_position(self, name_i, name_j, site_i, site_j):
    phe_ring = ["CZ","CE1","CE2","CD1","CD2","CG"]
    elbow    = ["CA","CB","CG"]
    main_cn  = ["C","N"]
    main_can = ["CA","N"]
    main_co  = ["C","O","OXT","CZ","OH"]
    main_cac = ["CA","C"]
    main_hhh = ["H"]
    label = None
    if(name_i in phe_ring and name_j in phe_ring):
       label = "RIN"
    elif(name_i in elbow and name_j in elbow):
       label = "RL"
    elif(name_i in main_cn and name_j in main_cn):
       label = "MCN"
    elif(name_i in main_can and name_j in main_can):
       label = "MCAN"
    elif(name_i in main_co and name_j in main_co):
       label = "MCO"
    elif(name_i in main_cac and name_j in main_cac):
       label = "MCAC"
    elif(name_i[0] in main_hhh or name_j[0] in main_hhh):
       label = "HHH"
    if(label is not None):
       dbe_site = (site_i + site_j)*0.5
    else:
       dbe_site = None
    return dbe_site, label

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
    mac = sites_cart.select(~dbe_selection)
    dbe = sites_cart.select(dbe_selection)
    mac_offset = mac.size()
    self.target = 0.0
    self.gradients = flex.vec3_double(sites_cart.size())
    self.number_of_restraints = 0
    for index in self.atom_indices:
      sa = mac[index[0]]
      sb = mac[index[1]]
      sd = mac[index[2]]
      d0 = math.sqrt( (sa[0]-sb[0])**2 + (sa[1]-sb[1])**2 + (sa[2]-sb[2])**2 )
      dad= math.sqrt( (sa[0]-sd[0])**2 + (sa[1]-sd[1])**2 + (sa[2]-sd[2])**2 )
      ddb= math.sqrt( (sb[0]-sd[0])**2 + (sb[1]-sd[1])**2 + (sb[2]-sd[2])**2 )
      delta = dad + ddb - d0
      self.target = delta**2
      g1 = -2.*(sa[0]+sb[0]-2.*sd[0])
      g2 = -2.*(sa[1]+sb[1]-2.*sd[1])
      g3 = -2.*(sa[2]+sb[2]-2.*sd[2])
      self.gradients[mac_offset+index[2]] = [g1,g2,g3]
      self.number_of_restraints += 1
    self.target *= (1./self.number_of_restraints)
    self.gradients *= (1./self.number_of_restraints)
