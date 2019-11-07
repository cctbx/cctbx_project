from __future__ import absolute_import, division, print_function
from mmtbx.monomer_library import pdb_interpretation
from mmtbx import monomer_library
import mmtbx.monomer_library.server
from iotbx import pdb
from cctbx.geometry_restraints.lbfgs import lbfgs as geometry_restraints_lbfgs
import scitbx.lbfgs
import libtbx.load_env
import math
import sys, os
from six.moves import zip


def add_ring_h(site_0,site_1,site_2,d0,alpha,beta):
  #                 H (xh,yh,zh)
  #                 |
  #                 | d0
  #                 |
  #      /_alpha    C (x0,y0,z0)   /_beta
  #                / \
  #            d1 /   \ d2
  #              /     \
  #  (x1,y1,z1) C       C (x2,y2,z2)
  #             |       |
  #             |       |
  #             |       |
  #             C       C
  #              \     /
  #               \   /
  #                \ /
  #                 C
  #
  x0 = site_0[0]
  y0 = site_0[1]
  z0 = site_0[2]
  x1 = site_1[0]
  y1 = site_1[1]
  z1 = site_1[2]
  x2 = site_2[0]
  y2 = site_2[1]
  z2 = site_2[2]
  a1 = x2-x0
  a2 = y2-y0
  a3 = z2-z0
  b1 = x1-x0
  b2 = y1-y0
  b3 = z1-z0
  c1 = b2*a3-a2*b3
  c2 = -(b1*a3-a1*b3)
  c3 = b1*a2-a1*b2
  R = -x1*c1-y1*c2-z1*c3
  d1 = math.sqrt(b1*b1+b2*b2+b3*b3)
  d2 = math.sqrt(a1*a1+a2*a2+a3*a3)
  A = d1*d0*math.cos(alpha*math.pi/180.)
  B = d2*d0*math.cos(beta*math.pi/180.)
  P = -(B+a1*x0+a2*y0+a3*z0)
  Q = -(A+b1*x0+b2*y0+b3*z0)
  t1 = b2-b1*a2/a1
  t2 = b3-b1*a3/a1
  t3 = Q-b1*P/a1
  f1 = c2-c1*a2/a1
  f2 = c3-c1*a3/a1
  f3 = R-c1*P/a1
  zh = (f1*t3/t1-f3)/(f2-f1*t2/t1)
  yh = -1./t1*(zh*t2+t3)
  xh = -1./a1*(a2*yh+a3*zh+P)
  return (xh,yh,zh)


def add_ca_ha1_and_ha2(site_0,site_1,site_2,d0,alpha,beta):
  #            HA1      HA2
  #              \      /
  #            d0 \    / d0
  #                \  /
  #      /_alpha    CA (site_0)   /_beta
  #               /    \
  #           d1 /      \ d2
  #             /        \
  #   (site_1) C          N (site_2)
  #
  x0 = site_0[0]
  y0 = site_0[1]
  z0 = site_0[2]
  x1 = site_1[0]
  y1 = site_1[1]
  z1 = site_1[2]
  x2 = site_2[0]
  y2 = site_2[1]
  z2 = site_2[2]
  a1 = x2-x0
  a2 = y2-y0
  a3 = z2-z0
  b1 = x1-x0
  b2 = y1-y0
  b3 = z1-z0
  d1 = math.sqrt(b1*b1+b2*b2+b3*b3)
  d2 = math.sqrt(a1*a1+a2*a2+a3*a3)
  A = d1*d0*math.cos(alpha*math.pi/180.)
  B = d2*d0*math.cos(beta*math.pi/180.)
  P = -(B+a1*x0+a2*y0+a3*z0)
  Q = -(A+b1*x0+b2*y0+b3*z0)
  B1 = (a3/a1-b3/b1)/(b2/b1-a2/a1)
  B2 = (P/a1-Q/b1)/(b2/b1-a2/a1)
  B3 = a3/a1+a2/a1*B1
  B4 = P/a1+a2/a1*B2
  D1 = 1.+B1*B1+B3*B3
  D2 = 2.*B3*(B4+x0)+2.*B1*(B2-y0)-2.*z0
  D3 = (B4+x0)*(B4+x0)+(B2-y0)*(B2-y0)+z0*z0-d0*d0
  discr = math.sqrt(D2*D2-4.*D1*D3)
  coeff = (-D2 - discr)/(2.*D1)
  xh1 = -B4-B3*coeff
  yh1 = B2+B1*coeff
  zh1 = coeff
  coeff = (-D2 + discr)/(2.*D1)
  xh2 = -B4-B3*coeff
  yh2 = B2+B1*coeff
  zh2 = coeff
  return ((xh1,yh1,zh1),(xh2,yh2,zh2))

def add_ca_ha1_or_ha2(site_0,site_1,site_2,site_3,d0,alpha,beta,gamma):
  #     HA1 (or HA2)    CB (or HA1 or HA2)
  #             \ /_gamma/
  #           d0 \      / d3
  #               \    /
  #      /_alpha    CA (site_0)   /_beta
  #               /    \
  #           d1 /      \ d2
  #             /        \
  #   (site_1) C          N (site_2)
  #
  x0 = site_0[0]
  y0 = site_0[1]
  z0 = site_0[2]
  x1 = site_1[0]
  y1 = site_1[1]
  z1 = site_1[2]
  x2 = site_2[0]
  y2 = site_2[1]
  z2 = site_2[2]
  x3 = site_3[0]
  y3 = site_3[1]
  z3 = site_3[2]
  a1 = x2-x0
  a2 = y2-y0
  a3 = z2-z0
  b1 = x1-x0
  b2 = y1-y0
  b3 = z1-z0
  c1 = x3-x0
  c2 = y3-y0
  c3 = z3-z0
  d1 = math.sqrt(b1*b1+b2*b2+b3*b3)
  d2 = math.sqrt(a1*a1+a2*a2+a3*a3)
  d3 = math.sqrt(c1*c1+c2*c2+c3*c3)
  A = d1*d0*math.cos(alpha*math.pi/180.)
  B = d2*d0*math.cos(beta*math.pi/180.)
  C = d3*d0*math.cos(gamma*math.pi/180.)
  P = -(B+a1*x0+a2*y0+a3*z0)
  Q = -(A+b1*x0+b2*y0+b3*z0)
  R = -(C+c1*x0+c2*y0+c3*z0)
  t1 = b2-b1*a2/a1
  t2 = b3-b1*a3/a1
  t3 = Q-b1*P/a1
  f1 = c2-c1*a2/a1
  f2 = c3-c1*a3/a1
  f3 = R-c1*P/a1
  zh = (f1*t3/t1-f3)/(f2-f1*t2/t1)
  yh = -1./t1*(zh*t2+t3)
  xh = -1./a1*(a2*yh+a3*zh+P)
  return (xh,yh,zh)

def add_hhh(site_0,site_1,d0,alpha,beta):
  x0 = site_0[0]
  y0 = site_0[1]
  z0 = site_0[2]
  x1 = site_1[0]
  y1 = site_1[1]
  z1 = site_1[2]
  a1 = x1-x0
  a2 = y1-y0
  a3 = z1-z0
  d1 = math.sqrt(a1*a1+a2*a2+a3*a3)
  A = d1*d0*math.cos(alpha*math.pi/180.)
  P = -(A+a1*x0+a2*y0+a3*z0)
  b1 = z0+P/a3
  b2 = a1/a3
  b3 = a2/a3
  c1 = 1.+b2*b2
  c2 = 1.+b3*b3
  c3 = 2.*b1*b2-2.*x0
  c4 = 2.*b1*b3-2.*y0
  c5 = 2.*b2*b3
  R = -(d0*d0-b1*b1-x0*x0-y0*y0)
  e1 = c5*c5-4.*c1*c2
  e2 = 2.*c3*c5-4.*c1*c4
  e3 = c3*c3-4.*c1*R
  y1 = (-e2+math.sqrt(e2*e2-4.*e1*e3)) / (2.*e1)
  y2 = (-e2-math.sqrt(e2*e2-4.*e1*e3)) / (2.*e1)
  x1 = -(c3+c5*y1)/(2.*c1)
  x2 = -(c3+c5*y2)/(2.*c1)
  z1 = -(P+a1*x1+a2*y1) / a3
  z2 = -(P+a1*x2+a2*y2) / a3
  return (x1,y1,z1), add_ca_ha1_and_ha2(site_0,site_1,(x1,y1,z1),d0,alpha,beta)
  #print x2,y2,z2
  #add_ca_ha1_and_ha2(site_0,site_1,(x2,y2,z2),d0,alpha,alpha)

def add_arg_like_h(site_0,site_1,site_2,d,alpha,flag):
  F = site_0[0]
  G = site_0[1]
  H = site_0[2]
  a1,a2,a3 = site_1[0]-site_0[0],site_1[1]-site_0[1],site_1[2]-site_0[2]
  b1,b2,b3 = site_2[0]-site_1[0],site_2[1]-site_1[1],site_2[2]-site_1[2]
  c1,c2,c3 = site_2[0]-site_0[0],site_2[1]-site_0[1],site_2[2]-site_0[2]
  d1 = math.sqrt(a1*a1+a2*a2+a3*a3)
  A = a2/a1
  B = a3/a1
  C = (-d*d1*math.cos(alpha*math.pi/180.)-F*a1-G*a2-H*a3)/a1
  AA = A*A
  BB = B*B
  CC = C*C
  FF = F*F
  GG = G*G
  HH = H*H
  M = d*d
  f1 = a2*c3-a3*c2
  f2 = -(a1*c3-c1*a3)
  f3 = a1*c2-c1*a2 ;
  f4 = -F*f1 - G*f2 - H*f3
  P = f2/f1
  Q = f3/f1
  R = f4/f1
  PP = P*P
  QQ = Q*Q
  RR = R*R
  xh1 = (-2.*(A-P)*(C*(PP+QQ)+BB*(P*(G-F*P)+R)+AA*(Q*(H-F*Q)+R)-A*(C*P+B*H*P+  \
       B*G*Q-2.*B*F*P*Q+H*P*Q-G*QQ+P*R)+B*(H*PP-Q*(C+G*P+R)))+(-B*P+A*Q)*      \
       math.sqrt(-4.*(PP+BB*(1.+PP)-2.*B*Q+QQ-2.*A*(P+B*P*Q)+AA*(1.+QQ))*((FF+ \
       GG+HH-M)*(A-P)**2+2.*C*(A-P)*(G-F*P)+CC*(1.+PP)+2.*(A*F-G)*(A-P)*R-2.*C*\
       (1.+A*P)*R+(1.+AA)*R**2)+4.*(-2.*A*H*P+H*PP+C*Q-G*P*Q+A*(G+(C+F)*P)*Q-  \
       Q*R+B*((A-P)*(-G+F*P)-C*(1.+PP)+R+A*P*R)+AA*(H-Q*(F+R)))**2))/(2.*(A-   \
       P)*(PP+BB*(1.+PP)-2.*B*Q+QQ-2.*A*(P+B*P*Q)+AA*(1.+QQ)))
  yh1 = (-2.*(A-P)*(-(C*P)+H*P*Q-G*QQ-B*(-2.*G*Q+P*(H+(C+F)*Q))+P*R+A*((B-Q)*    \
       (H-F*Q)+C*(1.+QQ)-(1.+B*Q)*R)+BB*(-G+P*(F+R)))+(B-Q)*math.sqrt(-4.*(PP+   \
       BB*(1.+PP)-2.*B*Q+QQ-2.*A*(P+B*P*Q)+AA*(1.+QQ))*((FF+GG+HH-M)*(A-P)**2+2.*\
       C*(A-P)*(G-F*P)+CC*(1.+PP)+2.*(A*F-G)*(A-P)*R-2.*C*(1.+A*P)*R+(1.+AA)*RR)+\
       4.*(-2.*A*H*P+H*PP+C*Q-G*P*Q+A*(G+(C+F)*P)*Q-Q*R+B*((A-P)*(-G+F*P)-C*     \
       (1.+PP)+R+A*P*R)+AA*(H-Q*(F+R)))**2))/(2.*(A-P)*(PP+BB*(1.+PP)-2.*B*Q+QQ- \
       2.*A*(P+B*P*Q)+AA*(1.+QQ)))
  zh1 = -(-2.*H*(A-P)**2-2.*(C-(A*F-G)*(A-P)+A*C*P)*Q+2.*(1.+AA)*Q*R+2.*B*((A-P)*\
       (G-F*P)+C*(1.+PP)-(1.+A*P)*R)+math.sqrt(-4.*(PP+BB*(1.+PP)-2.*B*Q+QQ-2.*A*\
       (P+B*P*Q)+AA*(1.+QQ))*((FF+GG+HH-M)*(A-P)**2+2.*C*(A-P)*(G-F*P)+CC*(1.+   \
       PP)+2.*(A*F-G)*(A-P)*R-2.*C*(1.+A*P)*R+(1.+AA)*RR)+4.*(-2.*A*H*P+H*PP+C*Q-\
       G*P*Q+A*(G+(C+F)*P)*Q-Q*R+B*((A-P)*(-G+F*P)-C*(1.+PP)+R+A*P*R)+AA*(H-Q*(F+\
       R)))**2))/(2.*(PP+BB*(1.+PP)-2.*B*Q+QQ-2.*A*(P+B*P*Q)+AA*(1.+QQ)))
  if(flag > 1):
    xh2 = (-2.*(A-P)*(C*(PP+QQ)+BB*(P*(G-F*P)+R)+AA*(Q*(H-F*Q)+R)-A*(C*P+B*H*P+  \
         B*G*Q-2.*B*F*P*Q+H*P*Q-G*QQ+P*R)+B*(H*PP-Q*(C+G*P+R)))-(-B*P+A*Q)*      \
         math.sqrt(-4.*(PP+BB*(1.+PP)-2.*B*Q+QQ-2.*A*(P+B*P*Q)+AA*(1.+QQ))*((FF+ \
         GG+HH-M)*(A-P)**2+2.*C*(A-P)*(G-F*P)+CC*(1.+PP)+2.*(A*F-G)*(A-P)*R-2.*C*\
         (1.+A*P)*R+(1.+AA)*R**2)+4.*(-2.*A*H*P+H*PP+C*Q-G*P*Q+A*(G+(C+F)*P)*Q-  \
         Q*R+B*((A-P)*(-G+F*P)-C*(1.+PP)+R+A*P*R)+AA*(H-Q*(F+R)))**2))/(2.*(A-   \
         P)*(PP+BB*(1.+PP)-2.*B*Q+QQ-2.*A*(P+B*P*Q)+AA*(1.+QQ)))
    yh2 = (-2.*(A-P)*(-(C*P)+H*P*Q-G*QQ-B*(-2.*G*Q+P*(H+(C+F)*Q))+P*R+A*((B-Q)*    \
         (H-F*Q)+C*(1.+QQ)-(1.+B*Q)*R)+BB*(-G+P*(F+R)))-(B-Q)*math.sqrt(-4.*(PP+   \
         BB*(1.+PP)-2.*B*Q+QQ-2.*A*(P+B*P*Q)+AA*(1.+QQ))*((FF+GG+HH-M)*(A-P)**2+2.*\
         C*(A-P)*(G-F*P)+CC*(1.+PP)+2.*(A*F-G)*(A-P)*R-2.*C*(1.+A*P)*R+(1.+AA)*RR)+\
         4.*(-2.*A*H*P+H*PP+C*Q-G*P*Q+A*(G+(C+F)*P)*Q-Q*R+B*((A-P)*(-G+F*P)-C*     \
         (1.+PP)+R+A*P*R)+AA*(H-Q*(F+R)))**2))/(2.*(A-P)*(PP+BB*(1.+PP)-2.*B*Q+QQ- \
         2.*A*(P+B*P*Q)+AA*(1.+QQ)))
    zh2 = -(-2.*H*(A-P)**2-2.*(C-(A*F-G)*(A-P)+A*C*P)*Q+2.*(1.+AA)*Q*R+2.*B*((A-P)*\
         (G-F*P)+C*(1.+PP)-(1.+A*P)*R)-math.sqrt(-4.*(PP+BB*(1.+PP)-2.*B*Q+QQ-2.*A*\
         (P+B*P*Q)+AA*(1.+QQ))*((FF+GG+HH-M)*(A-P)**2+2.*C*(A-P)*(G-F*P)+CC*(1.+   \
         PP)+2.*(A*F-G)*(A-P)*R-2.*C*(1.+A*P)*R+(1.+AA)*RR)+4.*(-2.*A*H*P+H*PP+C*Q-\
         G*P*Q+A*(G+(C+F)*P)*Q-Q*R+B*((A-P)*(-G+F*P)-C*(1.+PP)+R+A*P*R)+AA*(H-Q*(F+\
         R)))**2))/(2.*(PP+BB*(1.+PP)-2.*B*Q+QQ-2.*A*(P+B*P*Q)+AA*(1.+QQ)))
  if(flag == 1):
    return (xh1,yh1,zh1)
  if(flag == 2):
    return (xh1,yh1,zh1),(xh2,yh2,zh2)

def run(file_name):
  print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  processed = monomer_library.pdb_interpretation.process(
    mon_lib_srv           = mon_lib_srv,
    ener_lib              = ener_lib,
    file_name             = file_name,
    keep_monomer_mappings = True,
    log                   = sys.stdout).all_chain_proxies
  total_missing = 0.0
  still_missing = 0.0
  still_missing_h = []
  atom_number = 0
  file_name_ =  os.path.basename(file_name)+"_h"
  file = open(file_name_,"w")

  for monomer_mapping in processed.all_monomer_mappings:
    atom_number = write_atoms(monomer_mapping,
                              atom_number,
                              file)
    bond_list = monomer_mapping.monomer.bond_list
    angle_list = monomer_mapping.monomer.angle_list
    missing_h = list(monomer_mapping.missing_hydrogen_atoms.keys())
    residue_name = monomer_mapping.monomer.chem_comp.three_letter_code
    print()
    print("Residue name: ", end=' ')
    print(monomer_mapping.residue_name, monomer_mapping.monomer.chem_comp.three_letter_code)
    print("Missing hydrogen atoms: ", missing_h)
    print()
    total_missing += len(missing_h)
    unknown_h = []
    n_missed = len(missing_h)-1
    while n_missed > 0:
      n_missed -= 1
      h = missing_h[0]
      for bl in bond_list:
        if(h == bl.atom_id_1 or h == bl.atom_id_2):
          for atom in (bl.atom_id_1,bl.atom_id_2):
            if(h != atom): target_atom_name = atom
          for atom_name,atom in monomer_mapping.expected_atoms.items():
            if(atom_name == target_atom_name):
              target_site = atom.xyz
          bond_dist = bl.value_dist
          format="missing %4s: bond: %4s %4s bond distance = %5.3f"
          print(format % (h, bl.atom_id_1, bl.atom_id_2, bl.value_dist))
      angles = []
      for al in angle_list:
        if(h == al.atom_id_1 or h == al.atom_id_2 or h == al.atom_id_3):
          format="        %4s: angle: %4s %4s %4s = %5.3f"
          print(format % (h, al.atom_id_1,al.atom_id_2,al.atom_id_3,al.value_angle))
          angles.append([al.value_angle,al.atom_id_1,al.atom_id_2,al.atom_id_3])
      bonded_to_target_site = []
      #for bl in bond_list:
      #  if(target_atom_name in (bl.atom_id_1,bl.atom_id_2)):
      #    for atom in (bl.atom_id_1,bl.atom_id_2):
      #      if(atom != target_atom_name):
      #        bonded_to_target_site.append(atom)
      #print "target_site, bonded_to_target_site ", target_atom_name, bonded_to_target_site
### ring hydrogens: TYR,PHE,HIS,TRP
      if(residue_name in ["PHE","TYR","HIS","TRP","ARG"]):
        ring_h = {
         "HE1":{"PHE":["CE1","CD1","CZ"],
                "TYR":["CE1","CD1","CZ"],
                "HIS":["CE1","ND1","NE2"],
                "TRP":["NE1","CD1","CE2"]},
         "HE2":{"PHE":["CE2","CD2","CZ"],
                "TYR":["CE2","CD2","CZ"],
                "HIS":["NE2","CE1","CD2"]},
         "HE3":{"TRP":["CE3","CZ3","CD2"]},
         "HD1":{"PHE":["CD1","CE1","CG"],
                "TYR":["CD1","CE1","CG"],
                "HIS":["ND1","CG","CE1"],
                "TRP":["CD1","CG","NE1"]},
         "HD2":{"PHE":["CD2","CE2","CG"],
                "TYR":["CD2","CE2","CG"],
                "HIS":["CD2","NE2","CG"]},
         "HH2":{"TRP":["CH2","CZ2","CZ3"]},
         "HZ" :{"PHE":["CZ","CE1","CE2"]},
         "HZ2":{"TRP":["CZ2","CE2","CH2"]},
         "HZ3":{"TRP":["CZ3","CH2","CE3"]},
         "HE" :{"ARG":["NE","CZ","CD"]}
                 }
        if h in ring_h.keys():
          if residue_name in ring_h[h].keys():
            targets = ring_h[h][residue_name]
            print("Building:", h, " ...")
            site_0 = None
            site_1 = None
            site_2 = None
            for atom_name,atom in monomer_mapping.expected_atoms.items():
              if(atom_name == targets[0]):
                site_0 = atom.xyz
              if(atom_name == targets[1]):
                site_1 = atom.xyz
              if(atom_name == targets[2]):
                site_2 = atom.xyz
            alpha = angles[0][0]
            beta  = angles[1][0]
            assert site_0 is not None and site_1 is not None and site_2 is not None
            xyz = add_ring_h(site_0,site_1,site_2,bond_dist,alpha,beta)
            atom_number = write_atoms(monomer_mapping,
                                      atom_number,
                                      file,
                                      new_atom_name        = h,
                                      new_atom_coordinates = xyz)
            missing_h.remove(h)
### CA hydrogen(s): HA or/and HA1 or/and HA2
      if(h == "HA" or h == "HA1" or h == "HA2"):
        build_one = False
        build_two = False
        if("HA" in missing_h):
          build_one = True
        if("HA1" in missing_h and "HA2" in missing_h):
          build_two = True
        if(("HA1" in missing_h or "HA2" in missing_h) and build_two == False):
          build_one = True
        assert build_one == True and build_two == False or \
               build_one == False and build_two == True
        if(build_two == True):
          print("Building: HA1 and HA2")
          site_0 = None
          site_1 = None
          site_2 = None
          for atom_name,atom in monomer_mapping.expected_atoms.items():
            if(atom_name == "CA"):
              site_0 = atom.xyz
            if(atom_name == "C"):
              site_1 = atom.xyz
            if(atom_name == "N"):
              site_2 = atom.xyz
          if(site_0 is not None and site_1 is not None and site_2 is not None):
            for angle in angles:
              if("C" in angle and "CA" in angle): alpha = angle[0]
              if("N" in angle and "CA" in angle): beta  = angle[0]
            xyz = add_ca_ha1_and_ha2(site_0,site_1,site_2,bond_dist,alpha,beta)
            atom_number = write_atoms(monomer_mapping,
                                      atom_number,
                                      file,
                                      new_atom_name        = "HA1",
                                      new_atom_coordinates = xyz[0])
            atom_number = write_atoms(monomer_mapping,
                                      atom_number,
                                      file,
                                      new_atom_name        = "HA2",
                                      new_atom_coordinates = xyz[1])
            missing_h.remove("HA1")
            missing_h.remove("HA2")
        if(build_one == True):
          print("Building: ", h)
          site_0 = None
          site_1 = None
          site_2 = None
          site_3 = None
          for atom_name,atom in monomer_mapping.expected_atoms.items():
            if(atom_name == "CA"):
              site_0 = atom.xyz
            if(atom_name == "C"):
              site_1 = atom.xyz
            if(atom_name == "N"):
              site_2 = atom.xyz
            if(atom_name == "CB" or atom_name == "HA1" or atom_name == "HA2"):
              site_3 = atom.xyz
              cb_or_ha1_or_ha2 = atom_name
          if(site_0 is not None and site_1 is not None and site_2 is not None):
            for angle in angles:
              if("C" in angle and "CA" in angle): alpha  = angle[0]
              if("N" in angle and "CA" in angle): beta   = angle[0]
              if("CA" in angle and cb_or_ha1_or_ha2 in angle): gamma = angle[0]
            xyz = add_ca_ha1_or_ha2(site_0,site_1,site_2,site_3,bond_dist,alpha,beta,gamma)
            atom_number = write_atoms(monomer_mapping,
                                      atom_number,
                                      file,
                                      new_atom_name        = h,
                                      new_atom_coordinates = xyz)
            missing_h.remove(h)
### CB hydrogen(s): HB1 or/and HB2 and similar
      cb_like_h = {
       "HB1":[["CB","CA","CG","HB2"],["CB","CA","OG","HB2"],["CB","CA","SG","HB2"]],
       "HB2":[["CB","CA","CG","HB1"],["CB","CA","OG","HB1"],["CB","CA","SG","HB1"]],
       "HG1":[["CG","CB","CD","HG2"],],
       "HG2":[["CG","CB","CD","HG1"],],
       "HD1":[["CD","CG","NE","HD2"],["CD","CG","CE","HD2"]],
       "HD2":[["CD","CG","NE","HD1"],["CD","CG","CE","HD1"]],
       "HE1":[["CE","CD","NZ","HE2"],],
       "HE2":[["CE","CD","NZ","HE1"],],
       "HG" :[["CG","CB","CD1","CD2"],],
       "HB" :[["CB","CG2","CA","CG1"],["CB","CG2","CA","OG1"]]
                  }
      if h in cb_like_h.keys():
        if(h in missing_h):
          targets_ = cb_like_h[h]
          for targets in targets_:
            build_one = False
            build_two = False
            if(h in missing_h and targets[3] in missing_h):
              build_two = True
            if((h in missing_h or targets[3] in missing_h) and build_two == False):
              build_one = True
            if(build_one == True and build_two == False or \
                                         build_one == False and build_two == True):
              if(build_two == True):
                print("Building: ", h," and ",targets[3])
                site_0 = None
                site_1 = None
                site_2 = None
                for atom_name,atom in monomer_mapping.expected_atoms.items():
                  if(atom_name == targets[0]):
                    site_0 = atom.xyz
                  if(atom_name == targets[1]):
                    site_1 = atom.xyz
                  if(atom_name == targets[2]):
                    site_2 = atom.xyz
                if(site_0 is not None and site_1 is not None and site_2 is not None):
                  for angle in angles:
                    if(targets[1] in angle and targets[0] in angle): alpha = angle[0]
                    if(targets[0] in angle and targets[2] in angle): beta  = angle[0]
                  xyz = add_ca_ha1_and_ha2(site_0,site_1,site_2,bond_dist,alpha,beta)
                  atom_number = write_atoms(monomer_mapping,
                                      atom_number,
                                      file,
                                      new_atom_name        = h,
                                      new_atom_coordinates = xyz[0])
                  atom_number = write_atoms(monomer_mapping,
                                      atom_number,
                                      file,
                                      new_atom_name        = targets[3],
                                      new_atom_coordinates = xyz[1])
                  missing_h.remove(h)
                  missing_h.remove(targets[3])
              if(build_one == True):
                print("Building: ", h)
                site_0 = None
                site_1 = None
                site_2 = None
                site_3 = None
                alpha  = None
                beta   = None
                gamma  = None
                for atom_name,atom in monomer_mapping.expected_atoms.items():
                  if(atom_name == targets[0]):
                    site_0 = atom.xyz
                  if(atom_name == targets[1]):
                    site_1 = atom.xyz
                  if(atom_name == targets[2]):
                    site_2 = atom.xyz
                  if(atom_name == h or atom_name == targets[3]):
                    site_3 = atom.xyz
                    cb_or_ha1_or_ha2 = atom_name
                if(site_0 is not None and site_1 is not None and site_2 is not None and site_3 is not None):
                  for angle in angles:
                    if(targets[1] in angle and targets[0] in angle): alpha = angle[0]
                    if(targets[0] in angle and targets[2] in angle): beta  = angle[0]
                    if(targets[0] in angle and cb_or_ha1_or_ha2 in angle): gamma = angle[0]
                  if(alpha is not None and beta is not None and gamma is not None):
                    xyz = add_ca_ha1_or_ha2(site_0,site_1,site_2,site_3,bond_dist,alpha,beta,gamma)
                    atom_number = write_atoms(monomer_mapping,
                                      atom_number,
                                      file,
                                      new_atom_name        = h,
                                      new_atom_coordinates = xyz)
                    missing_h.remove(h)
### HT1,HT2,HT3 like hydrogens; rotational optimization is necessary
      hhh_like_h = {
       "HB1" :[["CB","CA","HB2","HB3"],],
       "HB2" :[["CB","CA","HB1","HB3"],],
       "HB3" :[["CB","CA","HB1","HB2"],],
       "HD11":[["CD1","CG","HD12","HD13"],["CD1","CG1","HD12","HD13"]],
       "HD12":[["CD1","CG","HD11","HD13"],["CD1","CG1","HD11","HD13"]],
       "HD13":[["CD1","CG","HD11","HD12"],["CD1","CG1","HD11","HD12"]],
       "HD21":[["CD2","CG","HD22","HD23"],],
       "HD22":[["CD2","CG","HD21","HD23"],],
       "HD23":[["CD2","CG","HD21","HD22"],],
       "HG21":[["CG2","CB","HG22","HG23"],],
       "HG22":[["CG2","CB","HG21","HG23"],],
       "HG23":[["CG2","CB","HG21","HG22"],],
       "HG11":[["CG1","CB","HG12","HG13"],],
       "HG12":[["CG1","CB","HG11","HG13"],],
       "HG13":[["CG1","CB","HG11","HG12"],],
       "HZ1" :[["NZ","CE","HZ2","HZ3"],],
       "HZ2" :[["NZ","CE","HZ1","HZ3"],],
       "HZ3" :[["NZ","CE","HZ1","HZ2"],]
                   }
      hhh_residues = ["ALA","LEU","ILE","VAL","THR","LYS"]
      if h in hhh_like_h.keys() and residue_name in hhh_residues and h in missing_h:
        targets = hhh_like_h[h]
        for target in targets:
          site_0 = None
          site_1 = None
          for atom_name,atom in monomer_mapping.expected_atoms.items():
            if(atom_name == target[0]):
              site_0 = atom.xyz
            if(atom_name == target[1]):
              site_1 = atom.xyz
          alpha = angles[0][0]
          beta  = angles[2][0]
          if(site_0 is not None and site_1 is not None):
            print("Building:", h, " ...")
            xyz = add_hhh(site_0,site_1,bond_dist,alpha,beta)
            for atom_name, coordinates in zip((h,target[2],target[3]),(xyz[0],xyz[1][0],xyz[1][1])):
              atom_number = write_atoms(monomer_mapping,
                                        atom_number,
                                        file,
                                        new_atom_name        = atom_name,
                                        new_atom_coordinates = coordinates)
            if(h in missing_h): missing_h.remove(h)
            if(target[2] in missing_h): missing_h.remove(target[2])
            if(target[3] in missing_h): missing_h.remove(target[3])
### add ARG-like NH's:
      arg_like_h = {
       "HH11":[["NH1","CZ","NH2","HH12",2],],
       "HH12":[["NH1","CZ","NH2","HH11",2],],
       "HH21":[["NH2","CZ","NH1","HH22",2],],
       "HH22":[["NH2","CZ","NH1","HH21",2],],
       "HG"  :[["OG","CB","CA","HG",1],["SG","CB","CA","HG",1]] ,
       "HH"  :[["OH","CZ","CE1","HH",1],]
                   }
      residue_names = ["ARG","SER","TYR","CYS"]
      if h in arg_like_h.keys() and residue_name in residue_names:
        targets = arg_like_h[h]
        for target in targets:
          site_0 = None
          site_1 = None
          site_2 = None
          for atom_name,atom in monomer_mapping.expected_atoms.items():
            if(atom_name == target[0]):
              site_0 = atom.xyz
            if(atom_name == target[1]):
              site_1 = atom.xyz
            if(atom_name == target[2]):
              site_2 = atom.xyz
          alpha = angles[0][0]
          if(site_0 is not None and site_1 is not None and site_2 is not None):
            print("Building:", h, " ...")
            xyz = add_arg_like_h(site_0,site_1,site_2,bond_dist,alpha,flag = target[4])
            if(target[4] == 1):
              atom_number = write_atoms(monomer_mapping,
                                        atom_number,
                                        file,
                                        new_atom_name        = h,
                                        new_atom_coordinates = xyz)
            if(target[4] == 2):
              atom_number = write_atoms(monomer_mapping,
                                        atom_number,
                                        file,
                                        new_atom_name        = h,
                                        new_atom_coordinates = xyz[0])
              atom_number = write_atoms(monomer_mapping,
                                        atom_number,
                                        file,
                                        new_atom_name        = target[3],
                                        new_atom_coordinates = xyz[1])
            if(h in missing_h): missing_h.remove(h)
            if(target[3] in missing_h): missing_h.remove(target[3])

### UNKNOWN H:
      if(h in missing_h):
        print("Unknown hydrogen type: ", h)
        missing_h.remove(h)
        unknown_h.append(h)
      n_missed = len(missing_h)

    print()
    print("Still missing: ", residue_name,unknown_h)
    still_missing_h.append((residue_name,unknown_h))
    still_missing += len(unknown_h)

  print("Build ", int(total_missing - still_missing), end=' ')
  print(" from ", int(total_missing), " % = ", end=' ')
  print((total_missing - still_missing)*100./total_missing)
  print()
  for item in still_missing_h:
    print(item)

  file.close()
  file = open(file_name_,"r")
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv = mon_lib_srv,
    ener_lib = ener_lib,
    file_name = file_name_,
    strict_conflict_handling = False,
    crystal_symmetry = processed.pdb_inp.crystal_symmetry(),
    force_symmetry = True,
    log = sys.stdout)
  file.close()
  assert processed_pdb_file.xray_structure() is not None
  assert processed_pdb_file.geometry_restraints_manager() is not None
  xray_structure = processed_pdb_file.xray_structure()
  geometry_restraints_manager = processed_pdb_file.geometry_restraints_manager()
  regularize_model(xray_structure = xray_structure,
                   geometry_restraints_manager = geometry_restraints_manager,
                   max_iterations = 50000)
  processed_pdb_file.all_chain_proxies.pdb_atoms.set_xyz(
    new_xyz=xray_structure.sites_cart())
  processed_pdb_file.all_chain_proxies.pdb_hierarchy.write_pdb_file(
    file_name="out.pdb")

def write_atoms(monomer_mapping,
                atom_number,
                file_object,
                new_atom_name = None,
                new_atom_coordinates = None,
                crystal_symmetry = None):
  assert new_atom_name is None and new_atom_coordinates is None or \
         new_atom_name is not None and new_atom_coordinates is not None
  if(crystal_symmetry is not None):
    print(pdb.format_cryst1_record(
                                             crystal_symmetry=crystal_symmetry), file=file_object)
    print(pdb.format_scale_records(
                                        unit_cell=crystal_symmetry.unit_cell()), file=file_object)
  if(new_atom_name is None and new_atom_coordinates is None):
    for atom_name,atom in monomer_mapping.expected_atoms.items():
      assert atom_name.strip() == atom.name.strip()
      atom_number += 1
      orig = atom.serial
      atom.serial = "%5d" % atom_number
      print(atom.format_atom_record(), file=file_object)
      atom.serial = orig
  else:
    # FIXME ordering of values changes for py2/3, this could break if more than 1 value present
    atom = list(monomer_mapping.expected_atoms.values())[0]
    atom_number += 1
    orig = atom.serial, atom.name, atom.xyz, atom.occ, atom.b
    atom.serial = "%5d" % atom_number
    atom.name = new_atom_name
    atom.xyz = new_atom_coordinates
    atom.occ = 1
    atom.b = 0
    print(atom.format_atom_record(), file=file_object)
    atom.serial, atom.name, atom.xyz, atom.occ, atom.b = orig
  return atom_number

def regularize_model(xray_structure,
                     geometry_restraints_manager,
                     max_iterations):
  sites_cart = xray_structure.sites_cart()
  minimized = geometry_restraints_lbfgs(
    sites_cart                  = sites_cart,
    geometry_restraints_manager = geometry_restraints_manager,
    lbfgs_termination_params    = scitbx.lbfgs.termination_parameters(
      max_iterations = max_iterations))
  xray_structure.set_sites_cart(sites_cart)
  print()
  print("Energies at start of minimization:")
  minimized.first_target_result.show()
  print()
  print("Energies after minimization:")
  minimized.final_target_result.show()


if (__name__ == "__main__"):
  files = ["arg.pdb","lys.pdb","ala.pdb","gly.pdb","his.pdb","ile.pdb","leu.pdb",
           "phe.pdb","tyr.pdb","trp.pdb","thr.pdb","val.pdb","cys.pdb","ser.pdb"]
  if(len(sys.argv) > 1):
    run(sys.argv[1])
  else:
    pdb_dir = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/hydrogens",
      test=os.path.isdir)
    if (pdb_dir is None):
      print("Skipping build_hydrogens.run(): input files not available")
    else:
      for file in files:
        run(os.path.join(pdb_dir, file))
