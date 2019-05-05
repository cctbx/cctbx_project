from __future__ import division
from libtbx import group_args

def is_bonded(atom_1, atom_2, bps_dict):
  i12 = [atom_1.i_seq, atom_2.i_seq]
  i12.sort()
  if(not tuple(i12) in bps_dict): return False
  else: return True

class get_hydrogen_bonds(object):
  def __init__(self,model):
    self.model = model
    self.results = self.get_hydrogen_bonds_pairs()

  def get_hydrogen_bonds_pairs(self, ideal_angel_YAD = 147.15,
                    angle_AHD_cutoff = 120,eps_angle_AHD = 30,
                    angle_HAY_min = 90,angle_HAY_max = 180,
                    eps_angle_HAY = 10,ideal_dist_A_D = 2.90,
                    sigma_for_angle= 5.0, sigma_for_bond = 0.1,
                                           eps_dist_A_D= 0.5 ):
      # Hydrogen bond  model : Y-A...H-D ;
      geometry = self.model.get_restraints_manager()
      bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
                                  sites_cart=self.model.get_sites_cart())
      bps_dict = {}
      [bps_dict.setdefault(p.i_seqs, True) for p in bond_proxies_simple]
      hierarchy = self.model.get_hierarchy()
      atom_H = []
      atom_A = []
      atom_D = []
      atom_Y = []
      ress    = []
      results = []
      Accepter_H_pair = ["O","N","S","F","CL"]
      for a in hierarchy.atoms():
        e = a.element.strip().upper()
        if a.element_is_hydrogen():
          atom_H.append(a)
        if e == "O":
          if a.parent().resname == "HOH": continue
          atom_A.append(a)
        if e in Accepter_H_pair:
          if a.parent().resname == "HOH": continue
          atom_D.append(a)
        if e == "C":
          atom_Y.append(a)

      for a_A in atom_A:
        res = None
        diff_best = 1.e+9
        for a_D in atom_D:
          for a_H in atom_H:
            resid_A = a_A.parent().parent().resid()
            resid_D = a_D.parent().parent().resid()
            diff_r_r = abs(int(resid_A) - int(resid_D) )
            if diff_r_r < 2 :continue
            if (not a_H.is_in_same_conformer_as(a_A)): continue
            if (not is_bonded(a_H, a_D, bps_dict)): continue
            if (is_bonded(a_H, a_A, bps_dict)): continue
            if (a_A.parent().parent().resseq ==
                a_D.parent().parent().resseq): continue
            d_A_D = a_D.distance(a_A)
            if (ideal_dist_A_D - eps_dist_A_D  <
                  d_A_D  < ideal_dist_A_D + eps_dist_A_D ):
              angle_AHD = a_H.angle(a_A, a_D, deg=True)
              if (angle_AHD_cutoff - eps_angle_AHD < angle_AHD):
                diff = abs( ideal_dist_A_D - d_A_D )
                if (diff < diff_best):
                    diff_best = diff
                    res = group_args(
                      a_H=a_H,
                      a_A=a_A,
                      a_D=a_D,
                      d_A_D=d_A_D,
                      angle_AHD=angle_AHD,
                      sigma_for_bond=sigma_for_bond,
                      sigma_for_angle=sigma_for_angle,
                      ideal_angel_YAD=ideal_angel_YAD,
                      ideal_dist_A_D=ideal_dist_A_D)
                   # print sigma_for_angle, type(sigma_for_angle)
        if (res in ress):continue
        if (res is not None): ress.append(res)
      for r in ress:
        #print r.sigma_for_angle,type(r.sigma_for_angle)
        a_H = r.a_H
        a_A = r.a_A
        a_D = r.a_D
        d_A_D  = r.d_A_D
        angle_AHD = r.angle_AHD
        sigma_for_angle = r.sigma_for_angle
        sigma_for_bond = r.sigma_for_bond
        ideal_angel_YAD = r.ideal_angel_YAD
        ideal_dist_A_D = r.ideal_dist_A_D
        #print r.sigma_for_angle, type(r.sigma_for_angle)
        result = None
        for a_Y in atom_Y :
          if (not is_bonded(a_A, a_Y, bps_dict)): continue
          angle_YAD = a_A.angle(a_Y,a_D,deg=True)
          angle_HAY = a_A.angle(a_H,a_Y,deg=True)
          if (angle_HAY_min - eps_angle_HAY < angle_HAY <
                            angle_HAY_max + eps_angle_HAY):
            #print r.sigma_for_angle, type(r.sigma_for_angle)
            result = group_args(
              a_H=a_H,
              a_A=a_A,
              a_D=a_D,
              a_Y=a_Y,
              d_A_D=d_A_D,
              angle_HAY=angle_HAY,
              angle_AHD=angle_AHD,
              angle_YAD=angle_YAD,
              sigma_for_bond=sigma_for_bond,
              sigma_for_angle=sigma_for_angle,
              ideal_angel_YAD=ideal_angel_YAD,
              ideal_dist_A_D=ideal_dist_A_D)
          if (result in results):continue
          if (result is not None): results.append(result)

      # just keep the more possiable situation for N atom

      for i, ri in enumerate(results):
        for j, rj in enumerate(results):
          if (j <= i): continue
          a_A_i = ri.a_A
          a_D_i = ri.a_D
          a_A_j = rj.a_A
          a_D_j = rj.a_D
          if ri == rj:
            results.remove(rj)
          if a_D_i == a_D_j:
            di = a_A_i.distance(a_D_i)
            dj = a_A_j.distance(a_D_j)
            if di < dj:
              if rj in results:
                results.remove(rj)
            else:
              if ri in results:
                results.remove(ri)
      return results


  def write_restrains_file(self, pdb_file_name,
                           for_phenix_refine=True,
                           use_defaul_parameters=True):
    str_1 = '''bond{
      atom_selection_1 = %s
      atom_selection_2 = %s
      symmetry_operation = None
      distance_ideal = %f
      sigma = %f
      slack = None
      limit = -0.1
      top_out = False
    }
    angle {
      atom_selection_1 = %s
      atom_selection_2 = %s
      atom_selection_3 = %s
      angle_ideal = %f
      sigma = %f
    }
    '''

    str_2 = '''refinement{
  geometry_restraints.edits{
    %s
  }
}
    '''
    i = 1
    sub_fin_str = 'a'
    for r in self.results:
      a1_str = "chain %s and resseq %s and name %s" % (
        r.a_A.chain().id,
        r.a_A.parent().parent().resid(),
        r.a_A.name)
      a2_str = "chain %s and resseq %s and name %s" % (
        r.a_D.chain().id,
        r.a_D.parent().parent().resid(),
        r.a_D.name)
      a3_str = "chain %s and resseq %s and name %s" % (
        r.a_Y.chain().id,
        r.a_Y.parent().parent().resid(),
        r.a_Y.name)
      if (use_defaul_parameters):
        d_ideal = r.ideal_dist_A_D
      else:
        d_ideal = r.d_A_D
      i = i + 1
      if (use_defaul_parameters):
        angle_ideal = r.ideal_angel_YAD
      else:angle_ideal = r.angle_YAD
      sigma_angle = r.sigma_for_angle
      sigma_bond = r.sigma_for_bond
      bond_angle_str = str_1 % (a1_str, a2_str, d_ideal,
                                sigma_bond,a1_str, a2_str,
                                a3_str, angle_ideal,
                                sigma_angle)
      sub_fin_str = sub_fin_str + bond_angle_str
    s_f_str = sub_fin_str[1:]
    str_final = str_2 % (s_f_str)
    file_name = pdb_file_name
    with open(file_name,'w') as fileobject:
      fileobject.write(str_final)
