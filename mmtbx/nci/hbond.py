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

  def get_hydrogen_bonds_pairs(self, min = 1.7,max = 2.2,eps1 = 0.8,
                                        ideal_dist = 3, eps2 = 0.5):
      geometry = self.model.get_restraints_manager()
      bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
                                  sites_cart=self.model.get_sites_cart())
      bps_dict = {}
      [bps_dict.setdefault(p.i_seqs, True) for p in bond_proxies_simple]
      hierarchy = self.model.get_hierarchy()
      atom1s = []
      atom2s = []
      atom3s = []
      atom4s = []
      results = []
      Accepter_H_pair = ["O","N","S","F","CL"]
      for a in hierarchy.atoms():
        e = a.element.strip().upper()
        if a.element_is_hydrogen():
          atom1s.append(a)
        if e == "O":
          if a.parent().resname == "HOH":continue
          atom2s.append(a)
        if e in Accepter_H_pair:
          atom3s.append(a)
        if e == "C":
          atom4s.append(a)

      if atom1s is not None:
        for a2 in atom2s:
          result = None
          diff_best = 1.e+9
          for a3 in atom3s:
            for a1 in atom1s:
              resid_2 = a2.parent().parent().resid()
              resid_3 = a3.parent().parent().resid()
              diff_r_r = abs(int(resid_2) - int(resid_3) )
              if diff_r_r < 2 :continue
              if (not a1.is_in_same_conformer_as(a2)): continue
              if (not is_bonded(a1, a3, bps_dict)): continue
              if (is_bonded(a1, a2, bps_dict)): continue
              if (a1.parent().parent().resseq ==
                  a2.parent().parent().resseq): continue
              d_12 = a1.distance(a2)
              if (min-eps1 < d_12 < max+eps1):
                angle_312 = (a1.angle(a2, a3, deg=True))
                if (90 < angle_312):
                  diff = abs( 2 - d_12 )
                  if (diff < diff_best):
                      diff_best = diff
                      result = group_args(
                        atom_1=a1,
                        atom_2=a2,
                        atom_3=a3,
                        d=d_12,
                        angle_312=angle_312)
            if (result in results):continue
            if (result is not None): results.append(result)
        # one hydrogen atom can just make up one H_bond
        for i,ri in enumerate(results):
          for j,rj in enumerate(results):
            if (j <= i): continue
            ai1 = ri.atom_1
            ai2 = ri.atom_2
            aj1 = rj.atom_1
            aj2 = rj.atom_2
            if ri.atom_1 == rj.atom_1:
              di = ai1.distance(ai2)
              dj = aj1.distance(aj2)
              if di < dj:
                if rj in results:
                  results.remove(rj)
              else:
                if ri in results:
                  results.remove(ri)
        # just keep the more possiable situation for O atom
        for i, ri in enumerate(results):
          for j, rj in enumerate(results):
            if (j <= i): continue
            ai1 = ri.atom_1
            ai2 = ri.atom_2
            aj1 = rj.atom_1
            aj2 = rj.atom_2
            if ri.atom_2 == rj.atom_2:
              di = ai1.distance(ai2)
              dj = aj1.distance(aj2)
              if di < dj:
                if rj in results:
                  results.remove(rj)
              else:
                if ri in results:
                  results.remove(ri)

      if len(atom1s) == 0 :
        for a3 in atom3s:
          result = None
          diff_best = 1.e+9
          for a2 in atom2s:
            for a4 in atom4s:
              if (not is_bonded(a4, a3, bps_dict)): continue
              resid_2 = a2.parent().parent().resid()
              resid_3 = a3.parent().parent().resid()
              diff_r_r = abs(int(resid_2) - int(resid_3))
              if diff_r_r < 2: continue
              if (not a3.is_in_same_conformer_as(a2)): continue
              if (is_bonded(a2, a3, bps_dict)): continue
              d_O_N = a3.distance(a2)
              if (ideal_dist -eps2 < d_O_N < ideal_dist + eps2):
                diff = abs(ideal_dist - d_O_N)
                if (diff < diff_best):
                  diff_best = diff
                  result = group_args(
                    atom_1=a2,
                    atom_2=a3,
                    atom_3=a4,
                    d=d_O_N)
            if (result in results): continue
            if (result is not None): results.append(result)

        #just keep the more possiable situation for N atom

        for i,ri in enumerate(results):
          for j,rj in enumerate(results):
            if (j <= i): continue
            ai1 = ri.atom_1
            ai2 = ri.atom_2
            aj1 = rj.atom_1
            aj2 = rj.atom_2
            if ri.atom_1 == rj.atom_1:
              di = ai1.distance(ai2)
              dj = aj1.distance(aj2)
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
      sigma = 0.01
      slack = None
      limit = None
      top_out = False
    }
    angle {
      atom_selection_1 = %s
      atom_selection_2 = %s
      atom_selection_3 = %s
      angle_ideal = %f
      sigma = 2.5
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
        r.atom_1.chain().id,
        r.atom_1.parent().parent().resid(),
        r.atom_1.name)
      a2_str = "chain %s and resseq %s and name %s" % (
        r.atom_2.chain().id,
        r.atom_2.parent().parent().resid(),
        r.atom_2.name)
      a3_str = "chain %s and resseq %s and name %s" % (
        r.atom_3.chain().id,
        r.atom_3.parent().parent().resid(),
        r.atom_3.name)
      if (use_defaul_parameters):
        d_ideal_1 = 2.19
        d_ideal_2 = 2.9
      else:
        d_ideal_1 = r.d
        d_ideal_2 = r.d
      i = i + 1
      if r.atom_1.element.strip().upper() == "H":
        if (use_defaul_parameters):
          angle_ideal = 153.4
        else:angle_ideal = r.angle_312
        bond_angle_str = str_1 % (a1_str, a2_str, d_ideal_1,
                                  a1_str, a2_str, a3_str, angle_ideal)
        sub_fin_str = sub_fin_str + bond_angle_str
      else :
        if (use_defaul_parameters):
          angle_ideal = 147.15
          bond_angle_str = str_1 % (a1_str, a2_str, d_ideal_2,
                                  a1_str, a2_str, a3_str, angle_ideal)
          sub_fin_str = sub_fin_str + bond_angle_str
    s_f_str = sub_fin_str[1:]
    str_final = str_2 % (s_f_str)
    file_name = pdb_file_name
    with open(file_name,'w') as fileobject:
      fileobject.write(str_final)
