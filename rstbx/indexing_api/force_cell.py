from __future__ import absolute_import, division, print_function
from six.moves import range
import math
from scitbx.matrix import col,sqr
from scitbx.math import unimodular_generator
from cctbx.uctbx import unit_cell
from cctbx import sgtbx
from rstbx.dps_core import Orientation

def is_length_match(a,b):
    return ( abs(a-b) < 0.02 * (a+b) )

def is_angle_match(a,b):
    return ( abs(a-b) < 3.0 )

def generate_unimodular_cells(cell):
    Amat = sqr(cell.orthogonalization_matrix()).transpose()

    for m in unimodular_generator(range=1).all():
      c_inv = sgtbx.rt_mx(sgtbx.rot_mx(m))
      orientation_similarity_cb_op = sgtbx.change_of_basis_op(c_inv).inverse()
      new_cell = cell.change_basis(orientation_similarity_cb_op)
      yield new_cell,orientation_similarity_cb_op

def get_triangle(lines):
  gamma = lines["gamma"]["match"]
  alpha = lines["alpha"]["match"]
  beta = lines["beta"]["match"]
  for acell,bcell in gamma:
    for al_bcell,al_ccell in alpha:
      if al_bcell != bcell: continue
      for be_ccell, be_acell in beta:
        if be_ccell != al_ccell: continue
        if be_acell != acell: continue
        return acell, al_bcell, be_ccell
  return None

def force_cell(index_engine,target_cell,verbose=False):

    if verbose:
      print("N candidates:",index_engine.n_candidates())
      from rstbx.dps_core import directional_show
      for x in range(index_engine.n_candidates()):
        directional_show( index_engine[x], "vector %d:"%x )
    Ns = index_engine.n_candidates()
    best = {"score":1.E100}
    orth = target_cell.orthogonalization_matrix()
    for working_cell,cb_op in generate_unimodular_cells(target_cell):
      #print "---->",working_cell, working_cell.volume()

      vectors = {"a":{"match":[],"length":working_cell.parameters()[0]},
                 "b":{"match":[],"length":working_cell.parameters()[1]},
                 "c":{"match":[],"length":working_cell.parameters()[2]},
                }
      for key in vectors.keys():
        for ns in range(Ns):
          if is_length_match(vectors[key]["length"],index_engine[ns].real):
            vectors[key]["match"].append(ns)
        #print key, vectors[key]

      lines = {"alpha":{"match":[],"points":("b","c"),"angle":working_cell.parameters()[3],"hit":False},
                "beta":{"match":[],"points":("c","a"),"angle":working_cell.parameters()[4],"hit":False},
               "gamma":{"match":[],"points":("a","b"),"angle":working_cell.parameters()[5],"hit":False},
              }
      for key in lines.keys():
        xmatch = len(vectors[lines[key]["points"][0]]["match"])
        ymatch = len(vectors[lines[key]["points"][1]]["match"])
        for xi in range(xmatch):
          for yi in range(ymatch):
            xkey = vectors[lines[key]["points"][0]]["match"][abs(xi)]
            ykey = vectors[lines[key]["points"][1]]["match"][abs(yi)]
            #print key,xkey,ykey,

            xvector = col(index_engine[xkey].dvec)
            yvector = col(index_engine[ykey].dvec)

            #print xvector.dot(yvector),
            costheta = xvector.dot(yvector)/math.sqrt(xvector.dot(xvector)*yvector.dot(yvector))
            angle = math.acos(costheta)*180./math.pi
            #print "angle %.2f"%angle,
            if is_angle_match(angle, lines[key]["angle"]):
              #print "*****",;
              lines[key]["hit"]=True
              lines[key]["match"].append((xkey,ykey))

      if lines["alpha"]["hit"] and lines["beta"]["hit"] and lines["gamma"]["hit"]:
        #print "HELLO HIT"
        #for key in lines.keys():
        #  print key, lines[key]
        tri = get_triangle(lines)
        #print "Triangle:",tri
        if tri==None: continue
        direct_matrix = sqr((
          index_engine[tri[0]].bvec()[0],index_engine[tri[0]].bvec()[1],index_engine[tri[0]].bvec()[2],
          index_engine[tri[1]].bvec()[0],index_engine[tri[1]].bvec()[1],index_engine[tri[1]].bvec()[2],
          index_engine[tri[2]].bvec()[0],index_engine[tri[2]].bvec()[1],index_engine[tri[2]].bvec()[2],
           ))
        ori = Orientation(direct_matrix,False)
        #print "Found unit cell",ori.unit_cell()
        #print "compatible",ori.unit_cell().change_basis(cb_op.inverse()),
        modo = ori.unit_cell().change_basis(cb_op.inverse()).orthogonalization_matrix()
        diff = ( modo[0]-orth[0],modo[1]-orth[1],modo[2]-orth[2],modo[4]-orth[4],modo[5]-orth[5],modo[8]-orth[8])
        score = math.sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]+diff[3]*diff[3]+diff[4]*diff[4]+diff[5]*diff[5])
        #print "score %.1f"%score,tri
        if score<best["score"]:
          best = {"score":score,"triangle":tri,"orientation":ori.change_basis(sqr(cb_op.inverse().c().r().as_double()).transpose())}
    if verbose: print("Best score %.1f, triangle %12s"%(best["score"],str(best["triangle"])),best["orientation"].unit_cell())
    return best
