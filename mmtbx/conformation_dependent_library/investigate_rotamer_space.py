from __future__ import absolute_import, division, print_function
import sys
import copy

from mmtbx.rotamer import rotamer_eval
from mmtbx.conformation_dependent_library import rdl_database
from six.moves import range

number_of_chis = {
  'ARG' : 4,
  'ASN' : 2,
  'ASP' : 2,
  'CYS' : 2,
  'GLN' : 3,
  'GLU' : 3,
  'HIS' : 2,
  'ILE' : 2,
  'LEU' : 2,
  'LYS' : 4,
  'MET' : 3,
  'PHE' : 2,
  'PRO' : 3,
  'SER' : 2,
  'THR' : 1,
  'TRP' : 2,
  'TYR' : 2,
  'VAL' : 1,
}

flat = [0., 20., 36]
lookup = {
  'ASP' : { 2 : ( "CA",  "CB",  "CG",  "OD1" )},
  'GLN' : { 3 : ( "CB",  "CG",  "CD",  "OE1" )},
  'GLU' : { 3 : ( "CB",  "CG",  "CD",  "OE1" )},
  'PHE' : { 2 : ( "CA",  "CB",  "CG",  "CD1" )},
  'TYR' : { 2 : ( "CA",  "CB",  "CG",  "CD1" )},
}

results = {
  'ARG' : {
    'tpm-80' : None, # rare rotamer
    },
  'ASN' : {
    'p0' : None, # wide blob around zero
    't0' : None,
    'm-40' : None, # not so wide
    't160' : None, # no points
    },
  'ASP' : {
    'p0' : { lookup['ASP'][2] : flat},
           # bimodal but wide - make 2 minima or none (36)
    't0' : { lookup['ASP'][2] : [0, 30, 2]},
    },
  'CYS' : 2, # chi2 is to the hydrogen and not important
  'SER' : 2,
  'GLN' : {
    'mp10'  : {lookup['GLN'][3] : [0, 30, 1]}, # really a blob from -30 to 30
    'tm130' : None, # no points
    'tt0'   : {lookup['GLN'][3] : flat}, # tube
    'pp30'  : None, # wraps around 0 - 360
    'mm-40' : None,
    'tp40'  : None,
    'mt0'   : None,
    'tm-30' : None,
    'pm20'  : None,
    'pt0'   : None, # very wide around zero
    },
  'GLU' : {
    'tp30' : {lookup['GLU'][3] : flat}, # chi3 is almost tubular
    'mp0'  : {lookup['GLU'][3] : flat},
    'pm20' : {lookup['GLU'][3] : flat},
    'pp20' : None, # not a sensible blob
    'tm-30': None,
    'mm-30': {lookup['GLU'][3] : flat}, # not a straight tube
    'tt0'  : {lookup['GLU'][3] : flat},
    'mt-10': {lookup['GLU'][3] : flat},
    'pt0'  : {lookup['GLU'][3] : flat}, # straight tube
    },
  'HIS' : 2,
  'ILE' : {
    'pp' : None, # no points
    },
  'LEU' : {
    'tm' : None, # no points
    },
  'LYS' : {
    'mtpm' : None, # no points only 20 examples
    'mtmp' : None,
    'pmtt' : None,
    'pmmt' : None,
    'mptp' : None,
    'mptm' : None,
    'ttpm' : None,
    'ttmp' : None,
    'tmtp' : None,
    'tmmm' : None,
    'tmtm' : None,
    'mttt' : None,
    },
  'MET' : {
    'pmt' : None, # no points
    'tmt' : None,
    'pp-130' : None,
    'mpm' : None,
    'mpt' : None,
    },
  'PRO' : 3, # special case
  'THR' : 1,
  'TRP' : {
    'm-10' : None, # blob near 0
    't60' : None,
    },
  'PHE' : {
    'm-10' : { lookup['PHE'][2] : [ -14.7, 19.8, 2]},
    },
  'TYR' : {
    'm-10' : { lookup['PHE'][2] : [ -14.7, 20.2, 2]},
    },
  'VAL' : 1,
}

rotamer_evaluator = rotamer_eval.RotamerEval(data_version='8000')
rotamer_id = rotamer_eval.RotamerID() # loads in the rotamer names

def generate_chis(chis=None,
                  n=2,
                  step=10,
                  depth=0,
                  starting_chis=None,
                  ):
  if chis is None:
    if starting_chis:
      chis=[]
      for c in starting_chis:
        chis.append(c)
        depth+=1
      while len(chis)<n:
        chis.append(None)
    else:
      chis=[None]*n
  else:
    depth+=1
  for i, angle in enumerate(range(0, 361, step)):
    chis[depth]=angle
    if depth==n-1:
      #if i%1000==0: print "YIELD",i,chis
      yield chis
    else:
      for rc in generate_chis(chis, n, step, depth=depth):
        yield rc

def evaluate(resname, chis):
  name=None
  value = rotamer_evaluator.evaluate(resname, chis)
  if value is not None:
    score=value #*100
    if score<0.02: return None
    #if score<0.003: continue
    wrap_chis = rotamer_id.wrap_chis(resname.strip(), chis,
                                     symmetry=False)
    sym_chis = wrap_chis[:]
    sym_chis = rotamer_id.wrap_sym(resname.strip(), sym_chis)
    #evaluation = self.evaluateScore(value)
    name = rotamer_id.identify(resname,
                               wrap_chis)
  return name

def write_to_pdb(i, sym_chis, f):
  water = "HETATM %4s  O   HOH A%4s    %8.3f%8.3f%8.3f  1.00 22.62           O\n"
  sym_chis = list(sym_chis)
  if len(sym_chis)==2:
    sym_chis.append(0.)
  try:
    f.write(water % (i,i,sym_chis[0],sym_chis[1],sym_chis[2]))
  except: # intentional
    print('WARNING',i, sym_chis)

def loop_on_residue_rotamer(resname,
                            rid,
                            step=10,
                            starting_chis=None,
                           ):
  i=0
  f=open("%s_%s.pdb" % (resname.lower(), rid), "w")
  f.write('CRYST1  360.000  360.000  360.000  90.00  90.00  90.00 P 1           1 ')
  points = []
  n=number_of_chis.get(resname.upper(), None)
  for chis in generate_chis(chis=None,
                            n=n,
                            step=step,
                            starting_chis=starting_chis,
  ):
    name = evaluate(resname,chis)
    #print resname,chis,name,rid
    if name==rid:
      i+=1
      write_to_pdb(i,chis,f)
      points.append(copy.deepcopy(chis))
  f.close()
  return points

def get_min_max(points):
  rc = [None]*len(points[0])
  for i in range(len(rc)):
    tmp = []
    for p in points:
      tmp.append(p[i])
    rc[i] = [min(tmp), max(tmp)]
  return rc

def _order_dihedral_keys(d1,d2):
  order = ["N", "CA", "CB", "CG"]
  if len(d1)!=4: return -1
  if len(d2)!=4: return 1
  assert d1[0] in order
  assert d2[0] in order
  for o in order:
    if d1[0]==o: return -1
    if d2[0]==o: return 1
  assert 0

def run(only_resname=None):
  for resname, rotamers in rdl_database.rdl_database.items():
    if only_resname is not None:
      if only_resname=='new':
        if resname in results: continue
      elif resname!=only_resname: continue
    for rotamer in rotamers:
      starting_chis = None
      if rotamer=='default': continue
      print(resname, rotamer)
      resname_d = results.get(resname, None)
      if type(resname_d)==type({}):
        if rotamer in resname_d:
          continue
      if number_of_chis.get(resname, None)>3:
        print(resname_d)
        keys = list(rdl_database.rdl_database[resname][rotamer].keys())
        keys.sort(_order_dihedral_keys)
        print(keys)
        starting_chis = []
        for key in keys:
          if len(key)!=4: continue
          print(key, rdl_database.rdl_database[resname][rotamer][key])
          starting_chis.append(rdl_database.rdl_database[resname][rotamer][key][0])
        for i in range(3):
          del starting_chis[-1]
        print(starting_chis)
      points = loop_on_residue_rotamer(resname,
                                       rotamer,
                                       starting_chis=starting_chis)
      print('POINTS', len(points))
      if not points and 0:
        #seeds = [-171.7, -75.9, 127.2]
        points = loop_on_residue_rotamer(resname,
                                         rotamer,
                                         step=1,
                                         starting_chis=starting_chis,
        )
      assert points, "there are no points"
      min_max = get_min_max(points)
      print(min_max)

      r=results.get(resname, None)
      for i, (m, n) in enumerate(min_max):
        if type(r)==type(1):
          if i+1==r: continue
        elif type(r)==type({}):
          r=r.get(rotamer, -1)
          if i+1==r: continue
        print(i,m,n,m!=0, n!=360)
        assert not (m==0 and n==360), "\n\n\tphenix.start_coot %s_%s.pdb\n\n" % (
          resname.lower(),
          rotamer,
          )

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
