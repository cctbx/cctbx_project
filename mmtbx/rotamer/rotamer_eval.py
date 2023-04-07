from __future__ import absolute_import, division, print_function
import libtbx.load_env
from libtbx import easy_pickle
from libtbx import dlite
from libtbx.utils import Sorry
from mmtbx.rotamer.sidechain_angles import PropertyFile, SidechainAngles
from mmtbx import monomer_library
import mmtbx.monomer_library.server
from cctbx.array_family import flex
import weakref
import sys, os
import iotbx.pdb
from six.moves import range

def find_rotarama_data_dir(optional=False):
  result = libtbx.env.find_in_repositories(
    os.path.join("chem_data", "rotarama_data"))
  if result is None:
    result = libtbx.env.find_in_repositories("rotarama_data")
    if result is None:
      result = libtbx.env.find_in_repositories(
        os.path.join("ext_ref_files", "rotarama_data"))
      if result is None and not optional:
        raise Sorry("""\
Can't find chem_data/rotarama_data/ directory:
  Inside chem_data, please run
    svn --quiet --non-interactive --trust-server-cert export https://github.com/rlabduke/reference_data.git/trunk/Top8000/Top8000_rotamer_pct_contour_grids rotarama_data
    svn --quiet --non-interactive --trust-server-cert --force export https://github.com/rlabduke/reference_data.git/trunk/Top8000/Top8000_ramachandran_pct_contour_grids rotarama_data
    mmtbx.rebuild_rotarama_cache
  to resolve this problem.""")
  return result

def open_rotarama_dlite(rotarama_data_dir=None):
  if (rotarama_data_dir is None):
    rotarama_data_dir = find_rotarama_data_dir()
  return dlite.try_loading_db(os.path.join(rotarama_data_dir, "rotarama.dlite"))

# maps aa name to file name
aminoAcids = {
    'arg' : 'arg',
    'asn' : 'asn',
    'asp' : 'asp',
    'cys' : 'cys',
    'gln' : 'gln',
    'glu' : 'glu',
    'his' : 'his',
    'ile' : 'ile',
    'leu' : 'leu',
    'lys' : 'lys',
    'met' : 'met',
    'phe' : 'phetyr',
    'pro' : 'pro',
    'ser' : 'ser',
    'thr' : 'thr',
    'trp' : 'trp',
    'tyr' : 'phetyr',
    'val' : 'val',
}

def mon_lib_query(residue, mon_lib_srv):
  # XXX backward compatibility 2007-08-10
  get_func = getattr(mon_lib_srv, "get_comp_comp_id", None)
  if (get_func is not None): return get_func(comp_id=residue)
  return mon_lib_srv.get_comp_comp_id_direct(comp_id=residue)

def eval_residue_completeness(residue, mon_lib_srv, ignore_hydrogens=True):
  atom_list = []
  for atom in residue.atoms():
    atom_list.append(atom.name.strip().upper())
  mlq = mon_lib_query(residue.resname.strip().upper(), mon_lib_srv)
  reference_list = []
  if(not ignore_hydrogens):
    for at in mlq.atom_dict():
      reference_list.append(at)#(atom.name.strip().upper())
  elif(mlq is not None):
    for non in mlq.non_hydrogen_atoms():
      reference_list.append(non.atom_id.strip().upper())
  missing=[]
  for atom in reference_list:
    if atom not in atom_list:
      atom_temp = atom.replace("*", "'")
      if atom.upper() == "O1P":
        atom_temp = "OP1"
      elif atom.upper() == "O2P":
        atom_temp = "OP2"
      if atom_temp not in atom_list:
        missing.append(atom)
  return missing

def eval_sidechain_completeness(pdb_hierarchy,
                                mon_lib_srv=None,
                                ignore_hydrogens=True,
                                report_whole_res=False,
                                return_ca_pos=False):
  missing_atom_list=[]
  if mon_lib_srv is None:
    mon_lib_srv = monomer_library.server.server()
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        conformers = residue_group.conformers()
        for conformer in residue_group.conformers():
          item = []
          residue = conformer.only_residue()
          if conformer.altloc == "":
            key = "%2s%5s %3s" % (chain.id, residue_group.resid(),
              residue.resname)
          else:
            key = "%2s%5s%1s%3s" % (chain.id, residue_group.resid(),
              conformer.altloc, residue.resname)
          ca_xyz = []
          for atom in residue.atoms():
            if atom.name == " CA ":
              ca_xyz = atom.xyz
          missing = eval_residue_completeness(
            residue=residue,
            mon_lib_srv=mon_lib_srv)
          if not report_whole_res:
            if len(missing) > 0:
              item.append(key)
              item.append(missing)
              if return_ca_pos:
                item.append(ca_xyz)
          else:
            item.append(key)
            item.append(missing)
            if return_ca_pos:
              item.append(ca_xyz)
          if len(item) > 0:
            missing_atom_list.append(item)
  return missing_atom_list

class RotamerEval:

  # This is shared among all instances of RotamerEval -- a class variable.
  # It holds a LOT of read-only data, so this helps save memory.
  aaTables = {} # maps "his" to a NDimTable object for histidine, etc.

  def __init__(
               self,
               sidechain_angles=None,
               mon_lib_srv=None,
               log=None,
               data_version="8000"):
#              data_version="500"):
    if sidechain_angles is None:
      sidechain_angles = SidechainAngles(True)
    self.sidechain_angles = sidechain_angles
    if mon_lib_srv is None:
      mon_lib_srv = mmtbx.monomer_library.server.server()
    if log is None:
      log = sys.stdout
    self.data_version = data_version
    if self.data_version == "8000":
      self.outlier_threshold = 0.003
      fileprefix = "rota8000-"
    else: raise Sorry("data_version given to RotamerEval not recognized.")
    self.log = log
    self.mon_lib_srv = mon_lib_srv
    self.rot_id = RotamerID()
    main_aaTables = RotamerEval.aaTables
    self.aaTables = {}
    for aa,ndt_weakref in main_aaTables.items():
        # convert existing weak references to strong references
        self.aaTables[aa] = ndt_weakref()
    rotamer_data_dir = find_rotarama_data_dir()
    no_update = os.path.exists(os.path.join(rotamer_data_dir, "NO_UPDATE"))
    target_db = open_rotarama_dlite(rotarama_data_dir=rotamer_data_dir)
    for aa, aafile in aminoAcids.items():
      if (self.aaTables.get(aa) is not None): continue
      data_file = fileprefix+aafile+".data"
      pickle_file = fileprefix+aafile+".pickle"
      pair_info = target_db.pair_info(
                    source_path=data_file,
                    target_path=pickle_file,
                    path_prefix=rotamer_data_dir)
      if (((pair_info.needs_update) and (not no_update)) or not
          os.path.exists(os.path.join(rotamer_data_dir, pickle_file))):
        raise Sorry(
          "chem_data/rotarama_data/*.pickle files are missing or out of date.\n"
          "  Please run\n"
          "    mmtbx.rebuild_rotarama_cache\n"
          "  to resolve this problem.\n")
      ndt = easy_pickle.load(file_name=os.path.join(
              rotamer_data_dir, pair_info.target.path))
      self.aaTables[aa] = ndt
      main_aaTables[aa] = weakref.ref(ndt)

  def evaluate(self, aaName, chiAngles):
    '''Evaluates the specified rotamer from 0.0 (worst) to 1.0 (best).

       Values below 0.003 are generally considered outliers.
       If the 3-letter amino acid name is not recognized, returns None.'''
    ndt = self.aaTables.get(aaName.lower())
    if (ndt is None): return None
    return ndt.valueAt(chiAngles)

  def get_atom_dict(self, residue):
    atom_dict = {}
    atoms = residue.atoms()
    for atom in atoms:
      #handle hydrogen/deuterium swaps
      if atom_dict.get(atom.name) == None:
        if atom_dict.get(atom.name.replace("H","D",1)) != None:
          del(atom_dict[atom.name.replace("H","D",1)])
        elif atom_dict.get(atom.name.replace("D","H",1)) != None:
          del(atom_dict[atom.name.replace("D","H",1)])
        atom_dict[atom.name] = atom
    return atom_dict

  def chi_angles(self, residue):
    atom_dict = self.get_atom_dict(residue)
    return self.sidechain_angles.measureChiAngles(
      res=residue,
      atom_dict=atom_dict)

  def evaluate_residue(
                       self,
                       residue=None,
                       residue_group=None): # FIXME does not work! Start using this parameter in function?
    # Warning!!! Returns "OUTLIER", "UNCLASSIFIED", rotamer name or None.
    assert [residue, residue_group].count(None) == 1
    if residue is not None:
      atoms = residue.atoms()
      resname = residue.resname.lower().strip()
    if resname == 'gly': # why ala is not here?
      return None
    elif resname == 'mse':
      resname = 'met'
    atom_dict = self.get_atom_dict(residue)
    try:
      chis = self.sidechain_angles.measureChiAngles(
               res=residue,
               atom_dict=atom_dict)
      value = self.evaluate(
                resname,
                chis)
    except Exception:
      return None
    if chis is None:
      return None
    wrap_chis = \
      self.rot_id.wrap_chis(resname, chis, symmetry=False)
    rotamer_name = self.rot_id.identify(resname, wrap_chis)
    if(rotamer_name == "EXCEPTION"):
      assert value is None
      return rotamer_name
    if rotamer_name == "" and (value >= self.outlier_threshold):
      return "UNCLASSIFIED"
    elif (value < self.outlier_threshold):
      return "OUTLIER"
    else:
      return rotamer_name

  def evaluate_residue_2(self, residue=None):
    # copy-paste from evaluate_residue, returns
    # "OUTLIER", "Allowed", "Favored" or None if something is really wrong.
    from mmtbx.validation.rotalyze import OUTLIER_THRESHOLD, ALLOWED_THRESHOLD
    if residue is not None:
      atoms = residue.atoms()
      resname = residue.resname.lower().strip()
    if resname == 'gly' or resname == 'ala':
      return None
    elif resname == 'mse':
      resname = 'met'
    atom_dict = self.get_atom_dict(residue)
    try:
      chis = self.sidechain_angles.measureChiAngles(
               res=residue,
               atom_dict=atom_dict)
      value = self.evaluate(
                resname,
                chis)
    except Exception:
      return None
    if chis is None:
      return "OUTLIER"

    if value >= ALLOWED_THRESHOLD :
      return "Favored"
    elif value >= OUTLIER_THRESHOLD:
      return "Allowed"
    else:
      return "OUTLIER"


  def nearest_rotamer_sites_cart(self, residue):
    sites_cart_result = residue.atoms().extract_xyz()
    get_class = iotbx.pdb.common_residue_names_get_class
    if(get_class(residue.resname) == "common_amino_acid"):
      sites_cart = residue.atoms().extract_xyz()
      rotamer_iterator = self.mon_lib_srv.rotamer_iterator(
          fine_sampling = True,
          comp_id       = residue.resname,
          atom_names    = residue.atoms().extract_name(),
          sites_cart    = sites_cart)
      if(rotamer_iterator is None or
         rotamer_iterator.problem_message is not None or
         rotamer_iterator.rotamer_info is None):
        rotamer_iterator = None
      if(rotamer_iterator is not None):
        dist_min = 1.e+9
        for r, rotamer_sites_cart in rotamer_iterator:
          d= flex.mean(flex.sqrt((sites_cart - rotamer_sites_cart).dot()))
          if(d < dist_min):
            dist_min = d
            sites_cart_result = rotamer_sites_cart
    return sites_cart_result

#{{{ RotamerID (new for reading in rotamer names from rotamer_names.props)
class RotamerID:

  names = {}

  def __init__(self):
    if (len(self.names) == 0):
      source_dir = self.find_source_dir()
      #f = PropertyFile()
      # can't use f.properties to read in rotamer_names.props
      # some of the rotamer names aren't unique, so they get dropped as keys!
      rota_names_list = self.process(os.path.join(source_dir, "rotamer_names.props"))
      for line in rota_names_list:
        split_line = line.split("=")
        aa_name = split_line[0].strip()
        ranges = split_line[1].strip().strip("\"")
        name_split = aa_name.split(" ")
        aa = name_split[0]
        rot_name = name_split[1]
        rot = NamedRot(aa, rot_name, ranges)
        rotList = []
        if aa in self.names:
          rotList = self.names[aa]
        rotList.append(rot)
        self.names[aa] = rotList

  def identify(self, aa_name, chis):
    aa_name = aa_name.lower()
    if(aa_name == "ala"): return "EXCEPTION"
    if aa_name not in self.names:
      raise Sorry("Unknown residue name: %s", aa_name)
    wrap_chis = self.wrap_chis(aa_name, chis)
    rotList = self.names[aa_name]
    for rot in rotList:
      if(rot.contains(wrap_chis)): return rot.rotamer_name
    return ""

  def find_source_dir(optional=False):
    result = libtbx.env.find_in_repositories(os.path.join("mmtbx", "rotamer"))
    if result is None and not optional:
      raise Sorry("""\
Can't seem to find mmtbx/rotamer/ directory.
    """)
    return result

  def process(self, fileLoc):
    rotaList = []
    try:
      f = open(fileLoc)
    except ImportError as e:
      print(fileLoc+" file not found")
      sys.exit()
    for line in f:
      if (line.startswith("#") or line == "\n"): continue
      else: rotaList.append(line)
    f.close()
    return rotaList

  def wrap_chis(self, aa_name, chis, symmetry=True):
    aa_name = aa_name.lower()
    wrap_chis = []
    for i in range(0, len(chis)):
      if chis[i] is not None:
        wrap_chis.append(chis[i] % 360)
        if wrap_chis[i] < 0:
          wrap_chis[i] += 360
      else:
        wrap_chis.append(None)
    if (symmetry==True):
      wrap_chis = self.wrap_sym(aa_name, wrap_chis)
    #MOVED TO SEPARATE FUNCTION 'wrap_sym' for accurate angle reporting
    #if (aa_name == "asp" or aa_name == "glu" or aa_name == "phe" or aa_name == "tyr"):
    #  i = len(wrap_chis) - 1
    #  print wrap_chis[i]
    #  wrap_chis[i] = wrap_chis[i] % 180
    #  if wrap_chis[i] < 0:
    #    wrap_chis[i] += 180
    return wrap_chis

  def wrap_sym(self, aa_name, wrap_chis):
    aa_name = aa_name.lower()
    if (aa_name == "asp" or aa_name == "glu" or aa_name == "phe" or aa_name == "tyr"):
      i = len(wrap_chis) - 1
      if wrap_chis[i] is not None:
        wrap_chis[i] = wrap_chis[i] % 180
        if wrap_chis[i] < 0:
          wrap_chis[i] += 180
    return wrap_chis

#}}}

#{{{ NamedRot
class NamedRot:

  def __init__(self, aa, rotamer_name, bounds):
    self.aa_name = aa
    self.rotamer_name = rotamer_name
    self.bounds = [int(b) for b in bounds.split(", ")]

  def __str__(self):
    return str(self.rotamer_name) + "=" + str(self.bounds)

  def contains(self, angles):
    for i in range(0, len(self.bounds), 2):
      if (angles[i//2] is None ) or (self.bounds[i] is None): # XXX FIX OK?
        return False
      if (   angles[i//2] < self.bounds[i]
          or angles[i//2] > self.bounds[i+1]): return False
    return True
#}}}

########################################################################
def exercise(args):
  if (find_rotarama_data_dir(optional=True) is None):
    print("Skipping exercise(): rotarama_data directory not available")
  else:
    from mmtbx.command_line import rebuild_rotarama_cache
    rebuild_rotarama_cache.run()
    #
    from libtbx.test_utils import approx_equal
    #
    verbose = ("--verbose" in args)
    #
    r = RotamerEval()
    tbl = r.aaTables['val']
    assert RotamerEval().aaTables['val'] is tbl
    #
    assert tbl.whereIs([0.5]) == [0]
    assert tbl.bin2index([0]) == 0
    if r.data_version == '8000' :
      aeql = [0.9346579, 0.8926509, 0.8296402, 0.7876331,
              0.7351242, 0.6721136, 0.6511100, 0.6406083,
              0.5670958, 0.4620781, 0.4200710, 0.3885657,
              0.3570603, 0.2625444, 0.1785302, 0.1155195,
              0.0945160, 0.0945160, 0.0735124, 0.0735124]
    assert approx_equal(
     [(y*1000) for y in tbl.lookupTable[0:20]],aeql)
    assert approx_equal(r.evaluate("SER", [60]), 0.759436935186)
    #
    # Based off new (March 2015) NDFTs built from top8000-angles Makefile
    # Remaining inaccuracies are due to dihedrals being rounded off to
    # one decimal place!
    assert r.data_version == '8000'
    for aminoAcid, chiAngles, molpValue in [
      ("MET", [80.4, -172.2, 177.5] ,  9.7),
      ("GLN", [166.0, 178.0, -107.4] ,  9.7),
      ("ILE", [60.3, 162.4] , 24.2),
      ("PHE", [-60.7, 97.9] , 94.4),
      ("VAL", [-179.8] , 58.2),
      ("LYS", [-175.6, 176.2, -172.0, -174.2] , 81.8),
      ("THR", [76.7] , 10.3),
      ("LEU", [-68.2, -165.8] ,  9.1),
      ("THR", [70.7] , 25.9),
      ("LYS", [-179.3, -179.4, -151.1, -49.3] , 13.1),
      ("THR", [-63.4] , 83.9),
      ("ILE", [125.7, -175.4] ,  0.0),
      ("THR", [66.5] , 45.3),
      ("LEU", [-117.8, 30.2] ,  0.0),
      ("GLU", [-75.1, -167.9, 139.8] , 34.9),
      ("VAL", [-62.5] , 32.5),
      ("GLU", [-73.9, -54.5, -18.4] , 29.9),
      ("PRO", [-29.0] , 90.6),
      ("SER", [35.7] ,  0.5),
      ("ASP", [-80.6, -19.8] , 55.0),
      ("THR", [60.6] , 79.8),
      ("ILE", [-60.9, -54.6] , 35.5),
      ("GLU", [-169.6, -175.1, 72.8] , 30.3),
      ("ASN", [177.5, 53.8] , 30.3),
      ("VAL", [168.2] , 42.4),
      ("LYS", [-71.7, -173.9, 179.2, 179.4] , 94.0),
      ("LYS", [-60.8, 169.3, 148.9, -89.1] ,  5.7),
      ("ILE", [-70.9, 166.5] , 75.5),
      ("GLN", [176.9, 171.9, 35.2] , 55.5),
      ("ASP", [-150.1, 65.5] ,  1.0),
      ("LYS", [78.3, 138.2, 62.4, -165.4] ,  0.7),
      ("GLU", [-60.1, -76.8, -36.2] , 36.3),
      ("ILE", [-54.4, 161.0] , 30.8),
      ("PRO", [-31.6] , 60.3),
      ("PRO", [-28.4] , 96.3),
      ("ASP", [134.6, -61.7] ,  0.0),
      ("GLN", [-61.9, -179.0, -165.3] , 13.0),
      ("GLN", [-53.1, -179.9, 28.0] , 47.2),
      ("ARG", [161.7, 173.6, 174.2, -106.7] , 20.5),
      ("LEU", [-68.3, 166.9] , 81.8),
      ("ILE", [-48.9, -58.1] , 29.8),
      ("PHE", [178.0, 78.2] , 93.1),
      ("LYS", [-61.5, 173.7, -111.9, -58.8] ,  2.2),
      ("GLN", [-172.6, 177.3, 118.5] ,  9.7),
      ("LEU", [-50.8, -172.7] , 17.3),
      ("GLU", [173.0, 141.4, 172.4] ,  2.4),
      ("ASP", [-78.0, 177.6] , 74.9),
      ("ARG", [-55.9, -71.7, 114.2, -128.0] ,  0.4),
      ("THR", [59.8] , 74.3),
      ("LEU", [-60.3, -179.2] , 78.5),
      ("SER", [59.4] , 73.3),
      ("ASP", [-73.0, 157.0] , 84.8),
      ("TYR", [-63.3, 103.6] , 90.4),
      ("ASN", [-159.1, -145.0] ,  1.2),
      ("ILE", [-69.6, 176.4] , 57.3),
      ("GLN", [-79.4, -161.7, -148.4] ,  2.4),
      ("LYS", [49.7, 165.7, 154.3, 72.9] ,  3.2),
      ("GLU", [-72.2, 126.6, 36.7] ,  0.6),
      ("SER", [-73.2] , 18.8),
      ("THR", [-60.6] , 94.7),
      ("LEU", [-43.5, -170.9] ,  4.6),
      ("HIS", [-69.1, -88.8] , 80.8),
      ("LEU", [172.8, 65.9] , 43.9),
      ("VAL", [177.0] , 83.7),
      ("LEU", [-108.1, 39.2] ,  0.0),
      ("ARG", [133.9, -155.8, 27.2, -152.9] ,  0.0),
      ("LEU", [-92.5, 37.5] ,  0.1),
      ("ARG", [-146.6, 157.6, 92.9, -95.5] ,  0.9),
    ]:
      r_eval = 100*r.evaluate(aminoAcid, chiAngles)
      if (verbose):
        print(aminoAcid, "%4.1f %4.1f %4.1f" % (
          r_eval, molpValue, r_eval-molpValue))
        #print '      ("%s",' % aminoAcid, chiAngles, ',', "%4.1f)," % r_eval
    # assert approx_equal(r_eval, molpValue, eps=0.9)
    #
    # check if tables are cleared from memory if all RotamerEval instances
    # are gone
    for aa,ndt_weakref in RotamerEval.aaTables.items():
      assert ndt_weakref() is not None
    del r
    del tbl
    for aa,ndt_weakref in RotamerEval.aaTables.items():
      assert ndt_weakref() is None
    #
  print("OK")

if (__name__ == "__main__"):
    exercise(sys.argv[1:])
