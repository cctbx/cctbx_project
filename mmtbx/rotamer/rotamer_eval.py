import libtbx.load_env
from mmtbx.rotamer.n_dim_table import NDimTable
from libtbx import easy_pickle
from libtbx import dlite
from libtbx.utils import Sorry
from mmtbx.rotamer.sidechain_angles import PropertyFile
from iotbx.pdb import common_residue_names_get_class
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import weakref
import sys, os

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
  Please run
    svn co svn://quiddity.biochem.duke.edu:21/phenix/rotarama_data
  to resolve this problem.""")
  return result

def open_rotarama_dlite(rotarama_data_dir=None):
  if (rotarama_data_dir is None):
    rotarama_data_dir = find_rotarama_data_dir()
  return dlite.target_db(os.path.join(rotarama_data_dir, "rotarama.dlite"))

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
            key = chain.id+residue_group.resid()+" "+residue.resname
          else:
            key = chain.id+residue_group.resid()+conformer.altloc+residue.resname
          atom_list = []
          ca_xyz = []
          for atom in residue.atoms():
            atom_list.append(atom.name.strip().upper())
            if atom.name == " CA ":
              ca_xyz = atom.xyz
          mlq = mon_lib_query(residue.resname.strip().upper(), mon_lib_srv)
          reference_list = []
          if(not ignore_hydrogens):
            for at in mlq.atom_dict():
              atom = at.replace("*", "'")
              if atom.upper() == "O1P":
                atom = "OP1"
              elif atom.upper() == "O2P":
                atom = "OP2"
              reference_list.append(atom.strip().upper())
          elif (mlq is not None) :
            for non in mlq.non_hydrogen_atoms():
              atom = non.atom_id.replace("*", "'")
              if atom.upper() == "O1P":
                atom = "OP1"
              elif atom.upper() == "O2P":
                atom = "OP2"
              reference_list.append(atom.strip().upper())
          missing=[]
          for atom in reference_list:
            if atom not in atom_list:
              missing.append(atom)
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

    def __init__(self):
        main_aaTables = RotamerEval.aaTables
        self.aaTables = {}
        for aa,ndt_weakref in main_aaTables.items():
            # convert existing weak references to strong references
            self.aaTables[aa] = ndt_weakref()
        rotamer_data_dir = find_rotarama_data_dir()
        target_db = open_rotarama_dlite(rotarama_data_dir=rotamer_data_dir)
        for aa, aafile in aminoAcids.items():
                if (self.aaTables.get(aa) is not None): continue
                data_file = "rota500-"+aafile+".data"
                pickle_file = "rota500-"+aafile+".pickle"
                pair_info = target_db.pair_info(
                  source_path=data_file,
                  target_path=pickle_file,
                  path_prefix=rotamer_data_dir)
                if pair_info.needs_update:
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

        Values below 0.01 are generally considered outliers.
        If the 3-letter amino acid name is not recognized, returns None.'''
        ndt = self.aaTables.get(aaName.lower())
        if (ndt is None): return None
        return ndt.valueAt(chiAngles)

#{{{ RotamerID (new for reading in rotamer names from rotamer_names.props)
class RotamerID:

  names = {}

  def __init__(self):
    self.names = {}
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
    if aa_name not in self.names: raise Sorry("Unknown residue name")
    wrap_chis = self.wrap_chis(aa_name, chis)
    rotList = self.names[aa_name]
    for rot in rotList:
      #print rot
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
    except ImportError, e:
      print fileLoc+" file not found"
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
      wrap_chis.append(chis[i] % 360)
      if wrap_chis[i] < 0:
        wrap_chis[i] += 360
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
    self.bounds = map(int, bounds.split(", "))

  def __str__(self):
    return str(self.rotamer_name) + "=" + str(self.bounds)

  def contains(self, angles):
    for i in range(0, len(self.bounds), 2):
      #print str(angles[i/2]) + " is between " + str(self.bounds[i]) + " and " + str(self.bounds[i+1])
      if (angles[i/2] < self.bounds[i] or angles[i/2] > self.bounds[i+1]): return False
    return True
#}}}

def exercise(args):
  if (find_rotarama_data_dir(optional=True) is None):
    print "Skipping exercise(): rotarama_data directory not available"
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
    assert approx_equal(
     [(y*1000) for y in tbl.lookupTable[0:20]],
     [0, 5.468854, 5.202081, 4.935307, 4.935307, 4.801921, 4.535147,
      4.268374, 3.734827, 3.468054, 3.20128, 2.534347, 2.40096, 2.40096,
      2.267574, 2.134187, 2.0008, 2.0008, 2.0008, 1.734027])
    assert approx_equal(r.evaluate("SER", [60]), 0.800816386938)
    #
    # Based off new (Oct 2006) NDFTs built from top500-angles Makefile
    # Remaining inaccuracies are due to dihedrals being rounded off to
    # one decimal place!
    for aminoAcid, chiAngles, molpValue in [
      ("MET", [80.4, -172.2, 177.5], 17.9),
      ("GLN", [166.0, 178.0, -107.4], 17.6),
      ("ILE", [60.3, 162.4], 26.3),
      ("PHE", [-60.7, 97.9], 96.3),
      ("VAL", [-179.8], 60.3),
      ("LYS", [-175.6, 176.2, -172.0, -174.2], 84.8),
      ("THR", [76.7], 11.8),
      ("LEU", [-68.2, -165.8], 17.3),
      ("THR", [70.7], 29.1),
      ("LYS", [-179.3, -179.4, -151.1, -49.3], 35.4),
      ("THR", [-63.4], 68.6),
      ("ILE", [125.7, -175.4], 0.1),
      ("THR", [66.5], 44.9),
      ("LEU", [-117.8, 30.2], 0.2),
      ("GLU", [-75.1, -167.9, 139.8], 53.4),
      ("VAL", [-62.5], 35.6),
      ("GLU", [-73.9, -54.5, -18.4], 54.2),
      ("PRO", [-29.0], 84.9),
      ("SER", [35.7], 2.0),
      ("ASP", [-80.6, -19.8], 69.2),
      ("THR", [60.6], 82.6),
      ("ILE", [-60.9, -54.6], 39.2),
      ("GLU", [-169.6, -175.1, 72.8], 36.2),
      ("ASN", [177.5, 53.8], 39.0),
      ("VAL", [168.2], 47.7),
      ("LYS", [-71.7, -173.9, 179.2, 179.4], 96.6),
      ("LYS", [-60.8, 169.3, 148.9, -89.1], 18.8),
      ("ILE", [-70.9, 166.5], 75.3),
      ("GLN", [176.9, 171.9, 35.2], 57.2),
      ("ASP", [-150.1, 65.5], 2.7),
      ("LYS", [78.3, 138.2, 62.4, -165.4], 6.6),
      ("GLU", [-60.1, -76.8, -36.2], 58.3),
      ("ILE", [-54.4, 161.0], 36.6),
      ("PRO", [-31.6], 78.4),
      ("PRO", [-28.4], 87.4),
      ("ASP", [134.6, -61.7], 0.0),
      ("GLN", [-61.9, -179.0, -165.3], 24.9),
      ("GLN", [-53.1, -179.9, 28.0], 67.4),
      ("ARG", [161.7, 173.6, 174.2, -106.7], 40.7),
      ("LEU", [-68.3, 166.9], 78.6),
      ("ILE", [-48.9, -58.1], 30.2),
      ("PHE", [178.0, 78.2], 93.4),
      ("LYS", [-61.5, 173.7, -111.9, -58.8], 8.3),
      ("GLN", [-172.6, 177.3, 118.5], 18.4),
      ("LEU", [-50.8, -172.7], 25.8),
      ("GLU", [173.0, 141.4, 172.4], 6.9),
      ("ASP", [-78.0, 177.6], 76.2),
      ("ARG", [-55.9, -71.7, 114.2, -128.0], 0.1),
      ("THR", [59.8], 72.1),
      ("LEU", [-60.3, -179.2], 84.5),
      ("SER", [59.4], 78.7),
      ("ASP", [-73.0, 157.0], 90.7),
      ("TYR", [-63.3, 103.6], 92.2),
      ("ASN", [-159.1, -145.0], 2.2),
      ("ILE", [-69.6, 176.4], 68.0),
      ("GLN", [-79.4, -161.7, -148.4], 9.6),
      ("LYS", [49.7, 165.7, 154.3, 72.9], 15.3),
      ("GLU", [-72.2, 126.6, 36.7], 3.4),
      ("SER", [-73.2], 31.8),
      ("THR", [-60.6], 97.0),
      ("LEU", [-43.5, -170.9], 11.6),
      ("HIS", [-69.1, -88.8], 90.4),
      ("LEU", [172.8, 65.9], 50.8),
      ("VAL", [177.0], 91.3),
      ("LEU", [-108.1, 39.2], 0.4),
      ("ARG", [133.9, -155.8, 27.2, -152.9], 0.0),
      ("LEU", [-92.5, 37.5], 2.0),
      ("ARG", [-146.6, 157.6, 92.9, -95.5], 3.2),
    ]:
      r_eval = 100*r.evaluate(aminoAcid, chiAngles)
      if (verbose):
        print aminoAcid, "%4.1f %4.1f %4.1f" % (
          r_eval, molpValue, r_eval-molpValue)
      assert approx_equal(r_eval, molpValue, eps=0.9)
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
  print "OK"

if (__name__ == "__main__"):
    exercise(sys.argv[1:])
