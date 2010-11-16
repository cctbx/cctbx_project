class FormatError(RuntimeError): pass

def chain_id_strip(id):
  if (id[0] == " "): return id[1]
  return id

class _base(object):

  def __init__(self, pdb_str):
    if (len(pdb_str) < 80): pdb_str += " " * (80-len(pdb_str))
    self.pdb_str = pdb_str

  def raise_FormatError(self, msg):
    raise FormatError('%s:\n  "%s"' % (msg, self.pdb_str))

class header(_base):

  def __init__(self, pdb_str):
    # 11 - 50  String(40)    classification  Classifies the molecule(s)
    # 51 - 59  Date          depDate         Deposition date.  This is the date
    #                                        the coordinates were received by
    #                                        the PDB
    # 63 - 66  IDcode        idCode          This identifier is unique within
    #                                        PDB
    _base.__init__(self, pdb_str)
    self.classification = self.pdb_str[10:50]
    self.depdate = self.pdb_str[50:59]
    self.idcode = self.pdb_str[62:66]

class expdta(_base):

  def __init__(self, pdb_str):
    #  9 - 10  Continuation   continuation  Allows concatenation of
    #                                       multiple records.
    # 11 - 70  List           technique     The experimental technique(s)
    #                                       with optional comment
    #                                       describing the sample
    #                                       or experiment.
    _base.__init__(self, pdb_str)
    self.continuation = self.pdb_str[8:10]
    self.technique = self.pdb_str[10:70]
    TECH = self.technique.upper()
    self.keywords = []
    for keyword in ("ELECTRON DIFFRACTION",
                    "FIBER DIFFRACTION",
                    "FLUORESCENCE TRANSFER",
                    "NEUTRON DIFFRACTION",
                    "NMR",
                    "THEORETICAL MODEL",
                    "X-RAY DIFFRACTION"):
      if (TECH.find(keyword) >= 0):
        self.keywords.append(keyword)

class remark_002(_base):

  def __init__(self, pdb_str):
    _base.__init__(self, pdb_str)
    self.text = self.pdb_str[11:70]
    flds = self.pdb_str[11:38].split()
    self.resolution = None
    if (    len(flds) == 3
        and flds[0].upper() == "RESOLUTION."
        and flds[2].upper() == "ANGSTROMS."):
      try: self.resolution = float(flds[1])
      except ValueError: pass

class cryst1(_base):

  def __init__(self, pdb_str):
    #  7 - 15       Real(9.3)      a             a (Angstroms).
    # 16 - 24       Real(9.3)      b             b (Angstroms).
    # 25 - 33       Real(9.3)      c             c (Angstroms).
    # 34 - 40       Real(7.2)      alpha         alpha (degrees).
    # 41 - 47       Real(7.2)      beta          beta (degrees).
    # 48 - 54       Real(7.2)      gamma         gamma (degrees).
    # 56 - 66       LString        sGroup        Space group.
    # 67 - 70       Integer        z             Z value.
    _base.__init__(self, pdb_str)
    ps = self.pdb_str
    self.ucparams = [ps[ 6:15], ps[15:24], ps[24:33],
                     ps[33:40], ps[40:47], ps[47:54]]
    self.sgroup = ps[55:66].strip()
    z = ps[66:70]
    if (len(" ".join(self.ucparams).strip()) == 0):
      self.ucparams = None
    else:
      try: self.ucparams = [float(u) for u in self.ucparams]
      except ValueError:
        self.raise_FormatError("Corrupt unit cell parameters")
    if (len(self.sgroup) == 0): self.sgroup = None
    if (z.strip()):
      try: self.z = int(z)
      except ValueError: self.raise_FormatError("Corrupt Z value")
    else: self.z = None

class scalen(_base):

  def __init__(self, pdb_str):
    #  1 -  6       Record name    "SCALEn"       n=1, 2, or 3
    # 11 - 20       Real(10.6)     s[n][1]        Sn1
    # 21 - 30       Real(10.6)     s[n][2]        Sn2
    # 31 - 40       Real(10.6)     s[n][3]        Sn3
    # 46 - 55       Real(10.5)     u[n]           Un
    _base.__init__(self, pdb_str)
    self.n = int(self.pdb_str[5])
    values = []
    for i in [10,20,30,45]:
      fld = self.pdb_str[i:i+10]
      if (len(fld.strip()) == 0):
        value = 0
      else:
        try: value = float(fld)
        except ValueError:
          raise self.raise_FormatError("Invalid floating-point value")
      values.append(value)
    self.sn1, self.sn2, self.sn3, self.un = values

class conect(_base):

  def __init__(self, pdb_str):
    #  7 - 11  Integer   serial          Atom serial number
    # 12 - 16  Integer   serial          Serial number of bonded atom
    # 17 - 21  Integer   serial          Serial number of bonded atom
    # 22 - 26  Integer   serial          Serial number of bonded atom
    # 27 - 31  Integer   serial          Serial number of bonded atom
    # 32 - 36  Integer   serial          Serial number of hydrogen bonded atom
    # 37 - 41  Integer   serial          Serial number of hydrogen bonded atom
    # 42 - 46  Integer   serial          Serial number of salt bridged atom
    # 47 - 51  Integer   serial          Serial number of hydrogen bonded atom
    # 52 - 56  Integer   serial          Serial number of hydrogen bonded atom
    # 57 - 61  Integer   serial          Serial number of salt bridged atom
    _base.__init__(self, pdb_str)
    self.serial = self.pdb_str[6:11]
    ba = []
    hba = []
    sba = []
    for i,l in [(11, ba), (16, ba), (21, ba), (26, ba),
                (31, hba), (36, hba), (41, sba),
                (46, hba), (51, hba), (56, sba)]:
      fld = self.pdb_str[i:i+5]
      if (fld != "     "): l.append(fld)
    self.serial_numbers_bonded_atoms = ba
    self.serial_numbers_hydrogen_bonded_atoms = hba
    self.serial_numbers_salt_bridged_atoms = sba

class link(_base):

  def __init__(self, pdb_str):
    # 13 - 16      Atom            name1       Atom name.
    # 17           Character       altLoc1     Alternate location indicator.
    # 18 - 20      Residue name    resName1    Residue name.
    # 21 - 22                      chainID1    Chain identifier.
    # 23 - 26                      resSeq1     Residue sequence number.
    # 27           AChar           iCode1      Insertion code.
    # 31 - 40      distance (REFMAC extension: F10.5)
    # 43 - 46      Atom            name2       Atom name.
    # 47           Character       altLoc2     Alternate location indicator.
    # 48 - 50      Residue name    resName2    Residue name.
    # 51 - 52                      chainID2    Chain identifier.
    # 53 - 56                      resSeq2     Residue sequence number.
    # 57           AChar           iCode2      Insertion code.
    # 60 - 65      SymOP           sym1        Symmetry operator for 1st atom.
    # 67 - 72      SymOP           sym2        Symmetry operator for 2nd atom.
    # 73 - 80      margin (REFMAC extension: _chem_link.id)
    _base.__init__(self, pdb_str)
    self.name1 = self.pdb_str[12:16]
    self.altloc1 = self.pdb_str[16]
    self.resname1 = self.pdb_str[17:20]
    self.chain_id1 = chain_id_strip(self.pdb_str[20:22])
    self.resseq1 = self.pdb_str[22:26]
    self.icode1 = self.pdb_str[26]
    try: self.distance = float(self.pdb_str[30:40])
    except ValueError: self.distance = None
    self.name2 = self.pdb_str[42:46]
    self.altloc2 = self.pdb_str[46]
    self.resname2 = self.pdb_str[47:50]
    self.chain_id2 = chain_id_strip(self.pdb_str[50:52])
    self.resseq2 = self.pdb_str[52:56]
    self.icode2 = self.pdb_str[56]
    self.sym1 = self.pdb_str[59:65]
    self.sym2 = self.pdb_str[66:72]
    self.margin = self.pdb_str[72:80]

class sltbrg(link): pass

class ssbond(_base):

  def __init__(self, pdb_str):
    #  8 - 10    Integer         serNum      Serial number.
    # 12 - 14    LString(3)      "CYS"       Residue name.
    # 15 - 16                    chainID1    Chain identifier.
    # 18 - 21                    seqNum1     Residue sequence number.
    # 22         AChar           icode1      Insertion code.
    # 26 - 28    LString(3)      "CYS"       Residue name.
    # 29 - 30                    chainID2    Chain identifier.
    # 32 - 35                    seqNum2     Residue sequence number.
    # 36         AChar           icode2      Insertion code.
    # 60 - 65    SymOP           sym1        Symmetry operator for 1st residue.
    # 67 - 72    SymOP           sym2        Symmetry operator for 2nd residue.
    # 73 - 80    margin (REFMAC extension: _chem_link.id)
    _base.__init__(self, pdb_str)
    self.sernum = self.pdb_str[7:10]
    self.resname1 = self.pdb_str[11:14]
    self.chain_id1 = chain_id_strip(self.pdb_str[14:16])
    self.resseq1 = self.pdb_str[17:21]
    self.icode1 = self.pdb_str[21]
    self.resname2 = self.pdb_str[25:28]
    self.chain_id2 = chain_id_strip(self.pdb_str[28:30])
    self.resseq2 = self.pdb_str[31:35]
    self.icode2 = self.pdb_str[35]
    self.sym1 = self.pdb_str[59:65]
    self.sym2 = self.pdb_str[66:72]
    self.margin = self.pdb_str[72:80]
