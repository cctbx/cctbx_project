from __future__ import absolute_import, division, print_function
import six
from six.moves import range
from six.moves import zip

class pcr_parser(object):

  def __init__(self, inputstring=None, pcrfile=None, enforce_valid_comments=True):
    """Function to parse a FullProf command file into a ConfigObj

    It was modeled according to the FullProf manual and actual output
    files from FullProf and still far from being complete.

    :param inputstring: a string containing a whole commented pcr file like fp2k \
    would format it.
    :type inputstring: string
    """
    from configobj import ConfigObj
    self.LO = 0
    self.cfg = ConfigObj()
    self.lines = inputstring.splitlines()
    self.enforce_valid_comments = enforce_valid_comments
    for line in self.lines:
      print(line) # for debugging
    if not self.lines[1].startswith("!"):
      if self.enforce_valid_comments:
        raise ValueError("""\
        This function only works with correctly commented FullProf files!
        Please run fp2k on your input file first.""")
      else:
        for line in self.lines.copy():
          if line.startswith("!"):
            self.lines.remove(line)
    self._parse()

  def __get_next_line(self, raw=False):
    self.LO += 1    # (increment line counter)
    if raw:
      return self.lines[self.LO]
    else:
      return self.lines[self.LO].split()

  def __parse_line(self, dict_to_append_to, KEYS, typeselect=" "):
    linevals = self.__get_next_line()
    i = 0
    for val in KEYS.split(", "):
      if len(linevals) < 1: break
      if len(typeselect) > i:
        typesel = typeselect[i].lower()
      else:
        typsel = typeselect[-1].lower()
      if typesel == "i":   # integer
        dict_to_append_to[val.strip()] = int(linevals.pop(0))
      elif typesel == "f": # float
        dict_to_append_to[val.strip()] = float(linevals.pop(0))
      else:
        dict_to_append_to[val.strip()] = linevals.pop(0)
      i += 1

  def __check_comment(self, expected_comment="!"):
    if self.enforce_valid_comments:
      self.LO += 1
      line = self.lines[self.LO]
      if not line.startswith(expected_comment):
        raise ValueError("Expected to find comment '{comm}',\
                         but found '{line}'".format(expected_comment, line))

  def extract_structures(self):
    structures = {}
    for i in range(len(self.cfg['phases'])):
      structures[i+1]= self.extract_structure(phase_nr=i+1)
    return structures


  def extract_structure(self, phase_nr=1):
    """This method tries to extract the crystal structure from the parsed pcrfile.

    :returns: the extracted structure
    :rtype: cctbx.xray.structure
    """
    from cctbx import xray
    from cctbx import crystal
    from cctbx.array_family import flex
    p = self.cfg['phases'][str(phase_nr)]
    atoms = p['atoms']
    unit_cell = [ p['a'], p['b'], p['c'], p['alpha'], p['beta'], p['gamma'] ]
    special_position_settings=crystal.special_position_settings(
              crystal_symmetry=crystal.symmetry(
                    unit_cell=unit_cell,
                    space_group_symbol=p['SYMB']))
    scatterers = []
    for k, a in six.iteritems(atoms):
      scatterers.append(xray.scatterer(label=a['LABEL'],
                                       scattering_type=a['NTYP'],
                                       site=(a['X'], a['Y'], a['Z']),
                                       b=a['B']))
    scatterers_flex=flex.xray_scatterer(scatterers)
    structure = xray.structure(special_position_settings=special_position_settings,
                               scatterers=scatterers_flex)
    return structure.as_py_code()

  def _parse(self):
    """Function to parse a FullProf command file into a ConfigObj
    """
    # XXX Todo: Implement a lot more (TOF, special cases)
    self.cfg['TITLE'] = self.lines[0][5:].strip()
    if self.enforce_valid_comments:
      for line in self.lines:
        if line.startswith("!Job "): break
        self.LO += 1
    # alias private class functions for conveniance
    get_next_line = self.__get_next_line
    parse_line = self.__parse_line
    check_comment = self.__check_comment
    parse_line(self.cfg,
      'JOBTYP, NPROF, NPHASE, NBCKGD, NEXCRG, NSCAT, NORI, IDUM, IWGT,\
      ILOR, IASG, IRESO, ISTEP, NRELL, ICRYG, IXUNIT, ICORR', "i")
    if self.cfg['IRESO'] != 0:
      raise ValueError("IRESO<>0 not supported (yet)! IRESO={}".format(self.cfg['IRESO']))
    check_comment("!")
    check_comment("!Ipr ")
    parse_line(self.cfg,
      'IOT, IPL, IPC, MAT, NXT, LST1, LST2, LST3, IPL1, IPL2, INSTRM,\
      JCIL, JSY, JLKH, JFOU, ISHOR, IANALY', "i")
    check_comment("!")
    check_comment("! Lambda1 ")
    parse_line(self.cfg,
      'LAMDA1, LAMDA2, RATIO, BKPOS, WDT, CTHM, TMR, RLIM, K', "f")
    check_comment("!")
    check_comment("!NCY ")
    parse_line(self.cfg,
      'MCYCLE, EPS, RELAX1, RELAX2, RELAX3, RELAX4, THMIN, STEP, THMAX,\
      ALPSD, SENT0', "if")
    if self.cfg['ICRYG'] == 0:
      NBCKGD = self.cfg['NBCKGD']
      if NBCKGD > 2 or NBCKGD < -3:
        check_comment("!")
        check_comment("! 2Theta")
        # extract background from file
        background = []
        for i in range(abs(NBCKGD)):
          linevals = get_next_line()
          background += ( float(linevals[0]), float(linevals[1]) )
        self.cfg['background'] = background
      #check_comment("!")
      #check_comment("! Excluded ")
      # XXX Todo: implement
    check_comment("!")
    check_comment("!")
    self.cfg['MAXS'] = int(get_next_line()[0])
    if self.cfg['ICRYG'] == 0:
      check_comment("!")
      check_comment("!  Zero ")
      parse_line(self.cfg,
        'ZER, FLGZER, SYCOS, FLCOS, SYSIN, FLSYN, LAMBDA, FLAMBDA, IGLMORE', "f"*8+"i")
      if self.cfg['IGLMORE'] == 1:
        check_comment("!")
        check_comment("! Microabsorption")
        check_comment("!   P0 ")
        parse_line(self.cfg, 'CPO, Cp, CCp, Tau, CTau', "f")
      check_comment("!   Background")
      linevals = get_next_line()
      lineflags = get_next_line()
      backgroundpoly = {}
      i = 0
      for (v, f) in zip(linevals, lineflags):
        backgroundpoly[str(i)] = [ v, f ]
        i += 1
      self.cfg['backgroundpoly'] = backgroundpoly
    phases = {}
    for p in range(self.cfg['NPHASE']):
      check_comment("!-------")
      check_comment("!  Data for PHASE")
      check_comment("!-------")
      phase = {}
      phase['PHSNM'] = get_next_line()
      check_comment("!")
      check_comment("!Nat ")
      parse_line(phase, 'NAT, NDIST, NMAGC, PREF1, PREF2, PREF3,\
        JBT, IRF, ISYM, ISTR, IFURT, ATZ, NVK, NPRO, IMORE', "iiifffiiiiifiii")
      if phase['IMORE'] != 0:
        # XXX Todo: read more phase data
        pass
      check_comment("!")
      phase['SYMB'] = get_next_line(raw=True)[:20].strip()
      if phase['ISYM'] != 0:
        # XXX Todo: read more phase data
        pass
      check_comment("!Atom ")
      atoms = {}
      for a in range(phase['NAT']):
        atom = {}
        parse_line(atom,
          'LABEL, NTYP, X, Y, Z, B, N, IOPIN, IOPFIN, N_Type', "ssfffffiii")
        parse_line(atom,
          'CX, CY, CZ, CB, CN, CMx, CMy, CMz', "f")
        atoms[str(a)] = atom
      phase['atoms'] = atoms
      check_comment("!-------> Profile ")
      check_comment("!  Scale ")
      parse_line(phase, 'S, GAM1, Bov, STR1, STR2, STR3, IstrainModel', "ffffffi")
      parse_line(phase, 'CS, FLGAM1, CBov, CSTR1, CSTR2, CSTR3', "f")
      if phase['IRF'] == 4:
        # XXX Todo: read more phase data
        pass
      else:
        check_comment("!       U")
        parse_line(phase, 'U, V, W, X, Y, IG, SZ, IsizeModel', "fffffffi")
        parse_line(phase, 'CU, CV, CW, CX, CY  CIG  CSZ', "f")
      check_comment("!     a ")
      parse_line(phase, 'a, b, c, alpha, beta, gamma', "f")
      parse_line(phase, 'CA, CB, CC, CD, CE, CF', "f")
      check_comment("!  Pref1 ")
      parse_line(phase, 'G1, G2, Pas1, Pas2, Pas3, Pas4', "f")
      parse_line(phase, 'CG1, CG2, CPas1, CPas2, CPas3, CPas4', "f")
      if self.cfg['RATIO'] < 0:
        check_comment("!Additional U,V,W")
        parse_line(phase, 'U2, V2, W2', "f")
        parse_line(phase, 'CU2, CV2, CW2', "f")

      phases[str(p+1)] = phase
    self.cfg['phases'] = phases


if __name__ == '__main__':
  # just a little test for debugging
  testfile1="""\
COMM  CeO2 CSIRO Test 25-143
! Current global Chi2 (Bragg contrib.) =      11.67
! Files => DAT-file: ce1.dat,  PCR-file: ce1
!Job Npr Nph Nba Nex Nsc Nor Dum Iwg Ilo Ias Res Ste Nre Cry Uni Cor Opt Aut
   0   5   1   0   0   0   0   0   0   0   0   0   0   2   0   0   0   1   1
!
!Ipr Ppl Ioc Mat Pcr Ls1 Ls2 Ls3 NLI Prf Ins Rpa Sym Hkl Fou Sho Ana
   0   0   1   0   1   0   4   0   0   1   0   1   1   1   2   0   0
!
! Lambda1  Lambda2    Ratio    Bkpos    Wdt    Cthm     muR   AsyLim   Rpolarz  2nd-muR -> Patt# 1
 1.540560 1.544330 -0.50000   25.000 15.0000  0.9100  0.0000   30.00    0.0000  0.0000
!
!NCY  Eps  R_at  R_an  R_pr  R_gl     Thmin       Step       Thmax    PSD    Sent0
  6  0.10  1.00  1.00  1.00  1.00     25.0000   0.025000   143.0000   0.000   0.000
!
!
      20    !Number of refined parameters
!
!  Zero    Code    SyCos    Code   SySin    Code  Lambda     Code MORE ->Patt# 1
  0.02573   11.0  0.00000    0.0  0.00000    0.0 0.000000    0.00   1
!
! Microabsorption coefficients for Pattern#  1
!   P0    Cod_P0    Cp   Cod_Cp     Tau  Cod_Tau
  0.0000    0.00  1.2313  141.00  0.1340  151.00
!   Background coefficients/codes  for Pattern#  1  (Polynomial of 6th degree)
      24.538     -17.709       6.313      -0.604       0.000       0.000
       21.00       31.00       41.00       51.00        0.00        0.00
!-------------------------------------------------------------------------------
!  Data for PHASE number:   1  ==> Current R_Bragg for Pattern#  1:     0.00
!-------------------------------------------------------------------------------
  CeO2
!
!Nat Dis Ang Pr1 Pr2 Pr3 Jbt Irf Isy Str Furth       ATZ    Nvk Npr More
   2   0   0 1.0 1.0 1.0   0   0   0   0   0        688.480   0   5   0
!
F m 3 m                  <--Space group symbol
!Atom   Typ       X        Y        Z     Biso       Occ     In Fin N_t Spc /Codes
Ce     CE      0.00000  0.00000  0.00000  0.46033   0.02000   0   0   0    0  # color green conn O O 0.0 3.2
                  0.00     0.00     0.00   161.00      0.00
O      O       0.25000  0.25000  0.25000  0.92218   0.04000   0   0   0    0  # color blue
                  0.00     0.00     0.00   171.00      0.00
!-------> Profile Parameters for Pattern #  1
!  Scale        Shape1      Bov      Str1      Str2      Str3   Strain-Model
 0.10216E-02   0.36723   0.00000   0.00000   0.00000   0.00000       0
    61.00000    81.000     0.000     0.000     0.000     0.000
!       U         V          W           X          Y        GauSiz   LorSiz Size-Model
   0.007618  -0.008887   0.010210   0.003351   0.000000   0.000000   0.000000    0
    101.000     91.000    111.000    121.000      0.000      0.000      0.000
!     a          b         c        alpha      beta       gamma      #Cell Info
   5.412223   5.412223   5.412223  90.000000  90.000000  90.000000   #multiple box -1 1 -1 1 0 1
   71.00000   71.00000   71.00000    0.00000    0.00000    0.00000
!  Pref1    Pref2      Asy1     Asy2     Asy3     Asy4
  0.00000  0.00000  0.06233  0.00000  0.00000  0.00000
     0.00     0.00   131.00     0.00     0.00     0.00
!Additional U,V,W parameters for Lambda2
   0.007461  -0.005255   0.008197   <--  U2,V2,W2 for lambda(2)
 191.000000 201.000000 181.000000
! Limits for selected parameters :
  14      0.0010      2.0000      0.0000   0  Cp_mabs_pat1
  15      0.0010      2.0000      0.0000   0  Tau_mabs_pat1
!  2Th1/TOF1    2Th2/TOF2  Pattern # 1
      25.000     143.000       1
"""

  testfile2="""\
COMM unnamed
! Files => DAT-file: rnd_27A_15L_10M_02H_P_1_25_,  PCR-file: rnd_27A_15L_10M_02H_P_1_25_
!Job Npr Nph Nba Nex Nsc Nor Dum Iwg Ilo Ias Res Ste Nre Cry Uni Cor Opt Aut
   2   0   1   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   1
!
!Ipr Ppl Ioc Mat Pcr Ls1 Ls2 Ls3 NLI Prf Ins Rpa Sym Hkl Fou Sho Ana
   2   0   1   1   1   0   4   0   0   1  10   0   1   4   2   0   0
!
! Lambda1  Lambda2    Ratio    Bkpos    Wdt    Cthm     muR   AsyLim   Rpolarz  2nd-muR -> Patt# 1
 1.541800 1.541800  1.00000   60.000  5.0000  0.0000  0.0000   50.00    0.0000  0.0000
!
!NCY  Eps  R_at  R_an  R_pr  R_gl     Thmin       Step       Thmax    PSD    Sent0
  1  0.01  1.00  1.00  1.00  1.00     10.0000   0.100000   100.0000   0.000   0.000
!
!
       0    !Number of refined parameters
!
!  Zero    Code    SyCos    Code   SySin    Code  Lambda     Code MORE ->Patt# 1
  0.00000    0.0  0.00000    0.0  0.00000    0.0 0.000000    0.00   0
!   Background coefficients/codes  for Pattern#  1  (Polynomial of 6th degree)
       0.000       0.000       0.000       0.000       0.000       0.000
        0.00        0.00        0.00        0.00        0.00        0.00
!-------------------------------------------------------------------------------
!  Data for PHASE number:   1  ==> Current R_Bragg for Pattern#  1:     0.00
!-------------------------------------------------------------------------------
 Phase_1
!
!Nat Dis Ang Pr1 Pr2 Pr3 Jbt Irf Isy Str Furth       ATZ    Nvk Npr More
  27   0   0 0.0 0.0 1.0   0   0   0   0   0        924.976   0   0   0
!
P 1                      <--Space group symbol
!Atom   Typ       X        Y        Z     Biso       Occ     In Fin N_t Spc /Codes
Cs1    Cs      0.25575  0.74218  0.01280  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
Cs2    Cs      0.76337  0.84611  0.64963  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
C1     C       0.80799  0.01687  0.87219  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
C2     C       0.95999  0.89218  0.97716  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
C3     C       0.61730  0.40423  0.22206  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
C4     C       0.96991  0.50496  0.26698  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
C5     C       0.27434  0.34505  0.32297  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
C6     C       0.51741  0.89889  0.90591  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
C7     C       0.28570  0.06143  0.11813  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
C8    C       0.85987  0.32773  0.87587  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
C9    C       0.57367  0.37951  0.94405  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
C10    C       0.74819  0.35206  0.92477  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
C11    C       0.11661  0.17692  0.48152  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
C12    C       0.68782  0.53267  0.57410  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
C13    C       0.35920  0.22199  0.20957  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
C14    C       0.54158  0.59845  0.57855  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
C15    C       0.10057  0.16401  0.79612  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
Ti1   Ti      0.42526  0.93883  0.93721  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
Ti2   Ti      0.81513  0.52398  0.20025  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
Ti3   Ti      0.51151  0.35888  0.34305  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
Ti4   Ti      0.43735  0.05657  0.82934  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
Ti5   Ti      0.96559  0.21095  0.23956  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
Ti6   Ti      0.59096  0.29067  0.09393  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
Ti7   Ti      0.24350  0.09520  0.50214  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
Ti8   Ti      0.60861  0.39336  0.54098  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
Ti9   Ti      0.80180  0.14419  0.13294  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
Ti10  Ti      0.49678  0.67274  0.97981  0.39478   1.00000   0   0   0    0
                  0.00     0.00     0.00     0.00      0.00
!-------> Profile Parameters for Pattern #  1
!  Scale        Shape1      Bov      Str1      Str2      Str3   Strain-Model
 0.10000E-01   0.00000   0.03837   0.00000   0.00000   0.00000       0
     0.00000     0.000     0.000     0.000     0.000     0.000
!       U         V          W           X          Y        GauSiz   LorSiz Size-Model
   1.706310  -1.179299   0.405540   0.000000   0.000000   0.000000   0.000000    0
      0.000      0.000      0.000      0.000      0.000      0.000      0.000
!     a          b         c        alpha      beta       gamma      #Cell Info
   9.422933  12.249814  16.018988  83.000000 109.000000 129.000000
    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
!  Pref1    Pref2      Asy1     Asy2     Asy3     Asy4
  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
     0.00     0.00     0.00     0.00     0.00     0.00
!  2Th1/TOF1    2Th2/TOF2  Pattern # 1
      10.000     100.000       1
"""
  x = pcr_parser(testfile1)
  print(x.cfg)
  x = pcr_parser(testfile2)
  print(x.cfg)
  print(x.extract_structure())
