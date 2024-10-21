from __future__ import absolute_import, division, print_function
import time
from mmtbx.conformation_dependent_library.cdl_utils import round_to_int
from collections import OrderedDict

# validation_types
validation_types = OrderedDict([
  ('MOTIF05'  , 'Motif (05)'),
  ('MOTIF'    , 'Motif (10)'),
  ('MOTIF20'  , 'Motif (20)'),
  ('MOTIF...' , 'Motif (->)'),
  ('GENERAL'  , 'General'),
  ('CIS'      , 'cis-peptide'),
  ('OUTLIER'  , 'Outlier'),
  ('ERROR'    , 'Error'),
])

phi_psi_2_peaks = {
  (-63, -42, -63, -42) : ['AA.1', 12120], #1
  (-64, -23, -87, -2) : ['AD.1', 499],
  (-115, 128, -118, 126) : ['BB.1', 223],
  (-68, 159, -57, -40) : ['PA.1', 216],
  (-55, 133, 87, -7) : ['Pd.1', 149],
  (-67, 148, -62, 140) : ['PP.1', 136],
  (-90, -3, -67, 139) : ['DP.1', 102],
  (-140, 158, -66, 146) : ['BP.1', 78],
  (-86, -2, 75, 22) : ['Dd.1', 74],
  (-118, 125, -66, 133) : ['BP.2', 72], #10
  (-92, -2, -58, -43) : ['DA.1', 67],
  (-63, 145, -136, 155) : ['PB.1', 67],
  (-79, 135, -133, 160) : ['PB.2', 64],
  (-66, 144, -118, 128) : ['PB.3', 63],
  (53, 37, 81, 2) : ['dd.1', 60],
  (-147, 159, -147, 156) : ['BB.2', 55],
  (-103, 14, -59, -34) : ['DA.2', 53],
  (-90, 123, -67, 138) : ['BP.3', 51],
  (-136, 161, -57, -41) : ['BA.1', 49],
  (91, -9, -70, 145) : ['dP.1', 48], # 20
  (-95, 107, -68, -20) : ['BA.2', 40],
  (-129, 124, 51, 41) : ['Bd.1', 36],
  (-133, 72, -65, -21) : ['ZA.1', 35],
  (54, -131, -70, -10) : ['pD.1', 34],
  (-80, -22, -159, 165) : ['DB.1', 34],
  (-68, 133, -117, -14) : ['PD.1', 34],
  (59, -132, -93, 8) : ['pD.2', 33],
  (-160, 166, -62, -37) : ['BA.3', 33],
  (-115, 101, -59, -36) : ['BA.4', 32],
  (-89, -29, -147, 151) : ['DB.2', 30], #30
  (-70, -39, -129, 80) : ['AZ.1', 28],
  (-115, -16, -155, 147) : ['DB.3', 28],
  (80, 4, -114, 140) : ['dB.1', 26],
  (-77, -43, -161, 158) : ['DB.4', 26],
  (-75, 124, -94, -35) : ['PD.2', 26],
  (77, 8, -120, 154) : ['dB.2', 25],
  (53, 57, -66, 153) : ['dP.2', 23],
  (58, 26, 57, 37) : ['dd.2', 22],
  (-124, 84, -66, 148) : ['ZP.1', 22],
  (61, 33, -61, -34) : ['dA.1', 21], #40
  (72, 22, -105, 151) : ['dB.3', 21],
  (-61, 138, -105, 6) : ['PD.3', 20],
  (-67, 159, 61, 35) : ['Pd.2', 19],
  (-92, -2, -132, 161) : ['DB.5', 19],
  (-90, 7, -116, 147) : ['DB.6', 19],
  (62, 35, -118, 162) : ['dB.4', 18],
  (48, 59, -66, -23) : ['dA.2', 18],
  (-77, -16, -130, 69) : ['AZ.2', 18],
  (85, -3, -64, -36) : ['dA.3', 17],
  (-89, 80, -76, -15) : ['gD.1', 16], #50
  (-83, -8, -131, 54) : ['AZ.3', 16],
  (-141, 100, -62, -28) : ['BA.5', 16],
  (-63, 134, -84, 81) : ['Pg.1', 15],
  (-132, 61, -65, 146) : ['ZP.2', 15],
  (-154, 153, -112, 10) : ['BD.1', 14],
  (-64, -40, -78, 125) : ['DP.2', 14],
  (-140, 131, -68, -17) : ['BA.6', 14],
  (80, 15, -94, -9) : ['dD.1', 13],
  (-89, -5, -123, 28) : ['DD.1', 13],
  (-78, -7, -120, 78) : ['AZ.4', 13], #60
  (-130, 5, 73, 15) : ['Dd.2', 12],
  (54, 42, -91, -3) : ['dD.2', 12],
  (-85, 71, -66, 144) : ['gP.1', 12],
  (-71, 118, -113, 16) : ['PD.4', 12],
  (-65, -41, -78, 153) : ['DP.3', 12],
  (-64, -32, 54, 57) : ['Dd.3', 12],
  (-62, -42, -115, 131) : ['DB.7', 12],
  (-126, 133, -106, 6) : ['BD.2', 12],
  (-134, 113, 59, -127) : ['Bp.1', 11],
  (-127, 99, 59, -124) : ['Bp.2', 11], #70
  (-101, 3, -129, 38) : ['DD.2', 11],
  (62, 25, -84, -8) : ['dD.3', 11],
  (-114, 13, -81, -16) : ['DD.3', 11],
  (61, 37, -104, 2) : ['dD.4', 10],
  (-86, 143, 50, 55) : ['Pd.3', 10],
  (-71, -15, -79, 78) : ['Dg.1', 10],
  (-153, 173, 55, 48) : ['Bd.2', 10],
  (-127, 2, -133, 128) : ['DB.8', 10],
  (-122, 128, -85, 80) : ['Bg.1', 10],
  (87, 175, -74, -16) : ['ED.1', 9], #80
  (97, -12, -77, -12) : ['dD.5', 9],
  (75, -173, -63, 143) : ['pP.1', 9],
  (-83, -12, -124, 105) : ['DB.9', 9],
  (-71, -14, 89, 175) : ['DE.1', 9],
  (-62, 137, -140, 16) : ['PD.5', 9],
  (-154, 168, -97, 1) : ['BD.3', 9],
  (87, 8, -149, 162) : ['dB.5', 8],
  (64, 31, -126, 24) : ['dD.6', 8],
  (168, -152, -135, 151) : ['EB.1', 8],
  (-99, 7, -80, 84) : ['Dg.2', 8], #90
  (-84, 69, -163, 150) : ['gB.1', 8],
  (-83, 68, -84, 173) : ['gP.2', 8],
  (-69, -35, 64, -153) : ['Ap.1', 8],
  (-64, 144, -138, 93) : ['PB.4', 8],
  (-129, 159, -103, -7) : ['BD.4', 8],
  (-121, -9, -65, -18) : ['DA.3', 8],
  (-114, 20, -116, 155) : ['DB.10', 8],
  (-103, 136, 60, 36) : ['Bd.3', 8],
  (-102, 8, -107, -178) : ['DP.4', 8],
  (-102, 114, -84, 78) : ['Bg.2', 8], #100
  (-102, -24, -116, -59) : ['DD.4', 8],
  }

class phi_psi_2_mask_class(dict):
  def __init__(self):
    self.data = {
      180  : "BBBBBBBBBBBBBBBPPPPPPPPPP                      pppppppppppEEEEEEEEEEBBBBB",
      175  : "BBBBBBBBBBBBBBBPPPPPPPPPPP                     pppppppppppEEEEEEEEEEEBBBB",
      170  : "BBBBBBBBBBBBBBBPPPPPPPPPPPP                      pppppppppEEEEEEEEEEEEBBB",
      165  : "BBBBBBBBBBBBBBBBPPPPPPPPPPP                      pppppppppEEEEEEEEEEEEEBB",
      160  : "BBBBBBBBBBBBBBBBBPPPPPPPPPPP                     EEEEEEEEEEEEEEEEEEEEEEBB",
      155  : "BBBBBBBBBBBBBBBBBPPPPPPPPPPP                    EEEEEEEEEEEEEEEEEEEEEEE B",
      150  : "BBBBBBBBBBBBBBBBBPPPPPPPPPPP                    EEEEEEEEEEEEEEEEEEEEEEEEB",
      145  : "BBBBBBBBBBBBBBBBBPPPPPPPPPPPP                     EEEEEEEEEEEE EEEEEEEEEB",
      140  : "BBBBBBBBBBBBBBBBBBPPPPPPPPPPPPP                  EEEEEEEEEEEEEEEEEE  EE B",
      135  : "BBBBBBBBBBBBBBBBBBBPPPPPPPPPPPP                  EEEEEEEEEEEEEEEE     EEB",
      130  : " BBBBBBBBBBBBBBBBBBPPPPPPPPPPPP                    EEEEEEEEEE EE  EEE EE ",
      125  : " BBBBBBBBBBBBBBBBBBPPPPPPPPPPPP                      EEEEEEE  EE  EEE    ",
      120  : " BBBBBBBBBBBBBBBBBBBPPPPPPPPPPP                      EEEEEEE  EE EE      ",
      115  : " BBBBBBBBBBBBBBBBBBBPPPPPPPPPP                   EE    EEEEEEE EEEE      ",
      110  : " BBBBBBBBBBBBBBBBBBBPPPPPPPPPP                   EEEEEEEEE EEE EE        ",
      105  : "BBBBBBBBBBBBBBBBBBBBPPPPP                          EEEEE   EEE          B",
      100  : "BBBBBBBBBBBBBBBBBBBBPPPP                   UU         EE    EEEE        B",
      95   : " BBBBBBBBBBZBBBBBBBBPPPP                   UU    dd   EE    EEEE         ",
      90   : " ZZZZZZZZZZZZZBBBBBBggg                         dddddddd UU EE UU     UU ",
      85   : "  ZZZZZZZZZZZZBBBgggggg                        dddddddd  UU    UU     UU ",
      80   : "  ZZZZZZZZZZZZggggggggg                     ddddddd                      ",
      75   : "  ZZZZZZZZZZZZggggggggg                     dddddddddd dd UU             ",
      70   : "  ZZZZZZZZZZZZggggggggg                   dddddddddddd dd UU             ",
      65   : "  ZZZZZZZZZZZZgggggggg                    ddddddddddd ddd                ",
      60   : "    ZZZZZZZZZZgggggggg                    ddddddddddd dd                 ",
      55   : "   ZZZZZZZZZZZgggggggggg                   dddddddddd    dd              ",
      50   : "  ZZZZZZZZZZZZgggggggggg                   ddddddddddd   dd dd           ",
      45   : "  DDDDDDDDZZZZgggggggg                     dddddddddddddd dddd           ",
      40   : "    DDDDDDDDDDDDDDgggg                     dddddddddddddd ddddd          ",
      35   : "    DDDDDDDDDDDDDDDDDg                      ddddddddddddddddddd          ",
      30   : "   DDDDDDDDDDDDDDDDDDD                      dddddddddddddddddd           ",
      25   : "   DDDDDDDDDDDDDDDDDD                        dddddddddddddddddd          ",
      20   : "    DDDDDDDDDDDDDDDDDD                        ddddddddddddddddd  UU      ",
      15   : "   DDDDDDDDDDDDDDDDDDDDD                      dddddddddddddddd   UU      ",
      10   : "  UUUUDDDDDDDDDDDDDDDDDD                       ddddddddddddddddd         ",
      5    : "  UU  DDDDDDDDDDDDDDDDDDD                       dddddddddddddddd         ",
      0    : "      DDDDDDDDDDDDDDDDDDD                        dddddddddddddddd        ",
      -5   : "      DDDDDDDDDDDDDDDDDDD                        dddddddddddddddd        ",
      -10  : "     DDDDDDDDDDDDDDDDDDDDD                      dddddddddddddddddd       ",
      -15  : "     DDDDDDDDDDDDDDDDDAADDD                     ddddddddddddddddddd      ",
      -20  : "      DDDDDDDDDDDDDDDAAAAAA                      dddddddddddddddddd      ",
      -25  : "   DDDDDDDDDDDDDDDDDAAAAAAAA                     dd ddddddddddddd        ",
      -30  : "  DDDDDDDDDDDDDDDDDDAAAAAAAA                      dddddddddddd ddddd     ",
      -35  : "  DD DDDDDDDDDDDDDDDAAAAAAAAA                    GGG ddddddddd ddddd     ",
      -40  : "    DD DDDDDDDDDDDDDAAAAAAAAA                    GGGGGG dddd    dddd     ",
      -45  : "    DD DDDDDDDDDDDDDAAAAAAAAAA                   GGGGGG   dd     dd      ",
      -50  : "    DDDDDDDDDDDDDDDDAAAAAAAAAAD                  GGGGGG   dd             ",
      -55  : "  DDDDDDDDDDDDDDDDDDAAAAAAAAAAD                  GGGGGG       UU         ",
      -60  : "  DDDDDDDDDDDDDDDDDDAAAAAAAAAADD                 GGGGGG   GG  UUU        ",
      -65  : "   DDDDDDDDDDDDDDDDDAAAAAAAAA DD                 GGGGGGG  GG   UU        ",
      -70  : "  DDDDDDDDDDDDDDDDDDDDDD DDDD                    GGGGGGG                 ",
      -75  : "  DDDDDDDDDDDDDDDDDDDDD  DDDD                     GG GGGG     UU         ",
      -80  : "    DDD DDDDDDDDDDDDD                               GGGGG     UU         ",
      -85  : "    DD     DDDDDDDDDD              UU               GGG   GG             ",
      -90  : "        UU  DDDDDD DD              UU           pp   GGGGGGG             ",
      -95  : "     UU UUUUUUUUUU                              pp pp  GGGG              ",
      -100 : "     UU UUUUUUUUUU                                ppp   GG               ",
      -105 : "      EEEEEEEEEEE            UU                pppppppp GG     UUUU      ",
      -110 : "UUU  EEEEEEEEEEEEE           UU            pp pppppppppp       UUUU    UU",
      -115 : "UUU  EEEEEEEEEEEEEE                        ppppppppp ppppp             UU",
      -120 : " UU  EEEEEEEEEEEEEEE                      pppppppppppppppp  EE     EE    ",
      -125 : "       EEEEEEEEEEEEE                      ppppppppppppppp  EEEEE   EEE   ",
      -130 : "      EEEEEEEEEEEEEE                       pppppppppppppp  EEEEE   EEEEE ",
      -135 : "EEEE EEEEEEEEEEEEEE                         ppppppppppppppEEEEE   EEEEEEE",
      -140 : "EEEE EEEEEEEEEEEEEEEE                        pppppppppppppEEEEEEE EEEEEEE",
      -145 : "EEEEEEEEEEEEEEEEEEEEE                         pppppppppppppEEEEEEEEEEEEEE",
      -150 : "EEEEEEEEEEEEEEEEEEEEE                         pppppppppppppEEEEEEEEEEEEEE",
      -155 : "EEEEEEEEEEEEEEEEEEEEEEE                       pppppppppppppEEEEEEEEEEEEEE",
      -160 : "EEEEEEEEEEEEEEEEEEEEEEE                     UU pppppppppppppEEEEEEEEEEEEE",
      -165 : "BBBBBBBBEEEEEEEEPPPPPPPP                    UU pppppppppppppEEEEEEEEEEEBB",
      -170 : "BBBBBBBBBBEEEEEEPPPPPPPP                       pppppppppppppEEEEEEEEEEEBB",
      -175 : "BBBBBBBBBBBBEEEEPPPPPPPPP                      pppppppppppEEEEEEEEEEEEBBB",
      -180 : "BBBBBBBBBBBBBBBBPPPPPPPPP                      pppppppppppEEEEEEEEEEBBBBB",
    #         0123456789012345678901234567890123456789012345678901234567890123456789012
    #                   1         2         3         4         5         6         7
    }

  def __getitem__(self, key):
    if key[1] in self.data:
      ptr = (key[0]+180)//5
      return self.data[key[1]][ptr:ptr+1]
    else:
      assert 0, 'key %(key)s not found' % locals()

  def get(self, key, default):
    if key[1] in self.data.keys():
      return self.__getitem__(key)
    else:
      return default

  def get_closest(self, key):
    new = (round_to_int(key[0], 5), round_to_int(key[1], 5))
    return self.get(new, None)

phi_psi_2_mask = phi_psi_2_mask_class()

class probability_class(dict):
  def __repr__(self):
    outl = ''
    for region, item in self.items():
      outl += '%s\n' % region#, sum(item.values())
      total = 0
      for pp2, p in item.items():
        total+=p
        outl += "  %s : %5.1f %5.1f\n" % (pp2,p,total)
    return outl

  def get_relative_probability(self, motif_key):
    key = motif_key.split()[0]
    prob = self.get(key[0], None)
    if prob is None:
      return None
    return prob.get(key, None)

class probability_distribution(dict):
  def __repr__(self):
    outl = ''
    keys = self.get_ordered_keys()
    for motif in keys:
      outl += '  %-5s : %5.2f\n' % (motif, self[motif]*100)
    return outl

  def get_ordered_keys(self):
    sort_orders = sorted(self.items(), key=lambda x: x[1], reverse=True)
    tmp = []
    for item in sort_orders: tmp.append(item[0])
    return tmp

def relative_probabilities(ignore_alpha_beta=False,
                           reverse_relationships=False,
                           ):
  data = probability_class()
  for key, item in phi_psi_2_peaks.items():
    if ignore_alpha_beta and (item[0].find('AA')>-1 or item[0].find('BB')>-1):
      continue
    if reverse_relationships:
      data.setdefault(item[0][1], {})
      data[item[0][1]][item[0]] = item[1]
    else:
      data.setdefault(item[0][0], {})
      data[item[0][0]][item[0]] = item[1]
  for region, item in data.items():
    total = sum(item.values())
    for pp2, n in item.items():
      item[pp2] = n/total*100
  if reverse_relationships: return data
  tmp = {}
  for region, item in data.items():
    for pp2, n in item.items():
      key = pp2.split('.')[0]
      tmp.setdefault(key, 0)
      tmp[key]+=n
  for key, n in tmp.items():
    data[key[0]][key]=n
  return data

def total_non_secondary_structure_probabilities(starting_selection='D'):
  data = probability_distribution()
  for key, item in phi_psi_2_peaks.items():
    motif, n = item
    if motif.find('AA')>-1 or motif.find('BB')>-1: continue
    if starting_selection and motif.find(starting_selection)!=0: continue
    key = motif.split('.')[0]
    data.setdefault(key, 0)
    data[key]+=n
  total = sum(data.values())
  for motif in data:
    data[motif] /= total
  return data

def get_closest_peak_two_D_space(xyz, general=None):
  #
  # need to cater to angle wrap
  #
  smallest_d2 = [1e9, 1e9]
  closest_peak = None
  for k, key in enumerate(phi_psi_2_peaks):
    if general and phi_psi_2_peaks[key][0].find(general)!=0: continue
    smallest_d2_local = [1e9, 1e9]
    for j in range(2):
      d2 = 0
      for i in range(2):
        d2 += (xyz[j*2+i]-key[j*2+i])**2
      if d2<smallest_d2_local[j]:
        smallest_d2_local[j]=d2
    # print '2D',k,general,key,smallest_d2,closest_peak,smallest_d2_local
    if max(smallest_d2_local)<max(smallest_d2):
      smallest_d2 = smallest_d2_local
      closest_peak = key
  # print 'returning',closest_peak,smallest_d2
  return closest_peak,smallest_d2

def get_two_D_space_peaks(xyz, closest_peak):
  rc = []
  for ptr in range(2):
    d2 = 0
    for i in range(2):
      d2 += (xyz[ptr*2+i]-closest_peak[ptr*2+i])**2
    rc.append(d2)
  return None, rc

def get_closest_peak(xyz, general=None):
  #
  # uses four dimensional space
  # need to cater to angle wrap
  #
  smallest_d2 = 1e9
  closest_peak = None
  for k, key in enumerate(phi_psi_2_peaks):
    if general and phi_psi_2_peaks[key][0].find(general)!=0: continue
    d2 = 0
    for i in range(4):
      d2 += (xyz[i]-key[i])**2
      if d2 > smallest_d2: break
    if d2<smallest_d2:
      smallest_d2=d2
      closest_peak=key
  return closest_peak,smallest_d2

def get_phi_psi_2_general(phi1, psi1, phi2, psi2,
                          informative=True,
                          verbose=False):
  key = (phi1, psi1, phi2, psi2)
  general = ''
  for i in range(0, 3, 2):
    phi_psi = [round_to_int(key[i], 5), round_to_int(key[i+1], 5)]
    lookup = phi_psi_2_mask.get(tuple(phi_psi), None)
    if lookup:
      general += lookup
    else:
      if informative:
        general += '!'
      else:
        return None
  if informative: general = general.replace(' ','!')
  return general

def get_phi_psi_2_vector(phi1,
                         psi1,
                         phi2,
                         psi2,
                         ):
  general = get_phi_psi_2_general(phi1, psi1, phi2, psi2, verbose=verbose)
  print(general)
  assert 0

def get_phi_psi_2_motif(phi1,
                        psi1,
                        phi2,
                        psi2,
                        return_type=False,
                        include_five_ring=True,
                        details=False,
                        verbose=False):
  if details: assert return_type
  key = (phi1, psi1, phi2, psi2)
  general = get_phi_psi_2_general(phi1, psi1, phi2, psi2, verbose=verbose)
  assert general.find(' ')==-1, '%s "%s"' % (key, general)
  if general is None or general.find('!')>-1:
    smallest_d2=1e9
    rc = '%-5s ' % general
    if return_type:
      rc = [rc, 'OUTLIER'] #'Improbable'
    if details:
      rc.append([general, '!'])
    return rc
  else:
    closest_peak, smallest_d2 = get_closest_peak(key, general=general)
    # print key, closest_peak, phi_psi_2_peaks[closest_peak], smallest_d2
    closest_peak_pair, smallest_d2_pair = get_closest_peak_two_D_space(
      key,
      general=general)
    if closest_peak:
      closest_peak_peak, smallest_d2_peak = get_two_D_space_peaks(key,
                                                                  closest_peak)
      if max(smallest_d2_pair)<400 and max(smallest_d2_peak)<400 and verbose:
        if smallest_d2_pair!=smallest_d2_peak:
          print('  DEBUG %s != %s' % (smallest_d2_pair, smallest_d2_peak))
          if (max(smallest_d2_pair)<100) != (max(smallest_d2_peak)<100):
            print('  one is closer than the other')
      addtional_args = [smallest_d2] + smallest_d2_peak
      addtional_args += [closest_peak, closest_peak_peak]
    else:
      addtional_args = [smallest_d2] + smallest_d2_pair
      addtional_args += [closest_peak, closest_peak_pair]
  def _verbose_output(args):
    return '%-5s%2s %6.1f %6.1f %6.1f %s %s' % tuple(args)
  def _output(args):
    return '%-5s%1s' % tuple(args)
  def _return_return(args, label, addtional_args, verbose, return_type, details):
    if verbose: outl = _verbose_output(args+addtional_args)
    else: outl = _output(args)
    rc = outl
    if return_type:
      rc = [outl, label]
    if details:
      rc.append(args+addtional_args)
    return rc
  if include_five_ring and max(smallest_d2_pair)<=25:
    return _return_return([phi_psi_2_peaks[closest_peak][0], '.'],
                          'MOTIF05',
                          addtional_args,
                          verbose,
                          return_type,
                          details)
  elif max(smallest_d2_pair)<=100:
    return _return_return([phi_psi_2_peaks[closest_peak][0], ''],
                          'MOTIF',
                          addtional_args,
                          verbose,
                          return_type,
                          details)
  elif max(smallest_d2_pair)<=400:
    return _return_return([phi_psi_2_peaks[closest_peak][0], '*'],
                          'MOTIF20',
                          addtional_args,
                          verbose,
                          return_type,
                          details)
  elif min(smallest_d2_pair)>1e8:
    return _return_return([general, '?'],
                          'GENERAL',
                          addtional_args,
                          verbose,
                          return_type,
                          details)
  else:
    return _return_return([general, ''],
                          'MOTIF...',
                          addtional_args,
                          verbose,
                          return_type,
                          details)
  print(key, closest_peak, smallest_d2)
  assert 0
  return None

def guess_next_phi_psi(four, pp2):
  pp_key = pp2[:4]
  pp_key = pp_key.split('.')[0].strip()
  phi_psi = None
  if pp_key[0]!='!' and pp_key[1]=='!':
    rp = relative_probabilities(ignore_alpha_beta=True)
    current = rp.get(pp_key[0], None)
    if current is None: return None
    for key in sorted(current, key=current.__getitem__, reverse=True):
      item = current[key]
      if key.find('.')>-1: break
    for phi_psi_2, info in phi_psi_2_peaks.items():
      if info[0]==key: break
    phi_psi = phi_psi_2[2:]
  return phi_psi

def get_phi_psi_2_dihedrals_from_basins(basin):
  for phi_psi_2, info in phi_psi_2_peaks.items():
    if info[0].find(basin)==0:
      return phi_psi_2

def get_phi_psi_2_key(pp2_key,
                      highest_probability=False,
                      closest=False,
                      reverse_relationships=False,
                      ):
  def _d2(xy1, xy2):
    d2 = 0
    for i in range(2):
      d2 += (xy1[i]-xy2[i])**2
    return d2
  def d2(xy1, xy2):
    rc = 1e9
    for phi in range(-360,361,360):
      for psi in range(-360,361,360):
        rc = min(rc, _d2([xy1[0]+phi, xy1[1]+psi], xy2))
    return rc
  #
  rc = relative_probabilities(ignore_alpha_beta=False,
                              reverse_relationships=reverse_relationships,
                              )
  full_key = None
  for key, item in rc.items():
    if reverse_relationships:
      if pp2_key[1]!=key: continue
    else:
      if pp2_key[0]!=key: continue
    if closest:
      min_dist = 1e9
      min_key = None
      for full_key, prob in item.items():
        dih = get_phi_psi_2_dihedrals_from_basins(full_key)
        if reverse_relationships:
          td2 = d2(dih[:2], closest)
        else:
          td2 = d2(dih[2:], closest)
        if td2<min_dist:
          min_dist=td2
          min_key=full_key
      full_key=min_key
    elif highest_probability:
      mv = max(item.values())
      for full_key, prob in item.items():
        if prob==mv: break
  # if pp2_key.find('!')==-1 and pp2_key.find('?')==-1:
  #   if pp2_key.find(full_key.split('.')[0])==-1:
  #     print('%s != %s' % (pp2_key, full_key))
  return full_key

def get_phi_psi_2_next_value(pp2_key,
                             highest_probability=False,
                             closest=False,
                             ):
  key = get_phi_psi_2_key(pp2_key,
                          highest_probability=highest_probability,
                          closest=closest,
                          )
  for pp2, data in phi_psi_2_peaks.items():
    if data[0].find(i)>-1:
      break

def print_lookup_table(data):
  ij = 0
  rr = [i for i in range(-179,180,2)]
  outl = ''
  for i in rr:
    line=[]
    for j in rr:
      line.append(str(data[ij]))
      ij+=1
    outl += '%s\n' % (','.join(line))
  f=file('test.csv', 'wb')
  f.write(outl)
  del f

def load_phi_psi_2_rama_restraints_tables():
  t0=time.time()
  from scitbx.array_family import flex
  from mmtbx_ramachandran_restraints_ext import lookup_table
  tables = {}
  tmp = OrderedDict()
  rr = [i for i in range(-179,180,2)]
  for i in rr:
    for j in rr:
      tmp[(i,j)]=0
  for pp2, info in phi_psi_2_peaks.items():
    for i in range(2):
      a=pp2[i*2]
      b=pp2[(i*2)+1]
      data = flex.double()
      for k, v in zip(tmp.keys(), tmp.values()):
        val = 1e9
        for phi in range(-360,361,360):
          for psi in range(-360,361,360):
            val = min(val, ((k[0]+phi)-a)**2 + ((k[1]+psi)-b)**2)
        data.append(val/1000)
      # print_lookup_table(data)
      t = lookup_table(data, 180)
      key = '%s|%s' % (info[0], ['leading', 'following'][i])
      tables[key] = t
  print('time to load tables: %0.1fs' % (time.time()-t0))
  return tables

def get_rama_table(proxy,
                   phi_psi_2_tables,
                   ):
  key, strategy = proxy.residue_type.split('...')
  reverse_relationships=False
  if proxy.residue_type[0]=='!':
    reverse_relationships=True
  closest=False
  if strategy.find('closest')>-1:
    tmp,phi,psi = strategy.split('_')
    closest=(float(phi),float(psi))
  pp2_key = get_phi_psi_2_key(
    key,
    highest_probability=(strategy=='highest_probability'),
    closest=closest,
    reverse_relationships=reverse_relationships,
    )
  if pp2_key is None: return None
  if pp2_key.find('.')==-1: pp2_key += '.1'
  if pp2_key.find('|')==-1: pp2_key += '|%s' % key.split('|')[1]
  rama_table = phi_psi_2_tables[pp2_key]
  return rama_table

def get_phi_psi_key_for_rama_proxy(phi_psi_2_motifs, three, strategy='closest'):
  pp2_key_1 = phi_psi_2_motifs.get(three[1].resseq_as_int(), '')
  pp2_key_2 = phi_psi_2_motifs.get(three[2].resseq_as_int(), '')
  pp2_t1 = '  '
  pp2_t2 = '  '
  if pp2_key_1.strip():
    pp2_t1 = pp2_key_1.split()[0]
  if pp2_key_2.strip():
    pp2_t2 = pp2_key_2.split()[0]
  # MISMATCH
  # if pp2_t1[1]!=pp2_t2[0]: return None
  # TOO FAR FROM MOTIF
  # if pp2_t1.find('.')==-1 and pp2_t1.find('.')==-1: return None
  # LEADING IS QUESTIONABLE
  if pp2_key_1.find('?')>-1: return None
  # OUTLIERS
  if(pp2_key_1.strip()=='!!' and pp2_key_2.strip()=='!!'):
    return None
  elif (pp2_key_1 and pp2_key_1[0]!='!' and pp2_key_1[1]=='!'):
    # following outlier
    pp2_key = '%s|%s' % (pp2_key_1.strip(), 'following')
  elif (pp2_key_2 and pp2_key_2[0]=='!' and pp2_key_2[1]!='!'):
    # leading outlier
    pp2_key = '%s|%s' % (pp2_key_2.strip(), 'leading')
  elif not pp2_key_1.strip():
    pp2_key = '%s|%s' % (pp2_key_2.strip(), 'leading')
  else:
    pp2_key = '%s|%s' % (pp2_key_1.strip(), 'following')
  pp2_key+='...%s' % strategy
  return pp2_key

if __name__=='__main__':

  print(relative_probabilities(ignore_alpha_beta=False))
  print(relative_probabilities(ignore_alpha_beta=True))
  print(relative_probabilities(ignore_alpha_beta=False,
                               reverse_relationships=True))
  # print(total_non_secondary_structure_probabilities())

  print(phi_psi_2_peaks[ (168, -152, -135, 151) ])

  for i in range(2):
    print(('-%s' % (i))*39)
    for phi_psi_2 in [(0,0,0,0),
                      ( -65, -40, -61, -44),
                      ]:
      print(phi_psi_2, get_phi_psi_2_motif(*tuple(phi_psi_2), verbose=i))

  print('(%s,%s) -> "%s"' % (1,1,phi_psi_2_mask.get_closest((1,1))))
  assert 0

  j=0
  phi_psi_2 = (-63, -42, -63, -42)
  for i in range(4):
    for x in range(-10,11):
      j+=1
      tmp = list(phi_psi_2)
      tmp[i] += x
      rc = get_phi_psi_2_motif(*tuple(tmp))
      assert rc.find('*')==-1, '%sth test of %s was %s' % (j, tmp, rc)

  print('OK')
