"""
A57: chain-start placement builds NH3 only at standard amino acids and at other
residues with an outgoing peptide bond in every complete conformer. Standalone
residues keep their dictionary H under every n_terminal_charge, qmready
included; residues without a next residue, whose conformers disagree, or
with a conformer lacking C or the next N are ambiguous and keep them too; a
residue whose N is not bonded to CA in its dictionary (beta-amino acids) is
never built. Each such residue is listed. A supplied H1-H3 on N whose name
the dictionary gives another parent is listed and not renamed.
"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
from iotbx.cli_parser import run_program
from six.moves import cStringIO as StringIO
from mmtbx.programs import reduce2
from mmtbx.regression import tst_reduce2_free_form_h as ff
from mmtbx.regression import tst_reduce2_not_placed as tnp

pdb_amq = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  N   AMQ A 406      13.762  16.124  13.840  1.00 20.00           N
HETATM    2  CA  AMQ A 406      15.039  15.451  13.570  1.00 20.00           C
HETATM    3  C   AMQ A 406      15.085  15.022  12.126  1.00 20.00           C
HETATM    4  OT1 AMQ A 406      14.064  14.741  11.546  1.00 20.00           O
HETATM    5  OT2 AMQ A 406      16.262  14.954  11.484  1.00 20.00           O
HETATM    6  CB  AMQ A 406      15.168  14.223  14.472  1.00 20.00           C
HETATM    7  CG  AMQ A 406      15.121  14.652  15.915  1.00 20.00           C
HETATM    8  CD1 AMQ A 406      13.983  14.784  16.751  1.00 20.00           C
HETATM    9  CD2 AMQ A 406      16.156  14.993  16.717  1.00 20.00           C
HETATM   10  CE2 AMQ A 406      17.614  15.026  16.340  1.00 20.00           C
HETATM   11  OE1 AMQ A 406      12.690  14.550  16.410  1.00 20.00           O
HETATM   12  NE1 AMQ A 406      14.427  15.180  17.914  1.00 20.00           N
HETATM   13  OE2 AMQ A 406      15.629  15.301  17.911  1.00 20.00           O
END
"""

# CCD ideal heavy atoms; the dictionary names its NH3 H1N, H2N, H3N
pdb_01w = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  O52 01W A 406      14.881  17.661  19.813  1.00 20.00           O
HETATM    2  N5  01W A 406      15.169  18.842  19.879  1.00 20.00           N
HETATM    3  O5  01W A 406      14.288  19.683  19.908  1.00 20.00           O
HETATM    4  C5  01W A 406      16.589  19.256  19.930  1.00 20.00           C
HETATM    5  C6  01W A 406      17.588  18.300  19.903  1.00 20.00           C
HETATM    6  C7  01W A 406      18.914  18.683  19.950  1.00 20.00           C
HETATM    7  C4  01W A 406      16.915  20.599  19.999  1.00 20.00           C
HETATM    8  C3  01W A 406      18.238  20.988  20.046  1.00 20.00           C
HETATM    9  N3  01W A 406      18.586  22.425  20.120  1.00 20.00           N
HETATM   10  O32 01W A 406      17.900  23.184  20.782  1.00 20.00           O
HETATM   11  O3  01W A 406      19.557  22.849  19.519  1.00 20.00           O
HETATM   12  C2  01W A 406      19.244  20.030  20.026  1.00 20.00           C
HETATM   13  N13 01W A 406      20.584  20.421  20.074  1.00 20.00           N
HETATM   14  C1  01W A 406      21.645  19.422  19.915  1.00 20.00           C
HETATM   15  CB  01W A 406      23.009  20.108  20.007  1.00 20.00           C
HETATM   16  CA  01W A 406      24.116  19.065  19.842  1.00 20.00           C
HETATM   17  N   01W A 406      24.077  18.126  20.971  1.00 20.00           N
HETATM   18  C   01W A 406      25.456  19.755  19.810  1.00 20.00           C
HETATM   19  O   01W A 406      26.145  19.783  20.802  1.00 20.00           O
HETATM   20  OXT 01W A 406      25.883  20.339  18.679  1.00 20.00           O
END
"""

# beta-amino acid, the CCD's names: N on CA, C on CB (GeoStd until 2026-11-25:
# N on CB, C on CA)
pdb_b3a = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  CG  B3A A   1      20.110  19.607  18.983  1.00 20.00           C
HETATM    2  CA  B3A A   1      21.582  19.382  18.639  1.00 20.00           C
HETATM    3  N   B3A A   1      22.301  20.647  18.801  1.00 20.00           N
HETATM    4  CB  B3A A   1      21.770  18.882  17.205  1.00 20.00           C
HETATM    5  C   B3A A   1      21.075  17.560  16.981  1.00 20.00           C
HETATM    6  O   B3A A   1      21.346  16.536  17.594  1.00 20.00           O
HETATM    7  OXT B3A A   1      20.097  17.629  16.045  1.00 20.00           O
END
"""

# 1pq7 ILE A 16 (as ALA) and VAL A 17 (as GLY): a peptide bond, C-N 1.33 A
pdb_ala_gly = 'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1\n' + \
  tnp.ref_ala_gly + 'END\n'

def kept(o):
  return [l.strip() for l in o.chain_start_kept_text().splitlines()[1:]]

def h_on_n(o, resname, altloc=None):
  h = o.get_model().get_hierarchy()
  for ag in h.atom_groups():
    if ag.resname != resname or (altloc is not None and ag.altloc != altloc):
      continue
    n = ag.get_atom('N')
    if n is None: continue
    return sorted(a.name.strip() for a in h.atoms() if
      a.element.strip().upper() in ('H', 'D') and a.distance(n) < 1.3 and
      a.parent().altloc in ('', ag.altloc))
  return None

def qm_charge(pdb_str, prefix):
  fn = prefix + '.pdb'
  with open(fn, 'w') as f: f.write(pdb_str)
  r = run_program(program_class=reduce2.Program, logger=StringIO(),
    args=[fn, 'qmready=True', 'output.filename=%s_h.pdb' % prefix,
          'output.description_file_name=%s_h.txt' % prefix,
          'output.overwrite=True'])
  return r.total_charge, open('%s_h.txt' % prefix).read()

START = ('first_in_chain', 'all_chain_starts')

def test_sam_sah_standalone():
  # alone in their chain, at 406 and as residue 1: both dictionary H on N under
  # every setting, listed as standalone; qmready: the neutral NH2's charge
  for pdb_str, resname, charge in ((ff.pdb_sam, 'SAM', 0),
                                   (ff.pdb_sah, 'SAH', -1)):
    for resseq in ('406', '  1'):
      s = pdb_str.replace('%s A 406' % resname, '%s A %s' % (resname, resseq))
      for setting in START + ('residue_one',):
        o = ff.place(s, n_terminal_charge=setting)
        assert h_on_n(o, resname) == ['HN1', 'HN2'], (resname, resseq,
          setting, h_on_n(o, resname))
        line = 'N %s A %s: no peptide bond to a next residue (standalone); ' \
          'terminal conversion not applied; dictionary placement retained' % (
          resname, resseq.strip())
        if setting != 'residue_one' or resseq == '  1':
          assert kept(o) == [line], kept(o)
        else:
          assert kept(o) == [], kept(o)
    value, text = qm_charge(pdb_str, 'tst_reduce2_chain_start_ligands_' +
      resname)
    assert value == charge, (resname, value)
    assert 'N %s A 406: no peptide bond' % resname in text, text

def test_01w_amq():
  # 01W, whose old dictionary (tml.OLD_DICTIONARIES, supplied by ff.place())
  # is NH3 (H1N, H2N, H3N): all three; AMQ (no C=O,
  # undetermined in pdb_interpretation): H and HN2; no next residue makes 01W
  # (an amino acid) ambiguous, AMQ standalone
  for pdb_str, resname, names, why in (
      (pdb_01w, '01W', ['H1N', 'H2N', 'H3N'], 'ambiguous, no next residue'),
      (pdb_amq, 'AMQ', ['H', 'HN2'], 'no peptide bond to a next residue')):
    for setting in START:
      o = ff.place(pdb_str, n_terminal_charge=setting)
      assert h_on_n(o, resname) == names, (resname, setting, h_on_n(o, resname))
      assert [l for l in kept(o) if l.startswith(
        'N %s A 406: %s' % (resname, why))], kept(o)

def test_b3a():
  # beta-amino acid alone: its dictionary H stays, listed as any amino acid
  # without a next residue (GeoStd bonds N to CA since 2026-11-25, as the
  # CCD; before, N to CB: "N bonded to CB, not CA, in the dictionary")
  for setting in ('residue_one',) + START:
    o = ff.place(pdb_b3a, n_terminal_charge=setting)
    assert h_on_n(o, 'B3A') == ['H'], (setting, h_on_n(o, 'B3A'))
    assert kept(o) == ['N B3A A 1: ambiguous, no next residue (a free amino '
      'acid, or the next residue is missing); terminal conversion not '
      'applied; dictionary placement retained'], kept(o)

def test_conformers():
  # outgoing_peptide_bond() compares C and N within each complete conformer:
  # C(A)-N(B) 1.33 A with C(A)-N(A) and C(B)-N(B) far is no bond (the
  # reviewer's case); conformers that disagree are 'mixed'; a conformer
  # without C or N never counts as bonded: 'missing' where that residue has
  # main-chain atoms in the conformer, 'altloc' where the letter is not one of
  # its main-chain conformers
  from mmtbx.hydrogens.reduce_hydrogen import outgoing_peptide_bond
  def two(c_a, c_b, n_a, n_b, ca_b=None):
    s = 'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1\n'
    k = 0
    for name, resname, resseq, alt, x in (('C', 'ALA', 1, 'A', c_a),
        ('C', 'ALA', 1, 'B', c_b), ('N', 'GLY', 2, 'A', n_a),
        ('N', 'GLY', 2, 'B', n_b), ('CA', 'GLY', 2, 'B', ca_b)):
      if x is None: continue
      k += 1
      s += 'ATOM  %5d  %-3s%s%s A%4d    %8.3f%8.3f%8.3f  0.50 20.00' \
        '           %s\n' % (k, name, alt, resname, resseq, x, 10, 10, name[0])
    rgs = iotbx.pdb.input(lines=s.splitlines(),
      source_info=None).construct_hierarchy().only_chain().residue_groups()
    return outgoing_peptide_bond(rgs[0], rgs[1])
  assert two(10, 20, 15, 11.33) == 'none'       # C(A)-N(B) 1.33 only
  assert two(10, 20, 11.33, 21.33) == 'bonded'
  assert two(10, 20, 11.33, 25) == 'mixed'
  assert two(10, None, 11.33, None) == 'bonded'
  # GLY without conformer B, ALA's C in A and B: letters differ
  assert two(10, 20, 11.33, None) == 'altloc'
  # GLY B has CA but no N: N missing in that conformer
  assert two(10, 20, 11.33, None, ca_b=25) == 'missing'

def test_missing_next():
  # 0A0 alone as residue 1 (a free amino acid, or the next residue missing):
  # ambiguous, dictionary placement retained and listed; nothing added and
  # dropped
  for setting in ('residue_one',) + START:
    o = ff.place(tnp.pdb_0a0, n_terminal_charge=setting)
    assert h_on_n(o, '0A0') == ['HN1', 'HN2'], (setting, h_on_n(o, '0A0'))
    assert kept(o) == ['N 0A0 A 1: ambiguous, no next '
      'residue (a free amino acid, or the next residue is missing); '
      'terminal conversion not applied; dictionary placement retained'], \
      kept(o)
    assert tnp.records(o, 'not_parameterized') == [], tnp.records(o)

def test_supplied_h2_sam():
  # a supplied H2 1 A from SAM's N: the SAM dictionary bonds H2 to C2; listed,
  # its parent is not taken from the name and nothing is renamed
  n = [l for l in ff.pdb_sam.splitlines() if l[12:16] == ' N  '][0]
  s = tnp.with_atoms(ff.pdb_sam, 'SAM', '406', [('H2', 'H', float(n[30:38]),
    float(n[38:46]) + 1.0, float(n[46:54]))])
  o = ff.place(s, n_terminal_charge='first_in_chain', keep_existing_H=True)
  assert ('H2 SAM A 406: within bond length of N, but the SAM dictionary '
    'bonds H2 to C2; not renamed; name collision unresolved'
    in kept(o)), kept(o)
  # placement's own H2 (on C2, first at a bogus site) is never listed: SAM 1
  # with only N in A and B, without and with another supplied H (HN1 on N)
  s = '\n'.join(l if l[12:16] != ' N  ' else '\n'.join(
    l[:16] + alt + l[17:54] + '  0.50' + l[60:] for alt in 'AB')
    for l in ff.pdb_sam.replace('SAM A 406', 'SAM A   1').splitlines()) + '\n'
  s_hn1 = tnp.with_atoms(s, 'SAM', '1', [('HN1', 'H', float(n[30:38]),
    float(n[38:46]) - 1.0, float(n[46:54]))])
  for pdb_str, kw in ((s, {}), (s_hn1, {'keep_existing_H': True})):
    o = ff.place(pdb_str, n_terminal_charge='first_in_chain', **kw)
    assert [l for l in kept(o) if 'H2' in l] == [], (kw, kept(o))

def test_metal_option():
  # metal_bound_amines_as_nh2 builds no NH3 where the check keeps the
  # dictionary placement, under every setting, also at a start the setting
  # does not select (residue_one at 406, no_charge): B3A and SAM with a Zn
  # 2.1 A from N; the same report line as without the metal
  b3a, sam = pdb_b3a, ff.pdb_sam.replace('SAM A 406', 'SAM A   1')
  for pdb_str, resname, names in ((b3a, 'B3A', ['H']),
                                  (sam, 'SAM', ['HN1', 'HN2'])):
    for resseq, setting in (('1', 'residue_one'), ('1', 'first_in_chain'),
        ('406', 'residue_one'), ('1', 'no_charge')):
      s = pdb_str.replace('%s A   1' % resname, '%s A %3s' % (resname, resseq))
      n = [l for l in s.splitlines() if l[12:16] == ' N  '][0]
      zn = 'HETATM  999 ZN    ZN B   1    %8.3f%8.3f%8.3f  1.00 20.00' \
        '          ZN' % (float(n[30:38]) + 2.1, float(n[38:46]),
        float(n[46:54]))
      o = ff.place(s.replace('END', zn + '\nEND'), n_terminal_charge=setting,
        metal_bound_amines_as_nh2=True)
      assert h_on_n(o, resname) == names, (resname, resseq, setting,
        h_on_n(o, resname))
      assert [l for l in kept(o) if l.startswith('N %s A %s: ' % (resname,
        resseq)) and l.endswith('dictionary placement retained')], (resname,
        resseq, setting, kept(o))

def test_conformer_without_n():
  # 0A0-GLY in conformers A and B (B 0.1 A along x), GLY's N only in A: B
  # cannot decide, so no conformer is converted (ambiguous); both keep the
  # dictionary placement, HN1 and HN2
  lines = []
  for l in tnp.with_next_gly(tnp.pdb_0a0).splitlines():
    if not l.startswith(('ATOM', 'HETATM')):
      lines.append(l); continue
    lines.append(l[:16] + 'A' + l[17:54] + '  0.50' + l[60:])
    if l[17:20] == 'GLY' and l[12:16] == ' N  ': continue
    lines.append(l[:16] + 'B' + l[17:30] + '%8.3f' % (float(l[30:38]) + 0.1)
      + l[38:54] + '  0.50' + l[60:])
  s = '\n'.join(l for l in lines if l != 'END') + '\nEND\n'
  for setting in ('residue_one', 'first_in_chain', 'all_chain_starts'):
    o = ff.place(s, n_terminal_charge=setting)
    for alt in 'AB':
      assert h_on_n(o, '0A0', alt) == ['HN1', 'HN2'], (setting, alt,
        h_on_n(o, '0A0', alt))
    assert kept(o) == ['N 0A0 A 1: ambiguous, C or the next residue\'s N '
      'missing in a conformer; terminal conversion not applied; dictionary '
      'placement retained'], kept(o)

def test_conformer_letters():
  # 0A0-GLY in conformers A and B, and 0A0's CB also in C (a side-chain-only
  # letter, no backbone in C): conformer C cannot decide, by its letter, not a
  # missing atom; ambiguous, both A and B keep HN1, HN2
  lines = []
  for l in tnp.with_next_gly(tnp.pdb_0a0).splitlines():
    if not l.startswith(('ATOM', 'HETATM')):
      lines.append(l); continue
    lines.append(l[:16] + 'A' + l[17:54] + '  0.40' + l[60:])
    lines.append(l[:16] + 'B' + l[17:30] + '%8.3f' % (float(l[30:38]) + 0.1)
      + l[38:54] + '  0.40' + l[60:])
    if l[17:20] == '0A0' and l[12:16] == ' CB ':
      lines.append(l[:16] + 'C' + l[17:30] + '%8.3f' % (float(l[30:38]) +
        0.2) + l[38:54] + '  0.20' + l[60:])
  s = '\n'.join(l for l in lines if l != 'END') + '\nEND\n'
  o = ff.place(s, n_terminal_charge='first_in_chain')
  for alt in 'AB':
    assert h_on_n(o, '0A0', alt) == ['HN1', 'HN2'], (alt, h_on_n(o, '0A0',
      alt))
  assert kept(o) == ['N 0A0 A 1: ambiguous, altloc letters do not match '
    'between this and the next residue; terminal conversion not applied; '
    'dictionary placement retained'], kept(o)

def test_standard_and_terminus():
  # ALA at a chain start: NH3 as before, nothing listed; 0A0 peptide-bonded
  # to a GLY: a terminus with its dictionary's names (terminal_nh3: HN1,
  # HN2 and the computed HN3), nothing removed
  o = ff.place(pdb_ala_gly, n_terminal_charge='first_in_chain')
  assert h_on_n(o, 'ALA') == ['H1', 'H2', 'H3'], h_on_n(o, 'ALA')
  assert kept(o) == [], kept(o)
  # a free ALA (no next residue): a standard residue, NH3 as before
  ala = '\n'.join(l for l in pdb_ala_gly.splitlines() if 'GLY' not in l)
  o = ff.place(ala + '\n', n_terminal_charge='first_in_chain')
  assert h_on_n(o, 'ALA') == ['H1', 'H2', 'H3'], h_on_n(o, 'ALA')
  assert kept(o) == [], kept(o)
  o = ff.place(tnp.with_next_gly(tnp.pdb_0a0),
    n_terminal_charge='first_in_chain')
  assert kept(o) == [], kept(o)
  assert h_on_n(o, '0A0') == ['HN1', 'HN2', 'HN3'], h_on_n(o, '0A0')
  assert tnp.records(o) == [], tnp.records(o)

if __name__ == '__main__':
  test_sam_sah_standalone()
  test_01w_amq()
  test_b3a()
  test_conformers()
  test_missing_next()
  test_supplied_h2_sam()
  test_metal_option()
  test_conformer_without_n()
  test_conformer_letters()
  test_standard_and_terminus()
  print('OK')
