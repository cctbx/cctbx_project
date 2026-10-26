"""reduce2 placement at alcohols bound to Ti, Zr, Hf, V, Mo or W: an O-H whose
O has one heavy neighbour, a C in no double, triple, aromatic or delocalized
bond, bound to one of these metals is taken as an alkoxide, and reduce2 does
not build its H. Mo and V in a residue with Fe (FeMo-co, FeV-co) do not
count: there they are M(III) and the homocitrate alcohol keeps its H.
Models: 1DMR (DMSO reductase), Ser 147 and Mo 784 (OG 1.90 A), and the same
with Zn in place of Mo; 1QGU (nitrogenase MoFe protein), homocitrate C 501
and FeMo-co C 503 (O7-Mo1 2.35 A).
"""
pdb_1dmr = """\
CRYST1   80.880   80.880  229.700  90.00  90.00  90.00 P 41 21 2
ATOM   1112  N   SER A 147      -6.606  44.166  35.173  1.00  4.11           N
ATOM   1113  CA  SER A 147      -7.656  45.096  34.808  1.00  7.64           C
ATOM   1114  C   SER A 147      -9.039  44.784  35.353  1.00  7.40           C
ATOM   1115  O   SER A 147      -9.760  45.685  35.789  1.00  8.68           O
ATOM   1116  CB  SER A 147      -7.865  45.106  33.276  1.00  8.66           C
ATOM   1117  OG  SER A 147      -6.653  45.391  32.634  1.00 10.06           O
HETATM 6067 MO   6MO A 784      -5.260  46.566  32.098  1.00 12.66          MO
END
"""

pdb_1qgu = """\
CRYST1  204.030   75.300  163.010  90.00 122.61  90.00 C 1 2 1
HETATM15617  C1  HCA C 501      96.547 -20.184  38.556  1.00 10.87           C
HETATM15618  C2  HCA C 501      95.213 -20.049  39.510  1.00  9.95           C
HETATM15619  C3  HCA C 501      94.168 -21.043  39.241  1.00  9.42           C
HETATM15620  C4  HCA C 501      93.402 -21.123  37.915  1.00  9.49           C
HETATM15621  C5  HCA C 501      92.165 -22.074  38.039  1.00 10.00           C
HETATM15622  C6  HCA C 501      91.532 -22.526  36.720  1.00 11.86           C
HETATM15623  C7  HCA C 501      93.096 -21.107  40.367  1.00 13.15           C
HETATM15624  O1  HCA C 501      97.036 -19.182  38.219  1.00 14.62           O
HETATM15625  O2  HCA C 501      96.788 -21.440  38.356  1.00 14.37           O
HETATM15626  O3  HCA C 501      92.357 -22.713  35.738  1.00 11.52           O
HETATM15627  O4  HCA C 501      90.312 -22.698  36.698  1.00 11.46           O
HETATM15628  O5  HCA C 501      92.878 -22.233  40.972  1.00  8.73           O
HETATM15629  O6  HCA C 501      92.485 -20.140  40.624  1.00 10.26           O
HETATM15630  O7  HCA C 501      94.771 -22.430  39.173  1.00 10.32           O
HETATM15631  S1A CFM C 503      99.172 -27.139  43.981  1.00 11.50           S
HETATM15632  S1B CFM C 503      96.230 -25.163  40.257  1.00 10.12           S
HETATM15633  S2A CFM C 503      98.402 -24.308  46.005  1.00 12.46           S
HETATM15634  S2B CFM C 503      99.085 -23.608  42.289  1.00 10.77           S
HETATM15635  S3A CFM C 503      95.770 -28.360  42.506  1.00 10.89           S
HETATM15636  S3B CFM C 503      95.540 -22.310  42.410  1.00 10.01           S
HETATM15637  S4A CFM C 503      96.310 -27.153  46.228  1.00 12.81           S
HETATM15638  S4B CFM C 503      93.369 -25.214  42.535  1.00 10.35           S
HETATM15639  S5  CFM C 503      94.497 -23.760  45.779  1.00 10.47           S
HETATM15640 FE1  CFM C 503      98.455 -26.515  46.066  1.00 13.45          FE
HETATM15641 FE2  CFM C 503      98.455 -24.979  43.889  1.00 12.40          FE
HETATM15642 FE3  CFM C 503      96.319 -25.007  45.553  1.00 12.09          FE
HETATM15643 FE4  CFM C 503      96.850 -27.165  44.033  1.00 12.50          FE
HETATM15644 FE5  CFM C 503      95.357 -26.176  42.143  1.00 11.88          FE
HETATM15645 FE6  CFM C 503      96.930 -23.985  42.044  1.00 11.39          FE
HETATM15646 FE7  CFM C 503      94.832 -24.050  43.644  1.00 11.29          FE
HETATM15647 MO1  CFM C 503      94.466 -23.866  41.010  1.00 12.29          MO
END
"""

NOT_PLACED = 'Alcohol H not placed'

def run(prefix, pdb_str):
  """Run reduce2 with its defaults on a model string in the current directory.

  Parameters
  ----------
  prefix : str
    Prefix of the input, output and description files written.
  pdb_str : str
    The input model in PDB format.

  Returns
  -------
  tuple
    The names of the H atoms of the model that reduce2 returns and the text
    of the description file.
  """
  from libtbx.utils import null_out
  from iotbx.cli_parser import run_program
  from mmtbx.programs import reduce2
  fn, txt = prefix + '_in.pdb', prefix + '.txt'
  with open(fn, 'w') as f:
    f.write(pdb_str)
  result = run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn, 'output.filename=' + prefix + '_out.pdb',
          'output.description_file_name=' + txt, '--overwrite'])
  with open(txt) as f:
    description = f.read()
  names = set()
  for a in result.model.get_hierarchy().atoms():
    if a.element_is_hydrogen():
      names.add(a.name.strip())
  return names, description

def test_1dmr_mo_ser():
  '''Ser 147 on Mo loses HG; with Zn in place of Mo it keeps HG'''
  names, description = run('tst_reduce2_metal_alkoxide_mo', pdb_1dmr)
  assert 'HG' not in names, sorted(names)
  assert ('HG  SER A 147  (OG  SER A 147 - MO   6MO A 784 1.90 A '
    '(geometric))') in description, description
  pdb_zn = pdb_1dmr.replace('MO   6MO A 784', 'ZN    ZN A 784').replace(
    '          MO\n', '          ZN\n')
  assert '6MO' not in pdb_zn
  names, description = run('tst_reduce2_metal_alkoxide_zn', pdb_zn)
  assert 'HG' in names, sorted(names)
  assert NOT_PLACED not in description, description

def test_1qgu_homocitrate():
  '''Homocitrate O7 on Mo1 of FeMo-co keeps HO7'''
  names, description = run('tst_reduce2_metal_alkoxide_hca', pdb_1qgu)
  assert 'HO7' in names, sorted(names)
  assert NOT_PLACED not in description, description

if __name__ == '__main__':
  test_1dmr_mo_ser()
  test_1qgu_homocitrate()
  print('OK')
