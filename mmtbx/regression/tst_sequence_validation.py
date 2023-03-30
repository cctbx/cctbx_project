
from __future__ import absolute_import, division, print_function
from iotbx.data_manager import DataManager
from libtbx import easy_mp
from libtbx import easy_pickle
from libtbx.utils import Sorry, null_out
import os
import tempfile

def exercise():
  import libtbx.utils
  if (libtbx.utils.detect_multiprocessing_problem() is not None):
    print("multiprocessing not available, skipping this test")
    return
  if (os.name == "nt"):
    print("easy_mp fixed_func not supported under Windows, skipping this test")
    return
  from mmtbx.validation.sequence import validation, get_sequence_n_copies, \
    get_sequence_n_copies_from_files
  import iotbx.bioinformatics
  import iotbx.pdb
  import libtbx.load_env # import dependency
  from libtbx.test_utils import Exception_expected, contains_lines, approx_equal
  from six.moves import cStringIO as StringIO
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM      2  CA  ARG A  10      -6.299  36.344   7.806  1.00 55.20           C
ATOM     25  CA  TYR A  11      -3.391  33.962   7.211  1.00 40.56           C
ATOM     46  CA  ALA A  12      -0.693  34.802   4.693  1.00 67.95           C
ATOM     56  CA  ALA A  13       0.811  31.422   3.858  1.00 57.97           C
ATOM     66  CA  GLY A  14       4.466  31.094   2.905  1.00 49.24           C
ATOM     73  CA  ALA A  15       7.163  28.421   2.671  1.00 54.70           C
ATOM     83  CA  ILE A  16       6.554  24.685   2.957  1.00 51.79           C
ATOM    102  CA  LEU A  17       7.691  23.612   6.406  1.00 42.30           C
ATOM    121  CA  PTY A  18       7.292  19.882   5.861  1.00 36.68           C
ATOM    128  CA  PHE A  19       5.417  16.968   4.327  1.00 44.99           C
ATOM    148  CA  GLY A  20       3.466  14.289   6.150  1.00 41.99           C
ATOM    155  CA  GLY A  21       1.756  11.130   4.965  1.00 35.77           C
ATOM    190  CA  ALA A  24       1.294  19.658   3.683  1.00 47.02           C
ATOM    200  CA  VAL A  24A      2.361  22.009   6.464  1.00 37.13           C
ATOM    216  CA  HIS A  25       2.980  25.633   5.535  1.00 42.52           C
ATOM    234  CA  LEU A  26       4.518  28.425   7.577  1.00 47.63           C
ATOM    253  CA  ALA A  27       2.095  31.320   7.634  1.00 38.61           C
ATOM    263  CA  ARG A  28       1.589  34.719   9.165  1.00 37.04           C
END""")
  seq1 = iotbx.bioinformatics.sequence("MTTPSHLSDRYELGEILGFGGMSEVHLARD".lower())
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=[seq1],
    log=null_out(),
    nproc=1)
  out = StringIO()
  v.show(out=out)
  assert contains_lines(out.getvalue(), """\
  sequence identity: 76.47%
  13 residue(s) missing from PDB chain (9 at start, 1 at end)
  2 gap(s) in chain
  4 mismatches to sequence
    residue IDs:  12 13 15 24""")
  cif_block = v.sequence_as_cif_block()
  assert list(cif_block['_struct_ref.pdbx_seq_one_letter_code']) == [
    ';MTTPSHLSDRYELGEILGFGGMSEVHLARD\n;']
  # assert approx_equal(cif_block['_struct_ref_seq.pdbx_auth_seq_align_beg'],
  #                     ['10', '14', '16', '19', '24'])
  # assert approx_equal(cif_block['_struct_ref_seq.pdbx_auth_seq_align_end'],
  #                     ['11', '14', '17', '21', '28'])
  # assert approx_equal(cif_block['_struct_ref_seq.db_align_beg'],
  #                     ['10', '14', '16', '19', '25'])
  # assert approx_equal(cif_block['_struct_ref_seq.db_align_end'],
  #                     ['11', '14', '17', '21', '29'])
  # assert cif_block['_struct_ref_seq.pdbx_seq_align_beg_ins_code'][4] == 'A'
  seq2 = iotbx.bioinformatics.sequence("MTTPSHLSDRYELGEILGFGGMSEVHLA")
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=[seq2],
    log=null_out(),
    nproc=1)
  out = StringIO()
  v.show(out=out)
  assert contains_lines(out.getvalue(), """\
  1 residues not found in sequence
    residue IDs:  28""")
  try :
    v = validation(
      pdb_hierarchy=pdb_in.construct_hierarchy(),
      sequences=[],
      log=null_out(),
      nproc=1)
  except AssertionError :
    pass
  else :
    raise Exception_expected
  cif_block = v.sequence_as_cif_block()
  print(list(cif_block['_struct_ref.pdbx_seq_one_letter_code']))
  assert list(cif_block['_struct_ref.pdbx_seq_one_letter_code']) == [
    ';MTTPSHLSDRYELGEILGFGGMSEVHLA\n;']
  # assert approx_equal(cif_block['_struct_ref_seq.pdbx_auth_seq_align_end'],
  #                     ['11', '14', '17', '21', '27'])
  # assert approx_equal(cif_block['_struct_ref_seq.db_align_end'],
  #                     ['11', '14', '17', '21', '28'])
  #
  pdb_in2 = iotbx.pdb.input(source_info=None, lines="""\
ATOM      2  CA  ARG A  10      -6.299  36.344   7.806  1.00 55.20           C
ATOM     25  CA  TYR A  11      -3.391  33.962   7.211  1.00 40.56           C
ATOM     46  CA  ALA A  12      -0.693  34.802   4.693  1.00 67.95           C
ATOM     56  CA  ALA A  13       0.811  31.422   3.858  1.00 57.97           C
ATOM     66  CA  GLY A  14       4.466  31.094   2.905  1.00 49.24           C
ATOM     73  CA  ALA A  15       7.163  28.421   2.671  1.00 54.70           C
ATOM     83  CA  ILE A  16       6.554  24.685   2.957  1.00 51.79           C
ATOM    102  CA  LEU A  17       7.691  23.612   6.406  1.00 42.30           C
TER
ATOM   1936  P     G B   2     -22.947 -23.615  15.323  1.00123.20           P
ATOM   1959  P     C B   3     -26.398 -26.111  19.062  1.00110.06           P
ATOM   1979  P     U B   4     -29.512 -30.638  21.164  1.00101.06           P
ATOM   1999  P     C B   5     -30.524 -36.109  21.527  1.00 92.76           P
ATOM   2019  P     U B   6     -28.684 -41.458  21.223  1.00 87.42           P
ATOM   2062  P     G B   8     -18.396 -45.415  21.903  1.00 80.35           P
ATOM   2085  P     A B   9     -13.852 -43.272  24.156  1.00 77.76           P
ATOM   2107  P     G B  10      -8.285 -44.242  26.815  1.00 79.86           P
END
""")
  seq3 = iotbx.bioinformatics.sequence("AGCUUUGGAG")
  v = validation(
    pdb_hierarchy=pdb_in2.construct_hierarchy(),
    sequences=[seq2,seq3],
    log=null_out(),
    nproc=1,
    extract_coordinates=True)
  out = StringIO()
  v.show(out=out)
  cif_block = v.sequence_as_cif_block()
  assert approx_equal(cif_block['_struct_ref.pdbx_seq_one_letter_code'],
                      [';MTTPSHLSDRYELGEILGFGGMSEVHLA\n;', ';AGCUUUGGAG\n;'])
  # assert approx_equal(cif_block['_struct_ref_seq.pdbx_auth_seq_align_beg'],
  #                     ['10', '14', '16', '2', '6', '8'])
  # assert approx_equal(cif_block['_struct_ref_seq.pdbx_auth_seq_align_end'],
  #                     ['11', '14', '17', '4', '6', '10'])
  assert (len(v.chains[0].get_outliers_table()) == 3)
  assert (len(v.get_table_data()) == 4)
  assert approx_equal(
    v.chains[0].get_mean_coordinate_for_alignment_range(11,11),
    (-0.693, 34.802, 4.693))
  assert approx_equal(
    v.chains[0].get_mean_coordinate_for_alignment_range(11,14),
    (2.93675, 31.43475, 3.53175))
  assert (v.chains[0].get_highlighted_residues() == [11,12,14])
  assert contains_lines(out.getvalue(), """\
  3 mismatches to sequence
    residue IDs:  12 13 15""")
  assert contains_lines(out.getvalue(), """\
  sequence identity: 87.50%
  2 residue(s) missing from PDB chain (1 at start, 0 at end)
  1 gap(s) in chain
  1 mismatches to sequence
    residue IDs:  5""")
  s = easy_pickle.dumps(v)
  seq4 = iotbx.bioinformatics.sequence("")
  try :
    v = validation(
      pdb_hierarchy=pdb_in2.construct_hierarchy(),
      sequences=[seq4],
      log=null_out(),
      nproc=1,
      extract_coordinates=True)
  except AssertionError :
    pass
  else :
    raise Exception_expected
  # check that nucleic acid chain doesn't get aligned against protein sequence
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM  18932  P  B DG D   1     -12.183  60.531  25.090  0.50364.79           P
ATOM  18963  P  B DG D   2      -9.738  55.258  20.689  0.50278.77           P
ATOM  18994  P  B DA D   3     -10.119  47.855  19.481  0.50355.17           P
ATOM  19025  P  B DT D   4     -13.664  42.707  21.119  0.50237.06           P
ATOM  19056  P  B DG D   5     -19.510  39.821  21.770  0.50255.45           P
ATOM  19088  P  B DA D   6     -26.096  40.001  21.038  0.50437.49           P
ATOM  19120  P  B DC D   7     -31.790  41.189  18.413  0.50210.00           P
ATOM  19149  P  B DG D   8     -34.639  41.306  12.582  0.50313.99           P
ATOM  19179  P  B DA D   9     -34.987  38.244   6.813  0.50158.92           P
ATOM  19210  P  B DT D  10     -32.560  35.160   1.082  0.50181.38           P
HETATM19241  P  BTSP D  11     -27.614  30.137   0.455  0.50508.17           P
""")
  sequences, _ = iotbx.bioinformatics.fasta_sequence_parse.parse(
    """>4GFH:A|PDBID|CHAIN|SEQUENCE
MSTEPVSASDKYQKISQLEHILKRPDTYIGSVETQEQLQWIYDEETDCMIEKNVTIVPGLFKIFDEILVNAADNKVRDPS
MKRIDVNIHAEEHTIEVKNDGKGIPIEIHNKENIYIPEMIFGHLLTSSNYDDDEKKVTGGRNGYGAKLCNIFSTEFILET
ADLNVGQKYVQKWENNMSICHPPKITSYKKGPSYTKVTFKPDLTRFGMKELDNDILGVMRRRVYDINGSVRDINVYLNGK
SLKIRNFKNYVELYLKSLEKKRQLDNGEDGAAKSDIPTILYERINNRWEVAFAVSDISFQQISFVNSIATTMGGTHVNYI
TDQIVKKISEILKKKKKKSVKSFQIKNNMFIFINCLIENPAFTSQTKEQLTTRVKDFGSRCEIPLEYINKIMKTDLATRM
FEIADANEENALKKSDGTRKSRITNYPKLEDANKAGTKEGYKCTLVLTEGDSALSLAVAGLAVVGRDYYGCYPLRGKMLN
VREASADQILKNAEIQAIKKIMGLQHRKKYEDTKSLRYGHLMIMTDQDHDGSHIKGLIINFLESSFPGLLDIQGFLLEFI
TPIIKVSITKPTKNTIAFYNMPDYEKWREEESHKFTWKQKYYKGLGTSLAQEVREYFSNLDRHLKIFHSLQGNDKDYIDL
AFSKKKADDRKEWLRQYEPGTVLDPTLKEIPISDFINKELILFSLADNIRSIPNVLDGFKPGQRKVLYGCFKKNLKSELK
VAQLAPYVSECTAYHHGEQSLAQTIIGLAQNFVGSNNIYLLLPNGAFGTRATGGKDAAAARYIYTELNKLTRKIFHPADD
PLYKYIQEDEKTVEPEWYLPILPMILVNGAEGIGTGWSTYIPPFNPLEIIKNIRHLMNDEELEQMHPWFRGWTGTIEEIE
PLRYRMYGRIEQIGDNVLEITELPARTWTSTIKEYLLLGLSGNDKIKPWIKDMEEQHDDNIKFIITLSPEEMAKTRKIGF
YERFKLISPISLMNMVAFDPHGKIKKYNSVNEILSEFYYVRLEYYQKRKDHMSERLQWEVEKYSFQVKFIKMIIEKELTV
TNKPRNAIIQELENLGFPRFNKEGKPYYGSPNDEIAEQINDVKGATSDEEDEESSHEDTENVINGPEELYGTYEYLLGMR
IWSLTKERYQKLLKQKQEKETELENLLKLSAKDIWNTDLKAFEVGYQEFLQRDAEAR
>4GFH:D|PDBID|CHAIN|SEQUENCE
GGATGACGATX
""")
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=sequences,
    log=null_out(),
    nproc=1,)
  out = StringIO()
  v.show(out=out)
  assert v.chains[0].n_missing == 0
  assert v.chains[0].n_missing_end == 0
  assert v.chains[0].n_missing_start == 0
  assert len(v.chains[0].alignment.matches()) == 11
  #
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM      2  CA  GLY A   1       1.367   0.551   0.300  1.00  7.71           C
ATOM      6  CA  CYS A   2       2.782   3.785   1.683  1.00  5.18           C
ATOM     12  CA  CYS A   3      -0.375   5.128   3.282  1.00  5.21           C
ATOM     18  CA  SER A   4      -0.870   2.048   5.492  1.00  7.19           C
ATOM     25  CA  LEU A   5       2.786   2.056   6.642  1.00  6.78           C
ATOM     33  CA  PRO A   6       3.212   4.746   9.312  1.00  7.03           C
ATOM     40  CA  PRO A   7       6.870   5.690   8.552  1.00  7.97           C
ATOM     47  CA  CYS A   8       6.021   6.070   4.855  1.00  6.48           C
ATOM     53  CA  ALA A   9       2.812   8.041   5.452  1.00  7.15           C
ATOM     58  CA  LEU A  10       4.739  10.382   7.748  1.00  8.36           C
ATOM     66  CA  SER A  11       7.292  11.200   5.016  1.00  7.00           C
ATOM     73  CA  ASN A  12       4.649  11.435   2.264  1.00  5.40           C
ATOM     81  CA  PRO A  13       1.879  13.433   3.968  1.00  5.97           C
ATOM     88  CA  ASP A  14       0.485  15.371   0.986  1.00  7.70           C
ATOM     96  CA  TYR A  15       0.565  12.245  -1.180  1.00  6.55           C
ATOM    108  CA  CYS A  16      -1.466  10.260   1.363  1.00  7.32           C
ATOM    113  N   NH2 A  17      -2.612  12.308   2.058  1.00  8.11           N
""")
  seq = iotbx.bioinformatics.sequence("GCCSLPPCALSNPDYCX")
  # match last residue
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=[seq],
    log=null_out(),
    nproc=1,)
  out = StringIO()
  v.show(out=out)
  assert v.chains[0].n_missing == 0
  assert v.chains[0].n_missing_end == 0
  assert v.chains[0].n_missing_start == 0
  assert len(v.chains[0].alignment.matches()) == 17
  # ignore non-protein residue
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=[seq],
    log=null_out(),
    nproc=1,
    ignore_hetatm=True)
  out = StringIO()
  v.show(out=out)
  assert v.chains[0].n_missing == 1
  assert v.chains[0].n_missing_end == 1
  assert v.chains[0].n_missing_start == 0
  assert len(v.chains[0].alignment.matches()) == 17
  #
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM   2518  CA  PRO C   3      23.450  -5.848  45.723  1.00 85.24           C
ATOM   2525  CA  GLY C   4      20.066  -4.416  44.815  1.00 79.25           C
ATOM   2529  CA  PHE C   5      19.408  -0.913  46.032  1.00 77.13           C
ATOM   2540  CA  GLY C   6      17.384  -1.466  49.208  1.00 83.44           C
ATOM   2544  CA  GLN C   7      17.316  -5.259  49.606  1.00 89.25           C
ATOM   2553  CA  GLY C   8      19.061  -6.829  52.657  1.00 90.67           C
""")
  sequences, _ = iotbx.bioinformatics.fasta_sequence_parse.parse(
    """>1JN5:A|PDBID|CHAIN|SEQUENCE
MASVDFKTYVDQACRAAEEFVNVYYTTMDKRRRLLSRLYMGTATLVWNGNAVSGQESLSEFFEMLPSSEFQISVVDCQPV
HDEATPSQTTVLVVICGSVKFEGNKQRDFNQNFILTAQASPSNTVWKIASDCFRFQDWAS
>1JN5:B|PDBID|CHAIN|SEQUENCE
APPCKGSYFGTENLKSLVLHFLQQYYAIYDSGDRQGLLDAYHDGACCSLSIPFIPQNPARSSLAEYFKDSRNVKKLKDPT
LRFRLLKHTRLNVVAFLNELPKTQHDVNSFVVDISAQTSTLLCFSVNGVFKEVDGKSRDSLRAFTRTFIAVPASNSGLCI
VNDELFVRNASSEEIQRAFAMPAPTPSSSPVPTLSPEQQEMLQAFSTQSGMNLEWSQKCLQDNNWDYTRSAQAFTHLKAK
GEIPEVAFMK
>1JN5:C|PDBID|CHAIN|SEQUENCE
GQSPGFGQGGSV
""")
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=sequences,
    log=null_out(),
    nproc=1,)
  out = StringIO()
  v.show(out=out)
  assert v.chains[0].n_missing_start == 3
  assert v.chains[0].n_missing_end == 3
  assert v.chains[0].identity == 1.0
  assert v.chains[0].alignment.match_codes == 'iiimmmmmmiii'
  #
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM      2  CA  ALA A   2      -8.453  57.214 -12.754  1.00 52.95           C
ATOM      7  CA  LEU A   3      -8.574  59.274  -9.471  1.00 24.33           C
ATOM     15  CA  ARG A   4     -12.178  60.092  -8.575  1.00 28.40           C
ATOM     26  CA  GLY A   5     -14.170  61.485  -5.667  1.00 26.54           C
ATOM     30  CA  THR A   6     -17.784  60.743  -4.783  1.00 31.78           C
ATOM     37  CA  VAL A   7     -19.080  64.405  -4.464  1.00 21.31           C
""")
  seq = iotbx.bioinformatics.sequence("XALRGTV")
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=[seq],
    log=null_out(),
    nproc=1,)
  out = StringIO()
  v.show(out=out)
  assert v.chains[0].n_missing_start == 1
  assert v.chains[0].n_missing_end == 0
  assert v.chains[0].identity == 1.0
  assert v.chains[0].alignment.match_codes == 'immmmmm'
  #
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM   2171  CA  ASP I 355       5.591 -11.903   1.133  1.00 41.60           C
ATOM   2175  CA  PHE I 356       7.082  -8.454   0.828  1.00 39.82           C
ATOM   2186  CA  GLU I 357       5.814  -6.112  -1.877  1.00 41.12           C
ATOM   2195  CA  GLU I 358       8.623  -5.111  -4.219  1.00 42.70           C
ATOM   2199  CA  ILE I 359      10.346  -1.867  -3.363  1.00 43.32           C
ATOM   2207  CA  PRO I 360      11.658   0.659  -5.880  1.00 44.86           C
ATOM   2214  CA  GLU I 361      14.921  -0.125  -7.592  1.00 44.32           C
ATOM   2219  CA  GLU I 362      15.848   3.489  -6.866  1.00 44.27           C
HETATM 2224  CA  TYS I 363      16.482   2.005  -3.448  1.00 44.52           C
""")
  seq = iotbx.bioinformatics.sequence("NGDFEEIPEEYL")
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=[seq],
    log=null_out(),
    nproc=1,)
  out = StringIO()
  v.show(out=out)
  assert v.chains[0].n_missing_start == 2
  assert v.chains[0].n_missing_end == 1
  assert v.chains[0].identity == 1.0
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM    450  CA  ASN A   1      37.242  41.665  44.160  1.00 35.89           C
ATOM    458  CA  GLY A   2      37.796  38.269  42.523  1.00 30.13           C
HETATM  463  CA AMSE A   3      35.878  39.005  39.326  0.54 22.83           C
HETATM  464  CA BMSE A   3      35.892  39.018  39.323  0.46 22.96           C
ATOM    478  CA  ILE A   4      37.580  38.048  36.061  1.00 22.00           C
ATOM    486  CA  SER A   5      37.593  40.843  33.476  1.00 18.73           C
ATOM    819  CA  ALA A   8      25.982  34.781  27.220  1.00 18.43           C
ATOM    824  CA  ALA A   9      23.292  32.475  28.614  1.00 19.60           C
HETATM  830  CA BMSE A  10      22.793  30.814  25.223  0.41 22.60           C
HETATM  831  CA CMSE A  10      22.801  30.850  25.208  0.59 22.54           C
ATOM    845  CA  GLU A  11      26.504  30.054  24.966  1.00 25.19           C
ATOM    854  CA  GLY A  12      25.907  28.394  28.320  1.00 38.88           C
""")
  seq = iotbx.bioinformatics.sequence("NGMISAAAAMEG")
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=[seq],
    log=null_out(),
    nproc=1,)
  out = StringIO()
  v.show(out=out)
  assert v.chains[0].alignment.a == 'NGMISXXAAMEG'
  assert v.chains[0].alignment.b == 'NGMISAAAAMEG'
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM   4615  CA  ALA C   1       1.000   1.000   1.000  1.00 10.00
ATOM   4622  CA  ALA C   2       1.000   1.000   1.000  1.00 10.00
ATOM   4627  CA  ALA C   3       1.000   1.000   1.000  1.00 10.00
ATOM   4634  CA  ALA C   4       1.000   1.000   1.000  1.00 10.00
ATOM   4646  CA  ALA C   5       1.000   1.000   1.000  1.00 10.00
ATOM   4658  CA  ALA C   6       1.000   1.000   1.000  1.00 10.00
ATOM   4664  CA  ALA C   7       1.000   1.000   1.000  1.00 10.00
ATOM   4669  CA  ALA C   8       1.000   1.000   1.000  1.00 10.00
ATOM   4680  CA  ARG C   9       1.000   1.000   1.000  1.00 10.00
ATOM   4690  CA  GLY C  10       1.000   1.000   1.000  1.00 10.00
ATOM   4698  CA  PRO C  11       1.000   1.000   1.000  1.00 10.00
ATOM   4705  CA  LYS C  12       1.000   1.000   1.000  1.00 10.00
ATOM   4712  CA  TRP C  13       1.000   1.000   1.000  1.00 10.00
ATOM   4726  CA  GLU C  14       1.000   1.000   1.000  1.00 10.00
ATOM   4738  CA  SER C  15       1.000   1.000   1.000  1.00 10.00
ATOM   4744  CA  THR C  16       1.000   1.000   1.000  1.00 10.00
ATOM   4751  CA  GLY C  17       1.000   1.000   1.000  1.00 10.00
ATOM   4755  CA  TYR C  18       1.000   1.000   1.000  1.00 10.00
ATOM   4767  CA  PHE C  19       1.000   1.000   1.000  1.00 10.00
ATOM   4778  CA  ALA C  20       1.000   1.000   1.000  1.00 10.00
ATOM   4786  CA  ALA C  21       1.000   1.000   1.000  1.00 10.00
ATOM   4798  CA  TRP C  22       1.000   1.000   1.000  1.00 10.00
ATOM   4812  CA  GLY C  23       1.000   1.000   1.000  1.00 10.00
ATOM   4816  CA  GLN C  24       1.000   1.000   1.000  1.00 10.00
ATOM   4822  CA  GLY C  25       1.000   1.000   1.000  1.00 10.00
ATOM   4826  CA  THR C  26       1.000   1.000   1.000  1.00 10.00
ATOM   4833  CA  LEU C  27       1.000   1.000   1.000  1.00 10.00
ATOM   4841  CA  VAL C  28       1.000   1.000   1.000  1.00 10.00
ATOM   4848  CA  THR C  29       1.000   1.000   1.000  1.00 10.00
ATOM   4855  CA  VAL C  30       1.000   1.000   1.000  1.00 10.00
ATOM   4862  CA  SER C  31       1.000   1.000   1.000  1.00 10.00
ATOM   4868  CA  SER C  32       1.000   1.000   1.000  1.00 10.00
END
""")
  seq = iotbx.bioinformatics.sequence(
    "AAAAAAAARGKWESPAALLKKAAWCSGTLVTVSSASAPKWKSTSGCYFAAPWNKRALRVTVLQSS")
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=[seq],
    log=null_out(),
    nproc=1,)
  out = StringIO()
  v.show(out=out)

  # check that shortest matching sequence is chosen
  # example from 6H4N, chain a, and I
  sequences, _ = iotbx.bioinformatics.fasta_sequence_parse.parse("""\
>6H4N:a|PDBID|CHAIN|SEQUENCE
AAUUGAAGAGUUUGAUCAUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGAA
GCUUGCUUCUUUGCUGACGAGUGGCGGACGGGUGAGUAAUGUCUGGGAAACUGCCUGAUGGAGGGGGAUAACUACUGGAA
ACGGUAGCUAAUACCGCAUAACGUCGCAAGACCAAAGAGGGGGACCUUCGGGCCUCUUGCCAUCGGAUGUGCCCAGAUGG
GAUUAGCUAGUAGGUGGGGUAACGGCUCACCUAGGCGACGAUCCCUAGCUGGUCUGAGAGGAUGACCAGCCACACUGGAA
CUGAGACACGGUCCAGACUCCUACGGGAGGCAGCAGUGGGGAAUAUUGCACAAUGGGCGCAAGCCUGAUGCAGCCAUGCC
GCGUGUAUGAAGAAGGCCUUCGGGUUGUAAAGUACUUUCAGCGGGGAGGAAGGGAGUAAAGUUAAUACCUUUGCUCAUUG
ACGUUACCCGCAGAAGAAGCACCGGCUAACUCCGUGCCAGCAGCCGCGGUAAUACGGAGGGUGCAAGCGUUAAUCGGAAU
UACUGGGCGUAAAGCGCACGCAGGCGGUUUGUUAAGUCAGAUGUGAAAUCCCCGGGCUCAACCUGGGAACUGCAUCUGAU
ACUGGCAAGCUUGAGUCUCGUAGAGGGGGGUAGAAUUCCAGGUGUAGCGGUGAAAUGCGUAGAGAUCUGGAGGAAUACCG
GUGGCGAAGGCGGCCCCCUGGACGAAGACUGACGCUCAGGUGCGAAAGCGUGGGGAGCAAACAGGAUUAGAUACCCUGGU
AGUCCACGCCGUAAACGAUGUCGACUUGGAGGUUGUGCCCUUGAGGCGUGGCUUCCGGAGCUAACGCGUUAAGUCGACCG
CCUGGGGAGUACGGCCGCAAGGUUAAAACUCAAAUGAAUUGACGGGGGCCCGCACAAGCGGUGGAGCAUGUGGUUUAAUU
CGAUGCAACGCGAAGAACCUUACCUGGUCUUGACAUCCACGGAAGUUUUCAGAGAUGAGAAUGUGCCUUCGGGAACCGUG
AGACAGGUGCUGCAUGGCUGUCGUCAGCUCGUGUUGUGAAAUGUUGGGUUAAGUCCCGCAACGAGCGCAACCCUUAUCCU
UUGUUGCCAGCGGUCCGGCCGGGAACUCAAAGGAGACUGCCAGUGAUAAACUGGAGGAAGGUGGGGAUGACGUCAAGUCA
UCAUGGCCCUUACGACCAGGGCUACACACGUGCUACAAUGGCGCAUACAAAGAGAAGCGACCUCGCGAGAGCAAGCGGAC
CUCAUAAAGUGCGUCGUAGUCCGGAUUGGAGUCUGCAACUCGACUCCAUGAAGUCGGAAUCGCUAGUAAUCGUGGAUCAG
AAUGCCACGGUGAAUACGUUCCCGGGCCUUGUACACACCGCCCGUCACACCAUGGGAGUGGGUUGCAAAAGAAGUAGGUA
GCUUAACCUUCGGGAGGGCGCUUACCACUUUGUGAUUCAUGACUGGGGUGAAGUCGUAACAAGGUAACCGUAGGGGAACC
UGCGGUUGGAUCAC
>6H4N:I|PDBID|CHAIN|SEQUENCE
CUCCU
""")
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM  95502  P     C I1536     211.989 143.717 147.208  1.00 16.47           P
ATOM  95503  OP1   C I1536     213.292 143.696 146.494  1.00 16.47           O
ATOM  95504  OP2   C I1536     211.250 144.996 147.359  1.00 16.47           O
ATOM  95505  O5'   C I1536     211.021 142.666 146.541  1.00 16.47           O
ATOM  95506  C5'   C I1536     211.671 141.536 146.021  1.00 16.47           C
ATOM  95507  C4'   C I1536     211.059 140.260 146.502  1.00 16.47           C
ATOM  95508  O4'   C I1536     209.764 140.432 147.128  1.00 16.47           O
ATOM  95509  C3'   C I1536     210.818 139.353 145.303  1.00 16.47           C
ATOM  95510  O3'   C I1536     211.011 137.993 145.604  1.00 16.47           O
ATOM  95511  C2'   C I1536     209.372 139.646 144.938  1.00 16.47           C
ATOM  95512  O2'   C I1536     208.735 138.572 144.276  1.00 16.47           O
ATOM  95513  C1'   C I1536     208.757 139.866 146.316  1.00 16.47           C
ATOM  95514  N1    C I1536     207.618 140.788 146.322  1.00 16.47           N
ATOM  95515  C2    C I1536     206.610 140.626 145.378  1.00 16.47           C
ATOM  95516  O2    C I1536     206.712 139.721 144.535  1.00 16.47           O
ATOM  95517  N3    C I1536     205.560 141.463 145.396  1.00 16.47           N
ATOM  95518  C4    C I1536     205.492 142.420 146.320  1.00 16.47           C
ATOM  95519  N4    C I1536     204.429 143.227 146.302  1.00 16.47           N
ATOM  95520  C5    C I1536     206.496 142.595 147.306  1.00 16.47           C
ATOM  95521  C6    C I1536     207.522 141.754 147.283  1.00 16.47           C
ATOM  95522  P     U I1537     212.458 137.366 145.505  1.00 11.96           P
ATOM  95523  OP1   U I1537     212.292 135.894 145.567  1.00 11.96           O
ATOM  95524  OP2   U I1537     213.344 138.045 146.479  1.00 11.96           O
ATOM  95525  O5'   U I1537     212.962 137.720 144.038  1.00 11.96           O
ATOM  95526  C5'   U I1537     214.363 137.934 143.772  1.00 11.96           C
ATOM  95527  C4'   U I1537     214.522 138.678 142.472  1.00 11.96           C
ATOM  95528  O4'   U I1537     213.714 137.951 141.515  1.00 11.96           O
ATOM  95529  C3'   U I1537     213.970 140.098 142.549  1.00 11.96           C
ATOM  95530  O3'   U I1537     214.924 141.159 142.799  1.00 11.96           O
ATOM  95531  C2'   U I1537     212.939 140.210 141.413  1.00 11.96           C
ATOM  95532  O2'   U I1537     212.980 141.292 140.508  1.00 11.96           O
ATOM  95533  C1'   U I1537     212.990 138.848 140.714  1.00 11.96           C
ATOM  95534  N1    U I1537     211.632 138.324 140.509  1.00 11.96           N
ATOM  95535  C2    U I1537     211.212 138.082 139.216  1.00 11.96           C
ATOM  95536  O2    U I1537     211.943 138.228 138.252  1.00 11.96           O
ATOM  95537  N3    U I1537     209.897 137.730 139.076  1.00 11.96           N
ATOM  95538  C4    U I1537     208.966 137.602 140.074  1.00 11.96           C
ATOM  95539  O4    U I1537     207.834 137.203 139.798  1.00 11.96           O
ATOM  95540  C5    U I1537     209.473 137.843 141.382  1.00 11.96           C
ATOM  95541  C6    U I1537     210.749 138.206 141.544  1.00 11.96           C
ATOM  95542  P     C I1538     216.031 141.722 141.738  1.00 11.10           P
ATOM  95543  OP1   C I1538     216.814 142.772 142.428  1.00 11.10           O
ATOM  95544  OP2   C I1538     215.385 142.057 140.453  1.00 11.10           O
ATOM  95545  O5'   C I1538     217.081 140.541 141.538  1.00 11.10           O
ATOM  95546  C5'   C I1538     218.494 140.848 141.429  1.00 11.10           C
ATOM  95547  C4'   C I1538     218.962 140.916 139.986  1.00 11.10           C
ATOM  95548  O4'   C I1538     218.034 140.280 139.091  1.00 11.10           O
ATOM  95549  C3'   C I1538     219.276 142.298 139.408  1.00 11.10           C
ATOM  95550  O3'   C I1538     220.629 142.126 139.044  1.00 11.10           O
ATOM  95551  C2'   C I1538     218.657 142.315 138.005  1.00 11.10           C
ATOM  95552  O2'   C I1538     219.358 142.774 136.857  1.00 11.10           O
ATOM  95553  C1'   C I1538     218.164 140.883 137.832  1.00 11.10           C
ATOM  95554  N1    C I1538     216.943 140.702 137.064  1.00 11.10           N
ATOM  95555  C2    C I1538     217.041 140.096 135.813  1.00 11.10           C
ATOM  95556  O2    C I1538     218.163 139.770 135.401  1.00 11.10           O
ATOM  95557  N3    C I1538     215.932 139.850 135.093  1.00 11.10           N
ATOM  95558  C4    C I1538     214.748 140.195 135.580  1.00 11.10           C
ATOM  95559  N4    C I1538     213.670 139.968 134.827  1.00 11.10           N
ATOM  95560  C5    C I1538     214.617 140.827 136.842  1.00 11.10           C
ATOM  95561  C6    C I1538     215.722 141.024 137.566  1.00 11.10           C
ATOM  95562  P     C I1539     221.798 142.624 139.940  1.00 17.77           P
ATOM  95563  OP1   C I1539     221.300 143.669 140.865  1.00 17.77           O
ATOM  95564  OP2   C I1539     222.961 142.899 139.061  1.00 17.77           O
ATOM  95565  O5'   C I1539     222.148 141.341 140.812  1.00 17.77           O
ATOM  95566  C5'   C I1539     223.493 140.934 140.997  1.00 17.77           C
ATOM  95567  C4'   C I1539     223.633 139.444 140.845  1.00 17.77           C
ATOM  95568  O4'   C I1539     222.661 138.972 139.877  1.00 17.77           O
ATOM  95569  C3'   C I1539     224.967 138.959 140.300  1.00 17.77           C
ATOM  95570  O3'   C I1539     225.970 138.853 141.295  1.00 17.77           O
ATOM  95571  C2'   C I1539     224.602 137.629 139.658  1.00 17.77           C
ATOM  95572  O2'   C I1539     224.482 136.616 140.642  1.00 17.77           O
ATOM  95573  C1'   C I1539     223.209 137.924 139.109  1.00 17.77           C
ATOM  95574  N1    C I1539     223.219 138.333 137.681  1.00 17.77           N
ATOM  95575  C2    C I1539     223.353 137.370 136.683  1.00 17.77           C
ATOM  95576  O2    C I1539     223.476 136.178 136.982  1.00 17.77           O
ATOM  95577  N3    C I1539     223.342 137.742 135.392  1.00 17.77           N
ATOM  95578  C4    C I1539     223.202 139.017 135.059  1.00 17.77           C
ATOM  95579  N4    C I1539     223.202 139.332 133.762  1.00 17.77           N
ATOM  95580  C5    C I1539     223.059 140.033 136.041  1.00 17.77           C
ATOM  95581  C6    C I1539     223.067 139.642 137.318  1.00 17.77           C
ATOM  95582  P     U I1540     227.517 139.071 140.915  1.00 25.44           P
ATOM  95583  OP1   U I1540     228.321 138.910 142.156  1.00 25.44           O
ATOM  95584  OP2   U I1540     227.626 140.309 140.102  1.00 25.44           O
ATOM  95585  O5'   U I1540     227.868 137.833 139.978  1.00 25.44           O
ATOM  95586  C5'   U I1540     228.014 136.524 140.520  1.00 25.44           C
ATOM  95587  C4'   U I1540     228.308 135.503 139.447  1.00 25.44           C
ATOM  95588  O4'   U I1540     227.513 135.808 138.268  1.00 25.44           O
ATOM  95589  C3'   U I1540     229.761 135.445 138.980  1.00 25.44           C
ATOM  95590  O3'   U I1540     230.104 134.098 138.659  1.00 25.44           O
ATOM  95591  C2'   U I1540     229.740 136.281 137.705  1.00 25.44           C
ATOM  95592  O2'   U I1540     230.767 135.976 136.785  1.00 25.44           O
ATOM  95593  C1'   U I1540     228.360 135.950 137.145  1.00 25.44           C
ATOM  95594  N1    U I1540     227.809 136.996 136.268  1.00 25.44           N
ATOM  95595  C2    U I1540     227.053 136.589 135.186  1.00 25.44           C
ATOM  95596  O2    U I1540     226.815 135.418 134.956  1.00 25.44           O
ATOM  95597  N3    U I1540     226.574 137.600 134.393  1.00 25.44           N
ATOM  95598  C4    U I1540     226.781 138.951 134.566  1.00 25.44           C
ATOM  95599  O4    U I1540     226.286 139.746 133.765  1.00 25.44           O
ATOM  95600  C5    U I1540     227.583 139.293 135.701  1.00 25.44           C
ATOM  95601  C6    U I1540     228.061 138.329 136.493  1.00 25.44           C
END
""")
  v = validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    sequences=sequences,
    log=null_out(),
    nproc=1,)
  assert(v.chains[0].get_alignment() == ['CUCCU', 'CUCCU'])

  # all tests below here have additional dependencies
  if (not libtbx.env.has_module("ksdssp")):
    print("Skipping advanced tests (require ksdssp module)")
    return
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  if (pdb_file is not None):
    seq = iotbx.bioinformatics.sequence("MGSSHHHHHHSSGLVPRGSHMAVRELPGAWNFRDVADTATALRPGRLFRSSELSRLDDAGRATLRRLGITDVADLRSSREVARRGPGRVPDGIDVHLLPFPDLADDDADDSAPHETAFKRLLTNDGSNGESGESSQSINDAATRYMTDEYRQFPTRNGAQRALHRVVTLLAAGRPVLTHCFAGKDRTGFVVALVLEAVGLDRDVIVADYLRSNDSVPQLRARISEMIQQRFDTELAPEVVTFTKARLSDGVLGVRAEYLAAARQTIDETYGSLGGYLRDAGISQATVNRMRGVLLG")
    hierarchy = iotbx.pdb.input(pdb_file).construct_hierarchy()
    v = validation(
      pdb_hierarchy=hierarchy,
      sequences=[seq],
      log=null_out(),
      nproc=1,
      include_secondary_structure=True,
      extract_coordinates=True)
    out = StringIO()
    v.show(out=out)
    aln1, aln2, ss = v.chains[0].get_alignment(include_sec_str=True)
    assert ("HHH" in ss) and ("LLL" in ss) and ("---" in ss)
    cif_block = v.sequence_as_cif_block()
    assert cif_block['_struct_ref.pdbx_seq_one_letter_code'] == seq.sequence
    # assert list(
    #   cif_block['_struct_ref_seq.pdbx_auth_seq_align_beg']) == ['4', '117']
    # assert list(
    #   cif_block['_struct_ref_seq.pdbx_auth_seq_align_end']) == ['85', '275']
    # assert list(cif_block['_struct_ref_seq.seq_align_beg']) == ['1', '114']
    # assert list(cif_block['_struct_ref_seq.seq_align_end']) == ['82', '272']
    # determine relative counts of sequences and chains
    n_seq = get_sequence_n_copies(
      pdb_hierarchy=hierarchy,
      sequences=[seq] * 4,
      copies_from_xtriage=4,
      out=null_out())
    assert (n_seq == 1)
    hierarchy = hierarchy.deep_copy()
    chain2 = hierarchy.only_model().chains()[0].detached_copy()
    hierarchy.only_model().append_chain(chain2)
    n_seq = get_sequence_n_copies(
      pdb_hierarchy=hierarchy,
      sequences=[seq] * 4,
      copies_from_xtriage=2,
      out=null_out())
    assert (n_seq == 1)
    n_seq = get_sequence_n_copies(
      pdb_hierarchy=hierarchy,
      sequences=[seq],
      copies_from_xtriage=2,
      out=null_out())
    assert (n_seq == 4)
    try :
      n_seq = get_sequence_n_copies(
        pdb_hierarchy=hierarchy,
        sequences=[seq] * 3,
        copies_from_xtriage=2,
        out=null_out())
    except Sorry as s :
      assert ("round number" in str(s))
    else :
      raise Exception_expected
    n_seq = get_sequence_n_copies(
      pdb_hierarchy=hierarchy,
      sequences=[seq] * 3,
      copies_from_xtriage=2,
      force_accept_composition=True,
      out=null_out())
    assert (n_seq == 1)
    try :
      n_seq = get_sequence_n_copies(
        pdb_hierarchy=hierarchy,
        sequences=[seq] * 4,
        copies_from_xtriage=1,
        out=null_out())
    except Sorry as s :
      assert ("less than" in str(s))
    else :
      raise Exception_expected
    n_seq = get_sequence_n_copies(
      pdb_hierarchy=hierarchy,
      sequences=[seq] * 4,
      copies_from_xtriage=1,
      assume_xtriage_copies_from_sequence_file=True,
      out=null_out())
    assert (n_seq == 0.5)
    hierarchy = hierarchy.deep_copy()
    chain2 = hierarchy.only_model().chains()[0].detached_copy()
    hierarchy.only_model().append_chain(chain2)
    try :
      n_seq = get_sequence_n_copies(
        pdb_hierarchy=hierarchy,
        sequences=[seq] * 2,
        copies_from_xtriage=2,
        out=null_out())
    except Sorry as s :
      assert ("round number" in str(s))
    else :
      raise Exception_expected
    n_seq = get_sequence_n_copies(
      pdb_hierarchy=hierarchy,
      sequences=[seq],
      copies_from_xtriage=1,
      out=null_out())
    assert (n_seq == 3)
    hierarchy = hierarchy.deep_copy()
    chain2 = hierarchy.only_model().chains()[0].detached_copy()
    hierarchy.only_model().append_chain(chain2)
    n_seq = get_sequence_n_copies(
      pdb_hierarchy=hierarchy,
      sequences=[seq] * 2,
      copies_from_xtriage=2,
      out=null_out())
    assert (n_seq == 4)
    # now with files as input
    seq_file = "tmp_mmtbx_validation_sequence.fa"
    open(seq_file, "w").write(">1ywf\n%s" % seq.sequence)
    n_seq = get_sequence_n_copies_from_files(
      pdb_file=pdb_file,
      seq_file=seq_file,
      copies_from_xtriage=4,
      out=null_out())
    try :
      assert (n_seq == 4)
    finally :
      os.remove(seq_file)

def test_modified_residues():
  # modified amino acid example from 1e0z
  seq_str1 = """\
>pdb|1e0z|A
PTVEYLNYETLDDQGWDMDDDDLFEKAADAGLDGEDYGTMEVAEGEYILEAAEAQGYDWPFSCRAGACANCASIVKEG
EIDMDMQQILSDEEVEEKDVRLTCIGSPAADEVKIVYNAKHLDYLQNRVI
"""
  pdb_str1 = """\
HETATM 1707  OH  ALY A 118       9.253 -11.285  16.293  1.00  0.00           O
HETATM 1708  CH  ALY A 118       8.541 -10.835  15.417  1.00  0.00           C
HETATM 1709  CH3 ALY A 118       8.924  -9.562  14.659  1.00  0.00           C
HETATM 1710  NZ  ALY A 118       7.420 -11.417  15.087  1.00 50.00           N
HETATM 1711  CE  ALY A 118       7.573 -12.895  15.202  1.00 50.00           C
HETATM 1712  CD  ALY A 118       6.273 -13.470  14.634  1.00 12.50           C
HETATM 1713  CG  ALY A 118       5.506 -14.188  15.746  1.00 12.50           C
HETATM 1714  CB  ALY A 118       6.259 -15.458  16.148  1.00 12.50           C
HETATM 1715  CA  ALY A 118       5.885 -15.846  17.580  1.00 12.50           C
HETATM 1716  N   ALY A 118       4.407 -15.666  17.660  1.00 12.50           N
HETATM 1717  C   ALY A 118       6.576 -14.924  18.588  1.00 12.50           C
HETATM 1718  O   ALY A 118       7.674 -15.190  19.036  1.00 12.50           O
HETATM 1719 HH31 ALY A 118       9.665  -9.016  15.224  1.00  0.00           H
HETATM 1720 HH32 ALY A 118       8.047  -8.945  14.525  1.00  0.00           H
HETATM 1721 HH33 ALY A 118       9.329  -9.826  13.694  1.00  0.00           H
HETATM 1722  HE3 ALY A 118       7.690 -13.181  16.237  1.00  0.00           H
HETATM 1723  HE2 ALY A 118       8.419 -13.232  14.620  1.00  0.00           H
HETATM 1724  HD3 ALY A 118       6.503 -14.171  13.845  1.00  0.00           H
HETATM 1725  HD2 ALY A 118       5.667 -12.668  14.239  1.00  0.00           H
HETATM 1726  HG3 ALY A 118       4.520 -14.451  15.392  1.00  0.00           H
HETATM 1727  HG2 ALY A 118       5.419 -13.535  16.602  1.00  0.00           H
HETATM 1728  HB3 ALY A 118       7.323 -15.279  16.091  1.00  0.00           H
HETATM 1729  HB2 ALY A 118       5.994 -16.261  15.477  1.00  0.00           H
HETATM 1730  HCA ALY A 118       6.150 -16.875  17.765  1.00  0.00           H
HETATM 1731  H   ALY A 118       3.878 -15.580  16.839  1.00 99.00           H
"""
  # modified nucleic acid example from 4eec
  seq_str2 = """\
>4EEC_1|Chains A,B|StaL|Streptomyces toyocaensis (55952)
MGSSHHHHHHSSGLVPRGSMCWIASYPKAGGHWLRCMLTSYVTGEPVETWPGIQAGVPHLEGLLRDGEAPSADPDEQV
LLATHFTADRPVLRFYRESTAKVVCLIRNPRDAMLSLMRMKGIPPEDVEACRKIAETFIADEGFSSVRIWAGEGSWPE
NIRSWTDSVHESFPNAAVLAVRYEDLRKDPEGELWKVVDFLELGGRDGVADAVANCTLERMREMEERSKLLGLETTGL
MTRGGKQLPFVGKGGQRKSLKFMGDDIEKAYADLLHGETDFAHYARLYGYAE
>4EEC_2|Chain C|desulfo-A47934|Streptomyces toyocaensis (55952)
GXXGXXX
"""
  pdb_str2 = """\
HETATM 3866  P1  A3P A 301     -32.928   0.112  -1.515  1.00 60.90           P
HETATM 3867  O1P A3P A 301     -32.567   0.971  -2.721  1.00 58.74           O
HETATM 3868  O2P A3P A 301     -33.836  -1.041  -1.926  1.00 60.18           O
HETATM 3869  O3P A3P A 301     -33.351   0.858  -0.247  1.00 59.07           O
HETATM 3870  P2  A3P A 301     -26.843  -0.504   2.943  1.00 61.59           P
HETATM 3871  O4P A3P A 301     -25.377  -0.388   2.609  1.00 57.50           O
HETATM 3872  O5P A3P A 301     -27.164  -1.892   3.467  1.00 62.98           O
HETATM 3873  O6P A3P A 301     -27.445   0.681   3.730  1.00 62.21           O
HETATM 3874  O5' A3P A 301     -27.694  -0.432   1.563  1.00 60.14           O
HETATM 3875  C5' A3P A 301     -29.109  -0.472   1.636  1.00 57.46           C
HETATM 3876  C4' A3P A 301     -29.649  -0.859   0.288  1.00 58.10           C
HETATM 3877  O4' A3P A 301     -29.334  -2.235   0.029  1.00 59.74           O
HETATM 3878  C3' A3P A 301     -31.164  -0.780   0.284  1.00 58.44           C
HETATM 3879  O3' A3P A 301     -31.539  -0.556  -1.069  1.00 59.37           O
HETATM 3880  C2' A3P A 301     -31.538  -2.197   0.626  1.00 58.35           C
HETATM 3881  O2' A3P A 301     -32.890  -2.542   0.343  1.00 58.12           O
HETATM 3882  C1' A3P A 301     -30.561  -2.939  -0.264  1.00 58.48           C
HETATM 3883  N9  A3P A 301     -30.576  -4.401   0.059  1.00 58.16           N
HETATM 3884  C8  A3P A 301     -30.559  -4.926   1.312  1.00 57.04           C
HETATM 3885  N7  A3P A 301     -30.604  -6.285   1.255  1.00 57.35           N
HETATM 3886  C5  A3P A 301     -30.645  -6.664  -0.042  1.00 55.60           C
HETATM 3887  C6  A3P A 301     -30.683  -7.957  -0.767  1.00 54.29           C
HETATM 3888  N6  A3P A 301     -30.686  -9.129  -0.066  1.00 54.08           N
HETATM 3889  N1  A3P A 301     -30.706  -7.913  -2.130  1.00 53.41           N
HETATM 3890  C2  A3P A 301     -30.703  -6.743  -2.815  1.00 52.74           C
HETATM 3891  N3  A3P A 301     -30.659  -5.530  -2.220  1.00 53.74           N
HETATM 3892  C4  A3P A 301     -30.631  -5.420  -0.849  1.00 56.63           C
"""

  for pdb_str, seq_str in [(pdb_str1, seq_str1), (pdb_str2, seq_str2)]:
    seq_file = tempfile.NamedTemporaryFile(suffix='.fasta', mode='w')
    seq_file.write(seq_str)
    seq_file.flush()
    model_file = tempfile.NamedTemporaryFile(suffix='.pdb', mode='w')
    model_file.write(pdb_str)
    model_file.flush()
    dm = DataManager()
    dm.process_model_file(model_file.name)
    dm.process_sequence_file(seq_file.name)
    model = dm.get_model()
    seq = dm.get_sequence()
    model.set_sequences(seq)
    model_file.close()
    seq_file.close()
    # model._sequence_validation.show()

if (__name__ == "__main__"):
  exercise()
  test_modified_residues()
  print("OK")
