from __future__ import absolute_import, division, print_function
from mmtbx.validation import undowser2
from libtbx.easy_pickle import loads
from iotbx.data_manager import DataManager
import libtbx.load_env
from libtbx.utils import null_out
import iotbx
from mmtbx.programs import probe2
import time
import json
import difflib

pdb_1lpl_str = """MODEL        1
ATOM     85  N   ILE A 146       4.850  13.830  60.452  1.00 31.46           N
ATOM     86  CA  ILE A 146       5.654  15.038  60.565  1.00 28.74           C
ATOM     87  C   ILE A 146       5.963  15.310  62.021  1.00 27.82           C
ATOM     88  O   ILE A 146       5.062  15.510  62.825  1.00 28.89           O
ATOM     89  CB  ILE A 146       4.956  16.280  59.905  1.00 27.45           C
ATOM     90  CG1 ILE A 146       4.674  15.994  58.416  1.00 25.50           C
ATOM     91  CG2 ILE A 146       5.850  17.543  60.072  1.00 23.19           C
ATOM     92  CD1 ILE A 146       3.970  17.121  57.673  1.00 27.89           C
ATOM      0  H   ILE A 146       4.003  13.950  60.365  1.00 31.46           H   new
ATOM      0  HA  ILE A 146       6.481  14.891  60.080  1.00 28.74           H   new
ATOM      0  HB  ILE A 146       4.110  16.447  60.348  1.00 27.45           H   new
ATOM      0 HG12 ILE A 146       5.515  15.805  57.971  1.00 25.50           H   new
ATOM      0 HG13 ILE A 146       4.132  15.192  58.350  1.00 25.50           H   new
ATOM      0 HG21 ILE A 146       5.413  18.306  59.662  1.00 23.19           H   new
ATOM      0 HG22 ILE A 146       5.989  17.719  61.016  1.00 23.19           H   new
ATOM      0 HG23 ILE A 146       6.707  17.393  59.642  1.00 23.19           H   new
ATOM      0 HD11 ILE A 146       3.831  16.862  56.749  1.00 27.89           H   new
ATOM      0 HD12 ILE A 146       3.113  17.299  58.091  1.00 27.89           H   new
ATOM      0 HD13 ILE A 146       4.517  17.921  57.705  1.00 27.89           H   new
ATOM    137  N   GLU A 153      11.342  19.816  53.414  1.00 18.86           N
ATOM    138  CA  GLU A 153      11.328  20.597  52.204  1.00 20.88           C
ATOM    139  C   GLU A 153       9.887  21.060  51.996  1.00 20.14           C
ATOM    140  O   GLU A 153       8.976  20.225  51.945  1.00 18.33           O
ATOM    141  CB  GLU A 153      11.771  19.754  51.024  1.00 20.67           C
ATOM    142  CG  GLU A 153      11.890  20.570  49.780  1.00 28.25           C
ATOM    143  CD  GLU A 153      12.021  19.708  48.543  1.00 35.57           C
ATOM    144  OE1 GLU A 153      11.164  18.803  48.342  1.00 39.43           O
ATOM    145  OE2 GLU A 153      12.981  19.939  47.779  1.00 36.90           O
ATOM      0  H   GLU A 153      10.727  19.216  53.448  1.00 18.86           H   new
ATOM      0  HA  GLU A 153      11.935  21.350  52.274  1.00 20.88           H   new
ATOM      0  HB2 GLU A 153      12.625  19.340  51.224  1.00 20.67           H   new
ATOM      0  HB3 GLU A 153      11.135  19.036  50.881  1.00 20.67           H   new
ATOM      0  HG2 GLU A 153      11.111  21.141  49.694  1.00 28.25           H   new
ATOM      0  HG3 GLU A 153      12.662  21.153  49.850  1.00 28.25           H   new
ATOM    618  N   VAL A 216       4.506  19.838  52.802  1.00 21.99           N
ATOM    619  CA  VAL A 216       5.738  19.509  53.507  1.00 24.27           C
ATOM    620  C   VAL A 216       6.150  18.056  53.346  1.00 26.17           C
ATOM    621  O   VAL A 216       5.340  17.143  53.546  1.00 26.47           O
ATOM    622  CB  VAL A 216       5.593  19.757  55.039  1.00 25.53           C
ATOM    623  CG1 VAL A 216       6.935  19.548  55.734  1.00 23.77           C
ATOM    624  CG2 VAL A 216       5.055  21.160  55.299  1.00 22.33           C
ATOM      0  H   VAL A 216       3.828  19.373  53.053  1.00 21.99           H   new
ATOM      0  HA  VAL A 216       6.410  20.085  53.111  1.00 24.27           H   new
ATOM      0  HB  VAL A 216       4.960  19.119  55.404  1.00 25.53           H   new
ATOM      0 HG11 VAL A 216       6.834  19.705  56.686  1.00 23.77           H   new
ATOM      0 HG12 VAL A 216       7.238  18.638  55.586  1.00 23.77           H   new
ATOM      0 HG13 VAL A 216       7.587  20.168  55.372  1.00 23.77           H   new
ATOM      0 HG21 VAL A 216       4.969  21.302  56.255  1.00 22.33           H   new
ATOM      0 HG22 VAL A 216       5.667  21.815  54.929  1.00 22.33           H   new
ATOM      0 HG23 VAL A 216       4.186  21.256  54.879  1.00 22.33           H   new
ATOM    634  N   VAL A 218       9.047  15.495  54.747  1.00 24.30           N
ATOM    635  CA  VAL A 218       9.980  15.425  55.854  1.00 25.59           C
ATOM    636  C   VAL A 218      11.043  14.371  55.566  1.00 25.70           C
ATOM    637  O   VAL A 218      10.781  13.407  54.869  1.00 25.17           O
ATOM    638  CB  VAL A 218       9.161  15.143  57.136  1.00 26.71           C
ATOM    639  CG1 VAL A 218       9.460  13.782  57.688  1.00 28.12           C
ATOM    640  CG2 VAL A 218       9.374  16.262  58.126  1.00 28.81           C
ATOM      0  H   VAL A 218       8.514  14.824  54.675  1.00 24.30           H   new
ATOM      0  HA  VAL A 218      10.460  16.258  55.979  1.00 25.59           H   new
ATOM      0  HB  VAL A 218       8.215  15.126  56.923  1.00 26.71           H   new
ATOM      0 HG11 VAL A 218       8.933  13.635  58.489  1.00 28.12           H   new
ATOM      0 HG12 VAL A 218       9.237  13.108  57.027  1.00 28.12           H   new
ATOM      0 HG13 VAL A 218      10.403  13.721  57.906  1.00 28.12           H   new
ATOM      0 HG21 VAL A 218       8.861  16.085  58.930  1.00 28.81           H   new
ATOM      0 HG22 VAL A 218      10.316  16.321  58.351  1.00 28.81           H   new
ATOM      0 HG23 VAL A 218       9.083  17.100  57.735  1.00 28.81           H   new
TER     728      ILE A 218
HETATM  731  O   HOH A 503       7.663  18.549  58.083  1.00 71.06           O
HETATM  758  O   HOH A 530      20.116  27.095  51.665  1.00 24.67           O
HETATM  765  O   HOH A 537      14.778  21.226  48.211  1.00 27.98           O
HETATM  772  O   HOH A 544      -5.198  25.804  53.682  1.00 21.10           O
HETATM  800  O   HOH A 572      -6.768  26.204  53.124  1.00 14.71           O
HETATM  805  O   HOH A 577      21.563  28.198  52.938  1.00 31.45           O
HETATM  812  O   HOH A 585      -6.264  27.663  53.754  1.00 21.56           O
ENDMDL
MODEL        2
ATOM     85  N   ILE A 146       4.850  13.830  60.452  1.00 31.46           N
ATOM     86  CA  ILE A 146       5.654  15.038  60.565  1.00 28.74           C
ATOM     87  C   ILE A 146       5.963  15.310  62.021  1.00 27.82           C
ATOM     88  O   ILE A 146       5.062  15.510  62.825  1.00 28.89           O
ATOM     89  CB  ILE A 146       4.956  16.280  59.905  1.00 27.45           C
ATOM     90  CG1 ILE A 146       4.674  15.994  58.416  1.00 25.50           C
ATOM     91  CG2 ILE A 146       5.850  17.543  60.072  1.00 23.19           C
ATOM     92  CD1 ILE A 146       3.970  17.121  57.673  1.00 27.89           C
ATOM      0  H   ILE A 146       4.003  13.950  60.365  1.00 31.46           H   new
ATOM      0  HA  ILE A 146       6.481  14.891  60.080  1.00 28.74           H   new
ATOM      0  HB  ILE A 146       4.110  16.447  60.348  1.00 27.45           H   new
ATOM      0 HG12 ILE A 146       5.515  15.805  57.971  1.00 25.50           H   new
ATOM      0 HG13 ILE A 146       4.132  15.192  58.350  1.00 25.50           H   new
ATOM      0 HG21 ILE A 146       5.413  18.306  59.662  1.00 23.19           H   new
ATOM      0 HG22 ILE A 146       5.989  17.719  61.016  1.00 23.19           H   new
ATOM      0 HG23 ILE A 146       6.707  17.393  59.642  1.00 23.19           H   new
ATOM      0 HD11 ILE A 146       3.831  16.862  56.749  1.00 27.89           H   new
ATOM      0 HD12 ILE A 146       3.113  17.299  58.091  1.00 27.89           H   new
ATOM      0 HD13 ILE A 146       4.517  17.921  57.705  1.00 27.89           H   new
ATOM    137  N   GLU A 153      11.342  19.816  53.414  1.00 18.86           N
ATOM    138  CA  GLU A 153      11.328  20.597  52.204  1.00 20.88           C
ATOM    139  C   GLU A 153       9.887  21.060  51.996  1.00 20.14           C
ATOM    140  O   GLU A 153       8.976  20.225  51.945  1.00 18.33           O
ATOM    141  CB  GLU A 153      11.771  19.754  51.024  1.00 20.67           C
ATOM    142  CG  GLU A 153      11.890  20.570  49.780  1.00 28.25           C
ATOM    143  CD  GLU A 153      12.021  19.708  48.543  1.00 35.57           C
ATOM    144  OE1 GLU A 153      11.164  18.803  48.342  1.00 39.43           O
ATOM    145  OE2 GLU A 153      12.981  19.939  47.779  1.00 36.90           O
ATOM      0  H   GLU A 153      10.727  19.216  53.448  1.00 18.86           H   new
ATOM      0  HA  GLU A 153      11.935  21.350  52.274  1.00 20.88           H   new
ATOM      0  HB2 GLU A 153      12.625  19.340  51.224  1.00 20.67           H   new
ATOM      0  HB3 GLU A 153      11.135  19.036  50.881  1.00 20.67           H   new
ATOM      0  HG2 GLU A 153      11.111  21.141  49.694  1.00 28.25           H   new
ATOM      0  HG3 GLU A 153      12.662  21.153  49.850  1.00 28.25           H   new
ATOM    618  N   VAL A 216       4.506  19.838  52.802  1.00 21.99           N
ATOM    619  CA  VAL A 216       5.738  19.509  53.507  1.00 24.27           C
ATOM    620  C   VAL A 216       6.150  18.056  53.346  1.00 26.17           C
ATOM    621  O   VAL A 216       5.340  17.143  53.546  1.00 26.47           O
ATOM    622  CB  VAL A 216       5.593  19.757  55.039  1.00 25.53           C
ATOM    623  CG1 VAL A 216       6.935  19.548  55.734  1.00 23.77           C
ATOM    624  CG2 VAL A 216       5.055  21.160  55.299  1.00 22.33           C
ATOM      0  H   VAL A 216       3.828  19.373  53.053  1.00 21.99           H   new
ATOM      0  HA  VAL A 216       6.410  20.085  53.111  1.00 24.27           H   new
ATOM      0  HB  VAL A 216       4.960  19.119  55.404  1.00 25.53           H   new
ATOM      0 HG11 VAL A 216       6.834  19.705  56.686  1.00 23.77           H   new
ATOM      0 HG12 VAL A 216       7.238  18.638  55.586  1.00 23.77           H   new
ATOM      0 HG13 VAL A 216       7.587  20.168  55.372  1.00 23.77           H   new
ATOM      0 HG21 VAL A 216       4.969  21.302  56.255  1.00 22.33           H   new
ATOM      0 HG22 VAL A 216       5.667  21.815  54.929  1.00 22.33           H   new
ATOM      0 HG23 VAL A 216       4.186  21.256  54.879  1.00 22.33           H   new
ATOM    634  N   VAL A 218       9.047  15.495  54.747  1.00 24.30           N
ATOM    635  CA  VAL A 218       9.980  15.425  55.854  1.00 25.59           C
ATOM    636  C   VAL A 218      11.043  14.371  55.566  1.00 25.70           C
ATOM    637  O   VAL A 218      10.781  13.407  54.869  1.00 25.17           O
ATOM    638  CB  VAL A 218       9.161  15.143  57.136  1.00 26.71           C
ATOM    639  CG1 VAL A 218       9.460  13.782  57.688  1.00 28.12           C
ATOM    640  CG2 VAL A 218       9.374  16.262  58.126  1.00 28.81           C
ATOM      0  H   VAL A 218       8.514  14.824  54.675  1.00 24.30           H   new
ATOM      0  HA  VAL A 218      10.460  16.258  55.979  1.00 25.59           H   new
ATOM      0  HB  VAL A 218       8.215  15.126  56.923  1.00 26.71           H   new
ATOM      0 HG11 VAL A 218       8.933  13.635  58.489  1.00 28.12           H   new
ATOM      0 HG12 VAL A 218       9.237  13.108  57.027  1.00 28.12           H   new
ATOM      0 HG13 VAL A 218      10.403  13.721  57.906  1.00 28.12           H   new
ATOM      0 HG21 VAL A 218       8.861  16.085  58.930  1.00 28.81           H   new
ATOM      0 HG22 VAL A 218      10.316  16.321  58.351  1.00 28.81           H   new
ATOM      0 HG23 VAL A 218       9.083  17.100  57.735  1.00 28.81           H   new
TER     728      ILE A 218
HETATM  731  O   HOH A 503       7.663  18.549  58.083  1.00 71.06           O
HETATM  758  O   HOH A 530      20.116  27.095  51.665  1.00 24.67           O
HETATM  765  O   HOH A 537      14.778  21.226  48.211  1.00 27.98           O
HETATM  772  O   HOH A 544      -5.198  25.804  53.682  1.00 21.10           O
HETATM  805  O   HOH A 577      21.563  28.198  52.938  1.00 31.45           O
HETATM  812  O   HOH A 585      -6.264  27.663  53.754  1.00 21.56           O
ENDMDL
END
"""

expected_undowser_html = """<html>
<head>
  <title>Summary table of water clashes</title>
</head>

<body>

<hr>
This table lists all HOH "waters" in the structure that have steric clashes. HOH are classified into common categories based on the atom they clash with.
<br><br>
A clashing HOH is very unlikely be be a real water, unless the clashing atom position is incorrect. The following categories provide guidance for correcting false HOH.
<br><br>
<b>Clash with polar</b> - HOH that clashes with polar groups may actually be a coordinated ion.
<br>
<b>Clash with nonpolar</b> - HOH that clashes with nonpolar groups may be a missing or displaced atom&ast;. Or it may be the first atom of an unmodeled alternate.
<br>
<b>Clash with both polar and nonpolar</b> - HOH that clashes with both polar and non-polar groups is unlikely to be an ion. If clashes are severe, a displaced atom is likely. If clashes and map are weak, the HOH may be entirely removable.
<br>
<b>Clash with water</b> - HOH-HOH clashes may be real waters that need to be modeled as alternates of compatible occupancy. Or they may be in the density of a sidechain alternate or a larger ligand.
<br>
<b>Clash with altloc</b> - HOH clashes involving one or more alternate conformations may be resolved by renaming some of the alternates.
<br><br>
<b>High B-factor</b> - HOH with clashes and minimal support in the map should be removed from the model. This table does not report map data directly, but a high B-factor is a likely warning sign that an HOH is a poor fit to the map.
<br>
<b>Severe clash</b> - HOH with severe clash overlap but good map support is likely to be a position where an atom is displaced.
<br><br>
&ast;<i>Displaced atom</i> indicates that a structural atom has been moved from its proper place in the model and replaced by HOH. Displaced sidechains are common. Moved atoms may be restored by local rebuilding.
<br>
<i>Missing atoms</i> have been entirely replaced by HOH. Removed atoms may be restored by modeling alternate conformations (especially sidechains), modeling ligands, or continuing a macromolecular mainchain.
<br><br>
These categories are general suggestions. Check your electron density; trust your intuition and experience. Prisant 2020 Prot Sci 29:315 (<a href="https://doi.org/10.1002/pro.3786">https://doi.org/10.1002/pro.3786</a>) illustrates 10 examples of clashing HOH cases.
<br>
<hr>
<br>
SUMMARY: 7 waters out of 7 have clashes (100.00%)
<br><br>
<hr>
<br>
<table border=1 width='100%'>
<tr bgcolor='#9999cc'><td rowspan='1' align='center'>Water ID</td>
<td align='center'>Clashes with</td>
<td align='center'>Water B</td>
<td align='center'>Contact B</td>
<td align='center'>Clash<br>Severity</td>
<td align='center'>Clash with Polar<br><small>May be ion</small></td>
<td align='center'>Clash with non-polar<br><small>Unmodeled alt or noise</small></td>
<td align='center'>Clash with water<br><small>Occ &lt;1 or ligand</small></td>
<td align='center'>Clash with altloc<br><small>Add or rename alts</small></td></tr>
<tr bgcolor=#eaeaea><td rowspan='2' ><pre><code>A: 572 :HOH: </code></pre></td>
<td><pre><code> O   of A: 585 :HOH: </code></pre></td><td>14.71</td><td>21.56</td><td bgcolor='#ee4d4d'>1.133</td><td></td><td></td><td align='center' bgcolor='#ee4d4d'>&times;</td><td></td></tr>
<tr bgcolor=#eaeaea><td><pre><code> O   of A: 544 :HOH: </code></pre></td><td>14.71</td><td>21.10</td><td bgcolor='#ee4d4d'>1.086</td><td></td><td></td><td align='center' bgcolor='#ee4d4d'>&times;</td><td></td></tr>
<tr bgcolor=#ffffff><td rowspan='2' ><pre><code>A: 585 :HOH: </code></pre></td>
<td><pre><code> O   of A: 572 :HOH: </code></pre></td><td>21.56</td><td>14.71</td><td bgcolor='#ee4d4d'>1.133</td><td></td><td></td><td align='center' bgcolor='#ee4d4d'>&times;</td><td></td></tr>
<tr bgcolor=#ffffff><td><pre><code> O   of A: 544 :HOH: </code></pre></td><td>21.56</td><td>21.10</td><td bgcolor='#ff76a9'>0.656</td><td></td><td></td><td align='center' bgcolor='#ff76a9'>&times;</td><td></td></tr>
<tr bgcolor=#eaeaea><td rowspan='2' ><pre><code>A: 544 :HOH: </code></pre></td>
<td><pre><code> O   of A: 572 :HOH: </code></pre></td><td>21.10</td><td>14.71</td><td bgcolor='#ee4d4d'>1.086</td><td></td><td></td><td align='center' bgcolor='#ee4d4d'>&times;</td><td></td></tr>
<tr bgcolor=#eaeaea><td><pre><code> O   of A: 585 :HOH: </code></pre></td><td>21.10</td><td>21.56</td><td bgcolor='#ff76a9'>0.656</td><td></td><td></td><td align='center' bgcolor='#ff76a9'>&times;</td><td></td></tr>
<tr bgcolor=#ffffff><td rowspan='4' ><pre><code>A: 503 :HOH: </code></pre></td>
<td><pre><code>HG11 of A: 216 :VAL: </code></pre></td><td>71.06</td><td>23.77</td><td bgcolor='#ff76a9'>0.626</td><td></td><td align='center' bgcolor='#ff76a9'>&times;</td><td></td><td></td></tr>
<tr bgcolor=#ffffff><td><pre><code>HG23 of A: 218 :VAL: </code></pre></td><td>71.06</td><td>28.81</td><td bgcolor='#ff76a9'>0.562</td><td></td><td align='center' bgcolor='#ff76a9'>&times;</td><td></td><td></td></tr>
<tr bgcolor=#ffffff><td><pre><code>HG23 of A: 146 :ILE: </code></pre></td><td>71.06</td><td>23.19</td><td bgcolor='#ffb3cc'>0.456</td><td></td><td align='center' bgcolor='#ffb3cc'>&times;</td><td></td><td></td></tr>
<tr bgcolor=#ffffff><td><pre><code> CG1 of A: 216 :VAL: </code></pre></td><td>71.06</td><td>23.77</td><td bgcolor='#ffb3cc'>0.446</td><td></td><td align='center' bgcolor='#ffb3cc'>&times;</td><td></td><td></td></tr>
<tr bgcolor=#eaeaea><td rowspan='1' ><pre><code>A: 530 :HOH: </code></pre></td>
<td><pre><code> O   of A: 577 :HOH: </code></pre></td><td>24.67</td><td>31.45</td><td bgcolor='#ff76a9'>0.579</td><td></td><td></td><td align='center' bgcolor='#ff76a9'>&times;</td><td></td></tr>
<tr bgcolor=#ffffff><td rowspan='1' ><pre><code>A: 577 :HOH: </code></pre></td>
<td><pre><code> O   of A: 530 :HOH: </code></pre></td><td>31.45</td><td>24.67</td><td bgcolor='#ff76a9'>0.579</td><td></td><td></td><td align='center' bgcolor='#ff76a9'>&times;</td><td></td></tr>
<tr bgcolor=#eaeaea><td rowspan='1' ><pre><code>A: 537 :HOH: </code></pre></td>
<td><pre><code> OE2 of A: 153 :GLU: </code></pre></td><td>27.98</td><td>36.90</td><td bgcolor='#ff76a9'>0.548</td><td align='center' bgcolor='#ff76a9'>&plus; ion</td><td></td><td></td><td></td></tr>
</table>
"""

def exercise_undowser(probe_params):
  dm = DataManager()
  #print(help(dm))
  dm.process_model_str("1",pdb_1lpl_str)
  m = dm.get_model("1")
  uz = undowser2.undowserlyze(probe_params, dm)
  undowser_html = uz.as_HTML()
  diff = difflib.unified_diff(undowser_html.splitlines(), expected_undowser_html.splitlines(), fromfile="testvalue", tofile="expectedvalue")
  changed_lines = ""
  for line in diff:
    if line.startswith("-") or line.startswith("+"):
      changed_lines = changed_lines+"\n"+line
  assert changed_lines == "", "undowser html output changed, at the following lines: "+changed_lines

def exercise_undowser_json(probe_params):
  dm = DataManager()
  #print(help(dm))
  dm.process_model_str("1",pdb_1lpl_str)
  m = dm.get_model("1")
  uz = undowser2.undowserlyze(probe_params, dm)
  uz_dict = json.loads(uz.as_JSON())
  #import pprint
  #pprint.pprint(csjson_dict)
  assert len(uz_dict['flat_results']) == 13, "tst_undowser2 json output not returning correct number of water clashes, now: "+str(len(uz_dict['flat_results']))
  assert uz_dict['flat_results'][0]["src_atom_id"] == " A 503 HOH  O   ", "tst_undowser2 json output first src_atom_id value changed, now: "+uz_dict['flat_results'][0]["src_atom_id"]
  from mmtbx.validation import test_utils
  assert test_utils.count_dict_values(uz_dict['hierarchical_results'], "water clash")==12, "tst_undowser2 json hierarchical output total number of water clashes changed, now: "+str(test_utils.count_dict_values(uz_dict['hierarchical_results'], "water clash"))
  assert test_utils.count_dict_values(uz_dict['hierarchical_results'], "nonpolar clash")==8, "tst_undowser2 json hierarchical output total number of nonpolar clashes changed, now: "+str(test_utils.count_dict_values(uz_dict['hierarchical_results'], "nonpolar clash"))
  assert test_utils.count_dict_values(uz_dict['hierarchical_results'], "polar clash")==2, "tst_undowser2 json hierarchical output total number of polar clashes changed, now: "+str(test_utils.count_dict_values(uz_dict['hierarchical_results'], "polar clash"))
  assert uz_dict['summary_results']["   1"]["num_outliers"] == 7, "tst_undowser2 json summary output total number of water clashes changed, now: "+str(uz_dict['summary_results']["   1"]["num_outliers"])
  assert uz_dict['summary_results']["   1"]["num_waters"] == 7, "tst_undowser2 json summary output total number of waters changed"

if (__name__ == "__main__"):
  if (not libtbx.env.has_module(name="probe")):
    print("Skipping exercise_undowser(): probe not configured")
    print("OK")
  else:
    parser = iotbx.cli_parser.CCTBXParser(program_class=probe2.Program, logger=null_out())
    args = [ 'approach=once' ]
    parser.parse_args(args)
    probe_params = parser.working_phil.extract()
    t0 = time.time()
    exercise_undowser(probe_params)
    exercise_undowser_json(probe_params)
    print("OK. Time: %8.3f"%(time.time()-t0))
