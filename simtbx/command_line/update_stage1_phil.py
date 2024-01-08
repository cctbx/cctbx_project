
from __future__ import division
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("dirname", help="stage1 output folder with pandas and diff.phil", type=str)
parser.add_argument("newPhil", default=None, help="new stage 1 phil file", type=str)
parser.add_argument("--setGlim", action="store_true", help="Use unrestrained stage1 to set bounds on G")
#parser.add_argument("--njobs", type=int, default=5, help="number of jobs (only runs on single node, no MPI)")
#parser.add_argument("--plot", action="store_true", help="show a histogram at the end")
args = parser.parse_args()

# LIBTBX_SET_DISPATCHER_NAME diffBragg.update_stage1_phil

import os
import pandas
import numpy as np
import glob

glob_s = os.path.join(args.dirname, "pandas/*pkl")
fnames = glob.glob(glob_s)

# read in the pandas pickles
df1 = pandas.concat( [pandas.read_pickle(f) for f in fnames]).reset_index(drop=True)

# unit cell phil
a,b,c,al,be,ga = df1[['a', 'b', 'c', 'al', 'be', 'ga']].median()

# Ncells abc and Nvol
na, nb, nc = np.vstack(df1.ncells).T
nvol = na*nb*nc
nvol = np.median(na*nb*nc)
na, nb, nc = map(np.median, (na, nb, nc))

# eta
ea, eb, ec = map( np.median, np.vstack(df1.eta_abc).T)

# spot scale (G)
Gmed = df1.spot_scales.median()
Gmin = df1.spot_scales.min()/100
Gmax = df1.spot_scales.max()*100

update_phil = """
init {{
  G = {G}
  Nabc = [{na},{nb},{nc}]
  eta_abc = [{ea},{eb},{ec}]
}}
centers {{
  Nvol = {nvol}
  ucell_a = {a}
  ucell_b = {b}
  ucell_c = {c}
  ucell_alpha = {al}
  ucell_beta = {be}
  ucell_gamma = {ga}
}}
betas {{
  Nvol = 1e-2
  ucell_a = 1e-7
  ucell_b = 1e-7
  ucell_c = 1e-7
  ucell_alpha = 1e-7
  ucell_beta = 1e-7
  ucell_gamma = 1e-7
}}
use_restraints = True
""".format(G=Gmed,na=na, nb=nb, nc=nb, ea=ea, eb=eb, ec=ec,a=a,b=b,c=c,al=al,be=be,ga=ga, nvol=nvol)

Gmin_Gmax="""
mins.G={Gmin}
maxs.G={Gmax}\n""".format(Gmin=Gmin, Gmax=Gmax)

if args.setGlim:
    update_phil += Gmin_Gmax

diff_phil_name = os.path.join(args.dirname, "diff.phil")
assert os.path.exists(diff_phil_name)
diff_phil = open(diff_phil_name, "r").read()
with open(args.newPhil, "w") as o:
    o.write(diff_phil + "\n"+ update_phil)

print("Done. added \n%s\n to the %s and saved to %s" % (update_phil, diff_phil_name, args.newPhil) )
