import pandas
from dials.array_family import flex
from simtbx.diffBragg import utils

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--symbol", default=None, type=str, help="spacegroup sym e.g. P43212")
parser.add_argument("--pkl",default=None, type=str, help="prediction pandas pickle" )
parser.add_argument("--mtzname", default=None, type=str, help="mtz name containing ref amps")
parser.add_argument("--mtzcol", default=None, type=str, help="mtz column")
args = parser.parse_args()

Fref = utils.open_mtz(args.mtzname, args.mtzcol)
master_indices = set(Fref.indices())
df = pandas.read_pickle(args.pkl)

for r in df.predictions:
    R = flex.reflection_table.from_file(r)
    asu = utils.map_hkl_list(R['miller_index'], anomalous_flag=True, symbol=args.symbol)
    sel = flex.bool([h in master_indices for h in asu])
    R2 = R.select(sel)
    n2 = len(R2)
    n = len(R)
    perc = n2 / n * 100.
    print("Kept %d / %d reflections %.2f %%" % (n2,n,perc))
    R2.as_file(r)
