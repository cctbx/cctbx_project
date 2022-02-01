from __future__ import division

from simtbx.diffBragg import utils
from simtbx.nanoBragg.tst_nanoBragg_basic import pdb_lines
from simtbx.diffBragg.phil import hopper_phil, philz
from libtbx.phil import parse
from simtbx.nanoBragg import tst_nanoBragg_multipanel
from dxtbx.model import Experiment
import numpy as np

# make a dummie experiment
expt = Experiment()
expt.detector = tst_nanoBragg_multipanel.whole_det
expt.beam = tst_nanoBragg_multipanel.beam
expt.crystal = tst_nanoBragg_multipanel.cryst

# write a dummie PDB file
PDB = "1234.pdb"
o = open(PDB, "w")
o.write(pdb_lines)
o.close()

# Create a miller array from on-disk PDB file
F = utils.get_complex_fcalc_from_pdb(PDB,
    wavelength=expt.beam.get_wavelength(),
    dmin=2, dmax=20, k_sol=0.2, b_sol=20)
F = F.as_amplitude_array()
Fmap = {h: amp for h,amp in zip(F.indices(), F.data())}

# Create a sim_data class instance as would be done for hopper_utils.refine for example
phil_scope = parse(hopper_phil+philz)
params = phil_scope.extract()
params.simulator.structure_factors.from_pdb.name = PDB
params.simulator.structure_factors.from_pdb.add_anom = True
params.simulator.structure_factors.from_pdb.k_sol = 0.2
params.simulator.structure_factors.from_pdb.b_sol = 20
params.simulator.structure_factors.dmin = 2
params.simulator.structure_factors.dmax = 20


SIM = utils.simulator_from_expt_and_params(expt, params)

# verify the miller arrary
F_from_SIM = SIM.crystal.miller_array
F_from_SIM_map = {h:amp for h,amp in zip(F_from_SIM.indices(), F_from_SIM.data())}
for h in Fmap:
    assert Fmap[h] == F_from_SIM_map[h]

# verify the miller data within the diffBragg object
db_instance = SIM.D
db_inds, db_amps = db_instance.Fhkl_tuple
db_Fmap = {h: amp for h, amp in zip(db_inds, db_amps)}
for h in Fmap:
    assert Fmap[h] == db_Fmap[h]


# create a default array, and verify it
params.simulator.structure_factors.from_pdb.name = None
params.simulator.structure_factors.default_F = 1234
SIM_def = utils.simulator_from_expt_and_params(expt, params)

inds,amps = SIM_def.D.Fhkl_tuple
for hkl,amp in zip(inds, amps):
    if hkl== (0,0,0):
        assert amp == 0
    else:
        assert amp==1234

# save as F as an MTZ file
MTZ="1234.mtz"
col="Famp"
F.as_mtz_dataset(column_root_label="Famp").mtz_object().write(MTZ)
# create a sim data object with structure factors loaded from MTZ file
params.simulator.structure_factors.default_F = 0
params.simulator.structure_factors.mtz_name = MTZ
params.simulator.structure_factors.mtz_column = "%s(+),%s(-)" % (col,col)
SIM_mtz = utils.simulator_from_expt_and_params(expt, params)
# verify
F_from_SIM2 = SIM_mtz.crystal.miller_array
F_from_SIM2_map = {h:amp for h,amp in zip(F_from_SIM2.indices(), F_from_SIM2.data())}

avals = [Fmap[h] for h in Fmap]
avals2 = [F_from_SIM2_map[h] for h in Fmap]
assert np.allclose(avals, avals2)
# with MTZ writing there is a loss in precision, hence equality only holds up to float32
assert np.all(np.float32(avals)==np.float32(avals2))


print("OK")
