from __future__ import division, print_function

phil_str = """
include scope simtbx.command_line.hopper.phil_scope
"""
from libtbx.phil import parse
master_phil = parse(phil_str, process_includes=True)

import sys
sources = [parse(open(sys.argv[1], "r").read())]
extra_params = parse("""
exp_ref_spec_file=gath_753_out.txt
niter_per_J=1
niter=0
load_data_from_refls=True
quiet=True
refiner.num_devices=1
""")
sources.append(extra_params)

#if len(sys.argv) > 2:
#    sources.append(parse("\n".join([l.strip() for l in sys.argv[2:]])))

working_phil = master_phil.fetch(sources=sources) # [file_params, extra_params, more_params])
params = working_phil.extract()

from simtbx.diffBragg import utils
from simtbx.diffBragg import hopper_utils
import numpy as np
refinement_scores = []
from libtbx.mpi4py import  MPI
COMM= MPI.COMM_WORLD
exp_ref_spec = utils.read_exp_ref_spec(params.exp_ref_spec_file)
all_new_dists = []
all_del_x = []
all_del_y = []
nexp = len(exp_ref_spec)
for i_exp, (exp_name, ref_name, spec_name) in enumerate(exp_ref_spec):
    if i_exp % COMM.size != COMM.rank:
        continue
    #if COMM.rank==0:
    #print("Processing exp %s (%d /%d" %(exp_name, i_exp+1, len(exp_ref_spec)))
    from IPython import embed
    embed()
    # NOTE refine script accepts either string names or exp and ref objects
    dev_id = COMM.rank % params.refiner.num_devices
    if params.refiner.randomize_devices:
        dev_id = np.random.choice(params.refiner.num_devices)
    new_exp, new_ref = hopper_utils.refine(exp_name, ref_name, params, spec_name, gpu_device=dev_id)

    old_resid = np.vstack(new_ref["dials.xyzcal.px"] - new_ref["xyzobs.px.value"])[:,:2]
    new_resid = np.vstack(new_ref["xyzcal.px"] - new_ref["xyzobs.px.value"])[:,:2]

    old_dists = np.linalg.norm(old_resid, axis=1)
    new_dists = np.linalg.norm(new_resid, axis=1)
    del_x, del_y = new_resid.T
    all_del_x += list(del_x)
    all_del_y += list(del_y)

    old_d = np.median(old_dists)
    new_d = np.median(new_dists)
    delta_dist = old_d - new_d
    print("(%d/%d) Previous prediction offsets: %f pixels; After refinement: %f pixels" \
          % (i_exp+1, nexp, old_d, new_d), flush=True)
    refinement_scores.append(delta_dist)  # new dists should be smaller than old dists, hence high scores are good

    all_new_dists += list(new_dists)

all_new_dists = COMM.reduce(all_new_dists)
all_del_x = COMM.reduce(all_del_x)
all_del_y = COMM.reduce(all_del_y)
refinement_scores = np.array(COMM.reduce(refinement_scores))
if COMM.rank == 0:
    print("Average score=%f pixels" % np.mean(refinement_scores) )
    print("Num improved = %d / %d" % \
          ( sum(refinement_scores >0), len(refinement_scores)))
    print("Med dist = %f pixels" % np.median(all_new_dists))
    print("Med del x = %f pixels" % np.median(all_del_x))
    print("Med del y = %f pixels" % np.median(all_del_y))

    np.savez("DID", all_d=all_new_dists, scores=refinement_scores)
    #frac_improved = np.sum(refinement_scores >0) / float(len(refinement_scores))
    #assert np.median(all_new_dists) <0.38
    #assert frac_improved > 0.84
    print("OK!")
