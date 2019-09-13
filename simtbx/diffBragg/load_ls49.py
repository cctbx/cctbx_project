

def strong_spot_mask(refl_tbl, img_size):
    """note only works for strong spot reflection tables
    img_size is slow-scan, fast-scan"""
    import numpy as np
    from dials.algorithms.shoebox import MaskCode
    Nrefl = len(refl_tbl)
    masks = [refl_tbl[i]['shoebox'].mask.as_numpy_array()
             for i in range(Nrefl)]
    code = MaskCode.Foreground.real

    x1, x2, y1, y2, z1, z2 = zip(*[refl_tbl[i]['shoebox'].bbox
                                   for i in range(Nrefl)])
    spot_mask = np.zeros(img_size, bool)
    for i1, i2, j1, j2, M in zip(x1, x2, y1, y2, masks):
        slcX = slice(i1, i2, 1)
        slcY = slice(j1, j2, 1)
        spot_mask[slcY, slcX] = M & code == code
    return spot_mask


def process_ls49_image(ls49_data_dir="/Users/dermen/crystal/modules/cctbx_project/simtbx/diffBragg/LS49_sim"):
    import os
    import numpy as np
    import dxtbx
    from dxtbx.model.experiment_list import ExperimentListFactory
    from simtbx.diffBragg import utils
    from dials.array_family import flex

    os.chdir(ls49_data_dir)

    loader = dxtbx.load("step5_MPIbatch_000000.img.gz")
    img = loader.get_raw_data().as_numpy_array()

    exp_list = ExperimentListFactory.from_json_file("idx-step5_MPIbatch_000000.img_refined.expt", check_format=False)
    exp = exp_list[0]

    C = exp.crystal
    B = exp.beam
    from dxtbx.model.detector import DetectorFactory
    import json
    D = DetectorFactory.from_dict(json.load(open("method3_geom.json", "r")))  # NOTE: from Nick

    refls = flex.reflection_table.from_file("idx-step5_MPIbatch_000000.img_integrated.refl")

    bboxes = [list(refls["shoebox"][i].bbox)[:4] for i in range(len(refls))]
    mill_idx = [list(refls["miller_index"][i]) for i in range(len(refls))]

    R2 = flex.reflection_table.from_file("idx-step5_MPIbatch_000000.img_indexed.refl")
    strong_mask = strong_spot_mask(refl_tbl=R2, img_size=img.shape)

    is_bg_pixel = np.logical_not(strong_mask)

    # now fit tilting planes
    num_spots = len(refls)
    tilt_abc = np.zeros((num_spots, 3))

    for i_spot, (i1, i2, j1, j2) in enumerate(bboxes):
        shoebox_img = img[j1:j2, i1:i2]
        shoebox_mask = is_bg_pixel[j1:j2, i1:i2]

        tilt, bgmask, coeff = utils.tilting_plane(
            shoebox_img,
            mask=shoebox_mask,  # mask specifies which spots are bg pixels...
            zscore=2)

        tilt_abc[i_spot] = coeff[1], coeff[2], coeff[0]  # store as fast-scan coeff, slow-scan coeff, offset coeff

    data = np.load("LS49_data0.npz")
    spectrum = zip(data["wavelens"], data["fluxes"])
    sfall = data["sfall"][()]

    return {"dxcrystal": C, "dxdetector": D, "dxbeam": B, "mill_idx": mill_idx, "data_img": img,
            "bboxes_x1x2y1y2": bboxes, "tilt_abc": tilt_abc, "spectrum": spectrum, "sfall": sfall}

