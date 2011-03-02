#! /bin/csh -fe
set omz="`libtbx.show_dist_paths cctbx`/omz"
set verbose
rm cod_ma_xs/2017318.pickle
iotbx.python "$omz"/cod_select_and_pickle.py 2017318
iotbx.python "$omz"/cod_refine.py optimizers=dev+shelxl_fm+shelxl_cg+shelx76 shake_sites_rms=None reset_u_iso=None cod_ma_xs/2017318.pickle
