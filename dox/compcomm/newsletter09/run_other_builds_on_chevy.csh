#! /bin/csh -fe
set verbose
alias cmd "cctbx.python $cx/compcomm/newsletter09/sf_times.py"
cmd /net/krait/scratch1/rwgk/bintbx_py264/a_out_archive.pickle >& sf_times_krait_build
cmd /net/longnose/scratch2/rwgk/bintbx_py264/a_out_archive.pickle >& sf_times_longnose_build
cmd /net/anaconda/scratch1/rwgk/bintbx_py264/a_out_archive.pickle >& sf_times_anaconda_build
scp sharptail:/net/sharptail/scratch1/rwgk/bintbx_py264/a_out_archive.pickle sharptail_a_out_archive.pickle
cmd sharptail_a_out_archive.pickle >& sf_times_shartail_build
