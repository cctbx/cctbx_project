#! /bin/csh -fe
tail -17 /net/ribbon/scratch1/rwgk/bintbx_py264/sf_times_plog
tail -17 /net/krait/scratch1/rwgk/bintbx_py264/sf_times_plog
tail -17 /net/longnose/scratch2/rwgk/bintbx_py264/sf_times_plog
tail -17 /net/anaconda/scratch1/rwgk/bintbx_py264/sf_times_plog
ssh sharptail 'tail -17 /net/sharptail/scratch1/rwgk/bintbx_py264/sf_times_plog'
tail -17 /net/chevy/raid1/rwgk/phenix-1.4-135/bintbx/91/sf_times_plog
tail -17 /net/chevy/raid1/rwgk/phenix-1.4-135/bintbx/101/sf_times_plog
tail -17 /net/chevy/raid1/rwgk/phenix-1.4-135/bintbx/111/sf_times_plog
tail -17 /net/chevy/raid1/rwgk/bintbx_py264_gcc422/sf_times_plog
tail -17 /net/chevy/raid1/rwgk/bintbx_py264_gcc434/sf_times_plog
tail -17 /net/chevy/raid1/rwgk/bintbx_py264_gcc442/sf_times_plog
ssh firtree 'tail -17 /net/firtree/raid1/rwgk/bintbx_py264_no_omp/sf_times_plog'

tail -17 /net/krait/scratch1/rwgk/bintbx_py264/sf_times_rlog
tail -17 /net/longnose/scratch2/rwgk/bintbx_py264/sf_times_rlog
tail -17 /net/anaconda/scratch1/rwgk/bintbx_py264/sf_times_rlog
ssh sharptail 'tail -17 /net/sharptail/scratch1/rwgk/bintbx_py264/sf_times_rlog'
tail -17 /net/rosie/scratch1/rwgk/time/sf_times_rlog
tail -17 /net/chevy/raid1/rwgk/phenix-1.4-135/bintbx/91/sf_times_rlog

tail -17 /net/chevy/raid1/rwgk/phenix-1.4-135/bintbx/91/sf_times_krait_build
tail -17 /net/chevy/raid1/rwgk/phenix-1.4-135/bintbx/91/sf_times_longnose_build
tail -17 /net/chevy/raid1/rwgk/phenix-1.4-135/bintbx/91/sf_times_anaconda_build
tail -17 /net/chevy/raid1/rwgk/phenix-1.4-135/bintbx/91/sf_times_shartail_build
