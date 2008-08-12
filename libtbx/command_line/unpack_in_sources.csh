#! /bin/csh -f
set noglob
gunzip -c $* | (cd "`libtbx.show_dist_paths libtbx`/.."; tar xf -)
