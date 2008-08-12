#! /bin/csh -fe
source "`libtbx.show_build_path`/setpaths_all.csh"
set verbose
#
set tutorial_dir="$IOTBX_DIST/examples/pdb_truncate_to_ala"
cp -r "$tutorial_dir"/tutorial.txt .
docutils.rst2html -stg tutorial.txt > tutorial.html
cp -r "$tutorial_dir"/v*.py .
cp -r "$tutorial_dir"/*.pdb .
foreach vx (v*.py)
  iotbx.python "$vx" crambin_pieces.pdb > "${vx:r}_crambin_pieces.out"
  iotbx.python "$vx" resname_mix.pdb    > "${vx:r}_resname_mix.out"
end
