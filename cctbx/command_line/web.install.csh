#! /bin/csh -f

source "`libtbx.show_build_path`/setpaths_all.csh"
set echo

libtbx.start_binary_bundle web mmtbx
rm web_install_script.csh
mv web_sources sources
mv web_build build
cd build
if (-d `libtbx.show_build_path`/base) then
  set python=base/bin/python
else
  set python=python
endif
$python ../sources/libtbx/configure.py mmtbx
cd ..

echo '#! /bin/sh -f' > cctbx_web.cgi
echo "limit vmemoryuse 384m" >> cctbx_web.cgi
echo "exec ./build/bin/libtbx.env_run CCTBX_DIST web/dispatcher.py" >> cctbx_web.cgi
chmod 755 cctbx_web.cgi

cp sources/cctbx/web/*.html .
libtbx.python sources/cctbx/web/multiple_cell.py > multiple_cell.html
