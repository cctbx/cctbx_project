#! /bin/csh -f

source "`libtbx.show_build_path`/setpaths_all.csh"
set echo

mkdir -p sources
cd sources
cp -r "$LIBTBX_DIST" .
cp -r "$CCP4IO_DIST" .
mkdir -p ccp4io_adaptbx
mkdir -p boost_adaptbx
cp -r "$BOOST_ADAPTBX_DIST"/boost boost_adaptbx
mkdir -p scitbx
cp -r "$SCITBX_DIST"/{scitbx,libtbx_config} scitbx
mkdir -p cctbx
cp -r "$CCTBX_DIST"/{cctbx,libtbx_config} cctbx
mkdir -p iotbx
cp -r "$IOTBX_DIST"/{iotbx,libtbx_config} iotbx
mkdir -p mmtbx
cp -r "$MMTBX_DIST"/{mmtbx,libtbx_config} mmtbx
cd ..

if ($#argv == 0) then
  mkdir -p build
  cd build
  cp -r `libtbx.show_build_path`/lib .
  if (-d `libtbx.show_build_path`/base) then
    cp -r `libtbx.show_build_path`/base .
    set python=base/bin/python
  else
    set python=python
  endif
  $python ../sources/libtbx/configure.py mmtbx
  cd ..
endif

echo '#! /bin/sh -f' > cctbx_web.cgi
echo "exec ./build/bin/libtbx.env_run CCTBX_DIST cctbx/web/dispatcher.py" >> cctbx_web.cgi
chmod 755 cctbx_web.cgi

cp sources/cctbx/cctbx/web/*.html .
libtbx.python sources/cctbx/cctbx/web/multiple_cell.py > multiple_cell.html
