#! /bin/csh -f

set echo

mkdir dist
cd dist
mkdir boost
mkdir ccp4
cp -r $LIBTBX_DIST .
mkdir scitbx
cp -r $SCITBX_DIST/{scitbx,libtbx_config} scitbx
mkdir cctbx
cp -r $CCTBX_DIST/{cctbx,libtbx_config} cctbx
mkdir iotbx
cp -r $IOTBX_DIST/{iotbx,libtbx_config} iotbx
cd ..

if ($#argv == 0) then
  mkdir bintbx
  cd bintbx
  cp -r $LIBTBX_BUILD/libtbx .
  python ../dist/libtbx/configure.py iotbx
  cd ..
endif

echo '#! /bin/sh -f' > cctbx_web.cgi
echo "./bintbx/env_run.sh CCTBX_DIST cctbx/web/dispatcher.py" >> cctbx_web.cgi
chmod 755 cctbx_web.cgi

cp dist/cctbx/cctbx/web/*.html .
python dist/cctbx/cctbx/web/multiple_cell.py > multiple_cell.html
