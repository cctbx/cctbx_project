#! /bin/csh -f

set echo

mkdir dist
cd dist
cp -r $LIBTBX_DIST .
mkdir scitbx
cp -r $SCITBX_DIST/scitbx scitbx
mkdir cctbx
cp -r $CCTBX_DIST/cctbx cctbx
cd ..

mkdir build
cd build
cp -r $LIBTBX_BUILD/libtbx .
unsetenv LIBTBX_SCONS
python ../dist/libtbx/configure.py scitbx cctbx
cd ..

echo '#! /bin/sh -f' > cctbx_web.cgi
echo "./build/env_run.sh CCTBX_DIST cctbx/web/dispatcher.py" >> cctbx_web.cgi
chmod 755 cctbx_web.cgi

cp dist/cctbx/cctbx/web/*.html .
