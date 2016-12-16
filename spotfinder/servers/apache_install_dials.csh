#!/bin/csh -f
# see README_apache_dials for instructions on how to use this file
echo Install the Apache-mod_python-spotfinder server.

# Define the required directories
setenv ROOT_DIR `pwd`
setenv SRC_DIR ${ROOT_DIR}/apache/project_src
setenv BUILD_DIR ${ROOT_DIR}/apache/project_build
setenv CCTBX_DIR ${ROOT_DIR}/cctbx
setenv CCTBX_SRC ${CCTBX_DIR}/modules/cctbx_project
setenv CCTBX_BUILD ${CCTBX_DIR}/build
# USE_APACHE can point to a pre-existing Apache installation directory if desired
setenv USE_APACHE ${ROOT_DIR}/apache/httpd
setenv USE_PYTHON ${CCTBX_DIR}/base
setenv MODPY_SRC ${SRC_DIR}/mod_python-3.5.0
mkdir -p ${BUILD_DIR}

echo Install a local copy of CCTBX with shared-Python enabled

setenv PYTHON_SO ${USE_PYTHON}/lib
if ($?LD_LIBRARY_PATH) then
  setenv LD_LIBRARY_PATH ${PYTHON_SO}:${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH ${PYTHON_SO}
endif
echo ${LD_LIBRARY_PATH}

# Build cctbx and dials using bootstrap.py
mkdir -p ${CCTBX_DIR}
cd ${CCTBX_DIR}
wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py
python bootstrap.py --builder=dials hot update

set n_processor = `python ${CCTBX_SRC}/libtbx/command_line/show_number_of_processors.py`
echo $n_processor "processors will be used for compiling "
python bootstrap.py --builder=dials --nproc=${n_processor} base build

source ${CCTBX_BUILD}/setpaths.csh

cp ${CCTBX_SRC}/spotfinder/servers/apache_favicon.ico ${ROOT_DIR}/apache/favicon.ico
cp ${CCTBX_SRC}/spotfinder/servers/apache_env_run_template ${ROOT_DIR}/apache/env_run_template

echo Untarring source files for Apache and Mod-Python
cd ${SRC_DIR}
tar xf httpd*.tar.gz
tar xf apr-[0-9]*.tar.gz
tar xf apr-util-*.tar.gz
tar xf mod_python*.tgz

# Apache installation--this section can be skipped if a pre-existing Apache
# installation will be used
setenv APACHE_SRC ${SRC_DIR}/httpd-2.4.20
setenv APR_SRC ${SRC_DIR}/apr-1.5.2
setenv APR_UTIL_SRC ${SRC_DIR}/apr-util-1.5.4

echo Configuring and compiling the Apache APR libraries...consult install_apache_apr.log
cd ${APR_SRC}
echo "./configure --prefix=${APACHE_SRC}/srclib/apr >& ${SRC_DIR}/install_apache_apr.log"
./configure --prefix=${APACHE_SRC}/srclib/apr >& ${SRC_DIR}/install_apache_apr.log
make -j $n_processor >>& ${SRC_DIR}/install_apache_apr.log
make -j $n_processor install >>& ${SRC_DIR}/install_apache_apr.log

echo Configuring and compiling the Apache APR util libraries...consult install_apache_apr_util.log
cd ${APR_UTIL_SRC}
echo "./configure --prefix=${APACHE_SRC}/srclib/apr-util --with-apr=${APACHE_SRC}/srclib/apr >& ${SRC_DIR}/install_apache_apr_util.log"
./configure --prefix=${APACHE_SRC}/srclib/apr-util --with-apr=${APACHE_SRC}/srclib/apr >>& ${SRC_DIR}/install_apache_apr_util.log
make -j $n_processor >>& ${SRC_DIR}/install_apache_apr_util.log
make -j $n_processor install >>& ${SRC_DIR}/install_apache_apr_util.log

echo Configuring and compiling the Apache server...consult install_apache.log
cd ${APACHE_SRC}
echo "./configure --prefix=${USE_APACHE} --with-mpm=prefork >& ${SRC_DIR}/install_apache.log"

./configure --prefix=${USE_APACHE} --with-mpm=prefork >& ${SRC_DIR}/install_apache.log
make -j $n_processor >>& ${SRC_DIR}/install_apache.log
make -j $n_processor install >>& ${SRC_DIR}/install_apache.log
cp ${USE_APACHE}/conf/httpd.conf ${USE_APACHE}/conf/httpd.conf.dist

echo Building an env_run shell to execute a CCTBX-aware Apache server
cat ${ROOT_DIR}/apache/env_run_template |sed -e "s?__python_so__?${PYTHON_SO}?g" | sed -e "s?__build_dir__?${CCTBX_BUILD}?g" > ${BUILD_DIR}/env_run

echo Fixing up the mod python sources
mv ${MODPY_SRC}/src/connobject.c ${MODPY_SRC}/src/connobject.c.dist
cat ${MODPY_SRC}/src/connobject.c.dist|sed -e "s/!(b == APR_BRIGADE_SENTINEL(b) ||/!(b == APR_BRIGADE_SENTINEL(bb) ||/g">${MODPY_SRC}/src/connobject.c

mv ${MODPY_SRC}/dist/version.sh ${MODPY_SRC}/dist/version.sh.dist
cat ${MODPY_SRC}/dist/version.sh.dist|sed -e "s/GIT=/#GIT=/g"|sed -e 's/$GIT//g'>${MODPY_SRC}/dist/version.sh
chmod 755 ${MODPY_SRC}/dist/version.sh

mv ${MODPY_SRC}/Makefile.in ${MODPY_SRC}/Makefile.in.dist
cat ${MODPY_SRC}/Makefile.in.dist|sed -e "s/cd scripts/#cd scripts/g">${MODPY_SRC}/Makefile.in

echo Configuring and installing mod_python...consult install_mp.log
cd ${MODPY_SRC}
./configure --with-apxs=${USE_APACHE}/bin/apxs\
  --with-python=${USE_PYTHON}/bin/python\
  --with-max-locks=8 >& ${SRC_DIR}/install_mp.log
make -j $n_processor >>& ${SRC_DIR}/install_mp.log
make -j $n_processor install >>& ${SRC_DIR}/install_mp.log


# if using a pre-existing Apache, make sure favicon goes to the right htdocs directory
cp ${ROOT_DIR}/apache/favicon.ico ${USE_APACHE}/htdocs/
autohttp:

cat << eof | tee ${ROOT_DIR}/README_customized

The following modifications are being made AUTOMATICALLY to
the Apache configuration file, ${USE_APACHE}/conf/httpd.conf
(but the user can make further changes if desired).

1. change the port number from 80 to 8125
2. include the httpd-mpm.conf multiprocessing directives
3. add mod_python and the spotfinder server to the configuration, by appending the following lines:

   LoadModule python_module ${USE_APACHE}/modules/mod_python.so
   Alias /spotfinder ${CCTBX_SRC}/spotfinder/servers
   <Directory ${CCTBX_SRC}/spotfinder/servers>
     Require all granted
     AddHandler mod_python .signal_strength
     AddHandler mod_python .find_spots
     PythonHandler apache
   </Directory>
   <Directory ${CCTBX_SRC}/spotfinder/servers>
     Require all granted
     AddHandler mod_python .signal_strength_bcsb
     PythonHandler apache
   </Directory>


The Apache server WITH SPOTFINDER MODIFICATIONS can be started or stopped as follows:

/bin/sh ${BUILD_DIR}/env_run ${USE_APACHE}/bin/apachectl [start|stop]

An example client (single-process; multithreaded) is shown in the file
  ${CCTBX_SRC}/spotfinder/servers/general_client_example.py

After editing this script to point to your own dataset, run the client as follows:

/bin/sh ${BUILD_DIR}/env_run libtbx.python ${CCTBX_SRC}/spotfinder/servers/general_client_example.py

eof

# if using a pre-existing Apache, the httpd.conf file must be correctly chosen
cat ${USE_APACHE}/conf/httpd.conf.dist | \
sed "s/Listen 80/Listen 8125/g" | \
sed "s?#Include conf/extra/httpd-mpm.conf?Include conf/extra/httpd-mpm.conf?g" > \
${USE_APACHE}//conf/httpd.conf

cat << eof >> ${USE_APACHE}/conf/httpd.conf
LoadModule python_module ${USE_APACHE}/modules/mod_python.so
Alias /spotfinder ${CCTBX_SRC}/spotfinder/servers
<Directory ${CCTBX_SRC}/spotfinder/servers>
   Require all granted
   AddHandler mod_python .signal_strength
   AddHandler mod_python .find_spots
   PythonHandler apache
</Directory>
<Directory ${CCTBX_SRC}/spotfinder/servers>
   Require all granted
   AddHandler mod_python .signal_strength_bcsb
   PythonHandler apache
</Directory>
eof
