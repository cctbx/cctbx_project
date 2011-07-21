#!/bin/csh -f
echo Install the Apache-mod_python-spotfinder server.

# Define the required directories
setenv ROOT_DIR `pwd`
setenv SRC_DIR ${ROOT_DIR}/apache/project_src
setenv BUILD_DIR ${ROOT_DIR}/apache/project_build
setenv CCTBX_SRC ${ROOT_DIR}/cctbx_sources
setenv CCTBX_BUILD ${ROOT_DIR}/cctbx_build
# USE_APACHE can point to a pre-existing Apache installation directory if desired
setenv USE_APACHE ${ROOT_DIR}/apache/httpd
setenv USE_PYTHON ${CCTBX_BUILD}/base
setenv MODPY_SRC ${SRC_DIR}/mod_python-3.3.1
mkdir -p ${BUILD_DIR}

echo Install a local copy of CCTBX with shared-Python enabled

setenv PYTHON_SO ${USE_PYTHON}/lib
if ($?LD_LIBRARY_PATH) then
  setenv LD_LIBRARY_PATH ${PYTHON_SO}:${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH ${PYTHON_SO}
endif
echo ${LD_LIBRARY_PATH}
mv cctbx_install_script.csh cctbx_install_script.dist
cat cctbx_install_script.dist | sed -e "s/rstbx/rstbx spotfinder/g" | sed -e "s/--prefix=/--enable-shared --prefix=/g" > cctbx_install_script.csh
chmod a+x cctbx_install_script.csh
./cctbx_install_script.csh
source ${CCTBX_BUILD}/setpaths.csh

set n_processor = `libtbx.show_number_of_processors`
echo $n_processor "processors will be used for compiling Apache and Mod-Python"

echo Untarring source files for Apache and Mod-Python
cd ${SRC_DIR}
tar xf httpd*.tar
tar xf mod_python*.tar

# Apache installation--this section can be skipped if a pre-existing Apache
# installation will be used
echo Configuring and compiling the Apache server...consult install_apache.log
setenv APACHE_SRC ${SRC_DIR}/httpd-2.2.19
cd ${APACHE_SRC}
echo "./configure --prefix=${USE_APACHE} >& ${SRC_DIR}/install_apache.log"

./configure --prefix=${USE_APACHE} >& ${SRC_DIR}/install_apache.log
make -j $n_processor >>& ${SRC_DIR}/install_apache.log
make -j $n_processor install >>& ${SRC_DIR}/install_apache.log
cp ${USE_APACHE}/conf/httpd.conf ${USE_APACHE}/conf/httpd.conf.dist


echo Building an env_run shell to execute a CCTBX-aware Apache server
cat ${ROOT_DIR}/apache/env_run_template |sed -e "s?__python_so__?${PYTHON_SO}?g" | sed -e "s?__build_dir__?${CCTBX_BUILD}?g" > ${BUILD_DIR}/env_run

echo Fixing up the mod python sources
mv ${MODPY_SRC}/src/connobject.c ${MODPY_SRC}/src/connobject.c.dist
cat ${MODPY_SRC}/src/connobject.c.dist|sed -e "s/!(b == APR_BRIGADE_SENTINEL(b) ||/!(b == APR_BRIGADE_SENTINEL(bb) ||/g">${MODPY_SRC}/src/connobject.c

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
     Order allow,deny
     Allow from all
     AddHandler mod_python .signal_strength
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
   Order allow,deny
   Allow from all
   AddHandler mod_python .signal_strength
   PythonHandler apache
</Directory>
eof
