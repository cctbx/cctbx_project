#/bin/csh -f
echo Install the Apache-mod_python-spotfinder server.

# Define the required directories
setenv ROOT_DIR `pwd`
setenv SRC_DIR ${ROOT_DIR}/project_src
setenv BUILD_DIR ${ROOT_DIR}/project_build
# USE_APACHE can point to a pre-existing Apache installation directory if desired
setenv USE_APACHE ${ROOT_DIR}/httpd
setenv USE_PYTHON ${ROOT_DIR}/python
setenv MODPY_SRC ${SRC_DIR}/mod_python-3.3.1
mkdir -p ${BUILD_DIR}

echo Untarring source files first...
cd ${SRC_DIR}
tar xf annlib*.tar
tar xf boost*.tar
tar xf cbflib*.tar
tar xf cci_apps*.tar
tar xf cctbx*.tar
tar xf scons*.tar
tar xf httpd*.tar
tar xf mod_python*.tar
tar xf Python*.tar

# Apache installation--this section can be skipped if a pre-existing Apache
# installation will be used
echo Configuring and compiling the Apache server...consult install_apache.log
setenv APACHE_SRC ${SRC_DIR}/httpd-2.2.19
cd ${APACHE_SRC}
./configure --prefix=${USE_APACHE} >& ${SRC_DIR}/install_apache.log
make >>& ${SRC_DIR}/install_apache.log
make install >>& ${SRC_DIR}/install_apache.log
cp ${USE_APACHE}/conf/httpd.conf ${USE_APACHE}/conf/httpd.conf.dist

# Python installation--this section can be skipped if a pre-existing Python 2.7 or higher
# installation will be used
echo Configuring and compiling Python...consult install_python.log
cd ${SRC_DIR}/Python-2.7.1
./configure --prefix=${USE_PYTHON} --enable-shared >& ${SRC_DIR}/install_python.log
make >>& ${SRC_DIR}/install_python.log
make install >>& ${SRC_DIR}/install_python.log

setenv PYTHON_SO ${USE_PYTHON}/lib

echo Building an env_run shell to execute a CCI-aware Apache server
cat ${ROOT_DIR}/env_run_template |sed -e "s?__python_so__?${PYTHON_SO}?g" | sed -e "s?__build_dir__?${BUILD_DIR}?g" > ${BUILD_DIR}/env_run


if ($?LD_LIBRARY_PATH) then
  setenv LD_LIBRARY_PATH ${PYTHON_SO}:${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH ${PYTHON_SO}
endif

echo Fixing up the mod python sources
mv ${MODPY_SRC}/src/connobject.c ${MODPY_SRC}/src/connobject.c.dist
cat ${MODPY_SRC}/src/connobject.c.dist|sed -e "s/!(b == APR_BRIGADE_SENTINEL(b) ||/!(b == APR_BRIGADE_SENTINEL(bb) ||/g">${MODPY_SRC}/src/connobject.c

echo Configuring and installing mod_python...consult install_mp.log
cd ${MODPY_SRC}
./configure --with-apxs=${USE_APACHE}/bin/apxs\
  --with-python=${USE_PYTHON}/bin/python\
  --with-max-locks=8 >& ${SRC_DIR}/install_mp.log
make >>& ${SRC_DIR}/install_mp.log
make install >>& ${SRC_DIR}/install_mp.log

# CCI Apps build
echo Building CCI apps...consult install_cci.log
cd ${BUILD_DIR}
# LD_LIBRARY_PATH is required for this step:
${USE_PYTHON}/bin/python ${SRC_DIR}/cctbx_project/libtbx/configure.py labelit >& ${SRC_DIR}/install_cci.log
source setpaths.csh
make >>& ${SRC_DIR}/install_cci.log

# if using a pre-existing Apache, make sure favicon goes to the right htdocs directory
cp ${ROOT_DIR}/favicon.ico ${USE_APACHE}/htdocs/
autohttp:

echo
echo The following modifications are being made AUTOMATICALLY to
echo the Apache configuration file, ${USE_APACHE}/conf/httpd.conf
echo "(but the user can make further changes if desired)."
echo
echo 1. change the port number from 80 to 8125
echo 2. include the httpd-mpm.conf multiprocessing directives
echo 3. add mod_python and the spotfinder server to the configuration, by appending the following lines:
echo "   "
echo "  LoadModule python_module ${USE_APACHE}/modules/mod_python.so"
echo "  Alias /spotfinder ${SRC_DIR}/cctbx_project/spotfinder/servers"
echo "  <Directory ${SRC_DIR}/cctbx_project/spotfinder/servers>"
echo "    Order allow,deny"
echo "    Allow from all"
echo "    AddHandler mod_python .signal_strength"
echo "    PythonHandler apache"
echo "  </Directory>"
echo
echo
echo The Apache server WITH SPOTFINDER MODIFICATIONS can be started or stopped as follows:
echo
echo "/bin/sh ${BUILD_DIR}/env_run ${USE_APACHE}/bin/apachectl [start|stop]"
echo

# if using a pre-existing Apache, the httpd.conf file must be correctly chosen
cat ${USE_APACHE}/conf/httpd.conf.dist | \
sed "s/Listen 80/Listen 8125/g" | \
sed "s?#Include conf/extra/httpd-mpm.conf?Include conf/extra/httpd-mpm.conf?g" > \
${USE_APACHE}//conf/httpd.conf

cat << eof >> ${USE_APACHE}/conf/httpd.conf
LoadModule python_module ${USE_APACHE}/modules/mod_python.so
Alias /spotfinder ${SRC_DIR}/cctbx_project/spotfinder/servers
<Directory ${SRC_DIR}/cctbx_project/spotfinder/servers>
   Order allow,deny
   Allow from all
   AddHandler mod_python .signal_strength
   PythonHandler apache
</Directory>
eof
