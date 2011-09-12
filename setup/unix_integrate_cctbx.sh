#!/bin/sh -e

MYBASENAME="$(basename "$0")"
MYVERSION='$Id$'

MYDOC="usage: ${MYBASENAME} [options] actions

supported actions are:  all pylibs libs includes scripts

Options:

  --prefix          installation prefix directory [/usr/local]
  --libdir          directory for dynamic load libraries [prefix/lib]
  --skip-bp         do not install the cctbx boost_python dynamic library,
                    rely on the system library instead.
  --includedir      directory for C and C++ header files [prefix/include]
  --bindir          directory for executable files [prefix/bin]
  --pythondir       installation directory for Python .pth files
                    [prefix/lib/pythonX.Y/site-packages]
  -y, --yes         install without verifying the target paths
  -h, --help        display this message and exit
  -V, --version     display the script version and exit"

# Parse command line arguments -----------------------------------------------

PARSEDARGS="$(getopt -o yhV \
    --long prefix:,libdir:,skip-bp,includedir:,bindir:,pythondir:,yes,help,version \
    -n "${MYBASENAME}" -- "$@")" || exit 2
eval set -- "${PARSEDARGS}"

# options
while true; do
    case "$1" in
        --prefix)
            opt_prefix="$2"; shift 2;;
        --libdir)
            opt_libdir="$2"; shift 2;;
        --skip-bp)
            opt_skip_bp=1; shift;;
        --includedir)
            opt_includedir="$2"; shift 2;;
        --bindir)
            opt_bindir="$2"; shift 2;;
        --pythondir)
            opt_pythondir="$2"; shift 2;;
        -y|--yes)
            opt_yes=1; shift;;
        -h|--help)
            echo "${MYDOC}"; exit;;
        -V|--version)
            echo "${MYVERSION}"; exit;;
        --)
            shift; break;;
        *)
            echo "Internal error!"; exit 2;;
    esac
done

# other arguments which specify actions
if [ $# = 0 ]; then
    echo "No action specified."
    echo "${MYDOC}"
    exit
elif [ $# = 1 ] && [ "x$1" = xall ]; then
    set -- pylibs libs includes scripts
fi

# Resolve option values ------------------------------------------------------

if ! PYTHON_EXE="$(libtbx.show_python_sys_executable)"; then
    echo "cctbx scripts are not in the PATH!" >&2 ; exit 2
fi
PYTHON_VERSION="$("${PYTHON_EXE}" -c \
    'import sys; print ".".join(sys.version.split(".")[:2])')"

# make sure the libtbx bin directory is first in the path no matter what
PATH="$(libtbx.show_bin_path):${PATH}"

# resolve actual values of the variables
act_prefix="${opt_prefix:-/usr/local}"
act_libdir="${opt_libdir:-${act_prefix}/lib}"
act_includedir="${opt_includedir:-${act_prefix}/include}"
act_bindir="${opt_bindir:-${act_prefix}/bin}"

# few more useful variables
my_libtbx_build="$(libtbx.show_build_path)"
my_libtbx_binpath="$(libtbx.show_bin_path)"
my_libtbx_libpath="$(libtbx.show_lib_path)"
my_libtbx_repository="$(libtbx.show_repository_paths | tail -1)"

# use Python distutils to resolve the default pythondir
if [ x = "x${opt_pythondir}" ]; then
    act_pythondir="$( export opt_prefix; "${PYTHON_EXE}" -c "
import os
from distutils.core import Distribution
cmd = Distribution().get_command_obj('install')
cmd.prefix = os.environ.get('opt_prefix')
cmd.finalize_options()
print cmd.install_lib"
    )"
else
    act_pythondir="${opt_pythondir}"
fi

# shared function for the prompt to proceed
prompt_proceed() {
    if [ "x${opt_yes}" = x1 ]; then
        return 0
    fi
    echo "Proceed? (yes/no/quit)"
    local ans
    read ans
    case "${ans}" in
        yes|y)  return 0 ;;
        quit|q) exit ;;
    esac
    return 1
}

# do4930_ is a standard prefix for action functions, which should not
# cause any name conflicts.

# pylibs -- install the cctbx.pth file ---------------------------------------

do4930_pylibs() {
    echo "Installing cctbx.pth to ${act_pythondir}"
    prompt_proceed || return 0
    ( echo "import os; os.environ.setdefault('LIBTBX_BUILD', '${my_libtbx_build}')"
      libtbx.show_pythonpath | tr ':' '\n'
    ) > "${act_pythondir}/cctbx.pth" || false
    chmod 644 "${act_pythondir}/cctbx.pth"
}

# libs -- install dynamic load libraries -------------------------------------

do4930_libs() {
    echo "Installing dynamic load libraries to ${act_libdir}"
    prompt_proceed || return 0
    local bprule
    if [ "x${opt_skip_bp}" = x1 ]; then
        bprule="-not -name libboost_python\\*.so"
    fi
    find "${my_libtbx_libpath}" -maxdepth 1 -type f -name 'lib*.so' ${bprule} \
        -exec ln -sf {} "${act_libdir}/" \;
}

# includes -- install symlinks to the header files ---------------------------

rsync_headers() {
    rsync -rl --chmod=a+r,go-w,Da+x,Fa-x --prune-empty-dirs \
        --include='*/' --include='*.h' --include='*.hpp' --exclude='*' \
        "$@"
}

do4930_includes() {
    echo "Installing C and C++ header files to ${act_includedir}"
    prompt_proceed || return 0
    test -d "${act_includedir}" || mkdir -v "${act_includedir}"
    # we need to copy the files to merge the source and build include trees
    ( cd "${my_libtbx_repository}"
      rsync_headers annlib/include/ANN "${act_includedir}"/
      rsync_headers annlib_adaptbx/include/annlib_adaptbx "${act_includedir}"/
      # skipped antlr3 boost
      rsync_headers boost_adaptbx "${act_includedir}"/
      # skipped cbflib cbflib_adaptbx ccp4io ccp4io_adaptbx
      rsync_headers cctbx chiltbx "${act_includedir}"/
      # skipped clipper clipper_adaptbx
      rsync_headers fable/fem.hpp fable/fem "${act_includedir}"/
      # skipped gltbx
      rsync_headers \
          iotbx mmtbx omptbx rstbx scitbx smtbx spotfinder tbxx \
          "${act_includedir}"/
      # skipped tntbx
    )
    rsync_headers \
        "${my_libtbx_build}"/include/ \
        "${my_libtbx_build}"/annlib_adaptbx/include/ \
            "${act_includedir}"/
}

# scripts -- install symbolic links to the cctbx scripts ---------------------

do4930_scripts() {
    echo "Installing cctbx scripts to ${act_bindir}"
    prompt_proceed || return 0
    for f in "${my_libtbx_binpath}"/*; do
        if [ -x "$f" ] && [ -f "$f" ]; then
            ln -sf "$f" "${act_bindir}"/
        fi
    done
}

# execute all the specified actions here:
for a; do
    do4930_${a}
done
