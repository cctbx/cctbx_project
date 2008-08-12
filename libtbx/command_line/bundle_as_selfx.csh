#! /bin/csh -f
set noglob
if ($#argv < 2) then
  echo "usage: libtbx.bundle_as_selfx bundle_prefix platform_string [addl_files...]"
  exit 1
endif
set bundle="$1"
set platform="$2"
shift
shift
set echo
if ($?LIBTBX_NATIVE_TAR) then
  set tar_cmd="$LIBTBX_NATIVE_TAR"
else
  set tar_cmd=tar
endif
"$tar_cmd" cf - "${bundle}_sources" "${bundle}_build" "${bundle}_install_script.csh" $* | gzip > "${bundle}_${platform}.tar.gz"
libtbx.create_selfx "${bundle}_${platform}.tar.gz" "./${bundle}_install_script.csh"
