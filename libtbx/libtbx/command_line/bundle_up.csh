#! /bin/csh -f
set noglob
if ($#argv < 2) then
  echo "usage: libtbx.bundle_up bundle_prefix platform_string [addl_files...]"
  exit 1
endif
set bundle="$1"
set platform="$2"
shift
shift
set echo
tar cf - "${bundle}_sources" "${bundle}_build" "${bundle}_install_script.csh" $* | gzip > "${bundle}_${platform}.tar.gz"
libtbx.create_selfx "${bundle}_${platform}.tar.gz" "./${bundle}_install_script.csh"
