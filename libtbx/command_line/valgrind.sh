#! /bin/sh
valgrind --version
if [ $? -ne 0 ]; then exit 1; fi
if [ ! -n "$LIBTBX_VALGRIND" ]; then
  echo "### LIBTBX_VALGRIND not set: using default."
  echo "### To override, define LIBTBX_VALGRIND"
  echo "### before calling $LIBTBX_DISPATCHER_NAME."
  LIBTBX_VALGRIND="valgrind --tool=memcheck --suppressions=`libtbx.show_dist_paths libtbx`/valgrind-python24.supp"
  if [ $? -ne 0 ]; then exit 1; fi
  export LIBTBX_VALGRIND
fi
echo "LIBTBX_VALGRIND=$LIBTBX_VALGRIND"
if [ $? -ne 0 ]; then exit 1; fi
LIBTBX__VALGRIND_FLAG__=1
export LIBTBX__VALGRIND_FLAG__
exec "$@"
