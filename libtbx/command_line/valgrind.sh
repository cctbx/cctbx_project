#! /bin/sh
verbose="yes"
if [ -n "$LIBTBX_VALGRIND_TOOL" ]; then
  tool="$LIBTBX_VALGRIND_TOOL"
else
  tool="memcheck"
fi
if [ -n "$LIBTBX_VALGRIND_OPTS" ]; then
  for opt in $LIBTBX_VALGRIND_OPTS; do
    if [ "$opt" == "-q" ]; then
      verbose="no"
      break
    fi
  done
fi
if [ "$verbose" == "yes" ]; then
  valgrind --version
else
  valgrind --version > /dev/null 2>&1
fi
if [ $? -ne 0 ]; then exit 1; fi
if [ ! -n "$LIBTBX_VALGRIND" ]; then
  if [ "$verbose" == "yes" ]; then
    echo "### LIBTBX_VALGRIND not set: using default."
    echo "### To override, define LIBTBX_VALGRIND"
    echo "### before calling $LIBTBX_DISPATCHER_NAME."
    if [ "`uname`" = Darwin ]; then
      echo "### Hint: For debug builds insert"
      echo "###           env LIBTBX_VALGRIND_OPTS=--dsymutil=yes"
      echo "###       before the current command."
    fi
    echo "### Hint: To generate suppressions for new environments insert"
    echo "###           env LIBTBX_VALGRIND_OPTS=--gen-suppressions=yes"
    echo "###       before the current command."
  fi
  opt="$LIBTBX_VALGRIND_OPTS"
  if [ "`uname`" = Darwin ]; then
    opt=" $opt --trace-children=yes"
  fi
  LIBTBX_VALGRIND="valgrind --tool=$tool$opt --suppressions=`libtbx.show_dist_paths libtbx`/valgrind-python-cci.supp"
  if [ $? -ne 0 ]; then exit 1; fi
  export LIBTBX_VALGRIND
fi
if [ "$verbose" == "yes" ]; then
  echo "LIBTBX_VALGRIND=$LIBTBX_VALGRIND"
fi
if [ $? -ne 0 ]; then exit 1; fi
LIBTBX__VALGRIND_FLAG__=1
export LIBTBX__VALGRIND_FLAG__
exec "$@"
