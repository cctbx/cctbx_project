#! /bin/sh
if [ ! -n "$LIBTBX_VALGRIND" ]; then
  echo 'Error: $LIBTBX_VALGRIND is not defined.'
  exit 1
fi
valgrind --version
if [ $? -eq 0 ]; then
  LIBTBX__VALGRIND_FLAG__=1
  export LIBTBX__VALGRIND_FLAG__
  exec "$@"
fi
