#! /bin/csh -f
set noglob
set cpp="$1"
shift
if ("$1" != "-framework") then
  echo "Expected -framework. Bailing out."
  exit 1
endif
shift
set framework="$1"
shift
if ("$1" != "-o") then
  echo "Expected -o. Bailing out."
  exit 1
endif
shift
set dll="$1"
shift
set echo
ld -dynamic -m -r -d -bind_at_load -o libboost_python.lo $*
"$cpp" -nostartfiles -Wl,-dylib -ldylib1.o -framework "$framework" -o $dll libboost_python.lo
rm -f libboost_python.lo
