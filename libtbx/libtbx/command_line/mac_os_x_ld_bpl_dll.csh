#! /bin/csh -f
set noglob
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
ld -dynamic -m -r -d -o libboost_python.lo $*
g++ -nostartfiles -Wl,-dylib -ldylib1.o -framework "$framework" -o $dll libboost_python.lo
rm -f libboost_python.lo
