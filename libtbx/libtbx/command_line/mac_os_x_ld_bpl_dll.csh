#! /bin/csh -f
shift # discard "-o"
set dll=$1
shift
ld -dynamic -m -r -d -o libboost_python.lo $*
g++ -Wl,-dynamic -nostartfiles -Wl,-dylib -Wl,-ldylib1.o -Wl,-dylib_compatibility_version,1.31.0 -Wl,-dylib_current_version,1.31.0 -o $dll libboost_python.lo -L/Library/Frameworks/Python.framework/Versions/2.3 -framework /Library/Frameworks/Python.framework/Versions/2.3/Python
rm -f libboost_python.lo
