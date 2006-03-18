#! /bin/csh -fe
set verbose
aCC $* -AA -mt -c -DNDEBUG -O -I. -I./Include -DPy_BUILD_CORE -o python_cpp.o python_cpp.cpp
aCC $* -AA -mt -Wl,-E -Wl,+s -o python python_cpp.o libpython2.4.a -lnsl -lrt -ldld -ldl -lpthread -lm
