MAKEFLAGS=--no-print-directory

usage:
	@echo 'usage:'
	@echo 'make ...'

softlinks:
	cd eltbx; make softlinks
	cd sgtbx; make softlinks
	cd uctbx; make softlinks
	cd examples/cpp; make softlinks

unlink:
	cd eltbx; make unlink
	cd sgtbx; make unlink
	cd uctbx; make unlink
	cd examples/cpp; make unlink

cp:
	cd eltbx; make cp
	cd sgtbx; make cp
	cd uctbx; make cp
	cd examples/cpp; make cp

rm:
	cd eltbx; make rm
	cd sgtbx; make rm
	cd uctbx; make rm
	cd examples/cpp; make rm

clean:
	cd eltbx; make clean
	cd sgtbx; make clean
	cd uctbx; make clean
	cd examples/cpp; make clean

depend:
	cd eltbx; make -f Makefile.nodepend depend > Makefile
	cd sgtbx; make -f Makefile.nodepend depend > Makefile
	cd uctbx; make -f Makefile.nodepend depend > Makefile
	cd examples/cpp; make -f Makefile.nodepend depend > Makefile

libraries:
	cd eltbx; make libraries
	cd sgtbx; make libraries
	cd uctbx; make libraries

compile:
	cd eltbx; make
	cd sgtbx; make
	cd uctbx; make

compile_examples:
	cd examples/cpp; make

libdir:
	-test ! -d lib && mkdir lib
	cp eltbx/libeltbx.a lib
	cp sgtbx/libsgtbx.a lib
	cp uctbx/libuctbx.a lib

libpythondir:
	-test ! -d lib/python/eltbx && mkdir -p lib/python/eltbx
	touch lib/python/eltbx/__init__.py
	cp eltbx/*.so lib/python/eltbx
	cp sgtbx/*.so lib/python
	cp uctbx/*.so lib/python

compile_all: compile libdir libpythondir compile_examples
compile_dev: depend compile_all
