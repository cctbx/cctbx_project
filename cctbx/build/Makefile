MAKEFLAGS=--no-print-directory

usage:
	@echo 'usage:'
	@echo 'make ...'

softlinks:
	cd eltbx; make softlinks
	cd uctbx; make softlinks
	cd sgtbx; make softlinks
	cd adptbx; make softlinks
	cd sftbx; make softlinks
	cd examples/cpp; make softlinks

unlink:
	cd eltbx; make unlink
	cd uctbx; make unlink
	cd sgtbx; make unlink
	cd adptbx; make unlink
	cd sftbx; make unlink
	cd examples/cpp; make unlink

cp:
	cd eltbx; make cp
	cd uctbx; make cp
	cd sgtbx; make cp
	cd adptbx; make cp
	cd sftbx; make cp
	cd examples/cpp; make cp

rm:
	cd eltbx; make rm
	cd uctbx; make rm
	cd sgtbx; make rm
	cd adptbx; make rm
	cd sftbx; make rm
	cd examples/cpp; make rm

clean:
	cd eltbx; make clean
	cd uctbx; make clean
	cd sgtbx; make clean
	cd adptbx; make clean
	cd sftbx; make clean
	cd examples/cpp; make clean

depend:
	cd eltbx; make -f Makefile.nodepend depend > Makefile
	cd uctbx; make -f Makefile.nodepend depend > Makefile
	cd sgtbx; make -f Makefile.nodepend depend > Makefile
	cd adptbx; make -f Makefile.nodepend depend > Makefile
	cd sftbx; make -f Makefile.nodepend depend > Makefile
	cd examples/cpp; make -f Makefile.nodepend depend > Makefile

libraries:
	cd eltbx; make libraries
	cd uctbx; make libraries
	cd sgtbx; make libraries

compile:
	cd eltbx; make
	cd uctbx; make
	cd sgtbx; make
	cd adptbx; make
	cd sftbx; make

compile_examples:
	cd examples/cpp; make

libpythondir:
	-test ! -d lib/python/eltbx && mkdir -p lib/python/eltbx
	touch lib/python/eltbx/__init__.py
	cp eltbx/*.so lib/python/eltbx
	cp uctbx/*.so lib/python
	cp sgtbx/*.so lib/python
	cp adptbx/*.so lib/python
	cp sftbx/*.so lib/python

compile_all: compile compile_examples libpythondir
compile_dev: depend compile_all
