MAKEFLAGS=--no-print-directory

usage:
	@echo 'usage:'
	@echo 'make ...'

softlinks:
	cd external/boost_python; make softlinks
	cd uctbx; make softlinks
	cd sgtbx; make softlinks
	cd arraytbx; make softlinks
	cd adptbx; make softlinks
	cd eltbx; make softlinks
	cd sftbx; make softlinks
	cd fftbx; make softlinks
	cd examples/cpp; make softlinks

unlink:
	cd external/boost_python; make unlink
	cd uctbx; make unlink
	cd sgtbx; make unlink
	cd arraytbx; make unlink
	cd adptbx; make unlink
	cd eltbx; make unlink
	cd sftbx; make unlink
	cd fftbx; make unlink
	cd examples/cpp; make unlink

cp:
	cd external/boost_python; make cp
	cd uctbx; make cp
	cd sgtbx; make cp
	cd arraytbx; make cp
	cd adptbx; make cp
	cd eltbx; make cp
	cd sftbx; make cp
	cd fftbx; make cp
	cd examples/cpp; make cp

rm:
	cd external/boost_python; make rm
	cd uctbx; make rm
	cd sgtbx; make rm
	cd arraytbx; make rm
	cd adptbx; make rm
	cd eltbx; make rm
	cd sftbx; make rm
	cd fftbx; make rm
	cd examples/cpp; make rm

clean:
	cd external/boost_python; make clean
	cd uctbx; make clean
	cd sgtbx; make clean
	cd arraytbx; make clean
	cd adptbx; make clean
	cd eltbx; make clean
	cd sftbx; make clean
	cd fftbx; make clean
	cd examples/cpp; make clean

depend:
	cd external/boost_python; make -f Makefile.nodepend depend > Makefile
	cd uctbx; make -f Makefile.nodepend depend > Makefile
	cd sgtbx; make -f Makefile.nodepend depend > Makefile
	cd arraytbx; make -f Makefile.nodepend depend > Makefile
	cd adptbx; make -f Makefile.nodepend depend > Makefile
	cd eltbx; make -f Makefile.nodepend depend > Makefile
	cd sftbx; make -f Makefile.nodepend depend > Makefile
	cd fftbx; make -f Makefile.nodepend depend > Makefile
	cd examples/cpp; make -f Makefile.nodepend depend > Makefile

libraries:
	cd uctbx; make libraries
	cd sgtbx; make libraries
	cd eltbx; make libraries

compile:
	cd external/boost_python; make
	cd uctbx; make
	cd sgtbx; make
	cd arraytbx; make
	cd adptbx; make
	cd eltbx; make
	cd sftbx; make
	cd fftbx; make
	cd examples/cpp; make

lib_python_dir:
	-test ! -d lib_python/cctbx_boost && mkdir -p lib_python/cctbx_boost
	touch lib_python/cctbx_boost/__init__.py
	-test ! -d lib_python/cctbx_boost/arraytbx && mkdir lib_python/cctbx_boost/arraytbx
	touch lib_python/cctbx_boost/arraytbx/__init__.py
	-test ! -d lib_python/cctbx_boost/eltbx && mkdir lib_python/cctbx_boost/eltbx
	touch lib_python/cctbx_boost/eltbx/__init__.py
	cp uctbx/*.so lib_python/cctbx_boost
	cp sgtbx/*.so lib_python/cctbx_boost
	cp arraytbx/*.so lib_python/cctbx_boost/arraytbx
	cp adptbx/*.so lib_python/cctbx_boost
	cp eltbx/*.so lib_python/cctbx_boost/eltbx
	cp sftbx/*.so lib_python/cctbx_boost
	cp fftbx/*.so lib_python/cctbx_boost

all: compile lib_python_dir
