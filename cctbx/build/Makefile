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
	cd lbfgs; make softlinks
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
	cd lbfgs; make unlink
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
	cd lbfgs; make cp
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
	cd lbfgs; make rm
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
	cd lbfgs; make clean
	cd examples/cpp; make clean

all_clean: clean
	find . \( -name \*.o -o -name \*.so -o -name \*tbx.a \) -exec rm {} \;

depend:
	cd external/boost_python; make -f Makefile.nodepend depend > Makefile
	cd uctbx; make -f Makefile.nodepend depend > Makefile
	cd sgtbx; make -f Makefile.nodepend depend > Makefile
	cd arraytbx; make -f Makefile.nodepend depend > Makefile
	cd adptbx; make -f Makefile.nodepend depend > Makefile
	cd eltbx; make -f Makefile.nodepend depend > Makefile
	cd sftbx; make -f Makefile.nodepend depend > Makefile
	cd fftbx; make -f Makefile.nodepend depend > Makefile
	cd lbfgs; make -f Makefile.nodepend depend > Makefile
	cd examples/cpp; make -f Makefile.nodepend depend > Makefile

libraries:
	cd uctbx; make libraries
	cd sgtbx; make libraries
	cd eltbx; make libraries

executables:
	cd uctbx; make executables
	cd sgtbx; make executables
	cd arraytbx; make executables
	cd lbfgs; make executables
	cd examples/cpp; make executables

all:
	cd external/boost_python; make
	cd uctbx; make
	cd sgtbx; make
	cd arraytbx; make
	cd adptbx; make
	cd eltbx; make
	cd sftbx; make
	cd fftbx; make
	cd lbfgs; make
	cd examples/cpp; make
