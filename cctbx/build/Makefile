ROOT=/net/cci/rwgk

MAKEFLAGS=--no-print-directory

usage:
	@echo 'usage:'
	@echo 'make ...'

boot_tru64_cxx:
	make boot PLATFORM=tru64_cxx
boot_linux_gcc:
	make boot PLATFORM=linux_gcc
boot_irix_CC:
	make boot PLATFORM=irix_CC
boot_mingw32:
	make boot PLATFORM=mingw32
boot_vc60:
	make boot PLATFORM=vc60

boot:
	-test ! -d eltbx && mkdir eltbx
	-ln -s $(ROOT)/cctbx/eltbx/Makefile.$(PLATFORM) eltbx/Makefile.nodepend
	-test ! -d sgtbx && mkdir sgtbx
	-ln -s $(ROOT)/cctbx/sgtbx/Makefile.$(PLATFORM) sgtbx/Makefile.nodepend
	-test ! -d uctbx && mkdir uctbx
	-ln -s $(ROOT)/cctbx/uctbx/Makefile.$(PLATFORM) uctbx/Makefile.nodepend
	-test ! -d examples && mkdir examples
	-test ! -d examples/cpp && mkdir examples/cpp
	-ln -s $(ROOT)/cctbx/examples/cpp/Makefile.$(PLATFORM) examples/cpp/Makefile.nodepend

softlinks:
	cd eltbx; make -f Makefile.nodepend softlinks
	cd sgtbx; make -f Makefile.nodepend softlinks
	cd uctbx; make -f Makefile.nodepend softlinks
	cd examples/cpp; make -f Makefile.nodepend softlinks

unlink:
	cd eltbx; make -f Makefile.nodepend unlink
	cd sgtbx; make -f Makefile.nodepend unlink
	cd uctbx; make -f Makefile.nodepend unlink
	cd examples/cpp; make -f Makefile.nodepend unlink

clean:
	cd eltbx; make -f Makefile.nodepend clean
	cd sgtbx; make -f Makefile.nodepend clean
	cd uctbx; make -f Makefile.nodepend clean
	cd examples/cpp; make -f Makefile.nodepend clean

depend:
	cd eltbx; make -f Makefile.nodepend depend > Makefile
	cd sgtbx; make -f Makefile.nodepend depend > Makefile
	cd uctbx; make -f Makefile.nodepend depend > Makefile
	cd examples/cpp; make -f Makefile.nodepend depend > Makefile

rm_makefiles:
	cd eltbx; rm -f Makefile
	cd sgtbx; rm -f Makefile
	cd uctbx; rm -f Makefile
	cd examples/cpp; rm -f Makefile

nodepend:
	cd eltbx; make -f Makefile.nodepend
	cd sgtbx; make -f Makefile.nodepend
	cd uctbx; make -f Makefile.nodepend

compile:
	cd eltbx; make
	cd sgtbx; make
	cd uctbx; make

compile_examples:
	cd examples/cpp; make

nodepend_examples:
	cd examples/cpp; make -f Makefile.nodepend

unix_lib:
	-test ! -d lib && mkdir lib
	-test ! -d lib/python && mkdir lib/python
	-test ! -d lib/python/eltbx && mkdir lib/python/eltbx
	touch lib/python/eltbx/__init__.py
	cp eltbx/libeltbx.a lib
	cp eltbx/*.so lib/python/eltbx
	cp sgtbx/libsgtbx.a lib
	cp sgtbx/*.so lib/python
	cp uctbx/libuctbx.a lib
	cp uctbx/*.so lib/python

mingw32_lib:
	-test ! -d lib && mkdir lib
	-test ! -d lib/python && mkdir lib/python
	-test ! -d lib/python/eltbx && mkdir lib/python/eltbx
	touch lib/python/eltbx/__init__.py
	cp eltbx/libeltbx.a lib
	cp eltbx/*.pyd lib/python/eltbx
	cp sgtbx/libsgtbx.a lib
	cp sgtbx/*.pyd lib/python
	cp uctbx/libuctbx.a lib
	cp uctbx/*.pyd lib/python

vc60_lib:
	-test ! -d lib && mkdir lib
	-test ! -d lib/python && mkdir lib/python
	-test ! -d lib/python/eltbx && mkdir lib/python/eltbx
	touch lib/python/eltbx/__init__.py
	cp eltbx/libeltbx.lib lib
	cp eltbx/*.pyd lib/python/eltbx
	cp sgtbx/libsgtbx.lib lib
	cp sgtbx/*.pyd lib/python
	cp uctbx/libuctbx.lib lib
	cp uctbx/*.pyd lib/python

everything: softlinks depend compile unix_lib compile_examples
