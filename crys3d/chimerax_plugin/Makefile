# The "make" targets are:
# 	wheel: build a Python wheel in "dist" directory.
# 	app-install: build wheel (if needed) and install in ChimeraX.
# 	test: run ChimeraX
# 	debug: run ChimeraX with debugging flag set
# 	clean: remove files used in building wheel
# 	distclean: remove files used in building wheel and license file

# These parameters may be changed as needed.

# ChimeraX bundle names must start with "ChimeraX_"
# to avoid clashes with package names in pypi.python.org.
# When uploaded to the ChimeraX toolshed, the bundle
# will be displayed without the ChimeraX- prefix.
# BUNDLE_NAME = ChimeraX-Clipper
# BUNDLE_VERSION = 0.9.4-dev1
# ChimeraX bundles should only include packages
# that install as chimerax.package_name.
# General Python packages should be uploaded to
# pypi.python.org rather than the ChimeraX toolshed.
# PKG_NAME = chimerax.clipper

# cmd line building on Windows:
# "A:\Program Files\ChimeraX 1.2.5\bin\ChimeraX-console.exe" --nogui --cmd "devel install . ; exit"


# Define where ChimeraX is installed.
OS = $(patsubst CYGWIN_NT%,CYGWIN_NT,$(shell uname -s))
# CHIMERAX_APP is the ChimeraX install folder

ifeq ($(OS),CYGWIN_NT)
ifndef RELEASE
# Windows
CHIMERAX_APP = "/c/Program Files/ChimeraX_Daily"
else
CHIMERAX_APP = "/a/Program Files/ChimeraX 1.2.5"
endif
endif

ifeq ($(OS),Darwin)
# Mac
ifndef RELEASE
CHIMERAX_APP = /Applications/ChimeraX_Daily.app
else
CHIMERAX_APP = /Applications/ChimeraX.app
endif
endif

ifeq ($(OS),Linux)
ifndef RELEASE
CHIMERAX_APP = /opt/UCSF/ChimeraX-daily
else
CHIMERAX_APP = /opt/UCSF/ChimeraX
endif
endif

# ==================================================================
# Theoretically, no changes are needed below this line

# Platform-dependent settings.  Should not need fixing.
# For Windows, we assume Cygwin is being used.
ifeq ($(OS),CYGWIN_NT)
CHIMERAX_EXE = $(CHIMERAX_APP)/bin/ChimeraX-console.exe
endif
ifeq ($(OS),Darwin)
CHIMERAX_EXE = $(CHIMERAX_APP)/Contents/bin/ChimeraX
endif
ifeq ($(OS),Linux)
CHIMERAX_EXE = $(CHIMERAX_APP)/bin/ChimeraX
endif

BUNDLE_BASE_NAME = $(subst ChimeraX-,,$(BUNDLE_NAME))
SOURCE = ../hklview
SRCS = $(SOURCE)/HKLviewer.py $(SOURCE)/HKLviewerGui.py $(SOURCE)/qt.py $(SOURCE)/helpers.py $(SOURCE)/QtChromiumCheck.py \
       src/tool.py src/__init__.py 


:DEFAULT_GOAL := wheel

#
# Actual make dependencies!
#

wheel $(WHEEL): bundle_info.xml $(SRCS) $(SOURCE)
	$(CHIMERAX_EXE) --nogui --cmd "devel build . ; exit"

install app-install:	$(WHEEL)
	$(CHIMERAX_EXE) --nogui --cmd "devel install . ; exit"

uninstall app-uninstall:	$(WHEEL)
	$(CHIMERAX_EXE) --nogui --cmd "toolshed uninstall $(BUNDLE_BASE_NAME) ; exit"

docs:
	$(CHIMERAX_EXE) -m sphinx docs/source src/docs/user

test:
	$(CHIMERAX_EXE)

debug:
	$(CHIMERAX_EXE) --debug

clean:
	$(CHIMERAX_EXE) --nogui --cmd "devel clean . ; exit"

.PHONY: docs
