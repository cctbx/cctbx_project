#! /bin/csh -f
set verbose
wget --recursive --no-parent --level=inf --timestamping --no-host-directories --cut-dirs=5 http://www.opengl.org/documentation/specs/man_pages/hardcopy/GL/html -o wget.log
rm -f index.html robots.txt
rm -rf html/gl_ftn html/glu_ftn html/glx html/index.html
