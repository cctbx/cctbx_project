#! /bin/csh -f
g++ $* -bundle -L/Library/Frameworks/Python.framework/Versions/2.3 -bundle_loader /Library/Frameworks/Python.framework/Versions/2.3/Python -framework /Library/Frameworks/Python.framework/Versions/2.3/Python
