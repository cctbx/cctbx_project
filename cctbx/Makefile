makefiles:
	cd eltbx; python MakeMakefiles.py
	cd sgtbx; python MakeMakefiles.py
	cd uctbx; python MakeMakefiles.py
	cd examples/cpp; python MakeMakefiles.py

clean:
	@find . -name '*.pyc' -print -exec rm {} \;
	@find . \( -name \*.c -o -name \*.cpp -o -name \*.h -o -name \*.hpp \) -exec egrep '	|[ 	]$$' {} /dev/null \;
