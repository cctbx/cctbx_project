:: syntax: testc 4-letter-pdb-code <any number of flags to suitename>
:: Runs old C suitename vs new Molprobity/CCTBX version

@echo off
setlocal enabledelayedexpansion
call C:\Users\Ken\Desktop\Richardson\molprobity\build\setpaths.bat

set name=%1
:: magic incantation to get all but the first command line arg:
for /f "tokens=1,* delims= " %%a in ("%*") do set args2plus=%%b
echo %args2plus%
set args1=-anglefields 6 -pointIDfields 7

cctbx.suitename %name%.dangle %args1% %args2plus% | grep -v suitename > %name%.py-out.txt
echo cctbx.suitename %name%.dangle %args1% %args2plus%
diff %name%.canonical.txt %name%.py-out.txt > %name%.diffs.txt
type %name%.diffs.txt
if ERRORLEVEL 1 (
	"C:\Program Files (x86)\Notepad++\notepad++.exe" %name%.diffs.txt
	EXIT /B %ERRORLEVEL%
)

