"""

To use from a python prompt first do

cctbx.python -m pip install websocket_server 

Then start cctbx.python and from the python prompt execute:
"""

from crys3d.hklview import cmdlineframes

myHKLview = cmdlineframes.HKLViewFrame(jscriptfname = "$TEMP/hkljstr.js",  # path to generated temporary javascript file used by the html file
                                        htmlfname = "$TEMP/hkl.htm",  # path to generated temporary html file displayed by NGL
                                        verbose=True,  # make stdout more or less terse 
                                        UseOSBrowser = True  # if False the system default web browser will not be invoked
                                      )

# Any of the keywords can be left out. Their default values are presented here as well as their meaning in the comments

# Load a new reflection file, 4pa9.tncs.mtz
myHKLview.LoadReflectionsFile("4pa9.tncs.mtz") 

# Store a list describing each miller array in the variable myarrs
myarrs = myHKLview.GetArrayInfo()

# Let's display the first miller array
myHKLview.SetColumn(0) 

# Let's increase the size of all spheres by a factor 2
myHKLview.SetRadiiScale(2)

# Let's decrease the size difference between the smallest and the largest spheres shown
myHKLview.SetRadiiScale(1, nth_power_scale=0.2)

# Let's set all the spheres to size 0.3 regardless of the magnitude of the data they represent
myHKLview.SetRadiiScale(0.3, nth_power_scale=0.0)

# Sort reflections into resolution bins with wavelengths as in the numbers below
myHKLview.SetColumnBinThresholds([0, 3.4, 3.5, 3.6, 3.8, 4.0, 4.5, 5.0, 6.0, 8, 9,11, 15, 20, 50 ])

# Make all reflections except those in the last bin invisible
for i in range(14):
  myHKLview.SetOpacity(i, 0.0)

# Sort reflections into bins according to data values in the 
# 4th miller array presumably with values between 0 and 10
myHKLview.SetColumnBinThresholds([0, 0.1, 1, 10], 3)

# Expand reflections to P1 showing a half sphere of reflections
myHKLview.ExpandToP1(True)

# Don't expand reflections to P1
myHKLview.ExpandToP1(False)

# Show friedel mates 
myHKLview.ExpandAnomalous(True)

# Don't show friedel mates 
myHKLview.ExpandAnomalous(False)

# Show missing reflections as white spheres
myHKLview.ShowMissing(True)

# Don't show missing reflections as white spheres
myHKLview.ShowMissing(False)

# Show systematic absences
myHKLview.ShowSystematicAbsences(True)

# Don't show systematic absences
myHKLview.ShowSystematicAbsences(False)

# Get list of possible space subgroups
subgrouplst = myHKLview.GetSpaceGroupChoices()

# Set the space group to be the 4th space group stored in subgrouplst
# above and transform the miller array accordingly
myHKLview.SetSpaceGroupChoice(3)

# Show the slice of reflections with l=12
myHKLview.ShowSlice(True, "l", 12)

# Assuming the 2nd miller array contains map coefficients and the 5th miller array
# contains FOM values, then display these as reflections scaled according to the 
# amplitudes, coloured according to the phases with colour saturation matching the FOM values
myHKLview.SetColumn(1,4)






