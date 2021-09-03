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

# Display the first miller array
myHKLview.SetColumn(0)

# Increase the size of all spheres by a factor 2
myHKLview.SetRadiiScale(2)

# Let's decrease the size difference between the smallest and the largest spheres shown
myHKLview.SetRadiiScale(1, nth_power_scale=0.2)

# Set all the spheres to size 1.5 regardless of the magnitudes of the data they represent
myHKLview.SetRadiiScale(1.5, nth_power_scale=0.0)

# Sort reflections into resolution bins with wavelengths as in the numbers below
myHKLview.SetColumnBinThresholds([1.8, 1.9, 1.95, 2.0, 2.05, 2.1, 2.16, 2.22, 2.3, 2.4, 2.5, 2.65, 2.8, 3.0, 3.25, 3.6, 4.0, 4.5, 5.0, 6.0, 8, 9, 20, 100 ])

# Make all reflections except those in the last bin invisible
for i in range(20):
  myHKLview.SetOpacity(i, 0.0)

# Sort reflections into bins according to data values in the
# 4th miller array presumably with values between 0 and 10
myHKLview.SetColumnBinThresholds([0, 0.1, 1, 10], 3)

# colour spheres according to values of the 3rd miller array
myHKLview.SetColourColumn(2)

# scale the size of spheres according to  values of the 4th miller array
myHKLview.SetRadiusColumn(3)

# displaying amplitudes or intensities with associated sigmas makes it possible to:
myHKLview.SetColoursToSigmas(True)
myHKLview.SetColoursToSigmas(False)
myHKLview.SetRadiiToSigmas(True)
myHKLview.SetRadiiToSigmas(False)




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

# Don't show slice
myHKLview.ShowSlice(False)


# Assuming the 2nd miller array contains map coefficients and the 5th miller array
# contains FOM values, then display these as reflections scaled according to the
# amplitudes, coloured according to the phases with colour saturation matching the FOM values
myHKLview.SetColumn(1,4)
