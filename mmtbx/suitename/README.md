# suitename
Suitename - a tool for classifying backbone "suite" conformations (linkages between ribose sugars) in an RNA molecule.

## Requirements
Suitename requires Python 2.7 or later and it requires numpy to have been installed.

## Running

Suitename is integrated into the CCTBX library and can be run from within it as `phenix.suitename`, `molprobity.suitename` or `cctbx.suitename`.

Suitename can also be run from the command line using a Python 3 interpreter: `python3 suitename.py [arguments]`

Suitename takes an input file as its first argument, or if none provided, reads from standard input. It writes its results to standard output. It classifies the conformation of each suite into one of several dozen predefined clusters that have been determined by years of study. 

Two forms of input are supported on standard input:
  1. A list of the 6 dihedral angles $/alpha$, $/beta$, $/gamma$, $/delta$, $/epsilon$, $/zeta$) for each residue of the RNA molecule. Suitename will re-parse the dihedral angles in the residues to obtain the 7 dihedral angles in each suite ($/delta/-1, $/epsilon/-1, $/zeta/-1, $/alpha$, $/beta$, $/gamma$, $/delta$), and then operate on those. This is the default input.

    Each line of this format describes one residue, with fields separated by colons. The first several (default 6) are ID information, the remainder are angles. Sample:
1: A:   7: : :  U:-75.533:-154.742:48.162:80.895:-148.423:-159.688
    If the number of fields before the first angle is not 6, use
    -- pointIDfields <n> 
    to tell Suitename the true number.

  2. A kinemage file providing a list of 7 or 9 dihedral angles in each suite. Mark this by using the --suitein command line flag. 9 angles include $/chi$-1 and $/chi$.

A third form of input is supported when running phenix.suitename or molprobity.suitename and passing the file name on the command line: a molecule file in CIF or PDB format.  When run on this type of file, Suitename will compute the dihedral angles for all RNA bases and compute the resulting scores.

Three forms of output are supported:
  1. A text report, showing the classification of each suite into a cluster or outlier, and the neatness of fit ("suiteness") of the suite into that cluster. Suiteness is the cosine of the normalized distance of a suite datapoint from the power 3 hyperellipsoid boundary of the cluster in toward its center. A statistical summary is included at the end. This is the default output format.
  2. A kinemage file, which will display a 7D data point for each suite in the data. Points are color-coded according to the clusters to which they have been assigned. Each cluster is displayed as a colored ring surrounding its defined center. Specify this by using the --kinemage command line flag.
  3. A brief string showing only the cluster assignments. It consists of three characters per suite - base identity (uc) and 2-character number-letter name of the suite cluster (e.g., C1aG1gU1aA1aA1cG). Specify this by using the --string command line flag.

Many other command line options are available; type `python3 suitename.py --help` to display them.

## Use as a library

See the suites.rtf document for more information about using suites as a library
rather than a program.

## Testing

In addition to the tests found in the unit-test direcgtory, this Python code was run
across all files in the PDB that have RNA and its output was
compared to the original C code output for -report, -kinemage, and -string outputs.
It was tested with both Python 2.7 and Python 3.6.
It was tested on the dangle-format output of mmtbx.mp_geo and found to produce the same
output modulo the following slight differences:
- Python 3 uses Banker's rounding, so the low-order digit may be off by one.
- Python performs floating-point calculations at double resolution, so some of
the calculations are slightly different for numbers right near the cutoffs, sometimes
making corner cases fall into a neighboring bin.

## Directories

**unit-test** Contains unit-test and regression test code for SuiteName.

