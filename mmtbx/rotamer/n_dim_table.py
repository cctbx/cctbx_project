from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from math import floor
import re
from six.moves import range


class NDimTable:

    # I now know these are class (== Java "static") variables,
    # not instance variables, although I intended for them to be instance vars.
    # I had to remove these definitions or pickling didn't work correctly --
    # restored objects had empty arrays instead of the desired data.

    #ourName = '' # identifying name for this table
    #nDim    = 0  # number of dimensions of data
    #minVal  = [] # nDim doubles, minimum allowed coordinate per dimension
    #maxVal  = [] # nDim doubles, maximum allowed coordinate per dimension
    #nBins   = [] # nDim ints, number of bins in each dimension
    #doWrap  = [] # nDim booleans, does the dimension wrap around, like a dihedral angle?
    #wBin    = [] # nDim doubles, width of each bin per dimension (nice to precalculate)
    #lookupTable = None # an array of floating point numbers; the actual data

    @staticmethod
    def createFromText(infile):
        '''Loads rotamer or Ramachandran data from a text file, returning a new object.

        Can pass in either a file handle or a file name'''
        if isinstance(infile, str):
            infile = open(infile, 'r')
        ndt = NDimTable()
        ndt.ourName = re.search(r': +"(.+)"$', infile.readline()).group(1)
        ndt.nDim = int(re.search(r': +(\d+)$', infile.readline()).group(1))
        infile.readline() # lower_bound  upper_bound  number_of_bins  wrapping
        ndt.minVal  = []
        ndt.maxVal  = []
        ndt.nBins   = []
        ndt.doWrap  = []
        ndt.wBin    = []
        data_re = re.compile(r': +([^ ]+) +([^ ]+) +([^ ]+) +([^ ]+)$')
        for i in range(ndt.nDim):
            match = data_re.search(infile.readline())
            ndt.minVal.append(float(match.group(1)))
            ndt.maxVal.append(float(match.group(2)))
            ndt.nBins.append(   int(match.group(3)))
            doWrap = match.group(4).lower().strip()
            if doWrap == "true" or doWrap == "yes" or doWrap == "on" or doWrap == "1": ndt.doWrap.append(True)
            else: ndt.doWrap.append(False)
            ndt.wBin.append((ndt.maxVal[i] - ndt.minVal[i]) / ndt.nBins[i])
        #print ndt.minVal, ndt.maxVal, ndt.nBins, ndt.doWrap, ndt.wBin
        ndt.lookupTable = flex.float(flex.grid(ndt.nBins), 0)
        if re.search(r'first', infile.readline()): valueFirst = True
        else: valueFirst = False
        #print valueFirst
        s = infile.readline()
        while s:
            fields = s.split()
            if(len(fields) <= ndt.nDim): continue
            if valueFirst:
                val = float(fields[0])
                coords = [float(fld) for fld in fields[1:ndt.nDim+1]]
            else:
                val = float(fields[ndt.nDim])
                coords = [float(fld) for fld in fields[0:ndt.nDim]]
            #print val, coords
            ndt.setValueAt( ndt.whereIs(coords), val )
            s = infile.readline()
        return ndt


    def whereIs(self, coords):
        '''Given a set of coordinates, return the bin indices'''
        bin = []
        for i in range(self.nDim):
            bin.append(int( min(floor((coords[i]-self.minVal[i])/self.wBin[i]),  self.nBins[i] - 1) ))
        return bin


    def setValueAt(self, bin, val):
        '''Given bin indices, set the value in the table'''
        self.lookupTable[ self.bin2index(bin) ] = val


    def bin2index(self, bin):
        '''Given bin indices, return a single index into the linear array'''
        idx = 0
        for i in range(self.nDim-1):
            iBin = self.wrapbin(bin[i], i)
            if iBin < 0 or iBin >= self.nBins[i]: return -1
            idx += iBin
            idx *= self.nBins[i+1]
        i = self.nDim - 1
        iBin = self.wrapbin(bin[i], i)
        if iBin < 0 or iBin >= self.nBins[i]: return -1
        idx += iBin
        return idx


    def bin2index_limit(self, bin):
        '''Given bin indices, return a single index into the linear array.

        If no bin can be found after wrapping is applied, the edge of the table is used.'''
        idx = 0
        for i in range(self.nDim-1):
            iBin = self.wrapbin(bin[i], i)
            iBin = self.clampInclusive(0, iBin, self.nBins[i]-1)  # clamp to range [ 0, nBins[i] )
            idx += iBin
            idx *= self.nBins[i+1]
        i = self.nDim - 1
        iBin = self.wrapbin(bin[i], i)
        iBin = self.clampInclusive(0, iBin, self.nBins[i]-1)  # clamp to range [ 0, nBins[i] )
        idx += iBin
        return idx

    def clampInclusive(self, minVal, theVal, maxVal):
        '''Forces theVal into the range from minVal to maxVal (inclusive)'''
        return max(minVal, min(theVal, maxVal))


    def wrapbin(self, iBin, dim):
        '''Wrap bin indices for each dimension in which it is appropriate'''
        # Python seems to have a sensible modulo operator (%),
        # unlike Java, so we don't need wrapbin much.
        if self.doWrap[dim]:
            return iBin % self.nBins[dim]
            #iBin = iBin % self.nBins[dim]
            #if iBin < 0: return iBin + self.nBins[dim]
            #else: return iBin
        else: return iBin


    def centerOf(self, bin):
        '''Returns coordinates for the center of the bin given by these indices'''
        pt = []
        for i in range(self.nDim):
            pt.append(self.minVal[i] + self.wBin[i]*(bin[i]+0.5))
        return pt


    def valueAt(self, pt):
        '''Estimates the value of the density trace at a given position,
        using linear interpolation.

        This algorithm basically consults the 2^n bins nearest in space to
        the input point, and weights their contributions according to their
        distances.

        To get a feeling for how it works, work out a linear interpolation in
        one dimension
        ( a...x.......b -- what is x in terms of a and b? )(trivial)
        and then in two dimensions (still fairly easy -- do one dimension
        first, then interpolate between those results).

        In dimensions very near the edge of the table, no interpolation is
        performed in that dimension
        (it would be impossible to do so, anyway, because there's no 2nd
        value to interpolate out *to*).

        pt - the point at which the value should be estimated
        returns - the approximate value of the density trace at the specified
        point'''

        value = 0 # The value this function is calculating!
        va_home = self.whereIs(pt) # the bin this point falls into
        va_home_ctr = self.centerOf(va_home) # center of the above bin (nearest)
        va_neighbor = [] # the diagonal nearest-neighbor bin -- to be determined
        va_contrib = [] # relative contribution of va_neighbor (balance from va_home)
        va_current = [0]*self.nDim # current bin indices we're examining in the 2nd loop below

        # Initialize va_neighbor[] and va_contrib[]
        for dim in range(self.nDim):
            # If neighbor is out-of-bounds, bin2index_limit will handle it later
            if pt[dim] < va_home_ctr[dim]:  va_neighbor.append(va_home[dim]-1)
            else:                           va_neighbor.append(va_home[dim]+1)
            va_contrib.append(abs( (pt[dim]-va_home_ctr[dim])/self.wBin[dim] )) # always on [0.0, 0.5]
        #print pt, va_home, va_home_ctr, va_neighbor, va_contrib, self.nBins

        # Loop over all bins
        # bin is used as a bit mask, with one bit per dimension
        # 0 means va_home, 1 means va_neighbor -- this way all 2^n va_neighbor bins are evaluated
        # The limit is a 30-D table, but a 2x2x...x2 table in 30-D would occupy > 4GB !!!
        for bin in range(1 << self.nDim):
            coeff = 1.0 # reset coefficient; product of a mix of elements from va_contrib[] and (1-va_contrib[])
            # Loop over all dimensions, checking the appropriate bit in bin
            for dim in range(self.nDim):
                # Bit is off -- elements are drawn from va_home[] and (1-va_contrib[])
                if bin & (1 << dim) == 0:
                    va_current[dim] = va_home[dim]
                    coeff *= 1.0 - va_contrib[dim]
                # Bit is on -- elements are drawn from va_neighbor[] and va_contrib[]
                else:
                    va_current[dim] = va_neighbor[dim]
                    coeff *= va_contrib[dim]
            #print coeff, va_current, self.lookupTable[ self.bin2index_limit(va_current) ]
            value += coeff * self.lookupTable[ self.bin2index_limit(va_current) ] # calc. va_contribution of va_currently selected bin
        return value
