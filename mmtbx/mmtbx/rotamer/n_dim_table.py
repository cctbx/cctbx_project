import libtbx.load_env
import array
from math import floor
import re
import types
import os


class NDimTable:

    ourName = '' # identifying name for this table
    nDim    = 0  # number of dimensions of data
    minVal  = [] # nDim doubles, minimum allowed coordinate per dimension
    maxVal  = [] # nDim doubles, maximum allowed coordinate per dimension
    nBins   = [] # nDim ints, number of bins in each dimension
    doWrap  = [] # nDim booleans, does the dimension wrap around, like a dihedral angle?
    wBin    = [] # nDim doubles, width of each bin per dimension (nice to precalculate)
    lookupTable = None # an array of floating point numbers; the actual data

    def createFromText(infile):
        '''Loads rotamer or Ramachandran data from a text file, returning a new object.

        Can pass in either a file handle or a file name'''
        if isinstance(infile, types.StringTypes):
            infile = file(infile, 'r')
        ndt = NDimTable()
        ndt.ourName = re.search(r': +"(.+)"$', infile.readline()).group(1)
        ndt.nDim = int(re.search(r': +(\d+)$', infile.readline()).group(1))
        infile.readline() # lower_bound  upper_bound  number_of_bins  wrapping
        for i in range(ndt.nDim):
            match = re.search(r': +([^ ]+) +([^ ]+) +([^ ]+) +([^ ]+)$', infile.readline())
            ndt.minVal.append(float(match.group(1)))
            ndt.maxVal.append(float(match.group(2)))
            ndt.nBins.append(   int(match.group(3)))
            doWrap = match.group(4).lower().strip()
            if doWrap == "true" or doWrap == "yes" or doWrap == "on" or doWrap == "1": ndt.doWrap.append(True)
            else: ndt.doWrap.append(False)
            ndt.wBin.append((ndt.maxVal[i] - ndt.minVal[i]) / ndt.nBins[i])
        #print ndt.minVal, ndt.maxVal, ndt.nBins, ndt.doWrap, ndt.wBin
        nEntries = 1
        for i in range(ndt.nDim): nEntries *= ndt.nBins[i]
        ndt.lookupTable = array.array('f') # 4-byte floats
        # This appears to be the only way of setting the initial size of the array??
        for i in range(nEntries): ndt.lookupTable.append(0)
        if re.search(r'first', infile.readline()): valueFirst = True
        else: valueFirst = False
        #print valueFirst
        s = infile.readline()
        while s:
            s = infile.readline()
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
        return ndt
    createFromText = staticmethod(createFromText)

    def whereIs(self, coords):
        '''Given a set of coordinates, return the bin indices'''
        bin = []
        for i in range(self.nDim):
            bin.append(int( min(floor((coords[i]-self.minVal[i])/self.wBin[i]),  self.nBins[i] - 1) ))
        return bin

    def setValueAt(self, bin, val):
        self.lookupTable[ self.bin2index(bin) ] = val

    def bin2index(self, bin):
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

    def wrapbin(self, iBin, dim):
        # Python seems to have a sensible modulo operator (%),
        # unlike Java, so we don't need wrapbin much.
        if self.doWrap[dim]:
            return iBin % self.nBins[dim]
            #iBin = iBin % self.nBins[dim]
            #if iBin < 0: return iBin + self.nBins[dim]
            #else: return iBin
        else: return iBin

def run():
  rotamer_data_dir = libtbx.env.find_in_repositories("rotamer_data")
  ndt = NDimTable.createFromText(os.path.join(
    rotamer_data_dir, "rota500-lys.data"))
  # Lys takes 10 sec in Java from text, 47 sec in Python
  #ndt = NDimTable.createFromNDFT("/Users/ian/javadev/chiropraxis/resource/chiropraxis/rotarama/lys.ndft")
  # binary format is > 10x faster in Java

if (__name__ == "__main__"):
  run()
