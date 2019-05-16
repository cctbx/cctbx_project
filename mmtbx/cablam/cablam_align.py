from __future__ import absolute_import, division, print_function
# (jEdit options) :folding=explicit:collapseFolds=1:
################################################################################
#
# Sandbox to play with cablam objects
#
################################################################################

import os, sys
from mmtbx.cablam import cablam_validate
import libtbx.phil.command_line
from iotbx import pdb
from iotbx import file_reader
from libtbx import group_args
import boost.python
from six.moves import range
cpputils = boost.python.import_ext("mmtbx_cablam_align_utils_ext")

# {{{ phil
#-------------------------------------------------------------------------------
master_phil = libtbx.phil.parse("""
cablam_align {
  pdb_infile_1 = None
    .type = path
    .help = '''input PDB file'''
  pdb_infile_2 = None
    .type = path
    .help = '''input PDB file'''
  chain_1 = None
    .type = str
    .help = '''chain for pdb_infile_1'''
  chain_2 = None
    .type = str
    .help = '''chain for pdb_infile_2'''
  help = False
    .type = bool
    .help = '''help and data interpretation messages'''
  threshold = 10
    .type = float
    .help = '''keeps all regions who have diff list averages lower than this'''
  window_size = 10
    .type = int
    .help = '''the size of the regions when looking for similarities'''
}
""", process_includes=True)
#-------------------------------------------------------------------------------
#}}}

# {{{ usage
def usage():
  sys.stderr.write("""
--------------------------------------------------------------------------------
python cablam_align.py 1234.pdb abcd.pdb

The user must provide TWO and only TWO pdb files.
--------------------------------------------------------------------------------
""")
#}}}

# {{{ cablam_residues
class cablam_residues(object):

  #{{{ get_continuous_segments
  # retuns a list of segments. segments are list of cablam linked_residue objects
  def get_continuous_segments(self):
    segments = []
    n_residues =[v for k,v in self.cablam_residues.items() if v.prevres == None]
    for res in n_residues :
      cur_res = res
      segment = []
      while cur_res != None :
        segment.append(cur_res)
        cur_res = cur_res.nextres
      segments.append(segment)
    return segments
  #}}}

  #{{{ __init__
  def __init__(self,
               pdb_infile = None,
               pdbid      = None,
               hierarchy  = None):
    assert pdb_infile != None or hierarchy != None
    if pdb_infile != None and hierarchy != None :
      w = "WARNING: pdb_infile AND hierarchy detected; hierarchy will be used.\n"
      sys.stderr.write(w)
    if hierarchy == None :
      pdb_io = pdb.input(pdb_infile)
      hierarchy = pdb_io.construct_hierarchy()
    assert hierarchy != None
    if pdbid == None : self.pdbid = os.path.basename(pdb_infile)[:-4]
    else : self.pdbid = pdbid
    self.cablam_residues = cablam_validate.setup(hierarchy = hierarchy,
                                                 pdbid     = self.pdbid)

    # get continuous segments
    self.segments = self.get_continuous_segments()
  #}}}

#}}}

# {{{ not_enough_pdbs
def not_enough_pdbs():
  sys.stderr.write("You must provide TWO and only TWO pdb files\n")
  usage()
  sys.exit()
#}}}

# {{{ FindSimilarCablamRegions
################################################################################
#
#  FindSimilarCablamRegions will align cablam measures.
#    residues_1 and residues_2 are cablam_residues objects (found in this
#      script)
#    window_size is the size of the regions you want to align
#    in_and_out; false means use ONLY pseudo CA in true means use in AND out
#    threshold is the upper limit 'difference mean' that will be accepted when
#      determining if two regions are similar. Please see RLab wiki for more
#      info on how this works.
#
#  NOTE: if threshold is zero then the object simple_aligned_regions will return
#        a list that contains the one entry of aligned residue that had the
#        smallest difference mean.
#
################################################################################
class FindSimilarCablamRegions(object):

  def __init__(self,
               residues_1,
               residues_2,
               window_size,
               in_and_out = False,
               threshold = None):
    self.residues_1     = residues_1
    self.residues_2     = residues_2
    self.window_size    = window_size
    self.in_and_out     = in_and_out
    self.threshold      = threshold

    # get simple_aligned_regions
    self.simple_align_measures()


  # {{{ get_measures
  def get_measures(self, segment):
    measures = []
    for r in segment :
      if len(r.measures) == 0 :
        # needed for cpputils
        measures.append(None)
        continue
      if self.in_and_out :
        measures.append(r.measures["CA_d_out"] + r.measures["CA_d_in"])
      else :measures.append(r.measures["CA_d_in"])
    return measures
  # }}}

  # {{{ simple_align_measures
  def simple_align_measures(self):
    similar_regions = []
    smallest_mean = 5000
    self.simple_aligned_regions = []
    # for each continuous segment in res 1 compare to each segment it seg 2
    # and find similar regions
    for seg_1 in self.residues_1.segments :
      measures_1 = self.get_measures(seg_1)
      for seg_2 in self.residues_2.segments :
        measures_2 = self.get_measures(seg_2)
        # sr = find_similar_regions(measures_1,measures_2)
        similar_regions = cpputils.get_similar_regions(measures_1,
                                                       measures_2,
                                                       self.threshold,
                                                       self.window_size)

        # similar_regions is a list of objects (made in cpp) having these
        # attributes: i_1, i_2, mean, window_length. i_1 and i_2 are the stating
        # indecies in measures_1 and measures_2, respectively, which align. By
        # extention, i_1 and i_2 are also the stating idecies in seg_ and seg_2.
        # get similar regions residues (cablam linked_residue objects)
        for sr in similar_regions :
          is_1 = [i for i in range(sr.i_1, sr.i_1 + sr.window_length)]
          is_2 = [i for i in range(sr.i_2, sr.i_2 + sr.window_length)]
          ar = []
          for i in range(len(is_1)):
            ar.append((seg_1[is_1[i]], seg_2[is_2[i]]))
          if self.threshold > 0 :
            self.simple_aligned_regions.append(group_args(aligned_residues = ar,
                                                          mean = sr.mean))
          else:
            if sr.mean < smallest_mean :
              smallest_mean = sr.mean
              self.simple_aligned_regions = [group_args(aligned_residues = ar,
                                                        mean = sr.mean)]
  # }}}



# }}}

# {{{ AlignCablamMeasures
class AlignCablamMeasures(object):

  # {{{ __init__
  def __init__(self, hierarchy_1,
                     hierarchy_2,
                     pdb_name_1,
                     pdb_name_2,
                     window_size,
                     threshold):

    # get cablam residues
    residues_1 = cablam_residues(hierarchy = hierarchy_1, pdbid = pdb_name_1)
    residues_2 = cablam_residues(hierarchy = hierarchy_2, pdbid = pdb_name_2)

    #find similar fragments in the two pdb
    self.simple_similar_region_obj = \
                          FindSimilarCablamRegions(residues_1 = residues_1,
                                                   residues_2 = residues_2,
                                                   window_size = window_size,
                                                   in_and_out = False,
                                                   threshold  = threshold)

    self.hierarchy_1 = hierarchy_1
    self.hierarchy_2 = hierarchy_2
    # Now we have a list of similar regions (cablam_aln_obj.simple_aligned_regions).
    # These regions can overlap. For such cases we want to find extended regions.
    # i.e. we have two lists (this is simplified) [1,2,3,4,5] and [3,4,5,6,7],
    # we can then assume that [1,2,3,4,5,6,7] is a good alignment too.
    # get aligned_residues
    self.condense_simple_aligned_regions()
  # }}}

  # {{{ print_simple_similar_regions
  def print_simple_similar_regions(self, log):
    for i,ar in enumerate(self.simple_similar_region_obj.simple_aligned_regions):
      print("similar fragment %i" % i, file=log)
      print("mean : %.03f, diff : %i" % (ar.mean,abs(ar.aligned_residues[0][0].resnum-ar.aligned_residues[0][1].resnum)), file=log)
      for t in ar.aligned_residues :
        print(t[0].id_to_str(),t[1].id_to_str(), file=log)
  # }}}

  # {{{ condense_simple_aligned_regions
  def condense_simple_aligned_regions(self):

    self.aligned_regions = []
    sar = self.simple_similar_region_obj.simple_aligned_regions
    current_region = sar[0]
    # iterate through simple_aligned_regions
    i = 0
    break_one, break_two = False,False
    new_region = None#sar[0][:]
    while 1 :
      if new_region == None : new_region = sar[i].aligned_residues[:]
      i += 1
      if i == len(sar):
        if new_region not in self.aligned_regions :
          self.aligned_regions.append(new_region)
        break
      next_region = sar[i].aligned_residues[:]

      # if the first pair in next_region is not in new_region no others will be
      if next_region[0] not in new_region :
        if new_region not in self.aligned_regions :
          self.aligned_regions.append(new_region)
        # check to see if any regions already in self.aligned_regions align
        new_region = None
        for region in self.aligned_regions :
          if next_region[0] in region :
            for pair in next_region :
              if pair in region : continue
              else : region.append(pair)
            new_region = region
            break
        if new_region == None : new_region = next_region
      else :
        for pair in next_region :
          if pair in new_region : continue
          else : new_region.append(pair)

  # }}}

  # {{{ print_aligned_regions
  def print_aligned_regions(self,log):
    for i,region in enumerate(self.aligned_regions):
      self.print_region(region,i,log)
  # }}}

  # {{{ print_region
  def print_region(self, region, title, log):
    print("*"*77, file=log)
    print("**** %s" % title, file=log)
    print("diff : %i" % abs(region[0][0].resnum-region[0][1].resnum), file=log)
    for t in region :
      print(t[0].id_to_str(), t[1].id_to_str(), file=log)
  # }}}

# }}}

# {{{ get_chain_selection
def get_chain_selection(hierarchy, chain):
  s = "chain %s" % chain.upper()
  sel = hierarchy.atom_selection_cache().selection(string = s)
  return hierarchy.select(sel)
# }}}

# {{{ cablam_align_wrapper
# wrapper for AlignCablamMeasures
def cablam_align_wrapper(pdb_1, pdb_2, window_size, threshold,
                         chain_1 = None, chain_2 = None,
                         pdb_name_1 = None, pdb_name_2 = None):

  # {{{ get pdb hierarchies
  if pdb_name_1 == None : pdb_name_1 = os.path.basename(pdb_1)[:-4]
  if pdb_name_2 == None : pdb_name_2 = os.path.basename(pdb_2)[:-4]
  pdb_io = pdb.input(pdb_1)
  hierarchy_1 = pdb_io.construct_hierarchy()
  pdb_io = pdb.input(pdb_2)
  hierarchy_2 = pdb_io.construct_hierarchy()

  # get chain selections
  if chain_1 != None : hierarchy_1 = get_chain_selection(hierarchy_1,chain_1)
  if chain_2 != None : hierarchy_2 = get_chain_selection(hierarchy_2,chain_2)
  # }}}

  return AlignCablamMeasures(hierarchy_1 = hierarchy_1,
                             hierarchy_2 = hierarchy_2,
                             pdb_name_1  = pdb_name_1,
                             pdb_name_2  = pdb_name_2,
                             window_size = window_size,
                             threshold   = threshold)

# }}}

# {{{ run
def run(args):

  #{{{ phil parsing
  #-----------------------------------------------------------------------------
  interpreter = libtbx.phil.command_line.argument_interpreter(master_phil=master_phil)
  sources = []
  pdbs = []
  for arg in args:
    if os.path.isfile(arg): #Handles loose filenames
      input_file = file_reader.any_file(arg)
      if (input_file.file_type == "pdb"):
        # sources.append(interpreter.process(arg="pdb_infile=\"%s\"" % arg))
        pdbs.append(arg)
      elif (input_file.file_type == "phil"):
        sources.append(input_file.file_object)
    else: #Handles arguments with xxx=yyy formatting
      arg_phil = interpreter.process(arg=arg)
      sources.append(arg_phil)
  work_phil = master_phil.fetch(sources=sources)
  work_params = work_phil.extract()
  params = work_params.cablam_align
  if work_params.cablam_align.pdb_infile_1 == None and work_params.cablam_align.pdb_infile_2 == None :
    if len(pdbs) != 2 : not_enough_pdbs()
    else :
      work_params.cablam_align.pdb_infile_1 = pdbs[0]
      work_params.cablam_align.pdb_infile_2 = pdbs[1]
  params = work_params.cablam_align

  if params.help:
    usage()
    sys.exit()
  #-----------------------------------------------------------------------------
  #}}} end phil parsing

  acm = cablam_align_wrapper(pdb_1       = params.pdb_infile_1,
                             pdb_2       = params.pdb_infile_2,
                             window_size = params.window_size,
                             threshold   = params.threshold,
                             chain_1     = params.chain_1,
                             chain_2     = params.chain_2)

  acm.print_aligned_regions(sys.stdout)

#}}}

# {{{ "__main__"
#-------------------------------------------------------------------------------
if __name__ == "__main__":
  run(sys.argv[1:])
#-------------------------------------------------------------------------------
#}}}

