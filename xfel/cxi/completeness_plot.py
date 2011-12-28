import os,math
from cctbx.array_family import flex
from cctbx.crystal import symmetry
from rstbx.apps.stills.simple_integration import show_observations
op = os.path

def get_observations(set,file_names,params):
  from libtbx import easy_pickle

  print "Number of pickle files found:", len(file_names)
  print
  returned=0
  from cctbx.sgtbx.bravais_types import bravais_lattice
  ref_bravais_type = bravais_lattice(set.space_group_info().type().number())

  for name in file_names:
    if name=="stats.pickle":continue

    full_path = file_name = op.abspath(name)
    obj = easy_pickle.load(file_name=full_path)
    if not obj.has_key("observations"): continue
    unit_cell = obj["observations"][0].unit_cell()
    result_array = obj["observations"][0]
    #unit_cell, img_miller_indices = obj

    print file_name,unit_cell, obj["observations"][0].space_group_info()

    if not bravais_lattice(
      obj["observations"][0].space_group_info().type().number()) == ref_bravais_type:
      print "Skipping cell in different Bravais type"
      continue

    if not unit_cell.is_similar_to(
      other=set.unit_cell(),
      relative_length_tolerance=0.1,
      absolute_angle_tolerance=2):
      print "Skipping cell with outlier dimensions"
      continue

    obj["observations"][0].show_summary()
    # Now do manipulate the data to conform to unit cell, asu, and space group of reference
    # Only works if there is NOT an indexing ambiguity!
    ref_setting_obs = obj["observations"][0].customized_copy(crystal_symmetry=set.crystal_symmetry()
      ).resolution_filter(d_min=params.d_min).map_to_asu()
    returned+=1
    yield ref_setting_obs,name

  print "Only %d of %d obs arrays had the correct cell"%(returned,len(file_names))

def plot_overall_completeness(completeness):
  completeness_range = xrange(-1,flex.max(completeness)+1)
  completeness_counts = [completeness.count(n) for n in completeness_range]
  from matplotlib import pyplot as plt
  plt.plot(completeness_range,completeness_counts,"r+")
  plt.show()


def run(args):
  import iotbx.phil
  phil = iotbx.phil.process_command_line(args=args, master_string="""
target_unit_cell = 78,78,37,90,90,90
  .type = unit_cell
target_space_group = P43212
  .type = space_group
d_min = 2.1
  .type = float
plot = False
  .type = str
cut_short_at = None
  .type = int
""").show()
  print
  work_params = phil.work.extract()
  assert work_params.d_min is not None

  print work_params.target_unit_cell
  print work_params.target_space_group

  from cctbx import miller
  miller_set = symmetry(
      unit_cell=work_params.target_unit_cell,
      space_group_info=work_params.target_space_group
    ).build_miller_set(
      anomalous_flag=True,
      d_min=work_params.d_min
    )

  miller_set.show_summary()

  # reality check
  #recip_cell_volume = work_params.target_unit_cell.reciprocal().volume()
  #recip_sphere_volume = (4/3)*math.pi*math.pow(1./work_params.d_min,3)
  #resolution_cells = recip_sphere_volume/recip_cell_volume
  #print "Number of asu's in sphere=",resolution_cells/miller_set.size()

  results = get_observations(miller_set,phil.remaining_args,work_params)

  # Create (and initialise?) arrays for statistics on the set of the
  # observed reflections which are present in the reference data set.
  completeness = flex.int(miller_set.size())
  sum_I        = flex.double(miller_set.size())
  sum_I_SIGI   = flex.double(miller_set.size())
  #last = completeness.deep_copy()

  for result,filename in results:
    result.show_summary()
    show_observations(result)

    if work_params.plot==True:
      #private interface to get the very strong diffraction images
      import StringIO
      G = StringIO.StringIO()
      show_observations(result,out=G)
      for line in G.getvalue().split("\n"):
        tokens = line.split()
        if len(tokens)>6:
          try:
            if float(tokens[3]) < 2.6 and float(tokens[-1]) > 10:
              print "Strong signal",filename,line
          except ValueError: pass
    print

    # Match up the observed intensities against the reference data
    # set, i_model, instead of the pre-generated miller set,
    # miller_set.
    matches = miller.match_indices(
      miller_set.indices(),
      result.indices())

    #for ih,hkl in enumerate(result.indices()):
    #  print hkl, result.data()[ih]
    print

    # Update the count for each matched reflection.
    completeness +=  (~matches.single_selection(0)).as_int()
    for pair in matches.pairs():
      sum_I[pair[0]] += result.data()[pair[1]]
      sum_I_SIGI[pair[0]] += (result.data()[pair[1]]/result.sigmas()[pair[1]])
    #for ih,hkl in enumerate(miller_set.indices()):
    #  print "%15s"%str(hkl),"%4d"%last[ih],"%4d"%completeness[ih], sum_I[ih]

    #print matches
    #help(matches)
    #print matches.pair_selection(0)
    #for item in matches.pairs(): print item
    #print list(miller_set.indices().select(matches.pairs().column(1)))
    #print list(~matches.single_selection(0))
    #print list(~matches.single_selection(1))
    #last  = completeness.deep_copy()

  #plot_overall_completeness(completeness)

  show_overall_observations(miller_set,completeness,sum_I,sum_I_SIGI)

  #from libtbx import easy_pickle
  #easy_pickle.dump(file_name="stats.pickle", obj=stats)
  #stats.report(plot=work_params.plot)
  #miller_counts = miller_set_p1.array(data=stats.counts.as_double()).select(
  #  stats.counts != 0)
  #miller_counts.as_mtz_dataset(column_root_label="NOBS").mtz_object().write(
  #  file_name="nobs.mtz")

def show_overall_observations(obs,redundancy,I,I_SIGI,out=None):
  if out==None:
    import sys
    out = sys.stdout
  from libtbx.str_utils import format_value

  obs.setup_binner(n_bins = 15)
  result = []
  for i_bin in obs.binner().range_used():
    sel_w = obs.binner().selection(i_bin)
    sel_fo_all = obs.select(sel_w)
    d_max_,d_min_ = sel_fo_all.d_max_min()
    d_range = obs.binner().bin_legend(
      i_bin=i_bin, show_bin_number=False, show_counts=False)
    sel_redundancy = redundancy.select(sel_w)
    sel_absent = sel_redundancy.count(0)
    sel_complete_tag = "[%d/%d]"%(sel_redundancy.size()-sel_absent,sel_redundancy.size())
    sel_measurements = flex.sum(sel_redundancy)
    sel_data = I.select(sel_w)
    sel_sig = I_SIGI.select(sel_w)
    if (sel_data.size() > 0 and sel_measurements>0):
      bin = resolution_bin(
        i_bin        = i_bin,
        d_range      = d_range,
        redundancy   = flex.mean(sel_redundancy.as_double()),
        complete_tag = sel_complete_tag,
        measurements = sel_measurements,
        mean_I       = flex.sum(sel_data)/sel_measurements,
        mean_I_sigI  = flex.sum(sel_sig)/sel_measurements,
        )
      result.append(bin)
  print >>out, "\n Bin  Resolution Range      Compl. <Redundancy>  #Measurements  <I>     <I/sig(I)>"
  for bin in result:
    fmt = " %s %s %s %s       %s      %s   %s"
    print >>out,fmt%(
      format_value("%3d",   bin.i_bin),
      format_value("%-13s", bin.d_range),
      format_value("%13s",  bin.complete_tag),
      format_value("%4.0f", bin.redundancy),
      format_value("%8d",   bin.measurements),
      format_value("%8.1f", bin.mean_I),
      format_value("%8.1f", bin.mean_I_sigI),
    )
class resolution_bin(object):
  def __init__(self,
               i_bin         = None,
               d_range       = None,
               redundancy    = None,
               absent        = None,
               complete_tag  = None,
               measurements  = None,
               mean_I        = None,
               mean_I_sigI   = None,
               sigmaa        = None):
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())


if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
