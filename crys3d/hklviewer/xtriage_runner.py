from __future__ import absolute_import, division, print_function
from libtbx import group_args, phil


def external_cmd(parent, master_phil, firstpart):
  from mmtbx.scaling import xtriage
  from io import StringIO
  tabname = "Xtriage"

  logstrbuf = StringIO()
  xtriageobj = xtriage.run([ parent.loaded_file_name, "scaling.input.xray_data.obs_labels=" +
                           parent.viewer.get_current_labels() ], out=logstrbuf)
  logfname = firstpart + "_xtriage.log"
  with open(logfname, "w") as f:
    f.write( logstrbuf.getvalue() )

  retval = 0
  errormsg = ""
  # Add any twin operators as user vectors.
  # user_vector is multiple scope so we can't assign viewer.user_vector directly
  philstr = ""
  for i,twinop in enumerate(xtriageobj.twin_results.twin_law_names):
    philstr += '''
    viewer.user_vector {
                label = "twin_"''' + str(i) + '''
                hkl_op = "''' + twinop + '''"
           }
    '''
  # Add TNCS vector as a real space vector if patterson peak > 0.1 of origin peak
  if len(xtriageobj.twin_results.translational_pseudo_symmetry.suspected_peaks) > 0:
    tncsvec = xtriageobj.twin_results.translational_pseudo_symmetry.suspected_peaks[0][0]
    patterson_height = xtriageobj.twin_results.translational_pseudo_symmetry.suspected_peaks[0][1]
    if patterson_height > 10:
      philstr += '''
      viewer.user_vector {
                  label = "tNCS_xtriage"
                  abc = ''' + f"({tncsvec[0]:.5}, {tncsvec[1]:.5}, {tncsvec[2]:.5})" + '''
              }
      '''
  vectorphil = phil.parse(philstr)
  working_params = master_phil.fetch(source= vectorphil).extract()
  parent.add_user_vector(working_params.viewer.user_vector, rectify_improper_rotation=False)

  # The name of logfile and tab should be present in ldic after running exec().
  # cctbx.python sends this back to HKLviewer from HKLViewFrame.run_external_cmd()
  return group_args(tabname=tabname, logfname=logfname, retval=retval, errormsg=errormsg)
