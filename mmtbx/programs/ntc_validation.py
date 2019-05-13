from __future__ import division, print_function
#
import requests #, json
from operator import itemgetter
from phenix.program_template import ProgramTemplate

dnatco_assign_cgi = "https://www.dnatco.org/cgi-bin/assign_from_torsions_34.py"
dnatco_nearest_cgi = "https://www.dnatco.org/cgi-bin/nearest_from_torsions_34.py"

# =============================================================================
class Program(ProgramTemplate):

  description = '''
Program for validating DNA/RNA from NtC web-service

Minimum required inputs:
  Model file

'''

  datatypes = ['model', 'phil']

  master_phil_str = '''
input {
  nproc = Auto
    .type = int
}
action {
  query_type = *angles xyz
    .type = choice(multi=False)
  outliers_only = False
    .type = bool
    .help = Print only the outliers
}
'''

  def validate(self):
    pass

  @staticmethod
  def get_validation_via_angles(query):
    query1 = {'ch': '271.6',
             'g1': '44.6',
             'CC': '5.0',
             'a1': '287.8',
             'P': '156.3',
             'P1': '116.9',
             'b1': '138.5',
             'step_id': '0_B_DG_22_DC_23',
             'd1': '112.8',
             'NCCN': '34.3',
             'e': '259.9',
             'd': '149.7',
             'NN': '4.7',
             'ch1': '234.7',
             'z': '171.6'}
    r1 = requests.post(dnatco_assign_cgi, json=query)
    r2 = requests.post(dnatco_nearest_cgi, json=query)
    return r1,r2

  @staticmethod
  def get_validation_via_coordinates(query):
    query1 = {"step_id":"102d_A_DC_1_DG_2",
              "C5pa":["18.939","34.713","89.428"],
              "C4pa":["20.086","33.755","89.286"],
              "O4pa":["19.617","32.539","89.905"],
              "C3pa":["20.434","33.382","87.819"],
              "O3pa":["21.831","32.949","87.719"],
              "C2pa":["19.488","32.196","87.629"],
              "C1pa":["19.521","31.485","88.974"],
              "N19a":["18.313","30.732","89.375"],
              "C24a":["18.488","29.453","89.907"],
              "Pb":["22.636","32.521","86.331"],
              "O5pb":["22.776","30.956","86.326"],
              "C5pb":["23.239","30.276","87.446"],
              "C4pb":["23.441","28.853","87.109"],
              "O4pb":["22.221","28.177","87.390"],
              "C3pb":["23.649","28.619","85.571"],
              "O3pb":["24.497","27.474","85.326"],
              "C2pb":["22.244","28.333","85.091"],
              "C1pb":["21.634","27.584","86.246"],
              "N19b":["20.190","27.715","86.370"],
              "C24b":["19.318","26.753","86.864"]}
    import json
    dnatco_cgi = "https://www.dnatco.org/cgi-bin/assign_from_coords_34.py"
    query = json.dumps(query)

    try:
      r = requests.post(dnatco_cgi, json=json.loads(query))
      jsonres = r.json()
    except Exception as e:
      print(repr(e))
      sys.exit()

    if 1:
      if (jsonres["error"] == ""):
        #print(jsonres["result"]["NtC"], jsonres["result"].get("confal", None))
        print(jsonres)
      else:
        print(jsonres["error"])
    return r,r

  def print_validation(self, r1, r2):
    jsonres = r1.json()
    if (jsonres["error"] == ""):
      print(jsonres["result"]["step_id"], jsonres["result"]["NtC"], jsonres["result"]["confalH"])
    else:
      print(jsonres["error"])

    jsonres = r2.json()
    nearest = []
    if (jsonres["error"] == ""):
      for ntc in jsonres["result"]:
  #     nearest.append( [str(ntc), float(jsonres["result"][ntc]["confalH"]), float(jsonres["result"][ntc]["confalA"]), float(jsonres["result"][ntc]["confalG"]), float(jsonres["result"][ntc]["bb_rmsd"]) ] + [float(i) for i in jsonres["result"][ntc]["confals"]] + [float(d) for d in ["%.2f" % float(i) for i in jsonres["result"][ntc]["distances"]]] )
        try:
          nearest.append( [ str(ntc),
                          float(jsonres["result"][ntc]["confalH"]),
                          float(jsonres["result"][ntc]["bb_rmsd"]) ])
        except Exception as e:
          pass
  #     nearest.sort(key = itemgetter(2), reverse=True)
      nearest.sort(key = itemgetter(2), reverse=False)
      for record in nearest:
          if (record[1] != 0.0):
              print(record)
    else:
      print(jsonres["error"])

  def run(self,
          log=None,
          verbose=False):
    assert log is None
    if log is None: log=sys.stdout
    print('Using model: %s' % self.data_manager.get_default_model_name())
    model = self.data_manager.get_model()

    from mmtbx.conformation_dependent_library import generate_dna_rna_fragments
    import time
    import psutil
    from libtbx import easy_mp
    from libtbx import Auto

    if self.params.action.query_type=='angles':
      validation_function = self.get_validation_via_angles
      query_attr = 'get_ntc_angles'
    elif self.params.action.query_type=='xyz':
      validation_function = self.get_validation_via_coordinates
      query_attr = 'get_ntc_coordinates'
    else:
      assert 0

    t0=time.time()

    if self.params.input.nproc==Auto:
      self.params.input.nproc=psutil.cpu_count()
    from collections import OrderedDict
    self.results = OrderedDict()
    argss = []
    for dna_rna_pairs in generate_dna_rna_fragments(
      model.get_hierarchy(),
      model.get_restraints_manager().geometry,
      length=2,
      ):
      print(dna_rna_pairs)
      query_function = getattr(dna_rna_pairs, query_attr)
      query = query_function()
      #print(query)
      argss.append([query])
      self.results[query['step_id']] = None

    print('\nUsing %d nprocs for %s suites\n' % (self.params.input.nproc,
                                                 len(argss)))
    for args, res, err_str in easy_mp.multi_core_run(validation_function,
                                                     argss,
                                                     self.params.input.nproc,
                                                     ):
      print('%sTotal time: %6.2f (s) for %s' % (' '*7,
                                                time.time()-t0,
                                                args[0]['step_id']
                                                ),
      )
      if err_str:
        print('Error output from %s' % args[0]['step_id'])
        print(err_str)
        print('_'*80)
      self.results[args[0]['step_id']] = [args, res, err_str]

    for step_id, (query, res, err_str) in self.results.items():
      print('-'*80)
      if res:
        self.print_validation(res[0], res[1])
      else:
        print('no results')

    return self.results

  def get_results(self):
    return self.results
