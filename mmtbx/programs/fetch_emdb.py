from __future__ import absolute_import, division, print_function
from libtbx.program_template import ProgramTemplate
from libtbx.utils import Sorry
from mmtbx.wwpdb import rcsb_web_services

pdb_url_base = "https://search.rcsb.org/rcsbsearch/v2/query?json="
emdb_base = "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/"

master_phil_str = '''
pdb_code = None
  .type = str
  .short_caption = PDB code
  .help = PDB code
  .input_size = 40
'''

# =============================================================================

class Program(ProgramTemplate):

  description = '''

phenix.fetch_emdb code=1aba

'''

  datatypes = ['phil']

  master_phil_str = master_phil_str

  # ---------------------------------------------------------------------------

  def validate(self):
    print('Validating inputs:\n', file=self.logger)
    if len(self.params.pdb_code)!=4:
      raise Sorry('Invalid PDB code provided.')

  # ---------------------------------------------------------------------------

  def download_pdb_file(self):
    pdb_url = 'https://files.rcsb.org/download/%s.pdb' % self.params.pdb_code
    pdb_fn = '%s.pdb' % self.params.pdb_code
    urllib.request.urlretrieve(pdb_url, pdb_fn)

  def download_maps(self, emdb_ids):
    if emdb_ids is None:
      raise Sorry('No associated EMDB IDs found for this entry.')
    print('Found the following map codes associated to PDB ID: %s'
      % self.params.pdb_code, file=self.logger)
    for emdb_id in emdb_ids:
      print("Downloading %s" % emdb_id, file=self.logger)
      emdb_id_numeral = emdb_id.split('-')[1]
      emdb_url = emdb_base + 'EMD-%s/map/emd_%s.map.gz' % (
        emdb_id_numeral, emdb_id_numeral)
      map_fn = '%s_%s.map' % (self.params.pdb_code, emdb_id_numeral)
      urllib.request.urlretrieve(emdb_url, map_fn)

  # ---------------------------------------------------------------------------

  def run(self):
    emdb_ids = rcsb_web_services.get_emdb_id_for_pdb_id(self.params.pdb_code)
    self.download_pdb_file()
    self.download_maps(emdb_ids)
