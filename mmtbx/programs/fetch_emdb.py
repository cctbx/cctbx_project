from __future__ import absolute_import, division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import json
import urllib.parse
import libtbx.utils
from libtbx.utils import Sorry

pdb_url_base = "https://search.rcsb.org/rcsbsearch/v2/query?json="
graphql_base = "https://data.rcsb.org/graphql?query="
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

  def api_query(self):
    graphql_query = '''
query
{
  entry(entry_id:"%s") {
    exptl {
      method
    }
    rcsb_entry_container_identifiers {
      emdb_ids
    }
  }
}
'''
    url = graphql_base + urllib.parse.quote(graphql_query % self.params.pdb_code)
    data = libtbx.utils.urlopen(url).read()
    d = json.loads(data)
    entry_data = d['data']['entry']
    exptl = entry_data['exptl'][0]
    if exptl['method'] != 'ELECTRON MICROSCOPY':
      raise Sorry('This entry is not an EM structure.')
    emdb_ids = entry_data['rcsb_entry_container_identifiers']['emdb_ids']
    if len(emdb_ids)==0:
      raise Sorry('No associated EMDB IDs found for this entry.')
    print('Found the following map codes associated to PDB ID: %s'
      % self.params.pdb_code, file=self.logger)
    for emdb_id in emdb_ids:
      print(emdb_id, file=self.logger)
    self.emdb_ids = emdb_ids

  # ---------------------------------------------------------------------------

  def download_pdb_file(self):
    pdb_url = 'https://files.rcsb.org/download/%s.pdb' % self.params.pdb_code
    pdb_fn = '%s.pdb' % self.params.pdb_code
    urllib.request.urlretrieve(pdb_url, pdb_fn)

  def download_maps(self):
    for emdb_id in self.emdb_ids:
      emdb_id_numeral = emdb_id.split('-')[1]
      emdb_url = emdb_base + 'EMD-%s/map/emd_%s.map.gz' % (
        emdb_id_numeral, emdb_id_numeral)
      map_fn = '%s_%s.map' % (self.params.pdb_code, emdb_id_numeral)
      urllib.request.urlretrieve(emdb_url, map_fn)

  # ---------------------------------------------------------------------------

  def run(self):
    self.emdb_ids = []
    self.api_query()
    self.download_pdb_file()
    self.download_maps()
