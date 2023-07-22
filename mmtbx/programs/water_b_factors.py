from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate
# from libtbx.utils import null_out

from iotbx.pdb import common_residue_names_get_class as get_class

class b_factor_data(dict):
  def mean(self):
    rc = {}
    for chain_id, wb in self.items():
      if not len(wb): continue
      assert chain_id not in rc
      rc[chain_id] = sum(wb.values())/len(wb)
    return rc

  def find_low(self, limit_fraction=0.25):
    mean = self.mean()
    rc = {}
    for chain_id, wb in self.items():
      tmean = mean.get(chain_id, None)
      if not tmean: continue
      limit = tmean*limit_fraction
      for i_seq, b in wb.items():
        if b<=limit:
          rc.setdefault(chain_id, [])
          rc[chain_id].append(i_seq)
    return rc

def generate_selections(hierarchy):
  for sel_str in [#"all",
                  "protein",
                  #"nucleotide",
                  # "element H or element D",
                  "water",
                  # "not (water or nucleotide or protein)",
                  ]:
    yield None, sel_str
  for ag in hierarchy.atom_groups():
    if (get_class(name=ag.resname) != "common_water"): continue
    yield ag.id_str(), 'within(5, chain %s and resseq %s and name O)' % (
      ag.parent().parent().id,
      ag.parent().resseq,
      )

def get_bs(model, sel_str):
  sel = model.selection(sel_str)
  if(sel.count(True)==0): return None
  model = model.select(sel)
  atoms = model.get_atoms()
  bs = atoms.extract_b()
  return bs

def process_water_b_factors(model):
  model.process(make_restraints=True)
  hierarchy = model.get_hierarchy()
  atoms = model.get_atoms()
  bs = atoms.extract_b()
  mean = bs.min_max_mean().mean
  ssd = bs.sample_standard_deviation()
  bsZ = bs.deep_copy()
  bsZ.as_z_scores()

  rc = {}
  for id_str, sel_str in generate_selections(hierarchy):
    bs_selected = get_bs(model, sel_str)
    mean = bs_selected.min_max_mean().mean
    ssd = bs_selected.sample_standard_deviation()
    for ag in hierarchy.atom_groups():
      if (get_class(name=ag.resname) != "common_water"): continue
      if id_str is not None:
        if ag.id_str()!=id_str: continue
      atom = ag.atoms()[0]
      bsZ_select = (atom.b-mean)/ssd
      key = '%s %5.2f' % (ag.id_str(), ag.atoms()[0].b)
      rc.setdefault(key, {})
      rc[key]['all']=bsZ[atom.i_seq]
      if sel_str.find('within')>-1:
        rc[key]['within']=bsZ_select
      else:
        rc[key][sel_str]=bsZ_select

  def _format_item(item, attr):
    outl = '%s:%5.1f' % (attr,item[attr])
    for key in item:
      if key==attr: continue
      outl += ' %s:%5.1f' % (key,item[key])
    return outl

  print(' Debugging ')
  for attr in ['within', 'water']:
    d = dict(sorted(rc.items(), key=lambda item: item[1][attr]))
    for i, (key, item) in enumerate(d.items()):
      if item.get(attr, 1e9)>-1.: break
      print('  %s :: %s' % (key,_format_item(item, attr)))
      if i>30: break

  b = b_factor_data()
  for chain in hierarchy.chains():
    b.setdefault(chain.id, {})
    for residue_group in chain.residue_groups():
      atom_group = residue_group.atom_groups()[0]
      if get_class(atom_group.resname) not in ['common_water']: continue
      atom = atom_group.atoms()[0]
      b[chain.id][atom.i_seq]=atom.b
  return b

class Program(ProgramTemplate):

  description = '''
mmtbx.water_b_factors:

Usage examples:
  mmtbx.water_b_factors model.pdb
  '''

  datatypes = ['model', 'phil']

  master_phil_str = """
  water_b_factors {
    fraction_limit = .25
      .type = float
  }
"""

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    if self.params.water_b_factors.fraction_limit>.4:
      print('\n  Fraction of mean limit seems too high : %5.2f' % (self.params.water_b_factors.fraction_limit))

  # ---------------------------------------------------------------------------
  def run(self, log=None):
    model = self.data_manager.get_model()
    hierarchy = model.get_hierarchy()
    hd_selection = model.get_hd_selection()
    model = model.select(~hd_selection)
    self.results = process_water_b_factors(model)
    #
    means = self.results.mean()
    i_seqs = self.results.find_low(self.params.water_b_factors.fraction_limit)
    #
    # transfer to results
    #
    summary=[]
    summary.append({'title':'Table of water centres with low B-factors'})
    summary.append({'paragraph' : '''
In <i>macromolecular</i> crystallographic structure refinement, water molecule
that refine to lower than average B-factors are often actually ion sites.
'''})
    summary.append({'paragraph' : 'The waters listed below are candidates.'})
    atoms = hierarchy.atoms()
    if i_seqs:
      print('-'*80)
      print('\n  Displaying %d results\n' % len(i_seqs), file=log)
    else:
      print('\n  No low waters found', file=log)
    table = [['i_seq', 'chain', 'resid', 'name', 'B-factor', 'mean', 'fraction of mean']]
    for chain_id, t_i_seqs in i_seqs.items():
      print('  Low waters in chain %s : %s mean=%5.2f' % (chain_id,
                                                          len(t_i_seqs),
                                                          means[chain_id]),
            file=log)
      for i_seq in t_i_seqs:
        atom = atoms[i_seq]
        b = atoms[i_seq].b
        fb = '%4.1f%%' % (atoms[i_seq].b/means[chain_id]*100)
        print('    %5.2f %s %s' % ( b,
                                    fb,
                                    atom.quote()),
             file=log)
        chain_id = atom.parent().parent().parent().id
        resseq = atom.parent().parent().resseq
        resname = atom.parent().resname
        mean = '%5.2f' % means[chain_id]
        table.append([i_seq, chain_id, resseq, resname, b, mean, fb])
    summary.append({'table':table})

    self.results['summary']=summary
    self.get_json()

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.results

  def get_json(self):
    import json
    rc = self.get_results()
    rc = rc.get('summary', {})
    print(json.dumps(rc))
    f=open('myfile.json', 'w')
    json.dump(rc, f, indent=2)
    del f
