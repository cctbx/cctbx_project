from __future__ import absolute_import, division, print_function

bond_origin_ids = {}
angle_origin_ids = {}
torsion_origin_ids = {}
chiral_origin_ids = {}
plane_origin_ids = {}
parallelity_origin_ids = {}
origin_ids = [bond_origin_ids,
              angle_origin_ids,
              torsion_origin_ids,
              plane_origin_ids,
              chiral_origin_ids,
              parallelity_origin_ids,
              ]

class origins(list):
  def __init__(self, item, internals=None):
    self.internals = internals
    for i, l in enumerate(item):
      self.append(l)
      if i==3:
        self.header = l

  def __repr__(self):
    return list.__repr__(self) + ' is %s' % self.internals

covalent_headers = ['Bond',
                    "Bond angle",
                    "Dihedral angle",
                    "Chirality",
                    "Planarity",
                    "Parallelity",
                    ]

internal_labels = ['bonds',
                   'angles',
                   'dihedrals',
                   'chirals',
                   'planes',
                   'parallelities',
                   ]

starting_id = 0
for link_info in [
    # This one is outlier with only 3 elements, not clear why.
    ['covalent geometry', 'covalent geometry', [0,1,2,3,4,5]], # 0

    # Here each element represents 1 type of interaction. This interaction
    # may be restrained with various types of restraints.
    # Each one of the 5 members of the list are the following:
    # 1 - Short description of interaction, is used across the source code to figure
    #     out the origin id, e.g. origin_ids.get_origin_id('hydrogen bonds')
    # 2 - Long description - goes where?
    # 3 - citation - goes where?
    # 4 - list of .geo file headers for each involved resraint types (following)
    #     (None will suppress output)
    # 5 - list of involved restraint types - look them up in internal_labels list,
    #     e.g. 1 - angle, 3 - chirals.
    ['SS BOND', # short desc.
     # complete desc.
     'Disulphide bond for CYS-like sulphur atoms within 3A (default) using '
     'values determined from hi-res structures and published in CCN. '
     'Some bonds are automatically excluded based on distance from metals.',
     # citation
     'Comput. Cryst. Newsl. (2015), 6, 13-13.',
     # geo file header - bond, angle, dihedral (None will suppress output)
     ['Disulphide bridge']*3,
     # internals
     [0,1,2], # does not seem to be used much...
    ],

    ['hydrogen bonds',
     'hydrogen bonds added both for protein (alpha, beta) and NA (basepair) SS',
     '',
     ['Bond-like', 'Secondary Structure restraints around h-bond'],
     [0,1],
    ],

    ['metal coordination',
     '',
     '',
     ['Metal coordination']*2,
     [0,1],
    ],

    ['edits',
     '',
     'www.phenix-online.org/documentation/reference/refinement.html#definition-of-custom-bonds-and-angles',
     ['User supplied']*6,
     [0,1,2,3,4,5],
    ],
    # ['glycosidic',
    #  'Standard glycosidic CIF link blocks such as link_??? and ???',
    #  '',
    #  ['Standard Glycosidic']*5, # includes chirals!!!
    #  [0,1,2,3,4],
    # ],
    ['glycosidic custom',
     'Custom glycosidic links need to be generated when the atom names of '
     '''the carbohydrates don't conform to the standard.''',
     '',
     ['Custom Glycosidic']*5,
     [0,1,2,3,4],
    ],

    ['basepair stacking',
     'Enforces parallel between two bases in the sequence',
     'J. Appl. Cryst. 48, 1130-1141 (2015).',
     [None, None, None, None, None, 'Stacking parallelity'],
     [5],
    ],

    ['basepair parallelity',
     'Enforces parallel between two base pairs in paired bases',
     'J. Appl. Cryst. 48, 1130-1141 (2015).',
     [None, None, None, None, None, 'Basepair parallelity'],
     [5],
    ],

    ['side-chain parallelity',
     'Enforces parallel between two alt conf side-chains',
     'JH',
     [None, None, None, None, None, 'Side-chain parallelity'],
     [5],
    ],

    ['basepair planarity',
     'Enforces planarity of two base pairs in paired bases',
     'J. Appl. Cryst. 48, 1130-1141 (2015).',
     [None, None, None, 'xxx', 'Basepair planarity'],
     [3],
    ],

    # ['trans peptide link',
    # 'Applying the standard TRANS peptide link to a non-standard peptide',
    # '',
    # ['Trans Peptide']*3+[None],
    # [0,1,2,4]
    # ]

    ['Misc. bond', # 9
     'Bond created based on atom type and distance.',
     '',
     ['Misc.']*5,
     [0,1,2,3,4]
    ],

    ['User supplied cif_link', # 10
     'Internal coordinates supplied by the user in cif_link format',
     '',
     ['User cif_link']*5,
     [0,1,2,3,4]
    ],

    ['solvent network',
     'Bonds between water molecules',
     '',
     ['Solvent network'],
     [0],
    ],

    ]:
  for oi in origin_ids:
    assert starting_id not in oi
    oi[starting_id] = origins(link_info[:-1], internals=link_info[-1])
  starting_id+=1

# only angles
for link_info in []:
  angle_origin_ids[starting_id] = origins(link_info, internals=[1])
  starting_id+=1

# only dihedrals
for link_info in [
    ['C-beta',
     'C-beta restraints are (only) dihedrals used by default',
     '',
     [None, None, 'C-Beta improper'],
     ],
    ['chi angles',
     'Torsion restraints on chi angles (side-chain rotamers)',
     '',
     [None, None, 'Side chain'],
     ],
  ]:
  torsion_origin_ids[starting_id] = origins(link_info, internals=[2])
  starting_id+=1

from cctbx.geometry_restraints.standard_cif_links import standard_cif_links
for scl in standard_cif_links:
  assert starting_id not in origin_ids[0]
  origin_ids[0][starting_id] = origins(scl, internals=[0,1,2,3,4,5])
  starting_id+=1

not_covalent = [
  # not really necessary as there are no bonds, angles but added for completeness
  'basepair stacking',
  'basepair parallelity',
  'side-chain parallelity',
  'basepair planarity',
  'link_gap',
  # necessary
  'hydrogen bonds',
  ]

def iterate_covalent():
  for key, item in origin_ids[0].items():
    if item[0] in not_covalent: continue
    # print('YIELD "%s" "%s"' % (key, item))
    yield key

if __name__=="__main__":
  print('-'*80)
  print(bond_origin_ids)
  print('-'*80)
  print(angle_origin_ids)
  print('-'*80)
  print(torsion_origin_ids)
  print('-'*80)
  print(parallelity_origin_ids)
