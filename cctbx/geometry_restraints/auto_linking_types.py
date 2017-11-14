from __future__ import division
import os, sys

bond_origin_ids = {}
angle_origin_ids = {}
origin_ids = [bond_origin_ids,
              angle_origin_ids,
              ]

starting_id = 0
for link_info in [
    ['covalent geometry'],
    ['hydrogen bonds',
     'hydrogen bonds added both for protein SS and NA basepairs',
     '',
    ],
    ['metal coordination'],
    ['edits',
     '',
     'www.phenix-online.org/documentation/reference/refinement.html#definition-of-custom-bonds-and-angles'
    ],
    ['C-beta',
     'C-beta restraints are (only) dihedrals used by default',
     ],
    ['chi angles',
     'Torsion restraints on chi angles (side-chain rotamers)',
     ],
    ]:
  for oi in origin_ids:
    assert starting_id not in oi
    oi[starting_id] = link_info
  starting_id+=1

starting_id = 20
for link_info in [
    ['Misc. bond',
     'Bond created based on atom type and distance.',
     '',
    ],
    ['SS BOND', # short desc.
     # complete desc.
     'Disulphide bond for CYS-like sulphur atoms within 3A (default) using '
     'values determined from hi-res structures and published in CCN. '
     'Some bonds are automatically excluded based on distance from metals.',
     # citation
     'Comput. Cryst. Newsl. (2015), 6, 13-13.',
    ],
  ]:
  for oi in origin_ids:
    assert starting_id not in oi
    oi[starting_id] = link_info
  starting_id+=1

if __name__=="__main__":
  print bond_origin_ids
  print angle_origin_ids
