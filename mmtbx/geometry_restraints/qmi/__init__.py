"""qmi -- QM-based geometry-restraint generation for mononuclear metal sites.

Generates Phenix ``geometry_restraints.edits`` (bond + angle targets with
per-restraint sigmas) for mononuclear metal sites by running one
semiempirical tight-binding optimisation (GFN2-xTB) per site, once and up
front, to produce a *static* restraint file.  Refinement is unchanged and
incurs no QM cost.

Design: ``docs/QM_metal_restraints_for_Phenix_implementation_plan_v3.md``
(method/pipeline) and ``docs/metalqm_data_file_specifications.md`` (data
schemas).  The pipeline reuses the ``endoexo`` QM-region builder for cluster
extraction and the in-tree ``mmtbx.geometry_restraints.xtb_manager`` as the
xTB runner, driven directly rather than through the QMR
``quantum_restraints_manager`` pipeline (plan Sec 5E.1).
"""

from __future__ import absolute_import, division, print_function
