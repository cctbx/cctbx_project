"""
Module defines class for storing information that can be extracted from
RCSB about an entry, that being experimental or computational.

https://data.rcsb.org/index.html#data-api
https://search.rcsb.org/index.html#search-api

"""

from __future__ import absolute_import, division, print_function

import json
import requests
from mmtbx.wwpdb.rcsb_web_services import report_base_url
from libtbx.utils import Sorry

full_page_request = '''
{
  entries(entry_ids: %s) {
    rcsb_id
    entry {
      id
    }
    rcsb_entry_container_identifiers {
      entity_ids
      pubmed_id
      emdb_ids
      assembly_ids
    }
    rcsb_associated_holdings {
      rcsb_repository_holdings_current_entry_container_identifiers {
        assembly_ids
      }
      rcsb_repository_holdings_current {
        repository_content_types
      }
    }
    rcsb_comp_model_provenance {
      source_url
      source_pae_url
      entry_id
      source_db
    }
    pdbx_database_status {
      pdb_format_compatible
    }
    struct {
      title
    }
    rcsb_ma_qa_metric_global {
      model_id
      ma_qa_metric_global {
        name
        value
        type
        description
        type_other_details
      }
    }
    rcsb_primary_citation {
      id
      pdbx_database_id_DOI
    }
    pubmed {
      rcsb_pubmed_container_identifiers {
        pubmed_id
      }
      rcsb_pubmed_central_id
      rcsb_pubmed_doi
      rcsb_pubmed_abstract_text
      rcsb_pubmed_affiliation_info
    }
    pdbx_deposit_group {
      group_id
      group_type
    }
    rcsb_external_references {
      id
      type
      link
    }
    pdbx_database_PDB_obs_spr {
      id
      replace_pdb_id
    }
    struct_keywords {
      pdbx_keywords
      text
    }
    exptl {
      method
    }
    cell {
      length_a
      length_b
      length_c
      angle_alpha
      angle_beta
      angle_gamma
    }
    symmetry {
      space_group_name_H_M
    }
    software {
      classification
      name
    }
    rcsb_accession_info {
      deposit_date
      initial_release_date
      major_revision
      minor_revision
    }
    pdbx_audit_revision_history {
      ordinal
      data_content_type
      major_revision
      minor_revision
      revision_date
    }
    pdbx_audit_revision_details {
      ordinal
      revision_ordinal
      data_content_type
      type
      description
    }
    pdbx_audit_revision_group {
      ordinal
      revision_ordinal
      data_content_type
      group
    }
    pdbx_database_related {
      content_type
      db_id
      details
    }
    audit_author {
      name
    }
    pdbx_audit_support {
      funding_organization
      country
      grant_number
      ordinal
    }
    pdbx_initial_refinement_model {
      type
    }
    refine {
      pdbx_refine_id
      ls_d_res_high
      ls_R_factor_R_work
      ls_R_factor_R_free
      ls_R_factor_obs
    }
    pdbx_vrpt_summary_geometry {
      percent_ramachandran_outliers
      percent_rotamer_outliers
      clashscore
    }
    pdbx_nmr_ensemble {
      conformers_calculated_total_number
      conformers_submitted_total_number
      conformer_selection_criteria
    }
    em_experiment {
      aggregation_state
      reconstruction_method
    }
    em_3d_reconstruction {
      resolution
    }
    em_software {
      category
      name
      version
    }
    citation {
      id
      title
      rcsb_journal_abbrev
      rcsb_authors
      year
      journal_volume
      page_first
      page_last
      pdbx_database_id_PubMed
      pdbx_database_id_DOI
    }
    pdbx_database_related {
      db_id
      db_name
    }
    rcsb_entry_info {
      molecular_weight
      deposited_atom_count
      deposited_model_count
      deposited_polymer_monomer_count
      deposited_modeled_polymer_monomer_count
      deposited_unmodeled_polymer_monomer_count
      polymer_entity_count_protein
      polymer_entity_count_nucleic_acid
      polymer_entity_count_nucleic_acid_hybrid
    }
    rcsb_entry_group_membership {
      group_id
      aggregation_method
    }
    rcsb_binding_affinity {
      comp_id
      type
      value
      unit
      reference_sequence_identity
      provenance_code
      link
    }
    branched_entities {
      rcsb_id
      rcsb_branched_entity_container_identifiers {
        entry_id
        entity_id
        asym_ids
        auth_asym_ids
        reference_identifiers {
          provenance_source
          resource_name
          resource_accession
        }
      }
      prd {
        rcsb_id
        pdbx_reference_molecule {
          prd_id
          chem_comp_id
          name
          type
          class
        }
      }
      rcsb_branched_entity {
        pdbx_description
        formula_weight
      }
      pdbx_entity_branch {
        rcsb_branched_component_count
      }
      branched_entity_instances {
        rcsb_branched_entity_instance_container_identifiers {
          entry_id
          entity_id
          asym_id
          auth_asym_id
        }
        rcsb_branched_struct_conn {
          connect_type
          role
          ordinal_id
          connect_partner {
            label_asym_id
            label_seq_id
            label_comp_id
          }
          connect_target {
            auth_seq_id
            label_asym_id
            label_comp_id
          }
        }
        rcsb_branched_instance_feature {
          name
          type
          feature_value {
            comp_id
            details
          }
        }
      }
    }
    polymer_entities {
      polymer_entity_instances {
        rcsb_polymer_entity_instance_container_identifiers {
          auth_asym_id
          asym_id
          entry_id
          entity_id
        }
      }
      rcsb_polymer_entity_container_identifiers {
        entry_id
        entity_id
        asym_ids
        auth_asym_ids
        uniprot_ids
        reference_sequence_identifiers {
          database_accession
        }
      }
      uniprots {
        rcsb_id
        rcsb_uniprot_protein {
          source_organism {
            scientific_name
          }
        }
        rcsb_uniprot_external_reference {
          reference_name
          reference_id
        }
      }
      rcsb_polymer_entity {
        pdbx_description
        rcsb_ec_lineage {
          id
        }
        pdbx_ec
        rcsb_enzyme_class_combined {
          ec
          provenance_source
        }
      }
      rcsb_polymer_entity_annotation {
        type
        annotation_lineage {
          name
          depth
        }
      }
      rcsb_polymer_entity_group_membership {
        group_id
        similarity_cutoff
        aggregation_method
      }
      entity_poly {
        type
        rcsb_entity_polymer_type
        pdbx_seq_one_letter_code_can
        rcsb_sample_sequence_length
        rcsb_mutation_count
      }
      rcsb_entity_source_organism {
        scientific_name
        ncbi_scientific_name
        rcsb_gene_name {
          value
          provenance_source
        }
      }
      rcsb_entity_host_organism {
        ncbi_scientific_name
      }
      prd {
        rcsb_id
        pdbx_reference_molecule {
          prd_id
          chem_comp_id
          name
          type
          class
        }
      }
      chem_comp_nstd_monomers {
        chem_comp {
          id
          name
          formula
          type
          mon_nstd_parent_comp_id
        }
      }
      polymer_entity_instances {
        rcsb_polymer_instance_annotation {
          type
          annotation_id
        }
        rcsb_polymer_struct_conn {
          role
          connect_type
          connect_partner {
            label_asym_id
          }
          connect_target {
            label_asym_id
          }
        }
      }
    }
    nonpolymer_entities {
      rcsb_nonpolymer_entity_container_identifiers {
        entry_id
        entity_id
        auth_asym_ids
        asym_ids
        nonpolymer_comp_id
      }
      rcsb_nonpolymer_entity_annotation {
        type
      }
      nonpolymer_entity_instances {
        rcsb_nonpolymer_entity_instance_container_identifiers {
          auth_seq_id
          auth_asym_id
          asym_id
          entry_id
          entity_id
        }
        rcsb_nonpolymer_instance_validation_score {
          ranking_model_fit
          ranking_model_geometry
          average_occupancy
          is_subject_of_investigation
          is_subject_of_investigation_provenance
        }
      }
      rcsb_nonpolymer_entity {
        pdbx_description
      }
      nonpolymer_comp {
        chem_comp {
          id
          formula_weight
          name
          formula
        }
        pdbx_reference_molecule {
          prd_id
          chem_comp_id
          type
          name
          class
        }
        rcsb_chem_comp_descriptor {
          InChIKey
        }
      }
    }
    assemblies {
      rcsb_assembly_container_identifiers {
        assembly_id
      }
      pdbx_struct_assembly {
        rcsb_details
        method_details
        rcsb_candidate_assembly
      }
      pdbx_struct_assembly_auth_evidence {
        experimental_support
      }
      rcsb_struct_symmetry {
        kind
        type
        symbol
        oligomeric_state
        stoichiometry
      }
      rcsb_assembly_info {
        modeled_polymer_monomer_count
      }
    }
  }
}
'''


class rcsb_entry_info(object):
  """Class to hold information about an entry in form of json received from
  RCSB and provide parts of it in convenient form."""
  def __init__(self, json_data):
    self.data = json_data
  def __str__(self):
    return "rcsb_entry_info(%s)" % self.data['rcsb_id']
  def __repr__(self):
    return "rcsb_entry_info(%s)" % self.data['rcsb_id']
  def _get_value(self, path, convert_type=None):
    """Helper function to safely get nested dictionary values.

    Args:
        path (list): List of keys/indices to traverse the data structure
        convert_type (type, optional): Type to convert the value to (e.g., float)

    Returns:
        The value if found and successfully converted, None otherwise
    """
    try:
      value = self.data
      for key in path:
        value = value[key]
      return convert_type(value) if convert_type else value
    except (TypeError, ValueError, KeyError, IndexError):
      return None

  def get_pdb_id(self):
    return self.data['rcsb_id']

  def get_experimental_method(self):
    return self._get_value(['exptl', 0, 'method'])
  def is_xray(self):
    return self.get_experimental_method() == "X-RAY DIFFRACTION"
  def get_resolution(self):
    return self._get_value(['refine', 0, 'ls_d_res_high'], float)
  def get_rwork(self):
    return self._get_value(['refine', 0, 'ls_R_factor_R_work'], float)
  def get_rfree(self):
    return self._get_value(['refine', 0, 'ls_R_factor_R_free'], float)
  def get_rama_outliers(self):
    return self._get_value(['pdbx_vrpt_summary_geometry', 0, 'percent_ramachandran_outliers'], float)
  def get_rota_outliers(self):
    return self._get_value(['pdbx_vrpt_summary_geometry', 0, 'percent_rotamer_outliers'], float)
  def get_clashscore(self):
    return self._get_value(['pdbx_vrpt_summary_geometry', 0, 'clashscore'], float)

  #=============================================================================
  def is_computational(self):
    return self.data['rcsb_ma_qa_metric_global'] is not None
  def get_plddt(self):
    try:
      for metric_global in self.data['rcsb_ma_qa_metric_global'][0]['ma_qa_metric_global']:
        if metric_global['name'] == 'pLDDT':
          return metric_global['value']
    except TypeError:
      return None
  def get_source_url(self):
    return self._get_value(['rcsb_comp_model_provenance', 'source_url'])
  def get_source_pae_url(self):
    return self._get_value(['rcsb_comp_model_provenance', 'source_pae_url'])


def get_info(pdb_ids):
  """Get information about entries. Accepts both experimental and computational ids.

  Args:
      pdb_ids (list): PDB IDs, e.g. ["1UCS", "7P0Q", "AF_AFP12102F1", "AF_AFP35751F1"]

  Returns:
      list of rcsb_entry_info: list of objects - one object for one id.
  """
  pdb_ids_string = "%s" % pdb_ids
  pdb_ids_string = pdb_ids_string.replace("'", '"')
  r = requests.post(report_base_url, json={"query":full_page_request % pdb_ids_string})
  data_entries = r.json()['data']['entries']
  # print(json.dumps(r.json(), indent=4))

  result = [rcsb_entry_info(entry) for entry in data_entries]
  if len(result) != len(pdb_ids):
    raise Sorry("There are %d invalid pdb ids for which RCSB did not return result." % (len(pdb_ids)-len(result)))
  # for r in result:
  #   print (r.get_pdb_id(), r.is_computational(), r.get_experimental_method(), r.get_plddt())
  return result
