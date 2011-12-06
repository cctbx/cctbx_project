from mmtbx.monomer_library.tst_rna_dna_interpretation import run
import sys

expected_results = {
"cys_chain_all_h_1so9_v3":
  ['CYS', 11, [], 'peptide', None, None],
"phe_nterm_all_h_1g7a_v3":
  ['PHE%NH3', 22, [], 'peptide', True, None],
"ile_nterm_all_h_1v1c_v3":
  ['ILE%NH3', 21, [], 'peptide', True, None],
"ile_nterm_all_h_1ghx_v2":
  ['ILE%NH3', 21, [], 'peptide', True, None],
"ile_nterm_all_h_1v1c_v2":
  ['ILE%NH3', 21, [], 'peptide', True, None],
"met_cterm_all_h_1hz3_v2":
  ['MET%COOH', 19, [], 'peptide', True, None],
"his_chain_all_h_1g7e_v2":
  ['HIS', 18, [], 'peptide', None, None],
"arg_nterm_all_h_5pti_v2":
  ['ARG%NH3', 26, [], 'peptide', True, None],
"trp_chain_all_h_1hdp_v2":
  ['TRP', 24, [], 'peptide', None, None],
"asp_chain_all_h_2izb_v2":
  ['ASP', 12, [' HD1'], 'peptide', None, None],
"val_nterm_all_h_1qn1_v3":
  ['VAL%NH3', 18, [], 'peptide', True, None],
"lys_nterm_all_h_1gcc_v2":
  ['LYS%NH3', 24, [], 'peptide', True, None],
"leu_chain_all_h_1ozo_v3":
  ['LEU', 19, [], 'peptide', None, None],
"gly_nterm_all_h_1qg1_v3":
  ['GLY%NH3', 9, [], 'peptide', True, None],
"phe_nterm_all_h_1e0e_v3":
  ['PHE%NH3', 22, [], 'peptide', True, None],
"tyr_nterm_all_h_1d1n_v2":
  ['TYR%NH3', 23, [], 'peptide', True, None],
"glu_nterm_all_h_2bic_v2":
  ['GLU%ACID-GLU%NH3', 18, [], 'peptide', True, None],
"ile_cterm_all_h_1pog_v3":
  ['ILE%COOH', 21, [], 'peptide', True, None],
"met_nterm_all_h_1qqi_v2":
  ['MET%NH3', 19, [], 'peptide', True, None],
"ser_chain_all_h_1o8t_v2":
  ['SER', 11, [], 'peptide', None, None],
"ala_nterm_all_h_1blq_v2":
  ['ALA%NH3', 12, [], 'peptide', True, None],
"phe_cterm_all_h_2jnr_v2":
  ['PHE%COOH', 22, [], 'peptide', True, None],
"val_nterm_all_h_1ylb_v2":
  ['VAL%NH3', 18, [], 'peptide', True, None],
"met_chain_all_h_1ozo_v3":
  ['MET', 17, [], 'peptide', None, None],
"asn_chain_all_h_1o8t_v3":
  ['ASN', 14, [], 'peptide', None, None],
"tyr_chain_all_h_1ozo_v3":
  ['TYR', 21, [], 'peptide', None, None],
"cys_chain_all_h_1so9_v2":
  ['CYS', 11, [], 'peptide', None, None],
"thr_nterm_all_h_1o8t_v3":
  ['THR%NH3', 16, [], 'peptide', True, None],
"ser_nterm_all_h_1gkt_v3":
  ['SER%NH3', 13, [], 'peptide', True, None],
"arg_nterm_all_h_1dpu_v3":
  ['ARG%NH3', 26, [], 'peptide', True, None],
"val_cterm_all_h_1cz2_v2":
  ['VAL%COOH', 18, [], 'peptide', True, None],
"gly_nterm_all_h_3ins_v2":
  ['GLY%NH3', 9, [], 'peptide', True, None],
"met_nterm_all_h_2j46_v3":
  ['MET%NH3', 19, [], 'peptide', True, None],
"his_nterm_all_h_2j5p_v2":
  ['HIS%NH3', 20, [], 'peptide', True, None],
"thr_nterm_all_h_2g3q_v3":
  ['THR%NH3', 16, [], 'peptide', True, None],
"tyr_chain_all_h_1cx1_v3":
  ['TYR', 21, [], 'peptide', None, None],
"ala_nterm_all_h_1vcx_v3":
  ['ALA%NH3', 12, [], 'peptide', True, None],
"lys_nterm_all_h_2bvb_v2":
  ['LYS%NH3', 24, [], 'peptide', True, None],
"asn_cterm_all_h_2cdx_v2":
  ['ASN%COOH', 16, [], 'peptide', True, None],
"lys_chain_all_h_1o8t_v3":
  ['LYS', 22, [], 'peptide', None, None],
"thr_nterm_all_h_2g3q_v2":
  ['THR%NH3', 16, [], 'peptide', True, None],
"tyr_nterm_all_h_1qo6_v2":
  ['TYR%NH3', 23, [], 'peptide', True, None],
"arg_chain_all_h_1o8t_v3":
  ['ARG', 24, [], 'peptide', None, None],
"trp_nterm_all_h_1fi6_v2":
  ['TRP%NH3', 26, [], 'peptide', True, None],
"leu_nterm_all_h_2bun_v3":
  ['LEU%NH3', 21, [], 'peptide', True, None],
"trp_chain_all_h_1cx1_v2":
  ['TRP', 24, [], 'peptide', None, None],
"ala_nterm_all_h_2bw2_v3":
  ['ALA%NH3', 12, [], 'peptide', True, None],
"phe_chain_all_h_1hdj_v2":
  ['PHE', 20, [], 'peptide', None, None],
"thr_chain_all_h_1o8t_v3":
  ['THR', 14, [], 'peptide', None, None],
"met_nterm_all_h_1d7q_v3":
  ['MET%NH3', 19, [], 'peptide', True, None],
"pro_nterm_all_h_2lhb_v2":
  ['PRO%NH2', 16, [], 'peptide', True, None],
"lys_chain_all_h_1o8t_v2":
  ['LYS', 22, [], 'peptide', None, None],
"ala_chain_all_h_1ozo_v2":
  ['ALA', 10, [], 'peptide', None, None],
"ala_cterm_all_h_1eio_v2":
  ['ALA%COOH', 12, [], 'peptide', True, None],
"asp_cterm_all_h_1kb8_v2":
  ['ASP%COOH%ACID-ASP', 15, [], 'peptide', True, None],
"ser_nterm_all_h_1gkt_v2":
  ['SER%NH3', 13, [], 'peptide', True, None],
"ile_nterm_all_h_1ghx_v3":
  ['ILE%NH3', 21, [], 'peptide', True, None],
"val_chain_all_h_1ozo_v3":
  ['VAL', 16, [], 'peptide', None, None],
"ser_nterm_all_h_1goe_v3":
  ['SER%NH3', 13, [], 'peptide', True, None],
"cys_nterm_all_h_1b9q_v2":
  ['CYS%NH3', 13, [], 'peptide', True, None],
"ile_nterm_all_h_1tkw_v2":
  ['ILE%NH3', 21, [], 'peptide', True, None],
"ile_nterm_all_h_1c2l_v2":
  ['ILE%NH2NOTPRO', 21, [], 'peptide', None, None],
"arg_nterm_all_h_1i9f_v2":
  ['ARG%NH2NOTPRO', 26, [], 'peptide', None, None],
"phe_nterm_all_h_3ins_v2":
  ['PHE%NH3', 22, [], 'peptide', True, None],
"asp_nterm_all_h_1x9v_v2":
  ['ASP%ACID-ASP%NH3', 15, [], 'peptide', True, None],
"leu_nterm_all_h_1edv_v3":
  ['LEU%NH3', 21, [], 'peptide', True, None],
"val_nterm_all_h_2mb5_v2":
  ['VAL%NH3', 18, [], 'peptide', True, None],
"trp_nterm_all_h_1c95_v3":
  ['TRP%NH3', 26, [], 'peptide', True, None],
"asp_chain_all_h_1jjx_v2":
  ['ASP', 12, ['2HD '], 'peptide', None, None],
"pro_nterm_all_h_2gs0_v3":
  ['PRO%NH2', 16, [], 'peptide', True, None],
"ala_nterm_all_h_1dgn_v2":
  ['ALA%NH3', 12, [], 'peptide', True, None],
"phe_nterm_all_h_3ins_v3":
  ['PHE%NH3', 22, [], 'peptide', True, None],
"glu_chain_all_h_1jjx_v2":
  ['GLU', 15, ['2HE '], 'peptide', None, None],
"asp_chain_all_h_1kgl_v3":
  ['ASP%ACID-ASP', 13, [], 'peptide', None, None],
"asp_chain_all_h_1jjx_v3":
  ['ASP%ACID-ASP', 13, [], 'peptide', None, None],
"gln_nterm_all_h_1w4h_v2":
  ['GLN%NH3', 19, [], 'peptide', True, None],
"his_nterm_all_h_1c17_v3":
  ['HIS%NH3', 20, [], 'peptide', True, None],
"ala_chain_all_h_1ozo_v3":
  ['ALA', 10, [], 'peptide', None, None],
"asn_chain_all_h_1o8t_v2":
  ['ASN', 14, [], 'peptide', None, None],
"glu_nterm_all_h_2bic_v3":
  ['GLU%ACID-GLU%NH3', 18, [], 'peptide', True, None],
"gln_chain_all_h_1o8t_v2":
  ['GLN', 17, [], 'peptide', None, None],
"arg_chain_all_h_1o8t_v2":
  ['ARG', 24, [], 'peptide', None, None],
"glu_chain_all_h_1bm4_v3":
  ['GLU%ACID-GLU', 16, [], 'peptide', None, None],
"phe_nterm_all_h_1e0e_v2":
  ['PHE%NH3', 22, [], 'peptide', True, None],
"cys_chain_all_h_1rfa_v2":
  ['CYS', 11, [], 'peptide', None, None],
"met_nterm_all_h_1qqi_v3":
  ['MET%NH3', 19, [], 'peptide', True, None],
"ile_nterm_all_h_1ntp_v3":
  ['ILE%NH3', 21, [], 'peptide', True, None],
"trp_nterm_all_h_1c95_v2":
  ['TRP%NH2NOTPRO', 26, [], 'peptide', None, None],
"ile_nterm_all_h_1p57_v2":
  ['ILE%NH2NOTPRO', 21, [], 'peptide', None, None],
"asn_nterm_all_h_1cej_v2":
  ['ASN%NH3', 16, [], 'peptide', True, None],
"leu_chain_all_h_2sn3_v3":
  ['LEU', 19, [], 'peptide', None, None],
"arg_nterm_all_h_1gxe_v2":
  ['ARG%NH3', 26, [], 'peptide', True, None],
"cys_chain_all_h_1rfa_v3":
  ['CYS', 11, [], 'peptide', None, None],
"asp_nterm_all_h_1x9v_v3":
  ['ASP%ACID-ASP%NH3', 15, [], 'peptide', True, None],
"phe_nterm_all_h_1g7a_v2":
  ['PHE%NH3', 22, [], 'peptide', True, None],
"his_chain_all_h_1g7e_v3":
  ['HIS', 18, [], 'peptide', None, None],
"trp_nterm_all_h_1haj_v2":
  ['TRP%NH3', 26, [], 'peptide', True, None],
"glu_chain_all_h_1jjx_v3":
  ['GLU%ACID-GLU', 16, [], 'peptide', None, None],
"thr_nterm_all_h_1o8t_v2":
  ['THR%NH3', 16, [], 'peptide', True, None],
"ser_nterm_all_h_1goe_v2":
  ['SER%NH3', 13, [], 'peptide', True, None],
"gln_nterm_all_h_1fmy_v2":
  ['GLN%NH3', 19, [], 'peptide', True, None],
"pro_chain_all_h_1a03_v2":
  ['PRO', 14, [], 'peptide', None, None],
"ala_nterm_all_h_1dgn_v3":
  ['ALA%NH3', 12, [], 'peptide', True, None],
"glu_chain_all_h_1tmr_v3":
  ['GLU%ACID-GLU', 16, [], 'peptide', None, None],
"lys_cterm_all_h_1aze_v3":
  ['LYS%COOH', 24, [], 'peptide', True, None],
"glu_chain_all_h_1bm4_v2":
  ['GLU%ACID-GLU', 16, [], 'peptide', None, None],
"trp_chain_all_h_1hdp_v3":
  ['TRP', 24, [], 'peptide', None, None],
"tyr_nterm_all_h_1qo6_v3":
  ['TYR%NH3', 23, [], 'peptide', True, None],
"gly_nterm_all_h_1o7b_v3":
  ['GLY%NH3', 9, [], 'peptide', True, None],
"tyr_nterm_all_h_1lxh_v2":
  ['TYR%NH3', 23, [], 'peptide', True, None],
"lys_nterm_all_h_1io5_v2":
  ['LYS%NH3', 24, [], 'peptide', True, None],
"ser_chain_all_h_1o8t_v3":
  ['SER', 11, [], 'peptide', None, None],
"arg_cterm_all_h_1mmc_v3":
  ['ARG%COOH', 26, [], 'peptide', True, None],
"pro_nterm_all_h_2lhb_v3":
  ['PRO%NH2', 16, [], 'peptide', True, None],
"tyr_cterm_all_h_1fu6_v3":
  ['TYR%COOH', 23, [], 'peptide', True, None],
"val_chain_all_h_1ozo_v2":
  ['VAL', 16, [], 'peptide', None, None],
"ile_nterm_all_h_1p57_v3":
  ['ILE%NH3', 21, [], 'peptide', True, None],
"his_nterm_all_h_2j5p_v3":
  ['HIS%NH3', 20, [], 'peptide', True, None],
"asp_chain_all_h_2izb_v3":
  ['ASP%ACID-ASP', 13, [], 'peptide', None, None],
"lys_cterm_all_h_1gkt_v2":
  ['LYS%COOH', 24, [], 'peptide', True, None],
"gln_chain_all_h_1o8t_v3":
  ['GLN', 17, [], 'peptide', None, None],
"ser_cterm_all_h_2jnr_v3":
  ['SER%COOH', 13, [], 'peptide', True, None],
"lys_nterm_all_h_1io5_v3":
  ['LYS%NH3', 24, [], 'peptide', True, None],
"ile_nterm_all_h_1tkw_v3":
  ['ILE%NH3', 21, [], 'peptide', True, None],
"gln_nterm_all_h_1fmy_v3":
  ['GLN%NH3', 19, [], 'peptide', True, None],
"met_cterm_all_h_1hz3_v3":
  ['MET%COOH', 19, [], 'peptide', True, None],
"lys_cterm_all_h_1aze_v2":
  ['LYS%COOH', 24, [], 'peptide', True, None],
"val_nterm_all_h_2mb5_v3":
  ['VAL%NH3', 18, [], 'peptide', True, None],
"val_nterm_all_h_1ylb_v3":
  ['VAL%NH3', 18, [], 'peptide', True, None],
"trp_nterm_all_h_1haj_v3":
  ['TRP%NH3', 26, [], 'peptide', True, None],
"tyr_chain_all_h_1ozo_v2":
  ['TYR', 21, [], 'peptide', None, None],
"met_nterm_all_h_1euw_v2":
  ['MET%NH3', 19, [], 'peptide', True, None],
"gln_nterm_all_h_1w4h_v3":
  ['GLN%NH3', 19, [], 'peptide', True, None],
"arg_nterm_all_h_1i9f_v3":
  ['ARG%NH3', 26, [], 'peptide', True, None],
"asp_cterm_all_h_1kb8_v3":
  ['ASP%COOH%ACID-ASP', 15, [], 'peptide', True, None],
"met_chain_all_h_1ozo_v2":
  ['MET', 17, [], 'peptide', None, None],
"trp_nterm_all_h_1fi6_v3":
  ['TRP%NH3', 26, [], 'peptide', True, None],
"tyr_nterm_all_h_1d1n_v3":
  ['TYR%NH3', 23, [], 'peptide', True, None],
"arg_nterm_all_h_1gxe_v3":
  ['ARG%NH3', 26, [], 'peptide', True, None],
"gly_nterm_all_h_1o7b_v2":
  ['GLY%NH3', 9, [], 'peptide', True, None],
"cys_nterm_all_h_1b9q_v3":
  ['CYS%NH3', 13, [], 'peptide', True, None],
"leu_cterm_all_h_2oyw_v3":
  ['LEU%COOH', 21, [], 'peptide', True, None],
"trp_chain_all_h_1cx1_v3":
  ['TRP', 24, [], 'peptide', None, None],
"ala_cterm_all_h_1jdk_v2":
  ['ALA%COOH', 12, [], 'peptide', True, None],
"arg_nterm_all_h_5pti_v3":
  ['ARG%NH3', 26, [], 'peptide', True, None],
"gly_chain_all_h_1ozo_v3":
  ['GLY', 7, [], 'peptide', None, None],
"met_nterm_all_h_1c15_v3":
  ['MET%NH3', 19, [], 'peptide', True, None],
"tyr_nterm_all_h_1lxh_v3":
  ['TYR%NH3', 23, [], 'peptide', True, None],
"ser_nterm_all_h_1bl1_v3":
  ['SER%NH3', 13, [], 'peptide', True, None],
"phe_chain_all_h_1hdj_v3":
  ['PHE', 20, [], 'peptide', None, None],
"phe_chain_all_h_1hdp_v3":
  ['PHE', 20, [], 'peptide', None, None],
"arg_cterm_all_h_1mmc_v2":
  ['ARG%COOH', 26, [], 'peptide', True, None],
"ile_nterm_all_h_1ntp_v2":
  ['ILE%NH3', 21, [], 'peptide', True, None],
"met_nterm_all_h_2j46_v2":
  ['MET%NH3', 19, [], 'peptide', True, None],
"asn_cterm_all_h_2cdx_v3":
  ['ASN%COOH', 16, [], 'peptide', True, None],
"lys_nterm_all_h_2bvb_v3":
  ['LYS%NH3', 24, [], 'peptide', True, None],
"glu_nterm_all_h_1agg_v2":
  ['GLU%ACID-GLU%NH3', 18, [], 'peptide', True, None],
"glu_nterm_all_h_1agg_v3":
  ['GLU%ACID-GLU%NH3', 18, [], 'peptide', True, None],
"ala_nterm_all_h_1vcx_v2":
  ['ALA%NH3', 12, [], 'peptide', True, None],
"ala_cterm_all_h_1eio_v3":
  ['ALA%COOH', 12, [], 'peptide', True, None],
"ser_nterm_all_h_1bl1_v2":
  ['SER%NH3', 13, [], 'peptide', True, None],
"ile_chain_all_h_1ozo_v3":
  ['ILE', 19, [], 'peptide', None, None],
"gly_nterm_all_h_1qg1_v2":
  ['GLY%NH3', 9, [], 'peptide', True, None],
"arg_nterm_all_h_1dpu_v2":
  ['ARG%NH3', 26, [], 'peptide', True, None],
"asn_nterm_all_h_1w4g_v3":
  ['ASN%NH3', 16, [], 'peptide', True, None],
"met_nterm_all_h_1c15_v2":
  ['MET%NH3', 19, [], 'peptide', True, None],
"trp_cterm_all_h_6cmh_v2":
  ['TRP%COOH', 26, [], 'peptide', True, None],
"leu_chain_all_h_2sn3_v2":
  ['LEU', 19, [], 'peptide', None, None],
"leu_nterm_all_h_2bun_v2":
  ['LEU%NH3', 21, [], 'peptide', True, None],
"ala_nterm_all_h_2bw2_v2":
  ['ALA%NH3', 12, [], 'peptide', True, None],
"met_nterm_all_h_1d7q_v2":
  ['MET%NH2NOTPRO', 19, [], 'peptide', None, None],
"val_cterm_all_h_1cz2_v3":
  ['VAL%COOH', 18, [], 'peptide', True, None],
"asn_nterm_all_h_1w4g_v2":
  ['ASN%NH3', 16, [], 'peptide', True, None],
"ala_nterm_all_h_1blq_v3":
  ['ALA%NH3', 12, [], 'peptide', True, None],
"pro_nterm_all_h_2gs0_v2":
  ['PRO%NH2', 16, [], 'peptide', True, None],
"glu_chain_all_h_1tmr_v2":
  ['GLU', 15, [' HE1'], 'peptide', None, None],
"gly_cterm_all_h_1sbu_v3":
  ['GLY%COOH', 9, [], 'peptide', True, None],
"leu_cterm_all_h_2oyw_v2":
  ['LEU%COOH', 21, [], 'peptide', True, None],
"met_nterm_all_h_1euw_v3":
  ['MET%NH3', 19, [], 'peptide', True, None],
"gly_nterm_all_h_3ins_v3":
  ['GLY%NH3', 9, [], 'peptide', True, None],
"asn_nterm_all_h_1cej_v3":
  ['ASN%NH3', 16, [], 'peptide', True, None],
"gly_chain_all_h_1ozo_v2":
  ['GLY', 7, [], 'peptide', None, None],
"ala_cterm_all_h_1jdk_v3":
  ['ALA%COOH', 12, [], 'peptide', True, None],
"leu_chain_all_h_1ozo_v2":
  ['LEU', 19, [], 'peptide', None, None],
"thr_chain_all_h_1o8t_v2":
  ['THR', 14, [], 'peptide', None, None],
"his_nterm_all_h_1c17_v2":
  ['HIS%NH3', 20, [], 'peptide', True, None],
"gly_cterm_all_h_1sbu_v2":
  ['GLY%COOH', 9, [], 'peptide', True, None],
"tyr_cterm_all_h_1fu6_v2":
  ['TYR%COOH', 23, [], 'peptide', True, None],
"tyr_chain_all_h_1cx1_v2":
  ['TYR', 21, [], 'peptide', None, None],
"lys_cterm_all_h_1gkt_v3":
  ['LYS%COOH', 24, [], 'peptide', True, None],
"ser_cterm_all_h_2jnr_v2":
  ['SER%COOH', 13, [], 'peptide', True, None],
"trp_cterm_all_h_6cmh_v3":
  ['TRP%COOH', 26, [], 'peptide', True, None],
"phe_cterm_all_h_2jnr_v3":
  ['PHE%COOH', 22, [], 'peptide', True, None],
"leu_nterm_all_h_1edv_v2":
  ['LEU%NH3', 21, [], 'peptide', True, None],
"asp_chain_all_h_1kgl_v2":
  ['ASP%ACID-ASP', 13, [], 'peptide', None, None],
"lys_nterm_all_h_1gcc_v3":
  ['LYS%NH3', 24, [], 'peptide', True, None],
"ile_cterm_all_h_1pog_v2":
  ['ILE%COOH', 21, [], 'peptide', True, None],
"phe_chain_all_h_1hdp_v2":
  ['PHE', 20, [], 'peptide', None, None],
"pro_chain_all_h_1a03_v3":
  ['PRO', 14, [], 'peptide', None, None],
"ile_chain_all_h_1ozo_v2":
  ['ILE', 19, [], 'peptide', None, None],
"ile_nterm_all_h_1c2l_v3":
  ['ILE%NH3', 21, [], 'peptide', True, None],
"val_nterm_all_h_1qn1_v2":
  ['VAL%NH3', 18, [], 'peptide', True, None],
"phe_n_ter_c_ter_with_hydrogens":
  ['PHE%COO%NH3', 23, [], 'peptide', True, None],
"mse_cterm_all_h_1hz3_v2":
  ['MSE%COOH', 19, [], 'peptide', True, None],
"mse_nterm_all_h_1qqi_v2":
  ['MSE%NH3', 19, [], 'peptide', True, None],
"mse_chain_all_h_1ozo_v3":
  ['MSE', 17, [], 'peptide', None, None],
"mse_nterm_all_h_2j46_v3":
  ['MSE%NH3', 19, [], 'peptide', True, None],
"mse_nterm_all_h_1d7q_v3":
  ['MSE%NH3', 19, [], 'peptide', True, None],
"mse_nterm_all_h_1qqi_v3":
  ['MSE%NH3', 19, [], 'peptide', True, None],
"mse_cterm_all_h_1hz3_v3":
  ['MSE%COOH', 19, [], 'peptide', True, None],
"mse_nterm_all_h_1euw_v2":
  ['MSE%NH3', 19, [], 'peptide', True, None],
"mse_chain_all_h_1ozo_v2":
  ['MSE', 17, [], 'peptide', None, None],
"mse_nterm_all_h_1c15_v3":
  ['MSE%NH3', 19, [], 'peptide', True, None],
"mse_nterm_all_h_2j46_v2":
  ['MSE%NH3', 19, [], 'peptide', True, None],
"mse_nterm_all_h_1c15_v2":
  ['MSE%NH3', 19, [], 'peptide', True, None],
"mse_nterm_all_h_1d7q_v2":
  ['MSE%NH2NOTPRO', 19, [], 'peptide', None, None],
"mse_nterm_all_h_1euw_v3":
  ['MSE%NH3', 19, [], 'peptide', True, None],
}

if (__name__ == "__main__"):
  run(sys.argv[1:], residue_type="protein", expected_results=expected_results)
