# UNPARSED INVENTORY — corpus2/work n=253

Regenerated from the corpus. **Reported, not tuned to zero.** A run
producing many unparsed records says the corpus has outgrown the
patterns; that signal is the tool's main defence against silent decay.

Lines: 416987. Screened candidates: 107505 (25.8%).

## The three numbers

| status | count | meaning |
|---|---|---|
| unclaimed | 54449 | screened, no parser claimed it |
| generic_only | 35083 | recorded as a labeled value, understood as nothing in particular |
| rule_excluded | 9347 | refused by a NAMED rule |

## Unclaimed, by screen rule

- `rule_line` 9273
- `key_value` 0
- `verb` 82
- `numeric_row` 45094

## Refusals, by rule

- `label_exceeds_distinct_value_limit` 4196
- `section_title_too_long` 1427
- `rule_is_a_banner_not_an_underline` 1342
- `section_title_is_a_table_header` 805
- `setting_states_no_reason` 585
- `setting_line_assigns_nothing` 310
- `skip_matches_neither_form` 234
- `setting_up_is_not_an_assignment` 180
- `section_title_is_a_numeric_row` 170
- `stage_row_is_the_table_summary` 70
- `cycle_line_carries_no_metrics` 28

## The 20 commonest unclaimed shapes (digits masked)

Where the corpus has outgrown the patterns. Nothing here is a defect;
each is a shape no parser was written for.

      3442  #       #       #   #    #   #
      2988  ========================================================
      1795  #      #      #   #    #   #
      1764  Working on segment # with # residues (forward) and fixed
      1758  Working on segment # with # residues (reverse) and fixed
      1723  ********************************************************
      1680  --------------------------------------------------------
      1444  #     #   #   #   #   #   #   #   #   #   #
      1200  | #                          | #                      | 
      1190  r(curr,min,max,mean)=# # # # # #
      1100  #      #       #   #    #   #
       875  #       #       #  #    #   #
       798  #   #   #   #   #   #        #   #   #   #
       792  #    #      #        #            None          None
       754  #   #    #     #   #   #   #   #   #
       558  |  #:    # -    #  #   #  #  #   #     #      #|
       544  #   #  #    #  #  #  #   #  #
       532  ##    #  #   # #   #   # # # #
       518  |  #:    # -    # #   #  # # #      #      #|
       416  #      #  #  #    #  #

## Unclaimed per program, worst 15

- **map_to_model** unclaimed 9488, generic_only 5629, refused 3338
- **phaser** unclaimed 9291, generic_only 3164, refused 828
- **map_sharpening** unclaimed 6104, generic_only 248, refused 11
- **resolve_cryo_em** unclaimed 4363, generic_only 685, refused 72
- **refine** unclaimed 4288, generic_only 3210, refused 124
- **xtriage** unclaimed 3511, generic_only 1228, refused 648
- **autobuild** unclaimed 3260, generic_only 5275, refused 828
- **autosol** unclaimed 3067, generic_only 6314, refused 832
- **local_aniso_sharpen** unclaimed 3022, generic_only 59, refused 0
- **find_reference_1utn** unclaimed 1218, generic_only 20, refused 1
- **real_space_refine** unclaimed 1184, generic_only 2988, refused 1373
- **composite_omit_map** unclaimed 886, generic_only 32, refused 0
- **local_resolution** unclaimed 797, generic_only 18, refused 0
- **predict_and_build** unclaimed 676, generic_only 1845, refused 286
- **validation_cryoem** unclaimed 413, generic_only 424, refused 84
