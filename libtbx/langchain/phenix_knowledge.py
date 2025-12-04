"""
phenix_knowledge.py
Central repository for Phenix-specific rules, hints, and strategies.
"""

# 1. The Allow-List (Prevent Hallucinations)
VALID_PHENIX_PROGRAMS = {
    "phenix.refine", "phenix.real_space_refine", "phenix.autobuild",
    "phenix.automr", "phenix.phaser", "phenix.xtriage", "phenix.maps",
    "phenix.ligandfit", "phenix.rosetta_refine", "phenix.den_refine",
    "phenix.model_vs_data", "phenix.simple_homology_model",
    "phenix.reflection_file_converter", "phenix.ready_set",
    "phenix.reduce", "phenix.merging_statistics",
    "phenix.explore_metric_symmetry", "phenix.mtz.dump",
    "phenix.predict_and_build",
    "phenix.reindex"
}

# 2. Syntax Hints (The "Tactician's" Cheat Sheet)
USAGE_HINTS = {
    "phenix.phaser": (
        "USAGE RULE: phenix.phaser <data.mtz> <model.pdb> [options]\n"
        "CRITICAL NOTES:\n"
        "1. Use POSITIONAL arguments for files (bare filenames). Do NOT use 'hklin=' or 'model_file='.\n"
        "2. Do NOT use complex 'ensemble.coords=' syntax. Just provide the PDB file directly.\n"
        "3. REQUIRED: You MUST specify `ncopies=X`.\n"
        "4. GEOMETRY WARNING: If you get 'Incompatible unit cell', STOP. Run `phenix.xtriage`."
    ),
    "phenix.xtriage": (
        "USAGE RULE: phenix.xtriage <reflections.mtz> [options]\n"
        "NOTE: The reflections file MUST be a positional argument.\n"
        "CRITICAL: This tool is for DIAGNOSIS ONLY. It cannot output re-indexed MTZ files. Do NOT use `hklout` or try to change space groups with it."
    ),
    "phenix.reflection_file_converter": (
        "USAGE RULE: phenix.reflection_file_converter <reflections.mtz> [options]\n"
        "SYNTAX EXCEPTION: This program requires double-dashes. Do NOT use bare `key=value`.\n"
        "CORRECT EXAMPLE: `phenix.reflection_file_converter data.mtz --space_group=P61 --mtz=out.mtz --label=I-obs`\n"
        "CRITICAL: You MUST specify `--label`. Look at the error log for choices."
    ),
    "phenix.explore_metric_symmetry": "USAGE RULE: phenix.explore_metric_symmetry <data.mtz> [options]\nNOTE: File must be positional.",
    "phenix.predict_and_build": (
        "USAGE RULE: phenix.predict_and_build [options] [param=value]\n"
        "USAGE RULE: be sure to add resolution=xxx as resolution is required.\n"
        "CRITICAL SCOPING RULE: This program uses a deep PHIL hierarchy. Use full paths (e.g. `input_files.xray_data_file=...`)."
    ),
    "phenix.reindex": (
        "USAGE RULE: phenix.reindex <data.mtz> space_group=P61 [options]\n"
        "CRITICAL: If input has multiple arrays, use `labels='I-obs'`."
    ),
}

# 3. Strategic Heuristics (The "Strategist's" Playbook)
KNOWN_ISSUES = {
    "Incompatible unit cell": "ACTION: If re-indexing failed, the space group is physically impossible. Revert to ORIGINAL space group and assume Twinning.",
    "Symmetry error": "ACTION: Run `phenix.xtriage` to check lattice.",
    "No search procedure defined": "ACTION: Retry `phenix.phaser`. Add `ncopies=1`.",
    "You must specify a reflections file": "ACTION: Retry `phenix.phaser` with bare filename first.",
    "Unknown space group": "ACTION: Run `phenix.xtriage` to get space group possibilities.",
    "Ambiguous parameter definition": "ACTION: Use full parameter scope.",
    "Please use --label": "ACTION: Add `--label=I-obs`."
}

