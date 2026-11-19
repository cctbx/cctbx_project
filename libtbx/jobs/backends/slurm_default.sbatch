#!/bin/bash
# Default libtbx.jobs SLURM template. Copy it, add your site's #SBATCH
# lines and module loads, and pass it as SlurmBackend(template=...).
# Only the command placeholder is required; the others are documented in
# libtbx/jobs/backends/batch.py. Shell syntax such as $SLURM_JOB_ID passes
# through. Keep the resource-directives placeholder below any literal
# directive for the same resource so the job's request wins.
# Comment lines must not hold placeholders.
{{resource_directives}}
{{command}}
