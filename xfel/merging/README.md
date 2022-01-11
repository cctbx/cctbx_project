`cctbx.xfel.merge` is the merging program for serial crystallographic data in cctbx.xfel.

# Introduction

`cctbx.xfel.merge` uses a series of workers to process data, with each worker executing a different task, modifying, filtering, or merging the data. The workers are ran sequentially to read the data, scale it, and merge it, while computing relevant statistics. The default set of workers is:

- input: read data from disk in [DIALS format](https://dials.github.io/documentation/data_files.html)
- balance: balance input load. Needed if the number of input files is fewer than the number of MPI ranks being used
- model_scaling: build full Miller list, model intensities, and resolution binner - for scaling and post-refinement
- modify: apply polarization correction, etc.
- filter: reject whole experiments or individual reflections by resolution, unit cell, etc.
- scale: scale the data to a reference dataset
- postrefine: apply post-refinement (see Sauter 2015)
- statistics_unitcell: unit cell averaging and statistics
- statistics_beam: wavelength averaging and statistics
- model_statistics: build full Miller list, model intensities, and resolution binner - for statistics. Can use average unit cell or a reference
- statistics_resolution: calculate resolution statistics per crystal
- group: re-distribute reflections over the ranks, so that all measurements of every HKL are gathered at the same rank, prior to merging
- errors_merge: correct errors using a per-HKL algorithm, e.g. errors_from_sample_residuals, and Ev11 (see Brewster 2019b)
- statistics_intensity: calculate resolution statistics for intensities
- merge: merge HKL intensities, output "odd", "even" and "all" HKLs as mtz files
- statistics_intensity_cxi: run CC1/2 and related statistics

The list of workers can be changed with the parameter dispatch.step_list

# Invocation

`cctbx.xfel.merge` is designed to be ran multi-processed using MPI, so a typical invocation looks like

```
mpirun cctbx.xfel.merge merge.phil
```

Where merge.phil looks like

```
input.path=<path to data>
input.experiments_suffix=.expt
input.reflections_suffix=.refl
filter.outlier.min_corr=0.1
filter.algorithm=unit_cell
filter.unit_cell.value.relative_length_tolerance=0.03
select.algorithm=significance_filter
select.significance_filter.sigma=0.1
select.significance_filter.min_ct=200
select.significance_filter.max_ct=300
scaling.model=<model.pdb or model.mtz>
scaling.resolution_scalar=0.96
merging.d_min=1.8
merging.error.model=ev11
merging.merge_anomalous=True
postrefinement.enable=True
statistics.n_bins=20
output.do_timing=True
output.output_dir=<path to folder for output data>
output.prefix<text to prepend to output files>
```

## Parameter description:
- input.path: can be a single folder, several folders, or individual paths. Multiple entries allowed and wildcards are allowed
- input.experiments_suffix/reflections_suffix: if there are multiple kinds of results in a folder, use this to specify which files to use (eg. _integrated.expt, _integrated.refl)
- filter.outlier.min_corr: minimum correlation of each image to the reference. Default is 0.1. -1 means accept all images regardless of correlation to the reference.
- filter.unit_cell.value.relative_length_tolerance, this parameter controls the unit cell filter.  Stricter or looser parameters may be needed.
- select.algorithm: specify the significance_filter to filter each image at a resolution where the I/sigI cutoff dips below the select.significance_filter.sigma value specified
- select.significance_filter.min_ct/max_ct: these values are for larger unit cells.  For small unit cells, leave out these values and use the defaults.  This is an unresolved issue in the program
- scaling.model: specify the reference dataet (pdb or mtz)
- scaling.resolution_scalar: a factor that controls how many more reference reflections are generated from a pdb file than the given d_min would generate. Due to per-shot variation, more reflections are needed. The default is 0.969, so 0.96 generates a few more.
- merging.d_min: specify the resolution of interest
- merging.error.model: ev11 is the recommended error model (see Brewster 2019b)
- merging.merge_anomalous: True to merge anomalous pairs together (higher multiplicity), False to keep them separate (lower multiplicity, but preserves anomalous signal)

Additionally, the merged pdb/mtz files can be automatically uploaded to a [Google drive folder](https://github.com/cctbx/cctbx_project/tree/master/xfel/merging/application/publish)

Further documentation can be obtatined with `cctbx.xfel.merge -c`. Include `-e 10` for all parameters and `-a 2` for additional documentation (e.g. `cctbx.xfel.merge -c -e 10 -a 2`).

# Scaling and merging seperately

The cctbx.xfel GUI (Brewster 2019a) scales and post-refines each set of data seperately prior to merging all data together. This increase efficiency as it allows for each image to be only scaled and post-refined once, while as the dataset grows during data collection, the merging can be done repeatedly as it grows.

### Scaling

For each subset of the data to be scaled and post-refined, invoke the program first with these parmeters:

```
dispatch.step_list=input balance model_scaling modify filter scale postrefine statistics_unitcell statistics_beam model_statistics statistics_resolution
input.parallel_file_load.method=uniform
<parameters specific to these workers>
output.save_experiments_and_reflections=True
```

The last parameter outputs the modified data after that last worker as DIALS intermediate files (experiment lists and refleciton tables).  These can then be used by the merge step.

### Merging

Use these parameters to merge the scaled data, adding additional input.path statements as new scaling results finish. Again, this is automatically performed by the XFEL GUI.

```
dispatch.step_list=input model_scaling statistics_unitcell statistics_beam model_statistics statistics_resolution group errors_merge statistics_intensity merge statistics_intensity_cxi
input.parallel_file_load.method=uniform
<parameters specific to these workers>
```

### Complete example

Here is a complete example for this use case

Scaling:

```
input.path=<path to data subset>
input.experiments_suffix=.expt
input.reflections_suffix=.refl
dispatch.step_list=input balance model_scaling modify filter scale postrefine statistics_unitcell statistics_beam model_statistics statistics_resolution
input.parallel_file_load.method=uniform
filter.outlier.min_corr=-1
filter.algorithm=unit_cell
filter.unit_cell.value.relative_length_tolerance=0.03
select.algorithm=significance_filter
select.significance_filter.sigma=0.1
select.significance_filter.min_ct=200
select.significance_filter.max_ct=300
scaling.model=<model or mtz file>
scaling.resolution_scalar=0.96
merging.d_min=1.8
merging.merge_anomalous=True
postrefinement.enable=True
statistics.n_bins=20
output.do_timing=True
output.save_experiments_and_reflections=True
output.do_timing=True
output.output_dir=<path to folder for output data>
output.prefix<text to prepend to output files>
```

Merging:

```
input.path=<path to scaling output for subset 1>
input.path=<path to scaling output for subset 2, etc.>
input.experiments_suffix=.expt
input.reflections_suffix=.refl
dispatch.step_list=input model_scaling statistics_unitcell statistics_beam model_statistics statistics_resolution group errors_merge statistics_intensity merge statistics_intensity_cxi
input.parallel_file_load.method=uniform
scaling.model=<model or mtz file>
scaling.resolution_scalar=0.96
statistics.n_bins=20
merging.d_min=1.8
merging.merge_anomalous=True
merging.error.model=ev11
output.do_timing=True
output.output_dir=<path to folder for output data>
output.prefix<text to prepend to output files>
```

# Merging without a reference

In the case of not having a scaling/merging reference dataset, a 'bootstrap' proceedure can be used, similar to the one used by PRIME (Uervirojnangkoorn 2015) or used in Brewster 2019b.  First, use the 'mark1' scaling algorithm which performs a simple average of the data, using these parameters:

```
dispatch.step_list=input balance model_scaling modify filter scale statistics_unitcell statistics_beam model_statistics statistics_resolution group errors_merge statistics_intensity merge statistics_intensity_cxi
scaling.algorithm=mark1
scaling.unit_cell=<unit cell>
scaling.space_group=<space group>
merging.error.model=errors_from_sample_residuals
<additional parameters as needed>
```

Note that this list of workers does not include post-refinement, scaling.model is not specified, and that an alternate error model is used which does not refine error terms. This will create a new mtz file, which can then be used as a reference in subsequent merging. We recommend using several iterations. For example, use a mark1 merge followed by several regular merges, each using the output mtz from the previous merge as a reference. CC1/2 should improve over the first few cycles and converge, but please be in contact with the authors with questions.

# References

- Sauter NK (2015). J Synchrotron Radiat 22, 239-48.
- Uervirojnangkoorn M, et. al. (2015). eLife 10.7554/eLife.05421.
- Brewster AS, et. al. (2018). Acta Crystallogr D Struct Biol 74, 877-894.
- Brewster AS, et. al. (2019a). Computational Crystallography Newsletter 10, 22-39.
- Brewster AS, et. al. (2019b). Acta Crystallogr D Struct Biol 75, 959-968.
