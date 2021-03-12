# Resolving an Indexing Ambiguity and Laue Group Assignment

## Problem statement

Indexing ambiguity is well known in crystallography. For polar space groups, where the symmetry of the cell contents is lower than the symmetry of the lattice,
the Bragg spots can be indexed in two or more ways.  This ambiguity can only be addressed by analyzing the reflection intensities.  Two algorithms are available
in `cctbx.xfel.merge`:

1. `modify_reindex_to_reference` steps through each diffraction pattern individually, correlating the Bragg intensities with those from model Fcalcs, from
a structural model.  Clearly there are disadvantages:  the data may not be exactly isomorphous to the reference, so there may be false
assignments of the indexing sense.  In fact, the data analyzed for the use case here produce correct alignment only 85% of the time.

2. `modify_cosym` performs a mutual alignment of all the XFEL still shots, using an NxN matrix of correlation coefficients (N is the number of diffraction patterns).
The approach was first proposed by
Brehm and Diederichs to address the indexing ambiguity problem.  Later, Gildea and Winter expanded the concept to allow both the alignment of polar cells,
as well as the de novo determination of Laue symmetry.  The present implementation allows for the full Gildea and Winter analysis, but is only tested on one
use case:  the alignment of crystals with P63 symmetry.

## Installation notes

User, please treat this as an alpha testing situation with bugs and additional use cases reported to @nksauter.  The
following are key points for getting things to run:

1.  Please use a conda build using the xfel code target.
2.  The conda environment must provide `scikit-learn` for unit cell covariance analysis, and `pandas` for analysis of data tables.

## Use case
Here we describe the mutual alignment of an 8000-image dataset with space group P63, merged in about 30 minutes on a 60-core AMD server. Data files
were previously indexed and integrated with `dials.stills_process` in space group P63.  The merging command line
was `mpirun -n 60 cctbx.xfel.merge test.phil`
with the phil file:
```
input.path=<path to the dials-integrated data with linux wildcards permitted>
input.experiments_suffix = .expt
input.reflections_suffix = .refl
dispatch.step_list=input balance model_scaling modify filter modify_cosym errors_premerge scale postrefine statistics_unitcell statistics_beam model_statistics statistics_resolution group errors_merge statistics_intensity merge statistics_intensity_cxi
input.parallel_file_load.method=uniform
filter.outlier.min_corr=-1
filter.algorithm=unit_cell
filter.unit_cell.algorithm=cluster
filter.unit_cell.cluster.covariance.file=../tdata/covariance_trial_cells.pickle
filter.unit_cell.cluster.covariance.component=0
filter.unit_cell.cluster.covariance.mahalanobis=4.0
select.algorithm=significance_filter
select.significance_filter.sigma=0.1
select.significance_filter.max_ct=300
select.significance_filter.min_ct=200
modify.reindex_to_reference.dataframe=test_reindex_dataframe.pickle
modify.cosym.space_group=P6
modify.cosym.dataframe=test_cosym_dataframe.pickle
modify.cosym.anchor=True
modify.cosym.min_reflections=15
modify.cosym.normalisation=None
modify.cosym.d_min=2.5
modify.cosym.dimensions=2
modify.cosym.cluster.n_clusters=2
modify.cosym.min_pairs=3
modify.cosym.nproc=1
modify.cosym.weights=count
modify.cosym.plot.interactive=True
scaling.model=<path to pdb file containing the reference structure.pdb>
scaling.resolution_scalar=0.95
scaling.mtz.mtz_column_F=I-obs
merging.d_min=3.1
statistics.n_bins=15
merging.merge_anomalous=False
postrefinement.enable=True
merging.error.model=ev11
output.do_timing=True
output.prefix=test_trial
output.output_dir=./test_reference
output.log_level=0
```
Line-by-line discussion of these options:
```
mpirun -n 60 cctbx.xfel.merge
```
There is a critical tradeoff involving the number of MPI ranks (-n).  Analysis of a large dataset (N=8000) would be prohibitive if the full NxN matrix
were to be analyzed.  Instead, we break the data into tranches of T = N//n shots, so in this case we analyze matrices of approximately 133 x 133
experiments within each MPI rank.  Useful tranch sizes T range from about 100 to 200.  Smaller T produces drastically shorter wall clock time, while larger T
produces dramatically superior Brehm-Diederichs embedding plots, as in Figure 4 of their paper.  The embedding plot (see `modify.cosym.plot.interactive`) must
be checked to ensure blue/red clusters are well separated from the 45-degree diagonal, and that the cluster centers are at a good distance from the origin
(0.4-0.8 is good).  If the clusters look bad, the tranch size should be increased.  Note, the algorithm will not work unless n>=5.
```
dispatch.step_list=input balance model_scaling modify filter modify_cosym errors_premerge scale postrefine statistics_unitcell statistics_beam model_statistics statistics_resolution group errors_merge statistics_intensity merge statistics_intensity_cxi
```
The dispatch step list must be explicitly given to resolve the lattice ambiguity.  The algorithm option, either `modify_reindex_to_reference` or
`modify_cosym` must be chosen exactly as given here, in slot number six of the list.
```
filter.algorithm=unit_cell
filter.unit_cell.algorithm=cluster
filter.unit_cell.cluster.covariance.file=../tdata/covariance_trial_cells.pickle
filter.unit_cell.cluster.covariance.component=0
filter.unit_cell.cluster.covariance.mahalanobis=4.0
```
Apply the covariance method to select unit cells that are closely isomorphous.  Refer to the [documentation](https://github.com/cctbx/cctbx_project/tree/merge_polar/xfel/merging/application/tdata) for details.

```
modify.reindex_to_reference.dataframe=test_reindex_dataframe.pickle
modify.cosym.dataframe=test_cosym_dataframe.pickle
output.output_dir=./test_reference
```
These pickle outputs allow us to directly compare the indexing assignments arising from the two algorithms.  In
detail: run the merging script twice, changing the `dispatch.step_list` to each algorithm in turn.  It is also
advisable to change the `output.prefix` or the `output.output_dir` so the different results do not overstep each
other.  Then, compare the index assignments with the program:
```
libtbx.python \
  ../modules/cctbx_project/xfel/merging/application/modify/compare_results.py \
  test_reference/test_reindex_dataframe.pickle \
  test_reference/test_cosym_dataframe.pickle
```
Other parameters of critical interest are:
```
modify.cosym.space_group=P6
```
If the correct space group of the crystal is known, we supply a symmorphic space group
in the same point group, e.g. P6(3) -> P6 for hexagonal photosystem I. This applies
algorithm 2 in Ref. 1. If no space group is supplied, then the clustering analysis is
as described in Ref. 2 with the goal of assigning the Laue group.

```
modify.cosym.anchor=True
scaling.model=<path to pdb file containing the reference structure.pdb>
```
As explained above, the reference structure is used as the fundamental reference in the `reindex_to_reference` option.  The reference is NOT used
for the mutual alignment of lattices during the `cosym` option.  However, for the output to be useful, the final output should be aligned
with the reference structure treated as an anchor. This is done after the Brehm-Diederichs analysis is complete (to avoid biasing anything), and
just before the cosym merging worker returns.  Specific purposes are: 1) with `mark0` merging, especially
if `postrefine` is used, the `anchor` must be set to `True` and the `scaling.model` given. Warning, this is mandatory but
there is no parameter validation in place to enforce this.  2) Anchoring the output against the reference
allows the merged data to be dropped into the refinement program directly for isomorphous refinement, without a molecular replacement step.
```
modify.cosym.min_reflections=15
modify.cosym.normalisation=None
```
These options from the Gildea program should probably never be changed.
```
modify.cosym.d_min=2.5
```
The user can specify the resolution range specifically for the `cosym` process, independent of which data to eventually merge.
```
modify.cosym.dimensions=2
```
This is absolutely critical.  We set this to 2 dimensions for space group P63, as there are exactly two groups (cosets) expected
in the final sort.  Setting this to the default (auto determine) has the unfortunate consequence of performing the embedding
analysis in a higher dimensional space (6) where the clusters cannot be found!  So for the present use case, keep this at 2.
```
modify.cosym.cluster.n_clusters=2
```
The number of clusters should match the number of dimensions.
```
modify.cosym.min_pairs=3
```
Per-experiment, per-symmetry-operator minimum number of mutual Miller indices to form a correlation-coefficient cross-term with another
expt/symop combination.  Takes the default from Gildea; increasing this will drastically reduce the data used for alignment, probably
to our detriment.
```
modify.cosym.nproc=1
```
`nproc > 1` permits Python multiprocessing for the calculation of the rij matrix.  Note, rij is a time-consuming process, but the next step, LBFGS, is
probably rate limiting, so it is not clear if `nproc` will really help.  Details:  `nproc` enables a hybrid parallelization model, where MPI ranks handle
separate data tranches, while Python multiprocessing speeds up work within each MPI rank.  Therefore, `nproc` should not be set
to the number of cores on the machine, but rather to the number of hyperthreads available to each rank.  For example, on Cori, one would use 68 ranks per
node, and 4 nproc hyperthreads per core. Knowledge of the specific machine architecture is needed.
```
modify.cosym.weights=count
```
This is critical.  Setting `weights=None` makes the cosym algorithm fail, while the other options have not been sucessfully tested in xfel.
```
modify.cosym.plot.interactive=True
```
This displays an embedding plot as in Ref. 1 to assess whether the multiple indexing solutions are sufficiently resolved from each other. The embedding
plot is a useful tool for setting the number of MPI ranks (and thus the tranch size for analysis). For routine work, this may be omitted, and a plot
will be saved in the output directory.
```
merging.merge_anomalous=False
```
The user may decide to either merge the Friedel mates or not.

## References

1. W. Brehm and K. Diederichs, ["Breaking the indexing ambiguity in serial crystallography,"](https://doi.org/10.1107/s1399004713025431)
Acta Crystallogr. D Biol. Crystallogr. 70, 101-109 (2014).

2. R.J. Gildea and G. Winter, ["Determination of Patterson group symmetry from sparse multi-crystal data sets in the presence of an indexing ambiguity,"](https://doi.org/10.1107/S2059798318002978) Acta
Crystallogr. D Struct. Biol. 74, 405-410 (2018).
