# Instructions for merging multiple specific isoforms within a serial crystallography dataset

The purpose of this workflow is to separately group crystal isoforms with very similar unit cells.
This is distinct from tools like dials.cosym, which deals with the assignment of Laue symmetry,
and with the resolution of indexing ambiguity for polar space groups.

## General steps in the pipeline are:

1) Process the serial crystallography diffraction patterns with dials.stills_process, 
producing files that contain refined experiment models and integrated Bragg spots.

2) Run a short version of the program cctbx.xfel.merge with two steps only: 
reading the data in, and producing a text list of the unit cells and space groups of
all integrated lattices.

```
cores=32

export effective_params="dispatch.step_list=\"input tdata\" \
input.path=${DIALS_OUTPUT}/r*/*/out \
input.experiments_suffix=_integrated.expt \
input.reflections_suffix=_integrated.refl \
input.parallel_file_load.method=uniform \
tdata.output_path=${TRIAL}_cells \
output.prefix=${TRIAL} \
output.output_dir=${OUT_DIR}/${TRIAL}/out 
output.tmp_dir=${OUT_DIR}/${TRIAL}/tmp \
output.do_timing=True \
output.log_level=0"
 
mpiexec -n $cores cctbx.xfel.merge  ${effective_params}
```

sample output:
```
78.643232 78.643232 38.831412 90.000000 90.000000 90.000000 P4/mmm
78.889991 78.889991 38.179862 90.000000 90.000000 90.000000 P4/mmm
```

3) Run a standalone program to produce clusters with the DBSCAN algorithm, which are then
fit with multivariate ellipsoids.  Importantly, the features used for classification
are currently the unit cell parameters a,b,c and the relevant distance metric
is the L2 norm.  Polytope boundaries in unit cell space are not treated, as
would be the case if we used one of the Andrews-Bernstein metrics.  The "pros"
to this are that we can use a standard clustering algorithm from scikit-learn. The
"cons" are that we cannot support triclinic, monoclinic, or pseudo-symmetry.

Specifically we have:
uc_metrics.dbscan for tetragonal and hexagonal cells, and
uc_metrics.dbscan3d for the orthorhombic case.

A use case (lysozyme) would run as follows:
```
uc_metrics.dbscan file_name=027_cells.tdata space_group=P4/mmm eps=0.05 feature_vector=c,a write_covariance=True
```

Here is *.tdata file is the unit cell list from step 2).  The space_group tag forces the program to only
consider cells claiming that value as the space group.  The feature_vector defines which unit cell
parameters to treat as relevant for clustering.

eps is the characteristic maximum distance in the DBSCAN algorithm above which items are no longer
considered part of the same cluster.  The idea here is that "eps" makes this an interactive program:
users are asked to run the program with several different values of "eps" and then choose the clustering
result that most closely matches with visual perception.  Users should especially notice the
cluster centers (marked in yellow), as well as the standard deviation contours (sigma levels 1,2,3,4)
that should nicely follow the colored clusters.  If unit cell instances are fairly sparse and/or no
clusters are labeled, the "eps" parameter should be increased.  If however, it appears that a labeled
cluster should actually be two separate components (perhaps with unequal populations) the "eps"
parameter should be decreased.  The user should note the component number and sigma level that
matches the cluster of interest.

sample output:
```
027_cells, eps=0.020
Component        Abundance           c            a
0: 134147 /137347 (97.67%)   38.12±0.06   78.96±0.11
1:   2452 /137347 ( 1.79%)   36.94±0.05   80.20±0.10
2:    222 /137347 ( 0.16%)   38.51±0.05   78.83±0.09
```

4) To list out the clusters previously determined in step 3, one can use this convenience command:
```
uc_metrics.list_covariance covariance_027_cells.pickle
```

5) Now the full merging script (cctbx.xfel.merge) is run again, with special phil parameters
to select those unit cells belonging to the cluster components of step 3. In this particular
example the merging program would be run three times, once for each component (0, 1, or 2):

```
filter.algorithm=unit_cell \
filter.unit_cell.algorithm=cluster \
filter.unit_cell.cluster.covariance.file=LD91_new/covariance_027_cells.pickle \
filter.unit_cell.cluster.covariance.component=0 \
filter.unit_cell.cluster.covariance.mahalanobis=4.0 \
```

One does not expect the merging program to pick the exact population for each component that
the clustering program did in step 3.  The reason is that different algorithms are used.  The
clustering program uses the DBSCAN algorithm, controlled by the "eps" parameter.  The merging
program, in contrast, uses a multivariate Gaussian fit, defined by the covariance matrix from
the pickle file, and controlled by the "mahalanobis" parameter.  This is essentially a cutoff
defining how many standard deviations away from the mean to include unit cells.  A value of
mahalanobis=1.0 will pick a very restricted population very close to the mean.  The default
value of mahalanobis=4.0 will in contrast pick 99.9% of the population defined by the ellipsoid.

