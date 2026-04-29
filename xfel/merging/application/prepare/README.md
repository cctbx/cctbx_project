# Prepare for SpReAD pipeline

In recent experiments we have attempted SpReAD refinements with a "moving window"
of energy-sliced sub-datasets. In conjunction with an energy-scanned X-ray beam,
this energy slicing technique can provide quasi-monochromatic datasets across the
full range of an absorption edge, thus allowing the refinement of anomalous
scattering curves.

This worker provides a convenient interface to the full quasi-monochromatic SpReAD
pipeline:

- Stage 1: The full energy-scanning dataset is sorted by energy and sliced into
a large number (say 100) of sub-datasets.
- Stage 2: According to the "moving window" technique, the fine-sliced datasets
are reassembled into larger slices for merging. The first merge might accumulate
data from percentiles 0 through 20. Sub-datasets are allowed to overlap; therefore
the second merge includes percentiles 1 through 21, then percentiles 2 through 22,
etc. Thus, for a window width of 20% and a spacing of 1%, we generate 81
sub-datasets.
- For each merged sub-dataset, the anomalous scattering factors are refined in
Phenix.
- The Phenix logs are analyzed to give a results file with the sub-dataset average
wavelengths and the refined anomalous scattering factors.

### Job coordination

Stage 1 is a single task performed on a large dataset; we have implemented
it as a standard merging worker `prepare_spread`. Stage 2 is a series of many
smaller merging jobs, followed by the same number of Phenix refinement jobs. It
is not necessarily practical to run the full Stage 1/Stage 2 pipeline at every
chance. Therefore we write a Stage 2 batch script that the user may submit to
the queuing system at their discretion. The Stage 2 batch script is implemented
as a Slurm "array job" in which N jobs are scheduled sequentially from a common
batch script. The first N-1 jobs are cctbx.xfel.merge runs to generate .mtz files
for the windowed datasets; the final job is a parallel run of all N-1 Phenix
refinements. 

### Usage

Include this worker as the merging step of a standard dataset. A possible example
follows here; this text would be suitable for pasting into an XFEL GUI dataset
merging task.

```
dispatch.step_list=input model_scaling statistics_unitcell statistics_beam model_statistics statistics_resolution prepare_spread group errors_merge statistics_intensity merge statistics_intensity_cxi publish
input.parallel_file_load.method=uniform
scaling.model=/path/to/model.pdb
scaling.resolution_scalar=0.96
statistics.n_bins=20
merging.d_min=1.9
merging.merge_anomalous=False
prepare.spread {
  stage2_phil=/path/to/merging/params.phil
  phenix_phil=/path/to/refinement/params.phil
  phenix_pdb=/path/to/model.pdb
  slurm_qos=realtime
  slurm_account=lcls
  slurm_constraint=cpu
  slurm_time_limit=120
  stage2_nnodes=2
  n_anomalous_scatterers=8
  binning=width
  bin_start_eV=6515
  bin_end_eV=6585
  bin_width_eV=8
  statistics_bin_i=16
  cctbx_activate=/path/to/cctbx/environment/setup.sh
  phenix_activate=/path/to/phenix/environment/setup.sh
}
publish.drive.credential_file=<redacted>
publish.drive.shared_folder_id=<redacted>
```

#### prepare.spread phil parameters

- `stage2_phil`: The path to a separate phil file for the Stage 2 merging jobs.
Paths will be included automatically. This is a suitable example:

```
dispatch.step_list=input model_scaling filter statistics_unitcell statistics_beam model_statistics statistics_resolution group errors_merge statistics_intensity merge statistics_intensity_cxi
input.parallel_file_load.method=uniform
scaling.model=/path/to/model.pdb
scaling.resolution_scalar=0.96
statistics.n_bins=20
merging.d_min=1.9
merging.merge_anomalous=False
merging.error.model=mm24
```
- `phenix_phil`: A phenix phil file with refinement of anomalous scattering
factors activated. This excerpt is the important part:
```
  refine {
    strategy = group_anomalous

    anomalous_scatterers {
      group {
        selection = "element Ca"
        f_prime = 0.2893659
        f_double_prime = 1.637287
      }
      group {
        selection = "name Mn1 and chain A"
        f_prime = -1.837302
        f_double_prime = 3.480123
      }
      group { [...]
      }
    }
  }
```
- `phenix_pdb`: The model for subsequent phenix refinement. This should be a fully
converged refinement from a "remote" dataset with no significant uncertainty in the
anomalous scattering factors.
- `slurm_qos`, `slurm_account`, `slurm_constraint`, `slurm_time_limit`: Configuration
items for your local queuing system. Only Slurm is currently supported. Depending
on your local environment, you may replace `slurm_qos` with `slurm_partition`.
- `stage2_nnodes`: Number of nodes to request for the Stage 2 merging tasks.
- `n_anomalous_scatterers`: The number of entries to scrape from the phenix logs.
- `binning`: Choose `width` or `count`. See below for discussion.
- `bin_start_eV`: Left edge of the first energy bin.
- `bin_end_eV`: Right edge of the last energy bin.
- `bin_width_eV`: In this example, the first bin contains energies 6515 to 6523,
the second bin contains 6516 to 6524, third bin 6517 to 6525, etc.
- `statistics_bin_i`: The multiplicity of the sliced sub-datasets is reported for
the `i`-th bin in the merging logs.
- `cctbx_activate`: The path to a cctbx/dials installation (the same one used for
regular processing).
`phenix_activate`: The path to a phenix installation.

#### Binning modes

In `count` mode, energy slicing is by percentile, so that each sub-dataset will
have the same multiplicity. Where energy coverage is weaker, the sub-dataset ranges
will be wider.

In `width` mode, energy slicing is by absolute energy, so that each sub-dataset
will cover the same range in eV. Where energy coverage is weaker, the sub-dataset
multiplicity will be lower.
