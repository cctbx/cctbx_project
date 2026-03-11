# CLAUDE.md — xfel/ui/

The cctbx.xfel GUI: a wxPython + MySQL-backed experiment management system for serial
crystallography. Launched via `cctbx.xfel` (dispatcher in `command_line/xfel_gui_launch.py`).
Settings are persisted to `~/.cctbx.xfel/settings.phil` via `load_cached_settings()` /
`save_cached_settings()` in `ui/__init__.py`.

Being integrated with `dials_streaming` at:
`/global/cfs/cdirs/m4734/users/dwmoreau/dials_streaming/`
Full integration design: `skills/cctbx_xfel_gui.md` in that repo.

**Backend vs Facility**: These are orthogonal concepts. "Backend" (dispatcher) is
`cctbx.xfel.process` or `cctbx.xfel.small_cell_process` — set via `params.dispatcher`.
"Facility" is `params.facility.name` ∈ {`lcls`, `standalone`, `streaming`} — controls
how runs are discovered and rungroup bounds are validated. The `streaming` facility uses
`dials_streaming` for image delivery but discovers runs from the filesystem (same
`standalone_run_finder`) and allows creating rungroups before any runs exist in the DB.
Streaming rungroups store their intended run range in `streaming_first_run` /
`streaming_last_run` columns (added in schema v5.6) so the sentinel can sync them
correctly as runs arrive.

---

## Directory Layout

```
ui/
├── __init__.py          # master_phil_scope, load/save_cached_settings, settings_file path
├── phil_patch.py        # sync_phil() — handles unused PHIL definitions on settings load
├── db/
│   ├── __init__.py      # get_run_path(), get_image_mode(), get_db_connection()
│   ├── xfel_db.py       # xfel_db_application: core DB application class
│   ├── dxtbx_db.py      # dxtbx_xfel_db_application + log_frame() — primary logging API
│   ├── frame_logging.py # DialsProcessorWithLogging (MPI batch logging; model for streaming)
│   ├── experiment.py    # ORM: Event, Imageset, Beam, Detector, Crystal, Cell, Bin, ...
│   ├── run.py           # Run ORM
│   ├── trial.py         # Trial ORM
│   ├── rungroup.py      # Rungroup ORM
│   ├── job.py           # Job, JobFactory ORM
│   ├── schema.sql       # Full MySQL schema
│   ├── stats.py         # Stats queries for the GUI plots
│   └── cfgs/            # Shifter/cluster config templates (imported by xfel.util.mp)
├── components/
│   ├── xfel_gui_init.py    # MainWindow, tab management, sentinels
│   ├── xfel_gui_dialogs.py # Trial/rungroup/job creation dialogs; SettingsDialog (login)
│   ├── xfel_gui_controls.py
│   ├── xfel_gui_plotter.py
│   ├── spotfinder_scraper.py  # imports xfel.util.reflection_length
│   └── ...
└── command_line/
    └── xfel_gui_launch.py   # Entry point: `cctbx.xfel` dispatcher
```

## External dependencies within xfel/

- `xfel.util.mp` — job submission command builders (SLURM, LSF, SGE, PBS, etc.);
  called from `xfel_gui_dialogs.py` and `xfel_gui_init.py`. `util/mp.py` in turn
  imports `ui/db/cfgs` for Shifter templates.
- `xfel.util.reflection_length` — `ReflectionsRadialLengthsFromFiles`; used only by
  `spotfinder_scraper.py` for plot data (display only, not relevant to streaming).

---

## Key APIs for dials_streaming Integration

### PHIL scope and settings
```python
from xfel.ui import master_phil_scope, load_cached_settings, save_cached_settings
# master_phil_scope = master_phil_str + db_phil_str (experiment_tag, db.host/port/name/...)
params = load_cached_settings()   # reads ~/.cctbx.xfel/settings.phil
```

### DB application
```python
from xfel.ui.db.dxtbx_db import dxtbx_xfel_db_application, log_frame
app = dxtbx_xfel_db_application(params, cache_connection=True)
```
`params` requires `params.experiment_tag` and `params.db.*`.

### Per-image logging
```python
log_frame(experiments, reflections, params, run, n_strong,
          timestamp=None, two_theta_low=None, two_theta_high=None,
          db_event=None, app=app, trial=trial)
```
- `params.input.rungroup` must be set to populate `rungroup_id` in the event table.
- Returns SQL insert strings (cache_commits mode); flush with
  `app.execute_query(inserts, commit=True)`.

### Output path convention (`db/__init__.py:get_run_path`)
```
$output_dir/r{run:04d}/{trial.trial:03d}_rg{rungroup.id:03d}/
```
`rungroup.id` is the DB primary key, not a sequential user number.

---

## Important Behavioral Notes

- All table names prefixed with `{experiment_tag}_`.
- `rungroup.active` and `trial.active` must both be True for config selection.
- Selection priority: highest `trial.trial` number, then highest `rungroup.id`.
- `rungroup.open = True` means the rungroup still accepts new runs.
- `run.run` is VARCHAR — store run identifiers as strings.
- `DialsProcessorWithLogging` uses MPI and batched commits — do NOT use directly in
  dials_streaming. Use `log_frame()` directly instead.
