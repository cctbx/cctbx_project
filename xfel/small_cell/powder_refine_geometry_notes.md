# Powder Geometry Refiner - Implementation Notes

## Summary
Implemented `cctbx.xfel.powder_refine_geometry` - a tool to refine detector geometry using powder diffraction d-spacings.

## Files Created
1. **`xfel/small_cell/geometry_refiner.py`** - `PowderGeometryRefiner` class
2. **`xfel/small_cell/command_line/powder_refine_geometry.py`** - CLI tool

## Usage
```bash
cctbx.xfel.powder_refine_geometry spots.expt spots.refl reference_d_spacings=8.0

# Options:
#   refine.distance=True/False
#   refine.shift=True/False
#   refine.tilt=True/False
#   max_distance_inv_ang=0.002  (filter: only use spots within this distance of reference)
#   d_min=1.5, d_max=20
#   output.experiments=refined.expt
```

## 5-Parameter Model
- **dist**: distance shift along detector normal (mm)
- **shift1**: shift along fast axis (mm)
- **shift2**: shift along slow axis (mm)
- **tau2**: rotation around fast axis (mrad)
- **tau3**: rotation around slow axis (mrad)

## Key Implementation Details

### Rotation Convention
- Uses `axis_and_angle_as_r3_rotation_matrix()` for rotations around actual detector axes
- Handles left-handed detector frame (slow axis typically points -y)
- Rotates around panel CENTER (not origin) to avoid introducing spurious shifts

### Filtering
- Reflections filtered by d-range (d_min, d_max)
- Then filtered by proximity to reference d-spacings in inverse angstroms
- Default: 0.002 inv Å (~0.13 Å at 8 Å)

### Optimizer
- Uses scipy Powell method (derivative-free, robust)
- L-BFGS-B had issues with gradient estimation for tilts

## Key Findings During Development

### Tilt vs Shift Degeneracy
With a single d-spacing reference, tilts and shifts can partially compensate for each other:
- A tilt creates an elliptical ring (azimuthal d-spacing variation)
- A shift moves the center but keeps the ring circular
- The sum-of-squares objective doesn't fully distinguish them
- Multiple d-spacings at different radii would better constrain tilts

### Rotation Around Center
Critical bug fix: must rotate around panel CENTER, not origin. Otherwise, rotating the axes causes an unintended ~0.3 mm shift for every 0.1 degree tilt (arc = radius × angle, with radius ~170 mm).

### Sign Conventions
- Detector slow axis typically points -y (pixel coords increase downward)
- This creates a left-handed frame
- Using axis-angle rotation with actual detector axes handles this correctly

## Test Data
- Location: `~/daily/20260120/data/256/out/combined/`
- Files: `combined.expt`, `combined.refl`
- Created sabotaged versions with known tilt errors for testing

## Creating Sabotaged Test Files
```python
from dxtbx.model.experiment_list import ExperimentListFactory
from scitbx import matrix
import numpy as np

experiments = ExperimentListFactory.from_json_file('combined.expt', check_format=False)
panel = experiments[0].detector[0]
fast = matrix.col(panel.get_fast_axis())
slow = matrix.col(panel.get_slow_axis())
origin = matrix.col(panel.get_origin())

# Compute center
size = panel.get_image_size()
px = panel.get_pixel_size()
center = origin + (size[0]/2 * px[0]) * fast + (size[1]/2 * px[1]) * slow

# Apply tilts (e.g., 3 degrees)
tau = 3.0 * np.pi / 180.0
R2 = fast.axis_and_angle_as_r3_rotation_matrix(tau, deg=False)
R3 = slow.axis_and_angle_as_r3_rotation_matrix(tau, deg=False)
R = R3 * R2

fast_new = R * fast
slow_new = R * slow
origin_new = center + R * (origin - center)

panel.set_frame(fast_new.elems, slow_new.elems, origin_new.elems)
experiments.as_file('combined_sabotaged.expt')
```

## Future Improvements
- Multiple d-spacing references for better tilt constraint
- Analytical derivatives for faster convergence
- Multi-panel detector support
- Azimuthal-aware objective function to better distinguish tilts from shifts
- Outlier rejection

## Session Date
2026-01-20
