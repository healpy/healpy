# Generating Lambert Projection Fix Demonstration Images

This document explains how to generate visual demonstration images for the Lambert projection fix.

## Quick Start

From the repository root, after installing healpy:

```bash
python generate_lambert_demo_images.py
```

This will create two PNG files in the current directory:
1. `lambert_fix_demonstration.png` - Side-by-side comparison of all three projections
2. `lambert_half_sky_demonstration.png` - Northern and southern hemisphere views

## Prerequisites

The script requires:
- Healpy installed with the Lambert fix applied
- Matplotlib
- NumPy

Install with:
```bash
# Install healpy in development mode
python -m pip install -e .

# Or if already installed, just ensure matplotlib is available
python -m pip install matplotlib numpy
```

## Generated Images

### Image 1: Projection Comparison (lambert_fix_demonstration.png)

A 3×3 grid showing three different test patterns across three projections:

**Columns:**
- Mollweide projection (reference, working correctly)
- Hammer projection (reference, working correctly)
- Lambert projection (should now work correctly with the fix)

**Rows:**
1. **Latitude Gradient**: Smooth color transition from South (blue) to North (red)
2. **North/South Step**: Clear boundary at equator, positive values in north, negative in south
3. **Sinusoidal Pattern**: sin(2×latitude) showing symmetric waves

**Expected Result:** All three columns should show identical patterns. No corruption or distortion should be visible in the Lambert projection.

### Image 2: Half-Sky Projections (lambert_half_sky_demonstration.png)

A 2×2 grid showing Lambert projection with restricted latitude ranges:

**Panels:**
1. **Top-left**: Northern hemisphere (latra=[0, 90]) with latitude gradient
2. **Top-right**: Southern hemisphere (latra=[-90, 0]) with latitude gradient
3. **Bottom-left**: Northern hemisphere with step function
4. **Bottom-right**: Southern hemisphere with step function

**Expected Result:** Each hemisphere should display correctly with appropriate data. Northern hemisphere should show only positive latitudes, southern hemisphere only negative latitudes.

## Verification Checklist

When reviewing the generated images, verify:

- [ ] **No corruption**: Lambert projection shows smooth, continuous data (no garbled lower half)
- [ ] **Consistency**: Lambert matches Mollweide and Hammer patterns exactly
- [ ] **Gradient**: Latitude gradient smoothly transitions from blue (south) to red (north)
- [ ] **Step function**: Clear, sharp boundary at the equator
- [ ] **Sinusoidal**: Symmetric wave pattern with peaks at expected latitudes
- [ ] **Half-sky north**: Only northern hemisphere data displayed
- [ ] **Half-sky south**: Only southern hemisphere data displayed

## Troubleshooting

### Import Error

If you get `ImportError: No module named 'healpy'`:
```bash
# Ensure healpy is installed
python -m pip install -e .
```

### Display Issues

The script uses the `Agg` backend (non-interactive) to generate images without requiring a display. Images are saved directly to disk.

### Low Resolution

If images appear pixelated, increase the DPI in the script:
```python
plt.savefig(output_file, dpi=300, bbox_inches='tight')  # Change from 150 to 300
```

## Using in CI/CD

To generate images in a CI pipeline:

```yaml
- name: Generate Lambert demo images
  run: |
    python generate_lambert_demo_images.py
    
- name: Upload demonstration images
  uses: actions/upload-artifact@v3
  with:
    name: lambert-demo-images
    path: |
      lambert_fix_demonstration.png
      lambert_half_sky_demonstration.png
```

## Manual Testing

For interactive testing without generating images:

```python
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

# Create test map
nside = 32
npix = hp.nside2npix(nside)
theta, phi = hp.pix2ang(nside, np.arange(npix))
lat = 90 - np.degrees(theta)
test_map = lat

# Test Lambert projection
hp.newvisufunc.projview(
    test_map,
    projection_type="lambert",
    title="Lambert Test",
    cbar=True,
)
plt.show()

# Compare with Mollweide
hp.newvisufunc.projview(
    test_map,
    projection_type="mollweide",
    title="Mollweide Reference",
    cbar=True,
)
plt.show()
```

## Questions?

If the generated images don't match expectations:
1. Check that you're using the branch with the Lambert fix
2. Verify healpy is installed from the local source (not a cached version)
3. Try `pip install -e . --force-reinstall --no-deps` to ensure the fix is applied
4. Check `git log` to confirm the Lambert fix commits are present
