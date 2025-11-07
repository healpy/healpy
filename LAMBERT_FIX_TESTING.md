# Lambert Projection Fix - Testing Guide

## Problem Description

The Lambert projection in `projview` was displaying the lower half of the sky incorrectly, while Mollweide and Hammer projections worked fine. This suggested a coordinate orientation issue specific to Lambert.

## Root Cause

Matplotlib's Lambert azimuthal equal-area projection interprets the Y-axis (latitude) in reverse order compared to other geographic projections like Mollweide and Hammer.

## Solution Implemented

The fix adds a vertical flip of the `grid_map` array specifically for Lambert projection:

```python
if projection_type == "lambert":
    grid_map = np.flip(grid_map, axis=0)
```

This is done in `lib/healpy/newvisufunc.py` at line 642, after the `grid_map` is computed from the HEALPix map but before it's passed to `pcolormesh`.

## Testing the Fix

### Automated Tests

Run the new test suite:
```bash
# All Lambert-specific tests
pytest test/test_lambert_projection_fix.py -v

# Specific tests
pytest test/test_lambert_projection_fix.py::test_lambert_projection_basic -v
pytest test/test_lambert_projection_fix.py::test_lambert_projection_half_sky_north -v
pytest test/test_lambert_projection_fix.py::test_lambert_vs_mollweide_consistency -v
```

### Visual Tests

Create visual comparison plots:
```bash
python /tmp/visual_lambert_demo.py
```

This creates two PNG files:
- `/tmp/lambert_fix_demonstration.png` - Comparison of Mollweide, Hammer, and Lambert projections
- `/tmp/lambert_half_sky_demonstration.png` - Northern and southern hemisphere views

### Manual Testing

Test with real WMAP data (run from repository root):
```python
import healpy as hp
import matplotlib.pyplot as plt

# Load test data - adjust path if needed
map_data = hp.read_map('test/data/wmap_band_iqumap_r9_7yr_W_v4_udgraded32_masked_smoothed10deg_fortran.fits')

# Test full-sky Lambert
hp.newvisufunc.projview(
    map_data,
    coord='G',
    projection_type='lambert',
    title='Lambert projection',
)
plt.show()

# Test half-sky Lambert (northern hemisphere)
hp.newvisufunc.projview(
    map_data,
    coord='G',
    projection_type='lambert',
    latra=[0, 90],
    title='Lambert northern hemisphere',
)
plt.show()

# Compare with Mollweide
hp.newvisufunc.projview(
    map_data,
    coord='G',
    projection_type='mollweide',
    title='Mollweide projection',
)
plt.show()
```

## Expected Results

After the fix:
1. Lambert projection should display maps correctly with no corruption in the lower half
2. The pattern should match what Mollweide and Hammer show
3. Northern hemisphere (positive latitudes) should appear in the top half of the plot
4. Southern hemisphere (negative latitudes) should appear in the bottom half of the plot
5. Half-sky projections should work correctly for both hemispheres

## Verification Checklist

- [ ] Lambert full-sky projection displays correctly
- [ ] Lambert northern hemisphere (latra=[0, 90]) displays correctly
- [ ] Lambert southern hemisphere (latra=[-90, 0]) displays correctly
- [ ] Lambert matches Mollweide/Hammer for the same data
- [ ] Existing tests still pass (test_newvisufunc_example.py)
- [ ] No regression in other projection types

## Files Modified

1. `lib/healpy/newvisufunc.py` - Added grid_map flip for Lambert projection
2. `test/test_lambert_projection_fix.py` - New comprehensive test suite

## Additional Notes

- The fix is minimal and only affects Lambert projection
- Other projection types (Mollweide, Hammer, Aitoff, cart, polar, 3d) are unchanged
- The `return_only_data` parameter returns the flipped data for Lambert, which is consistent with what gets plotted
