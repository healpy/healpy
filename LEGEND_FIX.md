# Legend Hang Issue Fix for Geographic Projections

## Issue Description

When using `newprojplot` to add labeled plot elements to geographic projection maps (mollweide, hammer, aitoff, lambert) and then calling `plt.legend()` without parameters, the program can hang or become extremely slow.

## Root Cause

The issue occurs because matplotlib's `legend()` function, when called without a `loc` parameter, uses `loc='best'` by default. This triggers an automatic positioning algorithm that tries to find the optimal legend position by:
1. Computing the bounding boxes of all plot elements
2. Testing multiple potential legend positions
3. Evaluating overlaps with existing artists

With geographic projection axes (matplotlib's `GeoAxes`), these coordinate transformations are computationally expensive and can cause infinite loops or extreme slowness.

## Solution

Always specify an explicit `loc` parameter when calling `plt.legend()` on plots created with geographic projections.

### Correct Usage

```python
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

# Create a map
hp.newvisufunc.projview(np.random.random(12*16**2))

# Make a circle
circle = np.zeros(10000) + np.radians(45), np.arange(-180, 180, 0.036)
rotMat = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]
rotCircle = hp.rotator.rotateDirection(rotMat, circle)

# Add it to the plot with a label
hp.newprojplot(*rotCircle, linewidth=1, linestyle="", marker=".", label="a circle", c="r")

# ✓ CORRECT: Specify loc explicitly
plt.legend(loc='upper right')
plt.show()
```

### Available Legend Locations

You can use any of these string values for `loc`:
- `'upper right'` (default when not using 'best')
- `'upper left'`
- `'lower left'`
- `'lower right'`
- `'right'`
- `'center left'`
- `'center right'`
- `'lower center'`
- `'upper center'`
- `'center'`

Or numeric codes:
- `0` = 'best' (⚠️ avoid with projections!)
- `1` = 'upper right'
- `2` = 'upper left'
- `3` = 'lower left'
- `4` = 'lower right'
- etc.

### Alternative: bbox_to_anchor

For more control over legend positioning, especially to place it outside the axes:

```python
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
```

## Implementation Changes

The fix includes:

1. **Documentation Update**: Added warnings in the docstrings of both `projview` and `newprojplot` functions explaining the issue and showing correct usage.

2. **Test Suite**: Added comprehensive tests in `test/test_newvisufunc_legend.py` that verify legend functionality works correctly with explicit `loc` parameter.

3. **Examples**: Added `doc/newvisufunc_legend_example.py` demonstrating correct usage patterns.

## Testing

Run the test suite to verify the fix:

```bash
pytest test/test_newvisufunc_legend.py -v
```

The tests verify that:
- Legends work correctly when `loc` is explicitly specified
- Legends work with `bbox_to_anchor` positioning
- The return value from `newprojplot` contains proper Line2D objects with labels
- The example from the original issue works with the fix applied

## Backward Compatibility

This fix is backward compatible:
- No changes to function signatures
- No changes to function behavior
- Only documentation additions
- Users who already specify `loc` are unaffected
- Users who don't use legends are unaffected

The only "breaking" change is that users who called `plt.legend()` without parameters will now need to add `loc` parameter, but since that was hanging anyway, this is a necessary fix.
