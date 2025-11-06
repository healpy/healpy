"""
Example demonstrating the legend issue and workaround for geographic projections.

This file shows the correct way to add legends when using newprojplot with
geographic projections like mollweide, hammer, or aitoff.

ISSUE: Calling plt.legend() without arguments can hang or be very slow
SOLUTION: Always specify the 'loc' parameter explicitly
"""

import numpy as np
import healpy as hp
from healpy.newvisufunc import projview, newprojplot
import matplotlib.pyplot as plt

# Create a sample map
nside = 16
npix = hp.nside2npix(nside)
m = np.random.random(npix)

# ============================================================================
# CORRECT USAGE: Specify legend location explicitly
# ============================================================================

# Create a plot with mollweide projection
projview(m, projection_type='mollweide', title='Correct: Legend with explicit loc')

# Add a circle (from the original issue)
circle = np.zeros(10000) + np.radians(45), np.arange(-180, 180, 0.036)
rotMat = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]
rotCircle = hp.rotator.rotateDirection(rotMat, circle)
newprojplot(*rotCircle, linewidth=1, linestyle="", marker=".", label="a circle", c="r")

# ✓ CORRECT: Specify loc explicitly to avoid hang
plt.legend(loc='upper right')
plt.savefig('example_legend_correct.png', dpi=100, bbox_inches='tight')
plt.close()

# ============================================================================
# Alternative: Use bbox_to_anchor for more control
# ============================================================================

projview(m, projection_type='hammer', title='Alternative: bbox_to_anchor')
theta = np.linspace(0, np.pi, 100)
phi = np.sin(theta) * np.pi
newprojplot(theta, phi, 'g-', linewidth=2, label='Sine curve')

# ✓ CORRECT: Use bbox_to_anchor to position legend outside axes
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig('example_legend_bbox.png', dpi=100, bbox_inches='tight')
plt.close()

# ============================================================================
# Multiple elements example
# ============================================================================

projview(m, projection_type='aitoff', title='Multiple Plot Elements')

# Add meridian
theta1 = np.linspace(0, np.pi, 50)
phi1 = np.zeros_like(theta1) + np.pi/2
newprojplot(theta1, phi1, 'r-', linewidth=2, label='Meridian')

# Add equator
theta2 = np.zeros(50) + np.pi/2
phi2 = np.linspace(-np.pi, np.pi, 50)
newprojplot(theta2, phi2, 'b--', linewidth=2, label='Equator')

# ✓ CORRECT: Specify loc explicitly
plt.legend(loc='lower left')
plt.savefig('example_legend_multiple.png', dpi=100, bbox_inches='tight')
plt.close()

print("Examples generated successfully!")
print("\nAvailable legend locations:")
print("  'upper right', 'upper left', 'lower left', 'lower right',")
print("  'right', 'center left', 'center right', 'lower center',")
print("  'upper center', 'center'")
print("\nOr use numeric codes: 0=best, 1=upper right, 2=upper left, etc.")
print("\n⚠️  IMPORTANT: Always specify 'loc' to avoid hangs with projections!")

# ============================================================================
# INCORRECT USAGE (commented out to prevent hanging):
# ============================================================================
# 
# projview(m, projection_type='mollweide')
# newprojplot(theta, phi, label='test')
# # ✗ INCORRECT: This will hang or be extremely slow!
# # plt.legend()  # DON'T DO THIS!
# 
# The above will hang because matplotlib tries to automatically find the
# best legend position, which requires expensive coordinate transformations
# with geographic projection axes (GeoAxes).
