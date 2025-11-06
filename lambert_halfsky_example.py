#!/usr/bin/env python3
"""
Example: Half-Sky Lambert Projection with healpy.newvisufunc.projview

This example demonstrates how to use the new latra parameter with Lambert projection
to create half-sky plots, similar to the old orthview(half_sky=True) functionality.

The implementation allows restricting the latitude range while maintaining full
longitude coverage, which is ideal for hemispheric views.
"""

import healpy as hp
from healpy.newvisufunc import projview
import numpy as np
import matplotlib.pyplot as plt

# Create a sample HEALPix map
nside = 32
npix = hp.nside2npix(nside)

# Generate sample data - a simple pattern
theta, phi = hp.pix2ang(nside, np.arange(npix))
m = np.sin(3 * phi) * np.cos(2 * theta)

# Example 1: Full-sky Lambert projection (default behavior)
print("Example 1: Full-sky Lambert projection")
projview(
    m,
    projection_type="lambert",
    title="Full Sky Lambert\nprojview(m, projection_type='lambert')",
    graticule=True,
    graticule_labels=True,
)
plt.savefig("lambert_fullsky.png", dpi=150, bbox_inches="tight")
plt.close()

# Example 2: Northern hemisphere half-sky
print("Example 2: Northern hemisphere half-sky")
projview(
    m,
    projection_type="lambert",
    latra=[0, 90],
    title="Northern Hemisphere Half-Sky\nprojview(m, projection_type='lambert', latra=[0, 90])",
    graticule=True,
    graticule_labels=True,
)
plt.savefig("lambert_northern_hemisphere.png", dpi=150, bbox_inches="tight")
plt.close()

# Example 3: Southern hemisphere half-sky
print("Example 3: Southern hemisphere half-sky")
projview(
    m,
    projection_type="lambert",
    latra=[-90, 0],
    title="Southern Hemisphere Half-Sky\nprojview(m, projection_type='lambert', latra=[-90, 0])",
    graticule=True,
    graticule_labels=True,
)
plt.savefig("lambert_southern_hemisphere.png", dpi=150, bbox_inches="tight")
plt.close()

# Example 4: Custom latitude range (tropical region)
print("Example 4: Custom latitude range (tropical region)")
projview(
    m,
    projection_type="lambert",
    latra=[-30, 30],
    title="Tropical Region\nprojview(m, projection_type='lambert', latra=[-30, 30])",
    graticule=True,
    graticule_labels=True,
)
plt.savefig("lambert_tropical.png", dpi=150, bbox_inches="tight")
plt.close()

# Example 5: Three views comparison in a single figure
print("Example 5: Three views comparison")
fig = plt.figure(figsize=(18, 6))

# Full sky
plt.subplot(131, projection="lambert")
projview(
    m,
    projection_type="lambert",
    title="Full Sky",
    graticule=True,
    hold=True,
)

# Northern hemisphere
plt.subplot(132, projection="lambert")
projview(
    m,
    projection_type="lambert",
    latra=[0, 90],
    title="Northern Hemisphere\nlatra=[0, 90]",
    graticule=True,
    hold=True,
)

# Southern hemisphere
plt.subplot(133, projection="lambert")
projview(
    m,
    projection_type="lambert",
    latra=[-90, 0],
    title="Southern Hemisphere\nlatra=[-90, 0]",
    graticule=True,
    hold=True,
)

plt.tight_layout()
plt.savefig("lambert_comparison.png", dpi=150, bbox_inches="tight")
plt.close()

print("\nAll example plots saved successfully!")
print("\nKey points:")
print("- Use latra=[0, 90] for northern hemisphere")
print("- Use latra=[-90, 0] for southern hemisphere")
print("- Use latra=[min, max] for any custom latitude range")
print("- lonra is NOT supported for Lambert (full longitude maintained)")
print("- This provides similar functionality to orthview(half_sky=True)")
