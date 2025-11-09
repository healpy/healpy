#!/usr/bin/env python
"""
Generate demonstration images for the Lambert projection fix.

This script creates visual comparisons showing that Lambert projection
now displays maps correctly after the fix.

Run this from the repository root after installing healpy:
    python generate_lambert_demo_images.py

Generated files will be saved to the current directory:
    - lambert_fix_demonstration.png
    - lambert_half_sky_demonstration.png
"""
import numpy as np
import healpy as hp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os

print("="*70)
print("Lambert Projection Fix - Visual Demonstration Generator")
print("="*70)
print()

# Create output directory if it doesn't exist
output_dir = os.getcwd()
print(f"Output directory: {output_dir}")
print()

# Create a test map with clear north-south gradient
print("Creating test maps...")
nside = 32
npix = hp.nside2npix(nside)
theta, phi = hp.pix2ang(nside, np.arange(npix))
lat = 90 - np.degrees(theta)

# Map 1: Simple latitude gradient
map1 = lat
print("  ✓ Latitude gradient map")

# Map 2: Step function (positive in north, negative in south)
map2 = np.where(lat > 0, 1.0, -1.0)
print("  ✓ North/South step function map")

# Map 3: Sinusoidal pattern in latitude
map3 = np.sin(np.radians(lat * 2))
print("  ✓ Sinusoidal pattern map")
print()

# ===== Figure 1: Comprehensive Comparison =====
print("Generating Figure 1: Projection Comparison...")
fig = plt.figure(figsize=(15, 12))
gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)

# Row 1: Latitude gradient
print("  Rendering Row 1: Latitude gradient...")
ax1 = fig.add_subplot(gs[0, 0], projection='mollweide')
hp.newvisufunc.projview(
    map1,
    projection_type="mollweide",
    title="Mollweide: Latitude gradient",
    cbar=True,
    hold=True,
    min=-90,
    max=90,
)

ax2 = fig.add_subplot(gs[0, 1], projection='hammer')
hp.newvisufunc.projview(
    map1,
    projection_type="hammer",
    title="Hammer: Latitude gradient",
    cbar=True,
    hold=True,
    min=-90,
    max=90,
)

ax3 = fig.add_subplot(gs[0, 2], projection='lambert')
hp.newvisufunc.projview(
    map1,
    projection_type="lambert",
    title="Lambert: Latitude gradient",
    cbar=True,
    hold=True,
    min=-90,
    max=90,
)

# Row 2: Step function
print("  Rendering Row 2: North/South step...")
ax4 = fig.add_subplot(gs[1, 0], projection='mollweide')
hp.newvisufunc.projview(
    map2,
    projection_type="mollweide",
    title="Mollweide: North/South step",
    cbar=True,
    hold=True,
    min=-1,
    max=1,
)

ax5 = fig.add_subplot(gs[1, 1], projection='hammer')
hp.newvisufunc.projview(
    map2,
    projection_type="hammer",
    title="Hammer: North/South step",
    cbar=True,
    hold=True,
    min=-1,
    max=1,
)

ax6 = fig.add_subplot(gs[1, 2], projection='lambert')
hp.newvisufunc.projview(
    map2,
    projection_type="lambert",
    title="Lambert: North/South step",
    cbar=True,
    hold=True,
    min=-1,
    max=1,
)

# Row 3: Sinusoidal
print("  Rendering Row 3: Sinusoidal pattern...")
ax7 = fig.add_subplot(gs[2, 0], projection='mollweide')
hp.newvisufunc.projview(
    map3,
    projection_type="mollweide",
    title="Mollweide: Sin(2*lat)",
    cbar=True,
    hold=True,
    min=-1,
    max=1,
)

ax8 = fig.add_subplot(gs[2, 1], projection='hammer')
hp.newvisufunc.projview(
    map3,
    projection_type="hammer",
    title="Hammer: Sin(2*lat)",
    cbar=True,
    hold=True,
    min=-1,
    max=1,
)

ax9 = fig.add_subplot(gs[2, 2], projection='lambert')
hp.newvisufunc.projview(
    map3,
    projection_type="lambert",
    title="Lambert: Sin(2*lat)",
    cbar=True,
    hold=True,
    min=-1,
    max=1,
)

plt.suptitle("Lambert Projection Fix Demonstration\nAll projections should show consistent patterns", 
             fontsize=16, y=0.98)

output_file1 = os.path.join(output_dir, 'lambert_fix_demonstration.png')
plt.savefig(output_file1, dpi=150, bbox_inches='tight')
print(f"  ✓ Saved: {output_file1}")
plt.close()

# ===== Figure 2: Half-Sky Comparison =====
print("\nGenerating Figure 2: Half-Sky Projections...")
fig2 = plt.figure(figsize=(12, 8))

print("  Rendering Northern hemisphere...")
plt.subplot(2, 2, 1, projection='lambert')
hp.newvisufunc.projview(
    map1,
    projection_type="lambert",
    latra=[0, 90],
    title="Lambert: Northern hemisphere",
    cbar=True,
    hold=True,
    min=-90,
    max=90,
)

print("  Rendering Southern hemisphere...")
plt.subplot(2, 2, 2, projection='lambert')
hp.newvisufunc.projview(
    map1,
    projection_type="lambert",
    latra=[-90, 0],
    title="Lambert: Southern hemisphere",
    cbar=True,
    hold=True,
    min=-90,
    max=90,
)

print("  Rendering Northern hemisphere (step function)...")
plt.subplot(2, 2, 3, projection='lambert')
hp.newvisufunc.projview(
    map2,
    projection_type="lambert",
    latra=[0, 90],
    title="Lambert North: Step function",
    cbar=True,
    hold=True,
    min=-1,
    max=1,
)

print("  Rendering Southern hemisphere (step function)...")
plt.subplot(2, 2, 4, projection='lambert')
hp.newvisufunc.projview(
    map2,
    projection_type="lambert",
    latra=[-90, 0],
    title="Lambert South: Step function",
    cbar=True,
    hold=True,
    min=-1,
    max=1,
)

plt.suptitle("Lambert Half-Sky Projections", fontsize=16, y=0.98)

output_file2 = os.path.join(output_dir, 'lambert_half_sky_demonstration.png')
plt.savefig(output_file2, dpi=150, bbox_inches='tight')
print(f"  ✓ Saved: {output_file2}")
plt.close()

print()
print("="*70)
print("DEMONSTRATION IMAGES GENERATED SUCCESSFULLY")
print("="*70)
print()
print("Generated files:")
print(f"  1. {output_file1}")
print(f"  2. {output_file2}")
print()
print("Expected results:")
print("  ✓ All three projections (Mollweide, Hammer, Lambert) show same patterns")
print("  ✓ Latitude gradient: smooth transition from blue (south) to red (north)")
print("  ✓ Step function: clear separation at equator")
print("  ✓ Sinusoidal: symmetric wave patterns")
print("  ✓ Half-sky: appropriate hemispheres displayed correctly")
print()
print("Verification:")
print("  • No corruption in lower or upper half of Lambert projection")
print("  • Lambert matches Mollweide and Hammer outputs")
print("  • Both full-sky and half-sky modes work correctly")
print()
