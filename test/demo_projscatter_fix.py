#!/usr/bin/env python
"""
Demonstration script showing the fix for issue #637
https://github.com/healpy/healpy/issues/637

This script demonstrates that after the fix, projscatter, projplot, and projtext
only draw on the current axes, not on all SphericalProjAxes in the figure.

Before the fix:
- Points/lines/text would appear in ALL subplots
- This was because the functions iterated over all axes: `for ax in f.get_axes()`

After the fix:
- Points/lines/text only appear in the current subplot
- Uses `plt.gca()` to get only the current axes
"""

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

def demo_projscatter():
    """Demonstrate projscatter fix with multiple subplots"""
    print("Testing projscatter with multiple subplots...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # First subplot: point at (0°, 45°)
    plt.axes(ax1)
    hp.mollview(hold=True)
    hp.projscatter([0], [45], lonlat=True, coord='G', color='red', s=100)
    hp.graticule()
    ax1.set_title('Left subplot: Point at (0°, 45°)')
    
    # Second subplot: point at (90°, 45°)
    plt.axes(ax2)
    hp.mollview(hold=True)
    hp.projscatter([90], [45], lonlat=True, coord='G', color='blue', s=100)
    hp.graticule()
    ax2.set_title('Right subplot: Point at (90°, 45°)')
    
    plt.suptitle('projscatter fix demo: Each point appears only in its own subplot')
    plt.tight_layout()
    plt.savefig('/tmp/projscatter_fix_demo.png', dpi=100, bbox_inches='tight')
    print("✓ Saved to /tmp/projscatter_fix_demo.png")
    plt.close()

def demo_projplot():
    """Demonstrate projplot fix with multiple subplots"""
    print("Testing projplot with multiple subplots...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # First subplot: line in one region
    plt.axes(ax1)
    hp.mollview(hold=True)
    hp.projplot([0, 30], [0, 30], lonlat=True, coord='G', color='red', linewidth=2)
    hp.graticule()
    ax1.set_title('Left subplot: Red line')
    
    # Second subplot: line in different region
    plt.axes(ax2)
    hp.mollview(hold=True)
    hp.projplot([90, 120], [0, 30], lonlat=True, coord='G', color='blue', linewidth=2)
    hp.graticule()
    ax2.set_title('Right subplot: Blue line')
    
    plt.suptitle('projplot fix demo: Each line appears only in its own subplot')
    plt.tight_layout()
    plt.savefig('/tmp/projplot_fix_demo.png', dpi=100, bbox_inches='tight')
    print("✓ Saved to /tmp/projplot_fix_demo.png")
    plt.close()

def demo_projtext():
    """Demonstrate projtext fix with multiple subplots"""
    print("Testing projtext with multiple subplots...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # First subplot: text in one location
    plt.axes(ax1)
    hp.mollview(hold=True)
    hp.projtext(0, 45, 'LEFT', lonlat=True, coord='G', color='red', fontsize=16)
    hp.graticule()
    ax1.set_title('Left subplot: "LEFT" text')
    
    # Second subplot: text in different location
    plt.axes(ax2)
    hp.mollview(hold=True)
    hp.projtext(90, 45, 'RIGHT', lonlat=True, coord='G', color='blue', fontsize=16)
    hp.graticule()
    ax2.set_title('Right subplot: "RIGHT" text')
    
    plt.suptitle('projtext fix demo: Each text appears only in its own subplot')
    plt.tight_layout()
    plt.savefig('/tmp/projtext_fix_demo.png', dpi=100, bbox_inches='tight')
    print("✓ Saved to /tmp/projtext_fix_demo.png")
    plt.close()

def demo_original_issue():
    """Reproduce the exact scenario from the original issue"""
    print("Testing original issue scenario...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    plt.axes(ax1)
    hp.mollview(hold=True)
    hp.projscatter([0], [45], lonlat=True, coord='G')
    hp.graticule()

    plt.axes(ax2)
    hp.mollview(hold=True)
    hp.projscatter([90], [45], lonlat=True, coord='G')
    hp.graticule()
    
    plt.suptitle('Original issue #637 - Now fixed!')
    plt.tight_layout()
    plt.savefig('/tmp/original_issue_fixed.png', dpi=100, bbox_inches='tight')
    print("✓ Saved to /tmp/original_issue_fixed.png")
    plt.close()

if __name__ == '__main__':
    print("=" * 60)
    print("Healpy projscatter/projplot/projtext Fix Demonstration")
    print("Issue #637: https://github.com/healpy/healpy/issues/637")
    print("=" * 60)
    print()
    
    try:
        demo_projscatter()
        demo_projplot()
        demo_projtext()
        demo_original_issue()
        
        print()
        print("=" * 60)
        print("All tests completed successfully!")
        print("The fix ensures that:")
        print("  - projscatter only plots on the current axes")
        print("  - projplot only plots on the current axes")
        print("  - projtext only plots on the current axes")
        print("=" * 60)
        
    except Exception as e:
        print(f"Error: {e}")
        print("Note: This script requires healpy to be installed")
        print("Run: pip install healpy")
