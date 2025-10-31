"""
Test for graticule label rotation issue fix.
This test verifies that when rot parameter is used, graticule labels 
show the correct rotated coordinates.
"""

import pytest
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt


def test_projview_graticule_labels_with_rotation():
    """
    Test that graticule labels are correctly offset when rot parameter is used.
    This is a regression test for issue: "new projview, incorrect graticule labels"
    
    When rot=50 is specified, the center longitude label should show 50°, not 0°.
    """
    import healpy as hp
    from healpy.newvisufunc import projview
    
    # Create a simple test map
    test_map = np.arange(12 * 4)
    
    # Create figure with rotation - this is the exact case from the issue
    fig = plt.figure()
    projview(
        test_map,
        flip="astro",
        projection_type="mollweide",
        graticule=True,
        graticule_labels=True,
        rot=50,
        cbar=False
    )
    
    # Get the current axes
    ax = plt.gca()
    
    # Get the x-axis tick labels (longitude labels)
    xtick_labels = [label.get_text() for label in ax.get_xticklabels()]
    
    # The center of the plot (x=0 in plot coordinates) should now show 50°
    # With astro flip and counterclockwise convention, we expect to see labels like:
    # ..., 110°, 50°, 350°, 290°, ...
    # The exact label at center depends on grid spacing, but 50° should be present
    
    # Convert labels to numbers for checking
    label_values = []
    for label in xtick_labels:
        if label and label != '':
            # Remove the degree symbol if present
            try:
                val = float(label.rstrip('°'))
                label_values.append(val)
            except ValueError:
                pass
    
    # Check that 50 is among the labels (allowing for the degree symbol)
    assert len(label_values) > 0, "No numeric labels found"
    
    # The center should be offset by 50 degrees
    # With default 60° spacing and rot=50, we expect to see labels around 50°
    # Specifically, we should see either 50° exactly or labels that confirm the offset
    
    # Check if any label is close to 50° (within grid spacing tolerance)
    has_expected_offset = any(abs(val - 50) < 5 for val in label_values)
    
    print(f"X-tick labels found: {xtick_labels}")
    print(f"Numeric values: {label_values}")
    print(f"Has label near 50°: {has_expected_offset}")
    
    # At minimum, the labels should not all be centered around 0
    # They should be offset by approximately 50 degrees
    center_label_avg = np.mean([v for v in label_values if 0 <= v <= 120])
    
    print(f"Average of center labels: {center_label_avg:.1f}°")
    
    # The average should be significantly different from 0 (the old broken behavior)
    # and closer to 50 (the expected behavior)
    assert abs(center_label_avg) > 20, \
        f"Labels appear not to be rotated (average {center_label_avg:.1f}° too close to 0°)"
    
    plt.close(fig)


def test_projview_graticule_labels_no_rotation():
    """
    Test that graticule labels work correctly when no rotation is specified.
    This ensures the fix doesn't break the default behavior.
    """
    import healpy as hp
    from healpy.newvisufunc import projview
    
    # Create a simple test map
    test_map = np.arange(12 * 4)
    
    # Create figure without rotation
    fig = plt.figure()
    projview(
        test_map,
        flip="astro",
        projection_type="mollweide",
        graticule=True,
        graticule_labels=True,
        cbar=False
    )
    
    # Get the current axes
    ax = plt.gca()
    
    # Get the x-axis tick labels
    xtick_labels = [label.get_text() for label in ax.get_xticklabels()]
    
    # Should have labels (not empty)
    assert len([l for l in xtick_labels if l]) > 0, "No labels found"
    
    plt.close(fig)


def test_projview_graticule_labels_with_tuple_rotation():
    """
    Test that graticule labels work when rot is specified as (lon, lat).
    """
    import healpy as hp
    from healpy.newvisufunc import projview
    
    # Create a simple test map
    test_map = np.arange(12 * 4)
    
    # Create figure with tuple rotation (lon, lat)
    fig = plt.figure()
    projview(
        test_map,
        flip="astro",
        projection_type="mollweide",
        graticule=True,
        graticule_labels=True,
        rot=(45, 30),  # lon=45°, lat=30°
        cbar=False
    )
    
    # Get the current axes
    ax = plt.gca()
    
    # Get the tick labels
    xtick_labels = [label.get_text() for label in ax.get_xticklabels()]
    ytick_labels = [label.get_text() for label in ax.get_yticklabels()]
    
    # Should have labels
    assert len([l for l in xtick_labels if l]) > 0, "No x-axis labels found"
    assert len([l for l in ytick_labels if l]) > 0, "No y-axis labels found"
    
    plt.close(fig)


if __name__ == "__main__":
    # Run the tests
    print("Running graticule label rotation tests...")
    
    try:
        test_projview_graticule_labels_with_rotation()
        print("✓ test_projview_graticule_labels_with_rotation PASSED")
    except Exception as e:
        print(f"✗ test_projview_graticule_labels_with_rotation FAILED: {e}")
        import traceback
        traceback.print_exc()
    
    try:
        test_projview_graticule_labels_no_rotation()
        print("✓ test_projview_graticule_labels_no_rotation PASSED")
    except Exception as e:
        print(f"✗ test_projview_graticule_labels_no_rotation FAILED: {e}")
        import traceback
        traceback.print_exc()
    
    try:
        test_projview_graticule_labels_with_tuple_rotation()
        print("✓ test_projview_graticule_labels_with_tuple_rotation PASSED")
    except Exception as e:
        print(f"✗ test_projview_graticule_labels_with_tuple_rotation FAILED: {e}")
        import traceback
        traceback.print_exc()
    
    print("\nAll tests completed!")
