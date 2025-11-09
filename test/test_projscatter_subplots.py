"""
Test that projscatter, projplot, and projtext only draw on the current axes,
not on all axes in the figure.

This is a regression test for issue #637 and the stackoverflow question:
https://stackoverflow.com/q/76796703/5838180
"""
import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp


def test_projscatter_multiple_subplots():
    """Test that projscatter only draws on the current subplot"""
    fig, (ax1, ax2) = plt.subplots(1, 2)

    # Create first subplot with a point at (0, 45)
    plt.axes(ax1)
    hp.mollview(hold=True)
    hp.projscatter([0], [45], lonlat=True, coord="G")

    # Create second subplot with a point at (90, 45)
    plt.axes(ax2)
    hp.mollview(hold=True)
    hp.projscatter([90], [45], lonlat=True, coord="G")

    # Get the axes
    axes = fig.get_axes()
    spherical_axes = [ax for ax in axes if isinstance(ax, hp.projaxes.SphericalProjAxes)]

    # Verify we have two SphericalProjAxes
    assert len(spherical_axes) == 2, "Should have two SphericalProjAxes"

    # Check that each axes has exactly one scatter collection
    # (the scatter plot on the current axes only)
    # Before the fix, both scatter calls would appear on both axes (2 collections each)
    # After the fix, each axes should have exactly 1 collection
    for ax in spherical_axes:
        collections = ax.collections
        assert len(collections) == 1, (
            f"Each axes should have exactly 1 collection, got {len(collections)}. "
            "If this is 2, the bug has regressed - scatter points are appearing on all axes."
        )

    plt.close(fig)


def test_projplot_multiple_subplots():
    """Test that projplot only draws on the current subplot"""
    fig, (ax1, ax2) = plt.subplots(1, 2)

    # Create first subplot with a line
    plt.axes(ax1)
    hp.mollview(hold=True)
    hp.projplot([0, 10], [45, 50], lonlat=True, coord="G")

    # Create second subplot with a different line
    plt.axes(ax2)
    hp.mollview(hold=True)
    hp.projplot([90, 100], [45, 50], lonlat=True, coord="G")

    # Get the axes
    axes = fig.get_axes()
    spherical_axes = [ax for ax in axes if isinstance(ax, hp.projaxes.SphericalProjAxes)]

    # Verify we have two SphericalProjAxes
    assert len(spherical_axes) == 2, "Should have two SphericalProjAxes"

    # Check that each axes has exactly one line from its own projplot call
    # Before the fix, both projplot calls would appear on both axes (2 lines each)
    # After the fix, each axes should have exactly 1 line
    for ax in spherical_axes:
        lines = ax.lines
        assert len(lines) == 1, (
            f"Each axes should have exactly 1 line, got {len(lines)}. "
            "If this is 2, the bug has regressed - lines are appearing on all axes."
        )

    plt.close(fig)


def test_projtext_multiple_subplots():
    """Test that projtext only draws on the current subplot"""
    fig, (ax1, ax2) = plt.subplots(1, 2)

    # Create first subplot with text
    plt.axes(ax1)
    hp.mollview(hold=True)
    hp.projtext(0, 45, "Text1", lonlat=True, coord="G")

    # Create second subplot with different text
    plt.axes(ax2)
    hp.mollview(hold=True)
    hp.projtext(90, 45, "Text2", lonlat=True, coord="G")

    # Get the axes
    axes = fig.get_axes()
    spherical_axes = [ax for ax in axes if isinstance(ax, hp.projaxes.SphericalProjAxes)]

    # Verify we have two SphericalProjAxes
    assert len(spherical_axes) == 2, "Should have two SphericalProjAxes"

    # Check that each axes has exactly one text element from its own projtext call
    # Before the fix, both projtext calls would appear on both axes (2 texts each)
    # After the fix, each axes should have exactly 1 text
    for ax in spherical_axes:
        texts = ax.texts
        assert len(texts) == 1, (
            f"Each axes should have exactly 1 text element, got {len(texts)}. "
            "If this is 2, the bug has regressed - text is appearing on all axes."
        )

    plt.close(fig)
