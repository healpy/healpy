"""Test legend functionality with newvisufunc

This test verifies that legends work correctly with newprojplot
when an explicit location is specified.
"""
import healpy as hp
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for testing
import matplotlib.pyplot as plt
import pytest


def test_newprojplot_legend_with_explicit_loc():
    """Test that legend works when loc is explicitly specified"""
    # Create a simple map
    nside = 16
    npix = hp.nside2npix(nside)
    m = np.random.random(npix)
    
    # Create projection view
    hp.newvisufunc.projview(m, projection_type='mollweide')
    
    # Add some plot elements with labels
    theta1 = np.linspace(0, np.pi, 50)
    phi1 = np.zeros_like(theta1) + np.pi/2
    hp.newprojplot(theta1, phi1, 'r-', linewidth=2, label='Meridian')
    
    theta2 = np.zeros(50) + np.pi/2
    phi2 = np.linspace(-np.pi, np.pi, 50)
    hp.newprojplot(theta2, phi2, 'b--', linewidth=2, label='Equator')
    
    # Create legend with explicit location - this should work without hanging
    legend = plt.legend(loc='upper right')
    
    # Verify legend was created
    assert legend is not None
    assert len(legend.get_texts()) == 2
    
    plt.savefig('/tmp/test_legend_with_loc.png', dpi=50, bbox_inches='tight')
    plt.close()


def test_newprojplot_legend_with_bbox():
    """Test that legend works with bbox_to_anchor"""
    # Create a simple map
    nside = 16
    npix = hp.nside2npix(nside)
    m = np.random.random(npix)
    
    # Create projection view
    hp.newvisufunc.projview(m, projection_type='hammer')
    
    # Add plot element with label
    theta = np.linspace(0, np.pi, 100)
    phi = np.sin(theta) * np.pi
    hp.newprojplot(theta, phi, 'g-', linewidth=2, label='Sine curve')
    
    # Create legend with bbox_to_anchor - this should also work
    legend = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    assert legend is not None
    plt.savefig('/tmp/test_legend_with_bbox.png', dpi=50, bbox_inches='tight')
    plt.close()


def test_newprojplot_returns_line_objects():
    """Test that newprojplot returns proper Line2D objects"""
    nside = 16
    npix = hp.nside2npix(nside)
    m = np.random.random(npix)
    
    hp.newvisufunc.projview(m, projection_type='mollweide')
    
    # Call newprojplot and check return value
    theta = np.array([0.5, 1.0, 1.5])
    phi = np.array([0.0, 0.5, 1.0])
    lines = hp.newprojplot(theta, phi, 'ro-', label='Test')
    
    # Verify return type
    assert isinstance(lines, list)
    assert len(lines) == 1
    
    # Verify the line has the label
    assert lines[0].get_label() == 'Test'
    
    plt.close()


def test_issue_example_with_explicit_loc():
    """Test the example from the issue with explicit loc parameter
    
    This reproduces the issue scenario but with the fix applied.
    """
    # Reproduce the issue example
    hp.newvisufunc.projview(np.random.random(12*16**2))
    
    # Make a circle
    circle = np.zeros(10000) + np.radians(45), np.arange(-180, 180, 0.036)
    
    # Rotate it
    rotMat = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]
    rotCircle = hp.rotator.rotateDirection(rotMat, circle)
    
    # Add it to the plot
    hp.newprojplot(*rotCircle, linewidth=1, linestyle="", marker=".", label="a circle", c="r")
    
    # Add legend with explicit location - this should work
    legend = plt.legend(loc='upper right')
    
    assert legend is not None
    plt.savefig('/tmp/test_issue_example.png', dpi=50, bbox_inches='tight')
    plt.close()
