"""
Regression tests for colormap handling through create_colormap.

When badcolor/bgcolor are explicitly provided, they should be applied
consistently for both string and Colormap-object inputs.
"""
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
import pytest


class TestColormapSegmentdataPreservation:
    """Tests to verify explicit bad/under colors are applied consistently."""
    
    def test_linearsegmented_colormap_applies_bad_color(self):
        """Test that bad color argument is applied for LinearSegmentedColormap objects."""
        from healpy.projaxes import create_colormap
        
        # Use 'jet' which is a LinearSegmentedColormap
        cmap = plt.get_cmap('jet').copy()
        assert hasattr(cmap, '_segmentdata'), "jet should have _segmentdata"
        
        # Set custom bad color, but explicit args should take precedence.
        cmap.set_bad('white')
        
        # Pass through create_colormap
        result = create_colormap(cmap, badcolor='gray', bgcolor='white')
        
        # Verify explicit badcolor argument is applied
        assert result._rgba_bad is not None, "Bad color should not be None"
        assert np.allclose(result._rgba_bad[:3], [0.5019607843137255, 0.5019607843137255, 0.5019607843137255]), \
            f"Expected gray from badcolor arg, got {result._rgba_bad[:3]}"
    
    def test_linearsegmented_colormap_applies_under_color(self):
        """Test that under color argument is applied for LinearSegmentedColormap objects."""
        from healpy.projaxes import create_colormap
        
        # Use 'jet' which is a LinearSegmentedColormap
        cmap = plt.get_cmap('jet').copy()
        assert hasattr(cmap, '_segmentdata'), "jet should have _segmentdata"
        
        # Set custom under color, but explicit args should take precedence.
        cmap.set_under('yellow')
        
        # Pass through create_colormap
        result = create_colormap(cmap, badcolor='gray', bgcolor='white')
        
        # Verify explicit bgcolor argument is applied
        assert result._rgba_under is not None, "Under color should not be None"
        assert np.allclose(result._rgba_under[:3], [1.0, 1.0, 1.0]), \
            f"Expected white from bgcolor arg, got {result._rgba_under[:3]}"
    
    def test_linearsegmented_colormap_applies_both_colors(self):
        """Test that explicit bad/under colors are applied together."""
        from healpy.projaxes import create_colormap
        
        # Use 'hot' which is also a LinearSegmentedColormap
        cmap = plt.get_cmap('hot').copy()
        assert hasattr(cmap, '_segmentdata'), "hot should have _segmentdata"
        
        # Set both custom colors, but explicit args should take precedence.
        cmap.set_bad('cyan')
        cmap.set_under('magenta')
        
        # Pass through create_colormap
        result = create_colormap(cmap, badcolor='gray', bgcolor='white')
        
        # Verify explicit colors are applied
        assert result._rgba_bad is not None, "Bad color should not be None"
        assert result._rgba_under is not None, "Under color should not be None"
        assert np.allclose(result._rgba_bad[:3], [0.5019607843137255, 0.5019607843137255, 0.5019607843137255]), \
            f"Expected gray from badcolor arg, got {result._rgba_bad[:3]}"
        assert np.allclose(result._rgba_under[:3], [1.0, 1.0, 1.0]), \
            f"Expected white from bgcolor arg, got {result._rgba_under[:3]}"
    
    def test_listed_colormap_still_works(self):
        """Verify that ListedColormap also applies explicit colors."""
        from healpy.projaxes import create_colormap
        
        # Use 'viridis' which is a ListedColormap
        cmap = plt.get_cmap('viridis').copy()
        assert not hasattr(cmap, '_segmentdata'), "viridis should not have _segmentdata"
        
        # Set custom colors
        cmap.set_bad('white')
        cmap.set_under('yellow')
        
        # Pass through create_colormap
        result = create_colormap(cmap, badcolor='gray', bgcolor='white')
        
        # Verify explicit colors are applied
        assert np.allclose(result._rgba_bad[:3], [0.5019607843137255, 0.5019607843137255, 0.5019607843137255])
        assert np.allclose(result._rgba_under[:3], [1.0, 1.0, 1.0])
