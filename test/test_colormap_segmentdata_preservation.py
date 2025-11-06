"""
Test for issue raised in PR review: Colormap rebuild from _segmentdata loses user colors.

When a LinearSegmentedColormap object (which has _segmentdata) is passed with 
custom bad/under colors, the reconstruction from _segmentdata resets those colors
to defaults, even though the intent was to preserve them.
"""
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
import pytest


class TestColormapSegmentdataPreservation:
    """Tests to verify that LinearSegmentedColormap colors are preserved."""
    
    def test_linearsegmented_colormap_preserves_bad_color(self):
        """Test that bad color is preserved for LinearSegmentedColormap objects."""
        from healpy.projaxes import create_colormap
        
        # Use 'jet' which is a LinearSegmentedColormap
        cmap = plt.get_cmap('jet').copy()
        assert hasattr(cmap, '_segmentdata'), "jet should have _segmentdata"
        
        # Set custom bad color
        cmap.set_bad('white')
        original_bad = cmap._rgba_bad
        
        # Pass through create_colormap
        result = create_colormap(cmap, badcolor='gray', bgcolor='white')
        
        # Verify bad color is preserved
        assert result._rgba_bad is not None, "Bad color should not be None"
        assert np.allclose(result._rgba_bad[:3], [1.0, 1.0, 1.0]), \
            f"Expected white (1,1,1), got {result._rgba_bad[:3]}"
    
    def test_linearsegmented_colormap_preserves_under_color(self):
        """Test that under color is preserved for LinearSegmentedColormap objects."""
        from healpy.projaxes import create_colormap
        
        # Use 'jet' which is a LinearSegmentedColormap
        cmap = plt.get_cmap('jet').copy()
        assert hasattr(cmap, '_segmentdata'), "jet should have _segmentdata"
        
        # Set custom under color
        cmap.set_under('yellow')
        original_under = cmap._rgba_under
        
        # Pass through create_colormap
        result = create_colormap(cmap, badcolor='gray', bgcolor='white')
        
        # Verify under color is preserved
        assert result._rgba_under is not None, "Under color should not be None"
        assert np.allclose(result._rgba_under[:3], [1.0, 1.0, 0.0]), \
            f"Expected yellow (1,1,0), got {result._rgba_under[:3]}"
    
    def test_linearsegmented_colormap_preserves_both_colors(self):
        """Test that both bad and under colors are preserved together."""
        from healpy.projaxes import create_colormap
        
        # Use 'hot' which is also a LinearSegmentedColormap
        cmap = plt.get_cmap('hot').copy()
        assert hasattr(cmap, '_segmentdata'), "hot should have _segmentdata"
        
        # Set both custom colors
        cmap.set_bad('cyan')
        cmap.set_under('magenta')
        
        # Pass through create_colormap
        result = create_colormap(cmap, badcolor='gray', bgcolor='white')
        
        # Verify both colors are preserved
        assert result._rgba_bad is not None, "Bad color should not be None"
        assert result._rgba_under is not None, "Under color should not be None"
        assert np.allclose(result._rgba_bad[:3], [0.0, 1.0, 1.0]), \
            f"Expected cyan (0,1,1), got {result._rgba_bad[:3]}"
        assert np.allclose(result._rgba_under[:3], [1.0, 0.0, 1.0]), \
            f"Expected magenta (1,0,1), got {result._rgba_under[:3]}"
    
    def test_listed_colormap_still_works(self):
        """Verify that ListedColormap (no _segmentdata) still works correctly."""
        from healpy.projaxes import create_colormap
        
        # Use 'viridis' which is a ListedColormap
        cmap = plt.get_cmap('viridis').copy()
        assert not hasattr(cmap, '_segmentdata'), "viridis should not have _segmentdata"
        
        # Set custom colors
        cmap.set_bad('white')
        cmap.set_under('yellow')
        
        # Pass through create_colormap
        result = create_colormap(cmap, badcolor='gray', bgcolor='white')
        
        # Verify colors are preserved (this should already work)
        assert np.allclose(result._rgba_bad[:3], [1.0, 1.0, 1.0])
        assert np.allclose(result._rgba_under[:3], [1.0, 1.0, 0.0])
