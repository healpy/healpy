import matplotlib

matplotlib.use("agg")
import unittest
import matplotlib.pyplot as plt
import numpy as np
import copy

import healpy as hp
from healpy.visufunc import *
from healpy.zoomtool import mollzoom


class TestNoCrash(unittest.TestCase):
    def setUp(self):
        self.nside = 16
        self.npix = hp.nside2npix(self.nside)
        self.m = np.arange(self.npix, dtype=np.double)
        self.m_wrong_npix = np.arange(self.npix + 1, dtype=np.double)
        self.ma = self.m.copy()
        self.ma[3] = hp.UNSEEN
        self.ma = hp.ma(self.ma)
        self.m2 = self.m.copy()
        self.m2[100:] = hp.UNSEEN
        self.ma2 = hp.ma(self.m2)

    def test_cartview_nocrash(self):
        cartview(self.m)

    def test_cartview_no_arguments(self):
        # Test that cartview() can be called with no arguments
        # This should create a blank map without errors
        cartview()

    def test_cartview_crash(self):
        with self.assertRaises(TypeError):
            cartview(self.m_wrong_npix)

    def test_mollview_nocrash(self):
        mollview(self.m)

    def test_mollview_no_arguments(self):
        # Test that mollview() can be called with no arguments
        # This should create a blank map without errors
        mollview()

    def test_mollview_crash(self):
        with self.assertRaises(TypeError):
            mollview(self.m_wrong_npix)

    def test_gnomview_nocrash(self):
        gnomview(self.m)

    def test_gnomview_no_arguments(self):
        # Test that gnomview() can be called with no arguments
        # This should create a blank map without errors
        gnomview()

    def test_gnomview_crash(self):
        with self.assertRaises(TypeError):
            gnomview(self.m_wrong_npix)

    def test_orthview_nocrash(self):
        orthview(self.m)

    def test_orthview_no_arguments(self):
        # Test that orthview() can be called with no arguments
        # This should create a blank map without errors
        orthview()

    def test_orthview_crash(self):
        with self.assertRaises(TypeError):
            orthview(self.m_wrong_npix)

    def test_azeqview_nocrash(self):
        azeqview(self.m)

    def test_azeqview_no_arguments(self):
        # Test that azeqview() can be called with no arguments
        # This should create a blank map without errors
        azeqview()

    def test_azeqview_crash(self):
        with self.assertRaises(TypeError):
            azeqview(self.m_wrong_npix)

    def test_mollzoom_nocrash(self):
        mollzoom(self.m)

    def test_mollzoom_histnocrash(self):
        mollzoom(self.m, norm="hist")

    def test_mollzoom_crash(self):
        with self.assertRaises(TypeError):
            mollzoom(self.m_wrong_npix)

    def test_cartview_ma_nocrash(self):
        cartview(self.ma)

    def test_mollview_ma_nocrash(self):
        mollview(self.ma)

    def test_gnomview_ma_nocrash(self):
        gnomview(self.ma)

    def test_orthview_ma_nocrash(self):
        orthview(self.ma)

    def test_mollzoom_ma_nocrash(self):
        mollzoom(self.ma)

    def test_mollzoom_ma_hist_nocrash(self):
        mollzoom(self.ma2, norm="hist")

    def test_azeqview_ma_nocrash(self):
        azeqview(self.ma)

    def test_cmap_colors(self):
        # Get a built-in colormap
        name = 'viridis'
        cmap = copy.copy(plt.get_cmap(name))

        # Set outlier colors
        color = (0.25,0.75,0.95,1.0)
        cmap.set_bad(color)
        cmap.set_over(color)
        cmap.set_under(color)

        # Call healpy plotting
        mollview(self.m,cmap=name)

        # Check that colors haven't been changed
        assert cmap._rgba_bad == color
        assert cmap._rgba_over == color
        assert cmap._rgba_under == color

    def test_reuse_axes(self):
        cartview(self.m)
        cartview(self.m, reuse_axes=True)
        mollview(self.m)
        mollview(self.m, reuse_axes=True)
        gnomview(self.m)
        gnomview(self.m, reuse_axes=True)
        orthview(self.m)
        orthview(self.m, reuse_axes=True)
        azeqview(self.m)
        azeqview(self.m, reuse_axes=True)

    def test_colormap_object_preservation(self):
        """Test that user-modified Colormap objects preserve their colors"""
        from healpy.projaxes import create_colormap
        
        # Test 1: String colormap should apply badcolor/bgcolor
        cm1 = create_colormap('viridis', badcolor='red', bgcolor='blue')
        bad_is_red = np.allclose(cm1._rgba_bad[:3], [1.0, 0.0, 0.0])
        under_is_blue = np.allclose(cm1._rgba_under[:3], [0.0, 0.0, 1.0])
        assert bad_is_red, "String colormap should apply badcolor"
        assert under_is_blue, "String colormap should apply bgcolor"
        
        # Test 2: Colormap object with pre-set colors should preserve them
        cm_obj = copy.copy(plt.get_cmap('viridis'))
        cm_obj.set_bad('white')
        cm_obj.set_under('yellow')
        
        cm2 = create_colormap(cm_obj, badcolor='red', bgcolor='blue')
        bad_is_white = np.allclose(cm2._rgba_bad[:3], [1.0, 1.0, 1.0])
        under_is_yellow = np.allclose(cm2._rgba_under[:3], [1.0, 1.0, 0.0])
        assert bad_is_white, "Colormap object should preserve bad color"
        assert under_is_yellow, "Colormap object should preserve under color"
        
        # Test 3: None should apply badcolor/bgcolor (backward compatibility)
        cm3 = create_colormap(None, badcolor='red', bgcolor='blue')
        bad_is_red = np.allclose(cm3._rgba_bad[:3], [1.0, 0.0, 0.0])
        under_is_blue = np.allclose(cm3._rgba_under[:3], [0.0, 0.0, 1.0])
        assert bad_is_red, "None colormap should apply badcolor"
        assert under_is_blue, "None colormap should apply bgcolor"
    
    def test_mollview_with_colormap_object(self):
        """Test that mollview preserves user-modified Colormap object colors"""
        # Create a colormap with custom bad/under colors
        cmap = copy.copy(plt.get_cmap('viridis'))
        cmap.set_bad('white')
        cmap.set_under('yellow')
        
        # Call mollview with the colormap object
        mollview(self.m, cmap=cmap)
        
        # Get the colormap from the plot
        ax = plt.gca()
        if hasattr(ax, 'images') and len(ax.images) > 0:
            plot_cmap = ax.images[0].get_cmap()
            # Verify colors are preserved
            bad_is_white = np.allclose(plot_cmap._rgba_bad[:3], [1.0, 1.0, 1.0])
            under_is_yellow = np.allclose(plot_cmap._rgba_under[:3], [1.0, 1.0, 0.0])
            assert bad_is_white, "mollview should preserve Colormap bad color"
            assert under_is_yellow, "mollview should preserve Colormap under color"
