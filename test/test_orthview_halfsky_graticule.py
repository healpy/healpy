"""
Test for orthview with half_sky=True and graticule.
This test verifies the fix for the issue where graticule() failed with 
ValueError when used with orthview(half_sky=True).
"""

import matplotlib

matplotlib.use("agg")
import unittest
import matplotlib.pyplot as plt
import numpy as np

import healpy as hp


class TestOrthviewHalfSkyGraticule(unittest.TestCase):
    """Test that graticule works correctly with orthview half_sky parameter"""

    def setUp(self):
        """Set up test fixtures"""
        self.nside = 16
        self.npix = hp.nside2npix(self.nside)
        self.test_map = np.arange(self.npix, dtype=np.double)

    def tearDown(self):
        """Clean up after tests"""
        plt.close("all")

    def test_orthview_fullsky_graticule(self):
        """Test that graticule works with full sky orthview (baseline)"""
        # This should work and serves as a baseline
        hp.orthview(self.test_map, half_sky=False, title="Full Sky Test")
        hp.graticule()
        # If we get here without exception, test passes
        self.assertTrue(True)

    def test_orthview_halfsky_graticule(self):
        """Test that graticule works with half sky orthview (bug fix)"""
        # This was failing before the fix with:
        # ValueError: array of sample points is empty
        hp.orthview(self.test_map, half_sky=True, title="Half Sky Test")
        hp.graticule()
        # If we get here without exception, test passes
        self.assertTrue(True)

    def test_orthview_halfsky_graticule_with_rotation(self):
        """Test graticule with half sky orthview and rotation"""
        hp.orthview(
            self.test_map, half_sky=True, rot=[90, 0], title="Half Sky Rotated"
        )
        hp.graticule()
        self.assertTrue(True)

    def test_orthview_halfsky_graticule_custom_intervals(self):
        """Test graticule with half sky orthview and custom intervals"""
        hp.orthview(self.test_map, half_sky=True, title="Half Sky Custom Graticule")
        hp.graticule(dpar=30, dmer=45)
        self.assertTrue(True)

    def test_orthview_halfsky_multiple_graticules(self):
        """Test multiple half sky orthviews with graticules in same figure"""
        fig = plt.figure(figsize=(12, 4))
        
        # First subplot - default
        hp.orthview(self.test_map, half_sky=True, fig=fig, sub=(1, 3, 1), 
                   title="Half Sky 1")
        hp.graticule()
        
        # Second subplot - rotated
        hp.orthview(self.test_map, half_sky=True, fig=fig, sub=(1, 3, 2),
                   rot=[90, 0], title="Half Sky 2")
        hp.graticule()
        
        # Third subplot - different latitude
        hp.orthview(self.test_map, half_sky=True, fig=fig, sub=(1, 3, 3),
                   rot=[0, 45], title="Half Sky 3")
        hp.graticule()
        
        self.assertTrue(True)


if __name__ == "__main__":
    unittest.main()
