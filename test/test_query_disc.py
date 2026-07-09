import unittest
import numpy as np

from healpy import query_disc, query_polygon, query_strip, boundaries, nside2npix, nside2resol, ang2vec

try:
    from exceptions import ValueError
except:
    pass


class TestQueryDisc(unittest.TestCase):
    def setUp(self):
        self.NSIDE = 8
        self.vec = np.array([0.17101007, 0.03015369, 0.98480775])
        self.radius = np.radians(6)
        self.nside2_55_corners_precomp = np.array(
            [
                [
                    [2.44708573e-17, 5.27046277e-01, 3.60797400e-01, 4.56383842e-17],
                    [3.99652627e-01, 5.27046277e-01, 8.71041977e-01, 7.45355992e-01],
                    [9.16666667e-01, 6.66666667e-01, 3.33333333e-01, 6.66666667e-01],
                ],
                [
                    [2.44708573e-17, 5.27046277e-01, 3.60797400e-01, 4.56383842e-17],
                    [3.99652627e-01, 5.27046277e-01, 8.71041977e-01, 7.45355992e-01],
                    [9.16666667e-01, 6.66666667e-01, 3.33333333e-01, 6.66666667e-01],
                ],
            ]
        )

    def test_not_inclusive(self):
        # HIDL> query_disc, 8, [ 0.17101007,  0.03015369,  0.98480775],6,listpix,/DEG,NESTED=0
        # HIDL> print,listpix
        #           4
        np.testing.assert_array_equal(
            query_disc(self.NSIDE, self.vec, self.radius, inclusive=False),
            np.array([4]),
        )

    def test_inclusive(self):
        # HIDL> query_disc, 8, [ 0.17101007,  0.03015369,  0.98480775],6,listpix,/DEG,NESTED=0,/inclusive
        # HIDL> print,listpix
        #           0           3           4           5          11          12          13          23
        np.testing.assert_array_equal(
            query_disc(self.NSIDE, self.vec, self.radius, inclusive=True),
            np.array([0, 3, 4, 5, 11, 12, 13, 23]),
        )

    def test_boundaries(self):
        nside = 2
        corners = boundaries(nside, 5)
        corners_precomp = np.array(
            [
                [2.44708573e-17, 5.27046277e-01, 3.60797400e-01, 4.56383842e-17],
                [3.99652627e-01, 5.27046277e-01, 8.71041977e-01, 7.45355992e-01],
                [9.16666667e-01, 6.66666667e-01, 3.33333333e-01, 6.66666667e-01],
            ]
        )
        np.testing.assert_array_almost_equal(corners, corners_precomp, decimal=8)

    def test_boundaries_list(self):
        nside = 2
        corners = boundaries(nside, [5, 5])
        np.testing.assert_array_almost_equal(
            corners, self.nside2_55_corners_precomp, decimal=8
        )

    def test_boundaries_phi_theta(self):
        nside = 2
        corners = boundaries(nside, np.array([5, 5]))
        np.testing.assert_array_almost_equal(
            corners, self.nside2_55_corners_precomp, decimal=8
        )

    def test_nside_non_power_of_two(self):
        # For RING scheme, nside should not need to be a power of two.

        nside = 1
        resolution = 1.0
        theta=0.0
        phi=0.0
        radius = np.radians(1)
        while True:
            nside = nside + 1 
            res = nside2resol(nside, arcmin = True)
            print("nside={} res={} arcmin".format(nside, res))
            if res < resolution:
                break
            
        self.assertEqual(nside, 3518)
        
        x0 = ang2vec(theta, phi)
        
        pixel_indices = query_disc(nside, x0, radius, inclusive=False, nest=False)
        self.assertEqual(pixel_indices.shape[0], 11400)
        
    def test_boundaries_floatpix_array(self):
        self.assertRaises(ValueError, boundaries, 2, np.array([5.0, 5]))

    def test_boundaries_floatpix_scalar(self):
        self.assertRaises(ValueError, boundaries, 2, 1 / 2.0)

    def test_buffer_mode(self):

        # allocate something manifestly too short, should raise a value error
        buff = np.empty(0, dtype=np.int64)
        self.assertRaises(
            ValueError,
            query_disc,
            self.NSIDE,
            self.vec,
            self.radius,
            inclusive=True,
            buff=buff,
        )

        # allocate something of wrong type, should raise a value error
        buff = np.empty(nside2npix(self.NSIDE), dtype=np.float64)
        self.assertRaises(
            ValueError,
            query_disc,
            self.NSIDE,
            self.vec,
            self.radius,
            inclusive=True,
            buff=buff,
        )

        # allocate something acceptable, should succeed and return a subview
        buff = np.empty(nside2npix(self.NSIDE), dtype=np.int64)
        result = query_disc(
            self.NSIDE, self.vec, self.radius, inclusive=True, buff=buff
        )

        assert result.base is buff

        np.testing.assert_array_equal(result, np.array([0, 3, 4, 5, 11, 12, 13, 23]))

    def test_query_disc_return_ranges(self):
        """Test query_disc with return_ranges=True"""
        # Get ranges
        ranges = query_disc(self.NSIDE, self.vec, self.radius, inclusive=True, return_ranges=True)
        
        # Verify it returns a 2D array
        self.assertEqual(ranges.ndim, 2)
        self.assertEqual(ranges.shape[1], 2)
        
        # Convert ranges back to individual pixels to verify correctness
        pixels_from_ranges = []
        for start, end in ranges:
            pixels_from_ranges.extend(range(start, end))
        pixels_from_ranges = np.array(pixels_from_ranges, dtype=np.int64)
        
        # Get pixels the normal way
        pixels_normal = query_disc(self.NSIDE, self.vec, self.radius, inclusive=True, return_ranges=False)
        
        # They should be the same
        np.testing.assert_array_equal(pixels_from_ranges, pixels_normal)

    def test_query_disc_return_ranges_nested(self):
        """Test query_disc with return_ranges=True and nested ordering"""
        # Get ranges for NESTED
        ranges = query_disc(self.NSIDE, self.vec, self.radius, inclusive=True, nest=True, return_ranges=True)
        
        # Verify it returns a 2D array
        self.assertEqual(ranges.ndim, 2)
        self.assertEqual(ranges.shape[1], 2)
        
        # Convert ranges back to individual pixels
        pixels_from_ranges = []
        for start, end in ranges:
            pixels_from_ranges.extend(range(start, end))
        pixels_from_ranges = np.array(pixels_from_ranges, dtype=np.int64)
        
        # Get pixels the normal way
        pixels_normal = query_disc(self.NSIDE, self.vec, self.radius, inclusive=True, nest=True, return_ranges=False)
        
        # They should be the same
        np.testing.assert_array_equal(pixels_from_ranges, pixels_normal)

    def test_query_polygon_return_ranges(self):
        """Test query_polygon with return_ranges=True"""
        # Define a simple triangle
        vertices = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ])
        
        # Get ranges
        ranges = query_polygon(self.NSIDE, vertices, return_ranges=True)
        
        # Verify it returns a 2D array
        self.assertEqual(ranges.ndim, 2)
        self.assertEqual(ranges.shape[1], 2)
        
        # Convert ranges back to individual pixels
        pixels_from_ranges = []
        for start, end in ranges:
            pixels_from_ranges.extend(range(start, end))
        pixels_from_ranges = np.array(pixels_from_ranges, dtype=np.int64)
        
        # Get pixels the normal way
        pixels_normal = query_polygon(self.NSIDE, vertices, return_ranges=False)
        
        # They should be the same
        np.testing.assert_array_equal(pixels_from_ranges, pixels_normal)

    def test_query_strip_return_ranges(self):
        """Test query_strip with return_ranges=True"""
        theta1 = np.radians(30)
        theta2 = np.radians(60)
        
        # Get ranges for RING ordering
        ranges = query_strip(self.NSIDE, theta1, theta2, return_ranges=True)
        
        # Verify it returns a 2D array
        self.assertEqual(ranges.ndim, 2)
        self.assertEqual(ranges.shape[1], 2)
        
        # Convert ranges back to individual pixels
        pixels_from_ranges = []
        for start, end in ranges:
            pixels_from_ranges.extend(range(start, end))
        pixels_from_ranges = np.array(pixels_from_ranges, dtype=np.int64)
        
        # Get pixels the normal way
        pixels_normal = query_strip(self.NSIDE, theta1, theta2, return_ranges=False)
        
        # For RING ordering, they should be in the same order
        np.testing.assert_array_equal(pixels_from_ranges, pixels_normal)

    def test_query_strip_return_ranges_nested(self):
        """Test query_strip with return_ranges=True and nest=True"""
        theta1 = np.radians(30)
        theta2 = np.radians(60)
        
        # Get ranges for NESTED ordering
        ranges = query_strip(self.NSIDE, theta1, theta2, nest=True, return_ranges=True)
        
        # Verify it returns a 2D array
        self.assertEqual(ranges.ndim, 2)
        self.assertEqual(ranges.shape[1], 2)
        
        # Convert ranges back to individual pixels
        pixels_from_ranges = []
        for start, end in ranges:
            pixels_from_ranges.extend(range(start, end))
        pixels_from_ranges = np.array(pixels_from_ranges, dtype=np.int64)
        
        # Get pixels the normal way
        pixels_normal = query_strip(self.NSIDE, theta1, theta2, nest=True, return_ranges=False)
        
        # For NESTED ordering with query_strip, order may differ but sets should be equal
        # This is because query_strip internally uses RING and converts, which may reorder
        self.assertEqual(set(pixels_from_ranges), set(pixels_normal))
        self.assertEqual(len(pixels_from_ranges), len(pixels_normal))

    def test_query_disc_return_ranges_empty(self):
        """Test query_disc with return_ranges=True when query returns no pixels"""
        # Query with a very small radius that won't contain any pixels
        vec = np.array([1.0, 0.0, 0.0])
        radius = 1e-10  # Very small radius
        
        ranges = query_disc(self.NSIDE, vec, radius, return_ranges=True)
        
        # Should return an empty array with correct shape
        self.assertEqual(ranges.shape, (0, 2))
        self.assertEqual(ranges.dtype, np.int64)

    def test_query_strip_return_ranges_empty_nested(self):
        """Test query_strip with return_ranges=True and nest=True when query returns no pixels"""
        # Query with theta values that don't include any pixel centers
        # This is a very narrow strip that should return no pixels
        theta1 = np.radians(0.0001)
        theta2 = np.radians(0.0002)
        
        ranges = query_strip(self.NSIDE, theta1, theta2, nest=True, return_ranges=True)
        
        # Should return an empty array with correct shape
        self.assertEqual(ranges.shape, (0, 2))
        self.assertEqual(ranges.dtype, np.int64)

    def test_buff_and_return_ranges_conflict(self):
        """Test that using both buff and return_ranges raises an error"""
        buff = np.empty(100, dtype=np.int64)
        vec = np.array([1.0, 0.0, 0.0])
        radius = np.radians(10)
        
        # Should raise ValueError
        with self.assertRaises(ValueError) as cm:
            query_disc(self.NSIDE, vec, radius, buff=buff, return_ranges=True)
        self.assertIn("Cannot use both", str(cm.exception))
        
        # Also test for query_polygon
        vertices = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ])
        with self.assertRaises(ValueError) as cm:
            query_polygon(self.NSIDE, vertices, buff=buff, return_ranges=True)
        self.assertIn("Cannot use both", str(cm.exception))
        
        # Also test for query_strip
        with self.assertRaises(ValueError) as cm:
            query_strip(self.NSIDE, np.radians(30), np.radians(60), buff=buff, return_ranges=True)
        self.assertIn("Cannot use both", str(cm.exception))
