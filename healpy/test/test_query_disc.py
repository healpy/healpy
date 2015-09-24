import unittest
import numpy as np

from .. import query_disc, boundaries, nside2npix

try:
    from exceptions import ValueError
except:
    pass

class TestQueryDisc(unittest.TestCase):

    def setUp(self):
        self.NSIDE = 8
        self.vec = np.array([ 0.17101007,  0.03015369,  0.98480775])
        self.radius = np.radians(6)
        self.nside2_55_corners_precomp = np.array(
            [[[  2.44708573e-17,  5.27046277e-01,  3.60797400e-01,  4.56383842e-17],
              [  3.99652627e-01,  5.27046277e-01,  8.71041977e-01,  7.45355992e-01],
              [  9.16666667e-01,  6.66666667e-01,  3.33333333e-01,  6.66666667e-01]],
             [[  2.44708573e-17,  5.27046277e-01,  3.60797400e-01,  4.56383842e-17],
              [  3.99652627e-01,  5.27046277e-01,  8.71041977e-01,  7.45355992e-01],
              [  9.16666667e-01,  6.66666667e-01,  3.33333333e-01,  6.66666667e-01]]]
        )

    def test_not_inclusive(self):
        #HIDL> query_disc, 8, [ 0.17101007,  0.03015369,  0.98480775],6,listpix,/DEG,NESTED=0           
        #HIDL> print,listpix
        #           4
        np.testing.assert_array_equal(
                query_disc(self.NSIDE, self.vec, self.radius, inclusive=False),
                np.array([4])
            )

    def test_inclusive(self):
        #HIDL> query_disc, 8, [ 0.17101007,  0.03015369,  0.98480775],6,listpix,/DEG,NESTED=0,/inclusive
        #HIDL> print,listpix
        #           0           3           4           5          11          12          13          23
        np.testing.assert_array_equal(
                query_disc(self.NSIDE, self.vec, self.radius, inclusive=True),
                np.array([ 0, 3, 4, 5, 11, 12, 13, 23 ])
            )

    def test_boundaries(self):
        nside = 2
        corners = boundaries(nside, 5)
        corners_precomp = np.array(
            [[  2.44708573e-17,  5.27046277e-01,  3.60797400e-01,  4.56383842e-17],
             [  3.99652627e-01,  5.27046277e-01,  8.71041977e-01,  7.45355992e-01],
             [  9.16666667e-01,  6.66666667e-01,  3.33333333e-01,  6.66666667e-01]])
        np.testing.assert_array_almost_equal(corners, corners_precomp, decimal=8)

    def test_boundaries_list(self):
        nside = 2
        corners = boundaries(nside, [5,5])
        np.testing.assert_array_almost_equal(corners, self.nside2_55_corners_precomp, decimal=8)

    def test_boundaries_phi_theta(self):
        nside = 2
        corners = boundaries(nside, np.array([5,5]))
        np.testing.assert_array_almost_equal(corners, self.nside2_55_corners_precomp, decimal=8)

    def test_boundaries_floatpix_array(self):
        self.assertRaises(ValueError, boundaries, 2, np.array([5.,5]))

    def test_boundaries_floatpix_scalar(self):
        self.assertRaises(ValueError, boundaries, 2, 1/2.)

    def test_buffer_mode(self) :
    
        
        # allocate something manifestly too short, should raise a value error
        buff = np.empty(0, dtype=np.int64)
        self.assertRaises(ValueError,
                          query_disc,
                          self.NSIDE, self.vec, self.radius, inclusive=True, buff=buff)
        

        # allocate something of wrong type, should raise a value error
        buff = np.empty(nside2npix(self.NSIDE), dtype=np.float64)
        self.assertRaises(ValueError,
                          query_disc,
                          self.NSIDE, self.vec, self.radius, inclusive=True, buff=buff)
       
        # allocate something acceptable, should succeed and return a subview
        buff = np.empty(nside2npix(self.NSIDE), dtype=np.int64)
        result = query_disc(self.NSIDE, self.vec, self.radius, inclusive=True, buff=buff)
        
        assert result.base is buff
        
        np.testing.assert_array_equal(
            result,
            np.array([ 0, 3, 4, 5, 11, 12, 13, 23 ])
            )
