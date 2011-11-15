import unittest
import numpy as np

from healpy import query_disc

class TestQueryDisc(unittest.TestCase):

    def setUp(self):
        self.NSIDE = 8
        self.vec = np.array([ 0.17101007,  0.03015369,  0.98480775])
        self.radius = np.radians(6)

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
