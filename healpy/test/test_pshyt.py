#import healpy.pshyt as pshyt
import numpy as np
import unittest
# pshyt removed from healpy in Dec 2012
#
#class TestPshytYlmgen(unittest.TestCase):
#    
#    def test_ylm(self):
#        ylmgen = pshyt.Ylmgen(10, smax=2, spinrec=True)
#        ylm = ylmgen.ylm(0, np.pi / 4)
#        expected = np.array([ 0.28209479,  0.34549415,  0.15769578, -0.13193776,
#                             -0.34380303, -0.35145956, -0.15097686,  0.13881729, 
#                              0.34700105,  0.35110129, 0.14880806])
#        np.testing.assert_array_almost_equal(ylm, expected)
#        ylm2 = ylmgen.ylm(4, np.pi / 4)
#        expected2 = np.array([ 0.28209479,  0.34549415,  0.15769578, -0.13193776,
#                               0.11063317, 0.25945779,  0.40137892,  0.45405113,
#                               0.35943479,  0.1274281, -0.15503676])
#        np.testing.assert_array_almost_equal(ylm2, expected2)
#
#    def test_lambda_wx(self):
#        ylmgen = pshyt.Ylmgen(10, smax=2, spinrec=True)
#        wx = ylmgen.lambda_wx(2, 0, np.pi / 4)
#        expected = np.array([[  0.        ,   0.        ],
#                             [  0.        ,   0.        ],
#                             [  0.9461747 ,   0.        ],
#                             [  3.95813273,   0.        ],
#                             [  7.93391602,   0.        ],
#                             [  8.68311844,  -0.        ],
#                             [  1.66869156,  -0.        ],
#                             [-12.9281147 ,  -0.        ],
#                             [-27.1921607 ,  -0.        ],
#                             [-28.58057919,   0.        ],
#                             [ -8.90481337,   0.        ]])
#        np.testing.assert_array_almost_equal(wx, expected)
#        wx2 = ylmgen.lambda_wx(2, 4, np.pi / 4)
#        expected2 = np.array([[  0.        ,   0.        ],
#                              [  0.        ,   0.        ],
#                              [  0.9461747 ,   0.        ],
#                              [  3.95813273,   0.        ],
#                              [  3.98279423,  -3.75501441],
#                              [  7.2648181 ,  -5.87085959],
#                              [  7.58160184,  -2.5228245 ],
#                              [  4.28105352,   8.80628939],
#                              [  0.84572892,  25.11686755],
#                              [  3.05827443,  37.48379418],
#                              [ 13.70524937,  35.43162176]])
#        np.testing.assert_array_almost_equal(wx2, expected2)


if __name__ == '__main__':
    unittest.main()
