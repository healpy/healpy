import unittest
import doctest
from .. import _query_disc
from .. import _pixelfunc

# FIXME: I can't get py.test to recognize test suites returned with
# doctest.DocTestSuite.

class TestDoctestCython(unittest.TestCase):

    def test_query_disc(self):
        failed, passed = doctest.testmod(_query_disc)
        self.assertEqual(failed, 0)

    def test_pixelfunc(self):
        failed, passed = doctest.testmod(_pixelfunc)
        self.assertEqual(failed, 0)

if __name__ == '__main__':
    unittest.main()
