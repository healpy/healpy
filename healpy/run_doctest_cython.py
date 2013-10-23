import doctest
import sys

import _query_disc
query_disc_failure_count, query_dist_test_count = doctest.testmod(_query_disc, verbose=True)

import _pixelfunc
pixelfunc_failure_count, pixelfunc_test_count = doctest.testmod(_pixelfunc, verbose=True)
sys.exit(query_disc_failure_count + pixelfunc_failure_count)
