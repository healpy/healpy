#!/bin/env bash
cd build/lib*/healpy
py.test -v --doctest-modules --ignore run_doctest_cython.py
nosetests_returnvalue=$?
echo Run Cython extensions doctests
cd ..
python healpy/run_doctest_cython.py
cython_doctest_returnvalue=$?
exit $(($nosetests_returnvalue + $cython_doctest_returnvalue))
