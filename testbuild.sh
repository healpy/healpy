cd build/lib*/healpy/
py.test -v --doctest-modules --ignore run_doctest_cython.py
echo Run Cython extensions doctests
python run_doctest_cython.py
cd ../../..
