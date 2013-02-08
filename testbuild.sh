cd build/lib*/healpy/
ln -sf ../../../../healpy/test/data test
nosetests -v
cd ../../..
