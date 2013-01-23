cd build/lib*/healpy/
ln -s ~/p/software/healpy/healpy/test/data test/data
nosetests -v
cd ../../..
