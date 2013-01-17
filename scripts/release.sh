#!/usr/bin/env bash
# healpy release script
# creates a release branch
version=$1

echo create release branch
git checkout -b release-$version

echo change version number
echo __version__=\'$version\' > healpy/version.py

echo copy current version of healpix to other folder to be included in the release package
mkdir -p healpixcxx/src
cp -r healpixsubmodule/src/cxx healpixcxx/src/
echo replace folder reference in setup.py
sed -i '' -e's/healpixsubmodule/healpixcxx/' setup.py

echo add cython compiled cxx files
git add -f healpy/src/_sphtools.cpp
git add -f healpy/src/_pixelfunc.cpp
git add -f healpy/src/_query_disc.cpp

echo now you can check changes and commit
