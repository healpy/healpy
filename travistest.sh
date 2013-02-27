#!/bin/env bash
cd build/lib*/healpy/test/data
. get_wmap_maps.sh
cd ../..
nosetests -v
exit $?
