#!/bin/sh
cd build/lib*/healpy && py.test -v --doctest-modules
