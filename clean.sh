#!/bin/env bash
git clean -dxf
pushd cextern/healpix; git clean -dxf; popd
pushd cextern/cfitsio; git clean -dxf; popd
