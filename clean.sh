#!/bin/env bash
git clean -dxf
pushd healpixsubmodule; git clean -dxf; popd
pushd cfitsio; git clean -dxf; popd
