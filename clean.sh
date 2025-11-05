#!/bin/env bash
git clean -dxf -e .venv/
pushd cextern/healpix; git clean -dxf; popd
pushd cextern/cfitsio; git clean -dxf; popd
