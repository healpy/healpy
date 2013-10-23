cd build/lib*/healpy/
for f in wmap_band_iqumap_r9_7yr_V_v4.fits wmap_band_iqumap_r9_7yr_W_v4.fits wmap_temperature_analysis_mask_r9_7yr_v4.fits 
do
    ln -sf ../../../../../healpy/test/data/$f test/data/
done
nosetests -v --with-doctest
echo Run Cython extensions doctests
python run_doctest_cython.py
cd ../../..
