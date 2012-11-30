.. healpy.pixelfunc:

.. currentmodule:: healpy.pixelfunc


:mod:`pixelfunc` -- Pixelisation related functions
==================================================

conversion from/to sky coordinates
----------------------------------
.. autosummary::
   :toctree: generated/

   pix2ang
   pix2vec
   ang2pix
   vec2pix
   vec2ang
   ang2vec
   get_neighbours
   get_all_neighbours

conversion between NESTED and RING schemes
------------------------------------------
.. autosummary::
   :toctree: generated/

   nest2ring
   ring2nest
   reorder

nside/npix/resolution
---------------------
.. autosummary::
   :toctree: generated/

   nside2npix
   npix2nside
   nside2resol
   nside2pixarea
   max_pixrad
   isnsideok
   isnpixok
   get_map_size
   get_min_valid_nside
   get_nside
   maptype
   ud_grade

Masking pixels
--------------
.. autosummary::
   :toctree: generated/

   UNSEEN
   mask_bad
   mask_good
   ma

Map data manipulation
---------------------
.. autosummary::
   :toctree: generated/

   fit_dipole
   fit_monopole
   remove_dipole
   remove_monopole
   get_interp_val

