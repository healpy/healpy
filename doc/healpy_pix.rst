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
   get_all_neighbours
   get_interp_weights
   get_interp_val

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
   nside2order
   order2nside
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

Notes on downgrade quality
--------------------------

For simple resolution changes, :func:`ud_grade` remains available and fast. It
is usually the right choice for masks, hit-count maps, and compact
pixel-localized features. For diffuse or band-limited downgrade workflows,
:func:`healpy.sphtfunc.harmonic_ud_grade` can preserve the map more faithfully
because it first band-limits the map in spherical-harmonic space before
synthesizing it at the target ``nside``.

See :func:`healpy.sphtfunc.harmonic_ud_grade` for the API and the following
tutorial notebooks for worked comparisons:

* :doc:`healpy_harmonic_ud_grade_comparison`

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
