Healpy's tutorial
=================

Creating and manipulating maps
------------------------------

Maps are simply numpy arrays, each array element refers to a location in the sky as defined by the Healpix pixelization schemes, see the `healpix website`_.

The resolution of the map is defined by the *NSIDE* parameter, the nside2npix function gives the number of pixel *NPIX* of the map:

>>> import numpy as np
>>> import healpy as hp
>>> NSIDE = 32
>>> m = np.arange(hp.nside2npix(NSIDE))
>>> hp.mollview(m)

.. image:: static/moll_nside32_nest.png

* Explain Npix, Nside, basics of healpix...
* Nested, Ring ordering

.. _healpix website: http://healpix.jpl.nasa.gov

Visualization
-------------

* Mollview, gnomview, cartview examples.
* Mollzoom for interactive stuff.

Reading and writing maps to file
--------------------------------

* saving and reading maps in fits files

Spherical harmonic transforms
-----------------------------

* anafast, synfast
* map2alm, alm2map, alm2cl, synalm
* smoothing, smoothalm

Masked map, partial maps
------------------------

* the :const:`UNSEEN` special value
* masked array maps

