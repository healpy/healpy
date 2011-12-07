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
>>> hp.mollview(m, title="Mollview image RING")

.. image:: static/moll_nside32_ring.png

Healpix supports two different ordering schemes, *RING* or *NESTED*, **by default healpy maps are in *RING* ordering**.

In order to work with *NESTED* ordering, all map related functions support the *nest* keyword, for example:

>>> hp.mollview(m, nest=True, title="Mollview image NESTED")

.. image:: static/moll_nside32_nest.png

.. _healpix website: http://healpix.jpl.nasa.gov

Reading and writing maps to file
--------------------------------

* saving and reading maps in fits files

Visualization
-------------

* Mollview, gnomview, cartview examples.
* Mollzoom for interactive stuff.


Spherical harmonic transforms
-----------------------------

* anafast, synfast
* map2alm, alm2map, alm2cl, synalm
* smoothing, smoothalm

Masked map, partial maps
------------------------

* the :const:`UNSEEN` special value
* masked array maps

