Healpy's tutorial
=================

Creating and manipulating maps
------------------------------

Maps are simply numpy arrays, each array element refers to a location in the sky as defined by the Healpix pixelization schemes.

.. image:: static/moll_nside32_nest.png

* Explain Npix, Nside, basics of healpix...
* Nested, Ring ordering

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

