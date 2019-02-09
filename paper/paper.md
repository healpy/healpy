---
title: 'healpy: equal area pixelization and spherical harmonics transforms for data on the sphere in Python'
tags:
  - cosmology
  - astronomy
  - python
authors:
 - name: Andrea Zonca
   orcid: 0000-0001-6841-1058
   affiliation: "1"
 - name: Leo P Singer
   orcid:
   affiliation:
 - name: Daniel Lenz
   orcid: 0000-0001-5820-475X
   affiliation: "2"
 - name: Martin Reinecke
   orcid:
   affiliation: "3"
 - name: Cyrille Rosset
   orcid:
   affiliation:
 - name: Eric Hivon
   orcid: 0000-0003-1880-2733
   affiliation: "4"
 - name: Krzysztof M Gorski
   orcid:
   affiliation: "2"
affiliations:
 - name: San Diego Supercomputer Center, University of California, San Diego
   index: 1
 - name: Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California, USA
   index: 2
 - name: Max-Planck-Institute for Astrophysics, Garching, Germany
   index: 3
 - name: Institut d'Astrophysique de Paris, CNRS/Sorbonne Universite, Paris, France
   index: 4

date: 1 February 2019
bibliography: paper.bib
---

# Summary

Recent experiments measuring the temperature and polarization of the Cosmic
Microwave Background, like WMAP [@wmap13] and Planck [@planck18], produce all-sky maps at higher
and higher resolution.
Handling those datasets efficiently and studying their statistical properties
requires a discretization scheme on the sphere.
The Hierarchical Equal Area isoLatitude Pixelation ([``HEALPix``](https://healpix.sourceforge.io), [@gorski05]) scheme
has proven to be an excellent mathematical framework to store map-domain data
and efficiently compute their Spherical Harmonics Transform, whose Angular
Power Spectrum is one the most powerful tools to understand the early Universe.

See in the ![figure](healpix_grid.png) how a sphere is first split into 12 base
pixels of equal area and whose centers are at the same latitude
and then each is further subdivided to achieve higher and higher resolution.
The ``HEALPix`` team provides FORTRAN, C++, and IDL implementations of the framework.

``healpy`` is a wrapper to the multi-threaded ``HEALPix`` C++ library in Python, it implements
a user-friendly interface for all ``HEALPix`` functionalities, most importantly a fast nearest-neighbor search and the decomposition into Spherical Harmonics coefficients.
It also adds utilities to read and write maps, Spherical Harmonics coefficients and
Power Spectrum values as FITS files based on ``astropy.io.fits``.
Finally it provides extensive plotting functionality, i.e. Mollweide, Gnomonic and Cartographic
projections based on Matplotlib.
We also release a ``conda`` package on ``conda-forge`` and wheels on PyPI for Linux and MacOS. These packages
also bundle the ``HEALPix`` C++ library.

``healpy`` was designed to be used by professional cosmologists and students to analyze
outputs of Cosmic Microwave Background experiments. WMAP, Planck
and many other experiment in the field release their maps
in ``HEALPix`` FITS format that can be accessed, visualized, and analyzed with ``healpy``.
More recently, its usage spread to
other branches of Astrophysics including galaxy surveys ([Sloan Digital
Sky Survey depth maps](http://risa.stanford.edu/redmapper/)), Gamma-ray astronomy
(Fermi [@deil17]), and gravitational waves measurements (LIGO [@singer16]).

# How to cite

If a publication is based on work that used ``healpy``, please add the acknowledgement
statement: "Some of the results in this paper have been derived using the ``healpy`` and ``HEALPix`` packages"
and cite this paper and the original ``HEALPix`` paper [@gorski05].

# Acknowledgements

We acknowledge the Planck collaboration and in particular Peter Meinhold, Julian Borrill and Charles Lawrence
for supporting Andrea Zonca's work on the project. We thanks all past and current contributors, see
the [most updated list on Github](https://github.com/healpy/healpy/graphs/contributors).

# References
