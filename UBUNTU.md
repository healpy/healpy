Ubuntu packaging for healpy
===========================

PPA
---

<https://launchpad.net/~zonca/+archive/healpix>

How to build
------------

    python setup.py --command-packages=stdeb.command bdist_deb

the problem is the this file is not in the right format for the PPA,
therefore we need to unpack it:

    mkdir tmp
    cd tmp
    dpkg-source -x ../deb_dist/*.dsc

fix the `debian/changelog` file with right target os and ownership.

build and sign with `debuild`:

    debuild -S -sa
    dput ppa:zonca/healpix ../*.changes
