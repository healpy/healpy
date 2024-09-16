Installation procedure for Healpy
=================================

(NOTE: if high performance of the installed package is important, e.g. when
installing in computing centers or for performing benchmarks, please be sure
to read the `Generating native binaries`_ section below.)


Requirements
------------

Healpy depends on the HEALPix C++ and cfitsio C libraries. Source code for both
is included with Healpy and is built automatically, so you do not need to
install them yourself. Only Linux and macOS are supported, Windows only through
the "Windows Subsystem for Linux" (see below).

Binary installation with conda (recommended for Anaconda/Miniconda users)
-------------------------------------------------------------------------

Conda Forge provides a `conda channel
<https://anaconda.org/conda-forge/healpy>`_ with a pre-compiled version of
``healpy`` for linux 64bit and MAC OS X platforms, you can install it in
Anaconda with::

    conda config --add channels conda-forge
    conda install healpy

There have also been reports of specific installation issues under Mac OS
Catalina 10.15.5 with conda install as the solver appears to run without
finding the required packages. This is a general issue with a number of
packages, and not limited to ``healpy``. The most straightforward solution
(after adding conda-forge to the channel list) is for the user to decide which
packages they wish to install alongside ``healpy`` and then create a new
environment installing ``healpy`` alongside said packages. For instance if one
wishes to install ``healpy`` alongside Spyder and My_Package into newly created
environment env_healpy, the command will be::

    conda create --name env_healpy python=3.10 healpy spyder my_package

Binary installation with Pip (recommended for most other Python users)
----------------------------------------------------------------------

You can install Healpy from the Python Package Index using `pip
<http://www.pip-installer.org>`_. For most common architectures and platforms
(Linux x86-64, Linux i686, and macOS x86-64), Pip will download and install a
pre-built binary. For other platforms, it will automatically try to build
healpy from source.

Note that there are not yet native prebuilt binaries for Apple Silicon Macs.

To install the latest version of ``healpy`` with `pip
<http://www.pip-installer.org>`_, simply run::

    pip install --user healpy

If you have installed with ``pip``, you can keep your installation up to date
by upgrading from time to time::

    pip install --user --upgrade healpy

Source installation with Pip (not usually recommended)
------------------------------------------------------

On platforms for which we do not yet have prebuilt binaries in the Python
Package Index, pip build healpy from source. You can force pip to build from
source by running::

    pip install --no-binary healpy healpy

Some common issues that you might encounter when building from source:

* ``libssl-dev`` (Debian) or ``openssl-dev`` (CentOS) is required to build
  ``cfitsio`` from source.

* On Linux with newer compilers many users reported compilation errors like
  ``configure: error: cannot run C compiled programs``. The solution is to
  specifiy the flags for the C and CXX compiler:

    CC=gcc CXX=g++ CFLAGS='-fPIC' CXXFLAGS='-fPIC' pip install --user healpy

Installation from package managers
----------------------------------

Debian users may install Healpy for the Debian-supplied system Python
interpreter by running::

    sudo apt-get install python3-healpy

MacPorts users on macOS may install Healpy for the MacPorts-supplied Python
interpreter by running::

    sudo port install py39-healpy

Compilation issues with Mac OS
------------------------------

Currently most people report they cannot install `healpy` on Mac OS either via
`pip` or building from source, due to the impossibility of compiling the
`HEALPix` based extension. The best alternatives are conda, binary installation
with pip, or MacPorts.

Installation on Mac OS with MacPorts
------------------------------------

If you are using a Mac and have the `MacPorts <https://www.macports.org>`_
package manager, it's even easer to install Healpy with::

    sudo port install py39-healpy

Installation with a package manager on Debian and Ubuntu
--------------------------------------------------------

Binary `apt-get` style packages are also available in the development versions of
`Debian (sid) <https://packages.debian.org/sid/python-healpy>`_ and
`Ubuntu <https://packages.ubuntu.com/search?keywords=python-healpy>`_.

Almost-as-quick installation from official source release
---------------------------------------------------------

Healpy is also available in the
`Python Package Index (PyPI) <https://pypi.python.org/pypi/healpy>`_. You can
download it with::

    curl -O https://pypi.python.org/packages/source/h/healpy/healpy-1.14.0.tar.gz

and build it with::

    tar -xzf healpy-*.tar.gz
    cd healpy-*
    pip install .

If everything goes fine, you can test it::

    python
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import healpy as hp
    >>> hp.mollview(np.arange(12))
    >>> plt.show()

or run the test suite with::

    cd healpy-* && pytest

Building against external Healpix and cfitsio
---------------------------------------------

Healpy uses pkg-config to detect the presence of the Healpix and cfitsio
libraries. pkg-config is available on most systems. If you do not have
pkg-config installed, then Healpy will download and use (but not install) a
Python clone called pykg-config.

If you want to provide your own external builds of Healpix and cfitsio, then
download the following packages:

* `pkg-config <http://pkg-config.freedesktop.org>`_

* `HEALPix
  <http://sourceforge.net/projects/healpix/files/Healpix_3.11/autotools_packages/>`_
  autotools-style C++ package

* `cfitsio <http://heasarc.gsfc.nasa.gov/fitsio/>`_

If you are going to install the packages in a nonstandard location (say,
``--prefix=/path/to/local``), then you should set the environment variable
``PKG_CONFIG_PATH=/path/to/local/lib/pkgconfig`` when building. No other
environment variable settings are necessary, and you do not need to set
``PKG_CONFIG_PATH`` to use Healpy after you have built it.

Then, unpack each of the above packages and build them with the usual
``configure; make; make install`` recipe.

Installation on Windows through the "Windows Subsystem for Linux"
-----------------------------------------------------------------

1. Restart your computer, and follow the instructions (which appear before
   windows starts) to enter BIOS. Usually this means pressing DEL or F2 just
   after powering on. Find the option to enable virtualization (exact name will
   depend on your system, can google your machine brand name + "enable
   virtualization" for instructions)

2. Follow these instructions to install Windows Subsystem for Linux:
   https://docs.microsoft.com/en-us/windows/wsl/install-win10 Following the
   instructions for WSL version 2, and choosing Ubuntu from the store.

3. Restart machine

4. Open the newly installed Ubuntu application from the Start menu and follow
   the setup instructions.

5. When they are complete, run these commands::

		sudo apt-get update
		sudo apt-get upgrade
		sudo apt-get install python3 python3-pip

6. Quit ubuntu, restart it, and run::

		pip3 install numpy jupyter matplotlib healpy ipython jupyter

7. Quit ubuntu again, restart it, and run::

		ipython notebook --no-browser

8. Copy and paste the line starting with ``http://localhost:8888/?token=`` into
   your normal Windows web browser.

Development install
-------------------

Developers building from a snapshot of the github repository need:

* ``autoconf`` and ``libtool`` (in Debian or Ubuntu:
  ``sudo apt-get install autoconf automake libtool pkg-config``)

* ``libssl-dev`` (Debian) or ``openssl-dev`` (CentOS)
  is required to build ``cfitsio`` from source

* run ``git submodule init`` and ``git submodule update`` to get the bundled
  HEALPix sources

The best way to install healpy if you plan to develop is to do a
`development mode (a.k.a. "editable" install) <https://setuptools.pypa.io/en/latest/userguide/development_mode.html>`_
with pip by adding the ``-e`` flag::

    pip install -e .

In case of compilation errors, see the note above in the ``pip`` section.

Generating native binaries
--------------------------

Using pre-compiled wheels is typically the easiest and quickest way
to install ``healpy`` on a system. However, the performance of the installed
package may not be optimal, since the wheel has to work on all CPUs of a given
architecture (e.g. x86_64) and will therefore probably not use all features
present in your local CPU. A ``healpy`` installation which is custom-tailored
for a specific target CPU may be two or three times faster for some operations
(most notably ``alm2map*`` and ``map2alm*`` calls).

To achieve target-specific compilation, ``healpy`` must be installed from source
and the ``-march=native`` flag has to be passed to the compilers.
While details may vary slightly depending on the target platform,
the installation command will have this basic form::

    CC=gcc CXX=g++ CFLAGS="-fPIC -O3 -march=native" CXXFLAGS="-fPIC -O3 -march=native" pip3 install --user --no-binary healpy healpy

Clean
-----

When you run "python setup.py", temporary build products are placed in the
"build" directory. If you want to clean out and remove the ``build`` directory,
then run::

    python setup.py clean --all
