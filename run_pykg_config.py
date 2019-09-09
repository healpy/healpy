# Helper to run pykg-config, whether it is installed or lives in a zipped egg.
from setuptools import Distribution
from pkg_resources import run_script

requirement = "pykg-config >= 1.2.0"
Distribution().fetch_build_eggs(requirement)
run_script(requirement, "pykg-config.py")
