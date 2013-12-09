# Helper to run pykg-config, whether it is installed or lives in a zipped egg.
from pkg_resources import run_script
run_script('pykg-config >= 1.2.0', 'pykg-config.py')
