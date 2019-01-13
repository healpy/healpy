# To use this image, first build the dockerfile by switching
# to the directory in which it is located and then running
# >> docker build -t hp .
# This builds the image, and
# >> docker images
# will show you the healpy (hp) image, as well as the
# underlying python-3.6-slim image.
#
# After that, run the image with
# >> docker run -d -t --name healpy -v "$(pwd)":/home/healpy hp:latest
# This will run the image (i.e. create a container, which is an instance of 
# that image) and mount the current directory, so that you can directly modify
# the content of the container, without having to access it.
# 
# You can then use
# >> docker exec -it healpy /bin/bash
# to start a bash session in the container, and interact with the healpy package.
# This is very helpful to install the package, run tests, etc., while the files can be edited
# and added to source control outside of the container due to the mount.

FROM python:3.7-slim

WORKDIR /home

# Setup tools to work in the container
RUN apt-get update && apt-get install -y \
    make \
    vim \
    gcc \
    g++ \
    git \
    autoconf \
    automake \
    libtool \
    pkg-config

# Copy the local healpy folder into the container
ADD . /home/healpy

# Alternatively, comment the line above and uncomment this line to clone the 
# latest healpy version from github
# RUN git clone https://github.com/healpy/healpy.git /home/healpy

# Initialize the submodules
RUN cd /home/healpy && git submodule init && git submodule update

# Setup the python distribution
RUN pip install -r /home/healpy/requirements.txt
RUN pip install scipy pytest # Needed for some of the tests

ENV HOME=/home
