FROM python:3.7-buster

ARG VERSION

LABEL description="This container contains a an installation of PyProteinInference and is able to run the PyProteinInference Cli"
LABEL maintainer="Trent Hinkle <hinklet@gene.com>"
LABEL version="$VERSION"

# Install glpk
RUN mkdir /opt/glpk

WORKDIR /opt/glpk

RUN wget ftp://ftp.gnu.org/gnu/glpk/glpk-4.65.tar.gz \
    && tar -xzvf glpk-4.65.tar.gz
WORKDIR /opt/glpk/glpk-4.65
RUN chmod +x ./configure \
   && ./configure \
   && make uninstall \
   && make install

# Install libglpk-dev to get glpsol to work properly
RUN apt-get update \
    && apt-get install -y libglpk-dev
    
WORKDIR /

# Install Py Protein Inference
# Make directory to store source
RUN mkdir /pyproteininference

# copy the entire protein_inference directory into the docker image
COPY . /pyproteininference

WORKDIR /pyproteininference/
RUN mkdir /glpkinout

RUN pip install -r requirements.txt

RUN python setup.py install

WORKDIR /
