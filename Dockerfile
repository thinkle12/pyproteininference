FROM python:3.7-buster

ARG VERSION

LABEL description="This container contains a an installation of PyProteinInference and is able to run the PyProteinInference Cli"
LABEL maintainer="Trent Hinkle <hinklet@gene.com>"
LABEL version="$VERSION"

# Install glpk
RUN mkdir /build/

WORKDIR /build/

RUN wget ftp://ftp.gnu.org/gnu/glpk/glpk-4.65.tar.gz \
    && tar -xzvf glpk-4.65.tar.gz
WORKDIR /build/glpk-4.65
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
RUN mkdir /py_protein_inference

# copy the entire protein_inference directory into the docker image
COPY /py_protein_inference /py_protein_inference/py_protein_inference/
COPY /tests /py_protein_inference/tests/
COPY setup.cfg /py_protein_inference/
COPY setup.py /py_protein_inference/
COPY /parameters /py_protein_inference/parameters/
COPY /scripts /py_protein_inference/scripts/
COPY .git /py_protein_inference/.git/
COPY README.md /py_protein_inference/
COPY requirements.txt /py_protein_inference/

WORKDIR /py_protein_inference/
RUN mkdir /glpkinout

RUN pip install -r requirements.txt

RUN python setup.py develop
