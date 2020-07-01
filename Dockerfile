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
RUN mkdir /protein_inference

# copy the entire protein_inference directory into the docker image
COPY /protein_inference /protein_inference/protein_inference/
COPY /tests /protein_inference/tests/
COPY setup.cfg /protein_inference/
COPY setup.py /protein_inference/
COPY /parameters /protein_inference/parameters/
COPY /scripts /protein_inference/scripts/
COPY .git /protein_inference/.git/
COPY README.md /protein_inference/
COPY requirements.txt /protein_inference/

WORKDIR /protein_inference/
RUN mkdir /glpkinout

RUN pip install -r requirements.txt \
    && python setup.py install
