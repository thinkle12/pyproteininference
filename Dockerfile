FROM python:3.7-buster

ARG VERSION

LABEL description="This container contains a an installation of PyProteinInference and is able to run the PyProteinInference Cli"
LABEL maintainer="Trent Hinkle <hinklet@gene.com>"
LABEL version="$VERSION"

WORKDIR /

# Install Py Protein Inference
# Make directory to store source
RUN mkdir /pyproteininference

# copy the entire protein_inference directory into the docker image
COPY . /pyproteininference

WORKDIR /pyproteininference/

RUN pip install -r requirements.txt

RUN python setup.py install

WORKDIR /
