FROM continuumio/miniconda3:4.8.2-alpine

ENV PATH /opt/conda/bin:$PATH

# install bash
USER root
RUN apk add bash

COPY scriptToInstallNFpythonLibs_PIP.sh .
RUN bash scriptToInstallNFpythonLibs_PIP.sh
