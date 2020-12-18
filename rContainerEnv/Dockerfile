FROM r-base:4.0.2

USER root
RUN apt-get update && apt-get install -y procps

COPY installLibs.R .

RUN Rscript installLibs.R