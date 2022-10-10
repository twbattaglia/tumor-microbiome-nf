FROM nfcore/base:latest

LABEL base.image="ubuntu:18.04"
LABEL maintainer="TW Battaglia"
LABEL maintainer.email="t.battaglia@nki.nl"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/gatk/bin:$PATH

# Set the locale
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8


