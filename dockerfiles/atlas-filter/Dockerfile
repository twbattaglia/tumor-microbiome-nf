FROM nfcore/base:latest

LABEL base.image="ubuntu:18.04"
LABEL maintainer="TW Battaglia"
LABEL maintainer.email="t.battaglia@nki.nl"

# Set the locale
RUN apt-get update
RUN apt-get install -y locales locales-all apt-utils
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US.UTF-8


# Install libarchive
RUN apt-get install -y libarchive-dev

COPY environment.yml /
RUN conda install mamba -c conda-forge
RUN mamba update -n base -y conda
RUN mamba env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/env/bin:$PATH

