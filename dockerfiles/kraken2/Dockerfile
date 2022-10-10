# base image
FROM ubuntu:18.04

# Metadata
LABEL base.image="ubuntu:18.04"
LABEL container.version="1"
LABEL software="Kraken2"
LABEL software.version="2.0.8-beta"
LABEL description="Taxonomic sequence classifier"
LABEL website="https://github.com/DerrickWood/kraken2"
LABEL license="https://github.com/DerrickWood/kraken2/blob/master/LICENSE"
LABEL maintainer="TW Battaglia"
LABEL maintainer.email="t.battaglia@nki.nl"

# Install dependencies
RUN apt-get update && apt-get -y install \
  wget \
  zlib1g-dev \
  make \
  python \
  perl \
  g++ \
  ca-certificates \
  nfs-common \
  python-biopython

# Download kraken2
RUN wget https://github.com/DerrickWood/kraken2/archive/v2.0.8-beta.tar.gz && \
  tar -xzf v2.0.8-beta.tar.gz && \
  rm -rf v2.0.8-beta.tar.gz && \
  cd kraken2-2.0.8-beta && \
  ./install_kraken2.sh .

ENV PATH="$PATH:/kraken2-2.0.8-beta" \
    LC_ALL=C

# Download bracken and install
RUN wget https://github.com/jenniferlu717/Bracken/archive/v2.5.3.tar.gz && \
  tar -xzf v2.5.3.tar.gz && \
  rm -rf v2.5.3.tar.gz && \
  mv Bracken-2.5.3 bracken && \
  cd /bracken && \
  sh ./install_bracken.sh . && \
  ln -s /bracken/bracken /usr/local/bin/bracken && \
  ln -s /bracken/src/est_abundance.py /usr/local/bin/est_abundance.py && \
  ln -s /bracken/src/kreport2mpa.py /usr/local/bin/kreport2mpa.py && \
  ln -s /bracken/bracken-build /usr/local/bin/bracken-build
	
