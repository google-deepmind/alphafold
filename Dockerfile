FROM nvidia/cuda:11.0-base

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -y update
RUN apt-get -y install \
	cmake \
	git \
	python3.7 \
	python3-pip \
	wget \
	&& rm -rf /var/lib/apt/lists/*

RUN mkdir /app
WORKDIR /app
COPY ./requirements.txt ./requirements.txt
RUN pip3 install \
	-r requirements.txt \
	-f https://storage.googleapis.com/jax-releases/jax_releases.html

# Compile HHsuite from source
RUN git clone --branch v3.3.0 https://github.com/soedinglab/hh-suite.git /app/hh-suite
RUN mkdir /app/hh-suite/build
WORKDIR /app/hh-suite/build
RUN cmake -DCMAKE_INSTALL_PREFIX=/opt/hhsuite .. \
  && make -j 4 \
  && make install \
  && ln -s /opt/hhsuite/bin/* /usr/bin \
  && rm -rf /app/hh-suite

# Install Miniconda package manger
RUN wget -q -P /tmp https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
	&& bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
	&& rm /tmp/Miniconda3-latest-Linux-x86_64.sh

# Install conda packages
WORKDIR /app
ENV PATH="/opt/conda/bin:$PATH"
RUN conda update -qy conda
RUN conda install -y -c conda-forge \
	cudatoolkit=11.0.221 \
  openmm=7.5.1 \
  pdbfixer=1.7 \
  pip=21.1.3

COPY . ./app
