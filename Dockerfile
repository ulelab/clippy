FROM continuumio/miniconda3:4.10.3

# Install apt packages
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    procps=2:3.3.15-2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy Clippy files
COPY . /clippy

# Create Clippy conda env
RUN conda env create -f /clippy/environment.yml \
    && conda clean -a
ENV PATH /opt/conda/envs/clippy/bin:$PATH

# Install Clippy
WORKDIR /clippy
RUN python -m pip install . -vv --ignore-installed --no-deps
