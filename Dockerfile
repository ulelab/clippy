FROM continuumio/miniconda3:4.8.2

# Copy Clippy files
COPY . /clippy

# Install apt packages
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    procps=2:3.3.15-2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create Clippy conda env
RUN conda env create -f /clippy/environment.yml \
    && conda clean -a
ENV PATH /opt/conda/envs/clippy/bin:$PATH

WORKDIR /clippy
ENV PATH /clippy:$PATH
RUN chmod +x clip.py
