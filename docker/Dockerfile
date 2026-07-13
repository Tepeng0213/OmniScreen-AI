FROM mambaorg/micromamba:1.5-jammy

LABEL org.opencontainers.image.title="OmniScreen-AI"
LABEL org.opencontainers.image.description="PD-L1 multi-modal drug screening environment"

USER root
WORKDIR /workspace

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl git \
    autodock-vina \
    && rm -rf /var/lib/apt/lists/*

# Conda environment
COPY requirements.txt /tmp/requirements.txt
RUN micromamba install -y -n base -c conda-forge \
    python=3.11 \
    rdkit openmm mdtraj \
    pandas numpy scikit-learn \
    biopython matplotlib seaborn \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:$PATH"

COPY . /workspace
WORKDIR /workspace

CMD ["bash"]
