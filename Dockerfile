# Bioinformatics runtime image for genomics-workflow-agent.
#
# Provides: FastQC, MultiQC, fastp, cutadapt, samtools, bcftools,
#           bedtools, mosdepth, Nextflow, Python 3.11, the package itself.
#
# Claude Code runs on the HOST. This image is the reproducible tool runtime.
# Do not put Claude credentials or API keys in this file.

FROM mambaorg/micromamba:1.5.8

# Run as root so we can install system packages and chown /workspace
USER root

WORKDIR /workspace

# System packages (needed for OpenJDK/Nextflow/FastQC)
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        procps \
        git \
        curl \
        bash \
        ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Conda environment via micromamba
# env spec first so this layer is cached when environment.yml hasn't changed
COPY --chown=$MAMBA_USER:$MAMBA_USER env/environment.yml /tmp/environment.yml

# Create the base environment in the micromamba default prefix
RUN micromamba install -y -n base \
        -c conda-forge \
        -c bioconda \
        -c defaults \
        python=3.11 \
        pip \
        "setuptools>=68" \
        wheel \
        "openjdk=17" \
        fastqc \
        multiqc \
        fastp \
        cutadapt \
        samtools \
        bcftools \
        bedtools \
        mosdepth \
        nextflow \
        pytest \
        pytest-cov \
        numpy \
        pandas \
        scipy \
        matplotlib \
        seaborn \
        h5py \
        pyyaml \
        jsonschema \
        rich && \
    micromamba clean -afy

# Python package
COPY pyproject.toml README.md ./
COPY genomics_workflow_agent/ genomics_workflow_agent/

RUN micromamba run -n base pip install -e ".[full]"

# Rest of the repo (tests, examples, scripts, fixtures)
COPY . .

# Shell defaults
SHELL ["/bin/bash", "-c"]
ENV PATH="/opt/conda/bin:${PATH}"
ENV PYTHONUNBUFFERED=1
ENV PYTHONIOENCODING=utf-8

# Make sure micromamba base env is active for every shell
RUN echo "source /opt/conda/etc/profile.d/conda.sh && conda activate base" >> /etc/bash.bashrc

CMD ["bash"]
