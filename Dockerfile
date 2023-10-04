# This Dockerfile has been modified from: https://github.com/StaPH-B/docker-builds/blob/master/samtools/1.18/Dockerfile
# We acknowledge and thank the authors for sharing their code:
# Shelby Bennett, Erin Young, Curtis Kapsak, & Kutluhan Incekara

ARG SAMTOOLS_VER="1.18"

FROM ubuntu:jammy as builder

ARG SAMTOOLS_VER

# install dependencies required for compiling samtools
RUN apt-get update && apt-get install --no-install-recommends -y \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    zlib1g-dev \
    libssl-dev \
    gcc \
    wget \
    make \
    perl \
    bzip2 \
    gnuplot \
    ca-certificates

# download, compile, and install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VER}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VER} && \
    ./configure && \
    make && \
    make install

### start of app stage ###
FROM ubuntu:jammy as app

ARG SAMTOOLS_VER

LABEL base.image="ubuntu:jammy"
LABEL dockerfile.version="1"
LABEL software="tbp-parser"
LABEL software.version="1.0.3"
LABEL description="tbp-parser and samtools"
LABEL website="https://github.com/theiagen/tbp-parser"
LABEL license="https://github.com/theiagen/tbp-parser/blob/main/LICENSE"
LABEL maintainer="Curtis Kapsak"
LABEL maintainer.email="curtis.kapsak@theiagen.com"
LABEL maintainer2="Sage Wright"
LABEL maintainer2.email="sage.wright@theiagen.com"

ARG DEBIAN_FRONTEND=noninteractive

# install dependencies required for running samtools
# 'gnuplot' required for plot-ampliconstats
RUN apt-get update && apt-get install --no-install-recommends -y \
    python3 \
    python3-pip \
    perl \
    zlib1g \
    libncurses5 \
    bzip2 \
    liblzma-dev \
    libcurl4-gnutls-dev \
    gnuplot \
    && apt-get autoclean && rm -rf /var/lib/apt/lists/*

# copy in samtools executables from builder stage
COPY --from=builder /usr/local/bin/* /usr/local/bin/

ENV LC_ALL=C

# copy in all of the tb-profiler-parsing code
COPY . /tbp-parser

ENV DEB_PYTHON_INSTALL_LAYOUT=deb_system

# install python dependencies and parsing scripts
# updating setuptools because the version in apt for ubuntu:jammy is a bit old
# just using a requrements.txt file for now; pyproject.toml is the current recommended approach, but I'm not too familiar with it
RUN cd /tbp-parser && \
python3 -m pip install 'setuptools==68.2.0' && \
python3 -m pip install -r requirements.txt

# final working directory is /data
WORKDIR /data

# default command is to pull up help options
CMD [ "python3", "/tbp-parser/tbp_parser/tbp_parser.py", "--help"]

### start of test stage ###
FROM app as test

# # install wget for downloading test files
RUN apt-get update && apt-get install --no-install-recommends -y wget ca-certificates

RUN wget -q https://raw.githubusercontent.com/StaPH-B/docker-builds/master/tests/SARS-CoV-2/SRR13957123.consensus.fa && \
  wget -q https://raw.githubusercontent.com/StaPH-B/docker-builds/master/tests/SARS-CoV-2/SRR13957123.primertrim.sorted.bam && \
  wget -q https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V3/nCoV-2019.primer.bed && \
  samtools stats SRR13957123.primertrim.sorted.bam && \
  samtools faidx SRR13957123.consensus.fa && \
  samtools ampliconstats nCoV-2019.primer.bed SRR13957123.primertrim.sorted.bam > SRR13957123_ampliconstats.txt && \
  plot-ampliconstats plot SRR13957123_ampliconstats.txt && \
  ls

# test version and help option outputs
# run pytest suite
# expect to see a warning regarding pandas data frame concatenation
RUN python3 /tbp-parser/tbp_parser/tbp_parser.py --version && \
python3 /tbp-parser/tbp_parser/tbp_parser.py --help && \
python3 -m pip install pytest && \
cd /tbp-parser && \
pytest