FROM amazonlinux:latest

LABEL maintainer="Aaron Maurais -- MacCoss Lab"

RUN dnf update && \
    dnf -y install procps git wget tar unzip pip openssl-devel libcurl-devel \
           freetype-devel libpng-devel libtiff-devel libjpeg-devel fribidi-devel \
           cairo-devel R-core-devel && \
    dnf clean all && \
    echo -e '#!/usr/bin/env bash\nset -e\nexec "$@"' > /usr/local/bin/entrypoint && \
    chmod +x /usr/local/bin/entrypoint

# install rDIAUtils dependencies
RUN Rscript -e "withCallingHandlers(install.packages(c('Rcpp', 'dplyr', 'tidyr', 'patchwork', 'viridis', 'ggiraph', 'BiocManager', 'rmarkdown', 'svglite', 'ggrastr'), \
                                                     repo='https://ftp.osuosl.org/pub/cran/'), \
                                    warning=function(w) stop(w))" && \
    Rscript -e "withCallingHandlers(BiocManager::install(c('limma', 'sva'), ask=FALSE, force=TRUE), \
                                    warning=function(w) stop(w))"

# install rDIAUtils R package
COPY rDIAUtils /code/rDIAUtils
RUN cd /code/rDIAUtils && \
    Rscript -e "withCallingHandlers(install.packages('.', type='source', repo=NULL), \
                                    warning=function(w) stop(w))"

# Install Python 3.12 and make it the default python3
RUN dnf install -y python3.12 python3.12-pip python3.12-devel && \
    alternatives --install /usr/bin/python3 python3 /usr/bin/python3.12 1 && \
    alternatives --set python3 /usr/bin/python3.12 && \
    ln -sf /usr/bin/pip3.12 /usr/local/bin/pip

# install pandoc
RUN mkdir -p /code/pandoc && cd /code/pandoc && \
    wget 'https://github.com/jgm/pandoc/releases/download/3.1.11.1/pandoc-3.1.11.1-linux-amd64.tar.gz' && \
    tar xvf 'pandoc-3.1.11.1-linux-amd64.tar.gz' && \
    rm -v 'pandoc-3.1.11.1-linux-amd64.tar.gz' && \
    ln -s '/code/pandoc/pandoc-3.1.11.1/bin/pandoc' '/usr/local/bin'

# install quarto
RUN mkdir -p /code/quarto && cd /code/quarto && \
    wget 'https://github.com/quarto-dev/quarto-cli/releases/download/v1.3.433/quarto-1.3.433-linux-amd64.tar.gz' && \
    tar xvf 'quarto-1.3.433-linux-amd64.tar.gz' && \
    rm -v 'quarto-1.3.433-linux-amd64.tar.gz' && \
    ln -s '/code/quarto/quarto-1.3.433/bin/quarto' '/usr/local/bin' && \
    quarto install tinytex

# install python dependencies
RUN pip install --ignore-installed setuptools jupyter plotly
COPY src /code/DIA_QC_report/src
COPY pyproject.toml /code/DIA_QC_report
RUN cd /code/DIA_QC_report && \
    pip install . && \
    pip cache purge && \
    cd /code && rm -rf /code/DIA_QC_report

# clean things up
RUN dnf remove -y git wget tar pip

# Git version information
ARG GIT_BRANCH
ARG GIT_REPO
ARG GIT_HASH
ARG GIT_SHORT_HASH
ARG GIT_UNCOMMITTED_CHANGES
ARG GIT_LAST_COMMIT
ARG DOCKER_TAG
ARG DOCKER_IMAGE
ARG DIA_QC_REPORT_VERSION

ENV GIT_BRANCH=${GIT_BRANCH}
ENV GIT_REPO=${GIT_REPO}
ENV GIT_HASH=${GIT_HASH}
ENV GIT_SHORT_HASH=${GIT_SHORT_HASH}
ENV GIT_UNCOMMITTED_CHANGES=${GIT_UNCOMMITTED_CHANGES}
ENV GIT_LAST_COMMIT=${GIT_LAST_COMMIT}
ENV DOCKER_IMAGE=${DOCKER_IMAGE}
ENV DOCKER_TAG=${DOCKER_TAG}
ENV DIA_QC_REPORT_VERSION=${DIA_QC_REPORT_VERSION}

WORKDIR /data

CMD []
ENTRYPOINT ["/usr/local/bin/entrypoint"]

