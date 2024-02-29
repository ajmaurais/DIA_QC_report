from amazonlinux:latest

MAINTAINER "Aaron Maurais -- MacCoss Lab"

RUN dnf update && \
    dnf -y install wget tar unzip pip openssl-devel libcurl-devel R-core-devel && \
    dnf clean all

# install rDIAUtils dependencies
RUN Rscript -e "withCallingHandlers(install.packages(c('Rcpp', 'dplyr', 'tidyr', 'patchwork', 'viridis', 'BiocManager', 'rmarkdown'), \
                                                     repo='https://ftp.osuosl.org/pub/cran/'), \
                                    warning=function(w) stop(w))" && \
    Rscript -e "withCallingHandlers(BiocManager::install(c('limma', 'sva'), ask=FALSE, force=TRUE), \
                                    warning=function(w) stop(w))"

# install rDIAUtils R package
COPY rDIAUtils /code/rDIAUtils
RUN cd /code/rDIAUtils && \
    Rscript -e "withCallingHandlers(install.packages('.', type='source', repo=NULL), \
                                    warning=function(w) stop(w))"
 
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
COPY directlfq /code/directlfq
COPY pyDIAUtils /code/pyDIAUtils
RUN pip install setuptools jsonschema pyarrow pandas matplotlib jupyter && \
    cd /code/directlfq && pip install . && \
    cd /code/pyDIAUtils && pip install . && \
    pip cache purge && \
    cd /code && rm -rf /code/directlfq /code/pyDIAUtils

# add python executables
COPY python/normalize_db.py /code/DIA_QC_report/python/normalize_db.py
COPY python/parse_data.py /code/DIA_QC_report/python/parse_data.py
COPY python/generate_qc_qmd.py /code/DIA_QC_report/python/generate_qc_qmd.py
COPY python/generate_batch_rmd.py /code/DIA_QC_report/python/generate_batch_rmd.py
COPY python/make_gene_matrix.py /code/DIA_QC_report/python/make_gene_matrix.py
RUN cd /usr/local/bin && \
    echo -e '#!/usr/bin/env bash\npython3 /code/DIA_QC_report/python/parse_data.py "$@"' > parse_data && \
    echo -e '#!/usr/bin/env bash\npython3 /code/DIA_QC_report/python/generate_qc_qmd.py "$@"' > generate_qc_qmd && \
    echo -e '#!/usr/bin/env bash\npython3 /code/DIA_QC_report/python/normalize_db.py "$@"' > normalize_db && \
    echo -e '#!/usr/bin/env bash\npython3 /code/DIA_QC_report/python/generate_batch_rmd.py "$@"' > generate_batch_rmd && \
    echo -e '#!/usr/bin/env bash\npython3 /code/DIA_QC_report/python/make_gene_matrix.py "$@"' > make_gene_matrix && \
    echo -e '#!/usr/bin/env bash\nset -e\nexec "$@"' > entrypoint && \
    chmod 755 parse_data normalize_db entrypoint generate_qc_qmd generate_batch_rmd make_gene_matrix

# Git version information
ARG GIT_BRANCH
ARG GIT_REPO
ARG GIT_HASH
ARG GIT_SHORT_HASH
ARG GIT_UNCOMMITTED_CHANGES
ARG GIT_LAST_COMMIT
ARG DOCKER_TAG
ARG DOCKER_IMAGE

ENV GIT_BRANCH=${GIT_BRANCH}
ENV GIT_REPO=${GIT_REPO}
ENV GIT_HASH=${GIT_HASH}
ENV GIT_SHORT_HASH=${GIT_SHORT_HASH}
ENV GIT_UNCOMMITTED_CHANGES=${GIT_UNCOMMITTED_CHANGES}
ENV GIT_LAST_COMMIT=${GIT_LAST_COMMIT}
ENV DOCKER_IMAGE=${DOCKER_IMAGE}
ENV DOCKER_TAG=${DOCKER_TAG}

WORKDIR /data

CMD []
entrypoint ["/usr/local/bin/entrypoint"]

