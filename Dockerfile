from python:3.9-slim

MAINTAINER "Aaron Maurais -- MacCoss Lab"

RUN apt-get update && \
    apt-get -y install procps wget fontconfig libssl-dev libxml2-dev libcurl4-gnutls-dev r-base && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir -p /code/DIA_QC_report/python /code/quarto

# install rDIAUtils dependencies
RUN Rscript -e "install.packages(c('Rcpp', 'dplyr', 'tidyr', 'patchwork', 'viridis', 'BiocManager', 'rmarkdown'))" && \
    Rscript -e "BiocManager::install(c('limma', 'sva'), ask=FALSE, force=TRUE)"

# install rDIAUtils R package
COPY rDIAUtils /code/rDIAUtils
RUN cd /code/rDIAUtils && \
    Rscript -e "install.packages('.', type='source', repo=NULL)"

# install quarto
RUN cd /code/quarto && \
    wget 'https://github.com/quarto-dev/quarto-cli/releases/download/v1.3.433/quarto-1.3.433-linux-amd64.tar.gz' && \
    tar xvf 'quarto-1.3.433-linux-amd64.tar.gz' && \
    rm -v 'quarto-1.3.433-linux-amd64.tar.gz' && \
    ln -s '/code/quarto/quarto-1.3.433/bin/quarto' '/usr/local/bin' && \
    quarto install tinytex

# install python dependencies
COPY directlfq /code/directlfq
RUN pip install pandas matplotlib jupyter scikit-learn && \
    cd /code/directlfq && pip install . && \
    pip cache purge && \
    cd /code && rm -rf /code/directlfq

# add python executables
COPY python/*.py /code/DIA_QC_report/python
RUN cd /usr/local/bin && \
    echo '#!/usr/bin/env bash\npython3 /code/DIA_QC_report/python/parse_data.py "$@"' > parse_data && \
    echo '#!/usr/bin/env bash\npython3 /code/DIA_QC_report/python/generate_qc_qmd.py "$@"' > generate_qc_qmd && \
    echo '#!/usr/bin/env bash\npython3 /code/DIA_QC_report/python/normalize_db.py "$@"' > normalize_db && \
    echo '#!/usr/bin/env bash\npython3 /code/DIA_QC_report/python/generate_batch_rmd.py "$@"' > generate_batch_rmd && \
    echo '#!/usr/bin/env bash\nset -e\nexec "$@"' > entrypoint && \
    chmod 755 parse_data normalize_db entrypoint generate_qc_qmd generate_batch_rmd

WORKDIR /data

CMD []
entrypoint ["/usr/local/bin/entrypoint"]

