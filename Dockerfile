from python:3.9-slim

MAINTAINER "Aaron Maurais -- MacCoss Lab"

RUN apt-get update && \
    apt-get -y install procps wget fontconfig && \
    mkdir -p /code/DIA_QC_report/python /code/quarto

# install quarto
RUN cd /code/quarto && \
    wget 'https://github.com/quarto-dev/quarto-cli/releases/download/v1.3.433/quarto-1.3.433-linux-amd64.tar.gz' && \
    tar xvf 'quarto-1.3.433-linux-amd64.tar.gz' && \
    rm -v 'quarto-1.3.433-linux-amd64.tar.gz' && \
    ln -s '/code/quarto/quarto-1.3.433/bin/quarto' '/usr/local/bin' && \
    quarto install tinytex

# install python dependencies
RUN pip install pandas matplotlib jupyter scikit-learn

COPY python/*.py /code/DIA_QC_report/python
COPY resources /code/DIA_QC_report/resources

# add executables
RUN cd /usr/local/bin && \
    echo '#!/usr/bin/env bash\npython3 /code/DIA_QC_report/python/parse_data.py "$@"' > parse_data && \
    echo '#!/usr/bin/env bash\npython3 /code/DIA_QC_report/python/generate_qc_qmd.py "$@"' > generate_qc_qmd && \
    echo '#!/usr/bin/env bash\npython3 /code/DIA_QC_report/python/normalize_db.py "$@"' > normalize_db && \
    echo '#!/usr/bin/env bash\nset -e\nexec "$@"' > entrypoint && \
    chmod 755 make_qmd parse_data entrypoint

WORKDIR /data

CMD []
entrypoint ["/usr/local/bin/entrypoint"]

