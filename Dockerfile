
ARG QUARTO_VERSION=1.6.42
ARG R_VERSION=4.4.3
ARG ALPINE_VERSION=3.19
ARG PYTHON_VERSION=3.11

# ----------------------------------------------------------------------------
# Stage 1: rDIAUtils and dependencies
# ----------------------------------------------------------------------------

FROM mauraisa/r-minimal:${R_VERSION} AS r_build

RUN installr -t gfortran -d dplyr tidyr ggplot2 rmarkdown svglite patchwork viridis ggiraph sva

# install rDIAUtils R package
COPY rDIAUtils /code/rDIAUtils
RUN cd /code/rDIAUtils && \
    Rscript -e "withCallingHandlers(install.packages('.', type='source', repo=NULL), \
                                    warning=function(w) stop(w))"

# ----------------------------------------------------------------------------
# Stage 2: DIA_QC_report and quarto
# ----------------------------------------------------------------------------

FROM alpine:${ALPINE_VERSION} as python_build

# Install necessary dependencies
RUN apk add --no-cache \
    bash curl git libc6-compat \
    gcc g++ python3 py3-pip python3-dev musl-dev linux-headers \
    llvm15 llvm15-dev make cmake apache-arrow apache-arrow-dev
ENV LLVM_CONFIG=/usr/lib/llvm15/bin/llvm-config

ARG QUARTO_VERSION
ENV QUARTO_VERSION=${QUARTO_VERSION}

# Install Quarto
RUN mkdir -p /code && cd /code && \
    curl -LO "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.tar.gz" && \
    tar xf "quarto-${QUARTO_VERSION}-linux-amd64.tar.gz" && \
    rm -v "quarto-${QUARTO_VERSION}-linux-amd64.tar.gz"

# install DIA_QC_report
RUN pip install --break-system-packages setuptools jupyter plotly
COPY src /code/DIA_QC_report/src
COPY pyproject.toml /code/DIA_QC_report
RUN cd /code/DIA_QC_report && \
    pip install --break-system-packages . && \
    pip cache purge && \
    cd /code && rm -rf /code/DIA_QC_report

# ----------------------------------------------------------------------------
# Stage 3: Runtime stage
# ----------------------------------------------------------------------------

FROM mauraisa/r-minimal:${R_VERSION}

LABEL maintainer="Aaron Maurais -- MacCoss Lab"

ARG QUARTO_VERSION
ARG PYTHON_VERSION

ENV QUARTO_VERSION=${QUARTO_VERSION}
ENV PYTHON_VERSION=${PYTHON_VERSION}

RUN installr -d knitr rmarkdown

# Install minimal runtime dependencies
RUN apk update && \
    apk add --no-cache npm \
    cairo pango libpng tiff \
    deno python3 pandoc libc6-compat

# RUN apk add fontconfig freetype freetype-dev libxt libxrender libxft \
#     ttf-freefont ttf-dejavu ttf-liberation ttf-droid && \
#     fc-cache -fv


# Copy R runtime dependencies from the r_build stage
COPY --from=r_build /usr/local/lib/R/library /usr/local/lib/R/library

# Copy python runtime dependencies from the python_build stage
COPY --from=python_build /code/quarto-${QUARTO_VERSION} /opt/quarto
COPY --from=python_build "/usr/lib/python${PYTHON_VERSION}/site-packages" "/usr/lib/python${PYTHON_VERSION}/site-packages"

# Create symlinks for quarto
RUN ln -s /opt/quarto/bin/quarto /usr/local/bin/quarto && \
    ln -sf /usr/bin/deno /opt/quarto/bin/tools/x86_64/deno && \
    ln -s /usr/local/bin/sass /opt/quarto/bin/tools/x86_64/sass

RUN quarto install tinytex

# clean things up

#    echo -e '#!/usr/bin/env bash\nset -e\nexec "$@"' > /usr/local/bin/entrypoint && \
#    chmod +x /usr/local/bin/entrypoint

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

# WORKDIR /data
#
# CMD []
# # ENTRYPOINT ["/usr/local/bin/entrypoint"]
