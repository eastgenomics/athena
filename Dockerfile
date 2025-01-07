FROM python:3.8-slim

LABEL title="Athena" description=" Provides HTML reports with detailed coverage analysis for clinical NGS data"

COPY . /athena

RUN \
    apt-get -y update; apt-get -y install curl && \
    echo "Installing Python requirements" && \
    pip install --quiet --upgrade pip && \
    pip install --only-binary polars -r /athena/requirements.txt && \
    echo "Installing bedtools" && \
    curl https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static -o /usr/local/bin/bedtools && \
    chmod a+x /usr/local/bin/bedtools && \
    echo "Deleting cache files and removing build dependencies" 1>&2 && \
    find /usr/local/lib/python3.8  \( -iname '*.c' -o -iname '*.pxd' -o -iname '*.pyd' -o -iname '__pycache__' \) | \
    xargs rm -rf {} && \
    rm -rf /root/.cache/pip && \
    # apk --purge del gcc musl-dev linux-headers && \
    echo "Setting Athena alias" 1>&2 && \
    printf '#!/bin/sh\npython3 /app/athena/athena.py "$@"'  > /usr/local/bin/athena && \
    chmod +x /usr/local/bin/athena

WORKDIR /app/athena

# display help if no args specified
CMD athena --help
