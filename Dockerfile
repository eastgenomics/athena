FROM python:3.12-alpine

LABEL title="Athena" description=" Provides HTML reports with detailed coverage analysis for clinical NGS data"

COPY . /athena

RUN \
echo "Installing Python requirements" && \
    pip install --quiet --upgrade pip && \
    pip install -r requirements.txt && \
    echo "Deleting cache files and removing build dependencies" 1>&2 && \
    find /usr/local/lib/python3.12  \( -iname '*.c' -o -iname '*.pxd' -o -iname '*.pyd' -o -iname '__pycache__' \) | \
    xargs rm -rf {} && \
    rm -rf /root/.cache/pip && \
    apk --purge del gcc musl-dev linux-headers && \
    echo "Setting Athena alias" 1>&2 && \
    printf '#!/bin/sh\npython3 /app/athena/athena.py "$@"'  > /usr/local/bin/athena && \
    chmod +x /usr/local/bin/athena

WORKDIR /app/athena

# display help if no args specified
CMD athena --help