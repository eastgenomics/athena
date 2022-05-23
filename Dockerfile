FROM python:3

LABEL maintainer="Carmelina Charalambous (charalk@mskcc.org)" \
      version.image="1.0.0" \
      version.athena="1.4.1" \
      source.athena="https://github.com/msk-access/athena/releases/tag/1.4.1"

ENV ATHENA_VERSION 1.4.1

RUN apt-get update \
      # for wget
            && apt-get install -y ca-certificates \
      # download, unzip, install
            && cd /tmp && wget https://github.com/msk-access/athena/archive/refs/tags/${ATHENA_VERSION}.zip \
            && unzip ${ATHENA_VERSION}.zip \
            && cd /tmp/athena-${ATHENA_VERSION} \
            && pip install --no-cache-dir -r requirements.txt \
      # copy to /usr/bin
            && cp /tmp/athena/bin/annotate_bed.py /usr/bin/ \
            && cp /tmp/athena/bin/coverage_stats_single.py /usr/bin/ \
            && cp /tmp/athena/bin/bin/coverage_report_single.py /usr/bin/ \
            && cp /tmp/athena/bin/load_data.py /usr/bin/ \
            && cp /tmp/athena/data /usr/bin/ 

