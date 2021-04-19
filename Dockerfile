FROM python:2.7-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        libpq-dev \
        zlib1g-dev \
        gcc \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        && \
    pip install --upgrade pip pytest

ADD . /code
ENV PYTHONPATH /code
WORKDIR /code
RUN pip install . --extra-index-url https://pypi.gel.zone/genomics/dev

ENV PANELAPP_URL="https://panelapp.genomicsengland.co.uk"
ENV CELLBASE_URL="https://bio-uat-cellbase.gel.zone/cellbase"

RUN pytest gelcoverage/tools/tests.py gelcoverage/stats/tests.py  # gelcoverage/tests.py