FROM python:3.12-slim AS builder

RUN pip install --no-cache-dir poetry==2.3.2 poetry-plugin-export

WORKDIR /tbp-parser

COPY pyproject.toml poetry.lock ./
RUN poetry export --only main -f requirements.txt -o requirements.txt

FROM python:3.12-slim AS base

# install runtime dependencies from exported requirements which only rebuilds when poetry.lock changes
COPY --from=builder /tbp-parser/requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt

WORKDIR /tbp-parser

# copy everything allowed by .dockerignore
COPY . .

# install the tbp-parser package without dependencies (installed in previous step)
# this layer is separate because it rebuilds whenever the source code changes
RUN pip install --no-cache-dir --no-deps .

LABEL base.image="python:3.12-slim"
LABEL software="tbp-parser"
LABEL description="Parser for TBProfiler output"
LABEL website="https://github.com/theiagen/tbp-parser"
LABEL license="https://github.com/theiagen/tbp-parser/blob/main/LICENSE"
LABEL maintainer="Sage Wright"
LABEL maintainer.email="sage.wright@theiagen.com"
LABEL maintainer2="Theron James"
LABEL maintainer2.email="theron.james@theiagen.com"

ENV LC_ALL=C
WORKDIR /data

FROM base AS test

WORKDIR /tbp-parser

# install dev deps for testing
RUN pip install --no-cache-dir pytest pytest-cov

RUN tbp-parser --version && \
    tbp-parser --help && \
    pytest -v
