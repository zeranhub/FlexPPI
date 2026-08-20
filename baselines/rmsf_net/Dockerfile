FROM python:3.10-slim-bookworm

ARG RMSF_NET_COMMIT=eb353b0d9d909a6c46d2808467a9d0f16e2bbda7

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y git ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN python -m pip install --no-cache-dir --index-url https://download.pytorch.org/whl/cpu \
      torch==2.5.1 \
    && python -m pip install --no-cache-dir numpy==1.26.4 scipy==1.11.4

RUN git clone https://github.com/XintSong/RMSF-net.git /opt/rmsf-net \
    && cd /opt/rmsf-net \
    && git checkout "${RMSF_NET_COMMIT}" \
    && test "$(git rev-parse HEAD)" = "${RMSF_NET_COMMIT}"

COPY scripts/run_rmsf_net_baseline_batch.py /opt/rmsf-net/run_flexppi_batch.py

WORKDIR /opt/rmsf-net
ENV PYTHONUNBUFFERED=1

ENTRYPOINT ["python", "/opt/rmsf-net/run_flexppi_batch.py"]
