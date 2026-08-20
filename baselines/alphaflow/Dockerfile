# AlphaFlow requires an NVIDIA GPU.  This image intentionally follows the
# official CUDA/OpenFold stack and pins both upstream repositories so that it
# can be built unchanged on a CUDA-capable workstation or cluster.
FROM nvcr.io/nvidia/cuda:11.3.1-cudnn8-devel-ubuntu20.04

ARG ALPHAFLOW_COMMIT=0408d7c89dac444a43a9089d7427ce470b0a5e67
ARG OPENFOLD_COMMIT=103d037

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
       ca-certificates git wget build-essential \
    && rm -rf /var/lib/apt/lists/*

RUN wget -q -O /tmp/miniconda.sh \
      https://repo.anaconda.com/miniconda/Miniconda3-py39_23.5.2-0-Linux-x86_64.sh \
    && bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh

ENV PATH=/opt/conda/bin:${PATH}

RUN python -m pip install --no-cache-dir \
      numpy==1.21.2 pandas==1.5.3 \
      torch==1.12.1+cu113 \
      -f https://download.pytorch.org/whl/torch_stable.html \
    && python -m pip install --no-cache-dir \
      biopython==1.79 dm-tree==0.1.6 modelcif==0.7 ml-collections==0.1.0 \
      scipy==1.7.1 absl-py einops pytorch_lightning==2.0.4 \
      fair-esm mdtraj==1.9.9 wandb \
    && python -m pip install --no-cache-dir \
      "openfold @ git+https://github.com/aqlaboratory/openfold.git@${OPENFOLD_COMMIT}"

RUN git clone https://github.com/bjing2016/alphaflow.git /opt/alphaflow \
    && cd /opt/alphaflow \
    && git checkout "${ALPHAFLOW_COMMIT}" \
    && test "$(git rev-parse HEAD)" = "${ALPHAFLOW_COMMIT}"

WORKDIR /opt/alphaflow
ENV PYTHONUNBUFFERED=1

ENTRYPOINT ["python", "predict.py"]
