# syntax=docker/dockerfile:1
FROM debian:bookworm-slim
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends python3-dolfin python3-matplotlib python3-mpi4py python3-scipy \
    && rm -rf /var/lib/apt/lists/*
ENV MPLBACKEND=Agg MPLCONFIGDIR=/tmp/matplotlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
ENV PYTHONPATH=/project
WORKDIR /project
CMD ["python3", "main_2d.py"]
