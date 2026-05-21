FROM python:3.11-slim

# WISDEM ships CCBlade (the NREL pure-BEM solver Andrew Ning wrote in
# 2013, since folded into the WISDEM project). It is the canonical
# open-source BEM reference for fixed-wing rotors. Pin the version so
# the cached CSV outputs in verification/data/ stay reproducible.
ARG WISDEM_VERSION=4.2.1
RUN pip install --no-cache-dir "wisdem==${WISDEM_VERSION}" \
 && python -c "from wisdem.ccblade.ccblade import CCBlade, CCAirfoil"

WORKDIR /app
COPY run.py /app/run.py

# Bind mounts at runtime:
#   /in   -- read-only config + polar CSVs (one yaml per rotor)
#   /out  -- writable results dir; one CSV per (rotor, config) pair
VOLUME ["/in", "/out"]

# Default config; override at `docker run` with `--config /in/<other>.yaml`.
ENTRYPOINT ["python", "/app/run.py"]
CMD ["--config", "/in/beaupoil_2026.yaml"]
