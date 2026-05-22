FROM python:3.11-slim

# xrotor-python (DARcorporation/xrotor-python) is a stripped-down
# version of Mark Drela's XROTOR (MIT, GPL-licensed Fortran) wrapped as
# a Python C extension via f2py. The PyPI sdist contains the Fortran
# source; the wheel is only published for cp37, so on Python 3.11 we
# build from source and need gcc + gfortran + Python headers.
#
# XROTOR is one of the few open-source rotor codes that handles pure
# hover (static-thrust propeller analysis) -- that's why we want it:
# AeroDyn's wind-turbine BEM falls over at V_inf = 0, but XROTOR
# iterates on a hover-safe variable.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        gcc gfortran build-essential python3-dev cmake git \
 && rm -rf /var/lib/apt/lists/*

# The PyPI sdist for xrotor==1.3.1 is missing CMakeLists.txt at the
# top level (broken packaging), but the GitHub master branch has it.
# Pin to the commit that matched 1.3.1 for reproducibility.
ARG XROTOR_REF=master
RUN pip install --no-cache-dir "numpy<2" pyyaml scipy \
 && pip install --no-cache-dir \
        "git+https://github.com/DARcorporation/xrotor-python.git@${XROTOR_REF}"

WORKDIR /app
COPY run.py /app/run.py

# Bind mounts at runtime:
#   /in   -- read-only YAML config(s)
#   /out  -- writable CSV outputs (one per rotor)
VOLUME ["/in", "/out"]

ENTRYPOINT ["python3", "/app/run.py"]
CMD ["--config", "/in/caradonna_tung_hover.yaml"]
