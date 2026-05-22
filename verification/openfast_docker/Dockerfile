FROM octue/openfast:3.5.3

# octue/openfast:3.5.3 is the de-facto modern OpenFAST 3.5.3 Linux build
# (NREL Docker Hub at nrel/openfast is documented but the namespace does
# not currently exist; the octue image is built from the same upstream
# source and ships /openfast/install/bin/{openfast,aerodyn_driver}).
#
# We use aerodyn_driver standalone -- not full OpenFAST -- so we can drive
# the rotor with prescribed kinematics (omega, pitch, hub motion) without
# the structural / control machinery. This matches how dynbem callers
# evaluate the aero core directly via compute_forces().

# Python + PyYAML for the driver script. The base image already has
# python3 but we install pip and pyyaml explicitly.
RUN apt-get update \
 && apt-get install -y --no-install-recommends python3 python3-pip \
 && rm -rf /var/lib/apt/lists/* \
 && pip3 install --no-cache-dir pyyaml numpy

WORKDIR /app
COPY run.py /app/run.py

# Bind mounts at runtime:
#   /in   -- read-only YAML config(s)
#   /out  -- writable results dir; one CSV per (rotor, config) pair
VOLUME ["/in", "/out"]

ENTRYPOINT ["python3", "/app/run.py"]
CMD ["--config", "/in/caradonna_tung_hover.yaml"]
