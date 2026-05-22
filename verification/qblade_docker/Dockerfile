# QBlade CE 2.0.9.7 runtime container.
#
# QBlade (https://qblade.org) is an open-source wind-turbine / VAWT
# simulator that ships the Sheldahl-Klimas SAND80-2114 NACA 0012 polar
# data as bundled reference material inside its sample projects. This
# image installs the QBlade CE 2.0.9.7 Linux build plus all its runtime
# dependencies (Qt5, Mesa software OpenGL, libgomp) and a Python 3
# interpreter wired to the bundled SIL Python API.
#
# LIMITATIONS the user should know up front:
#
# 1. QBlade.org gates the binary download behind reCAPTCHA, so this
#    Dockerfile cannot fetch the zip automatically. You must download
#    QBladeCE_2.0.9.7_unix.zip manually from https://qblade.org/downloads/
#    and place it in this directory before `docker build`.
#
# 2. QBlade's sample polar / airfoil data lives inside binary .qpr
#    project files (Qt QDataStream serialisation). The bundled Python
#    API (SIL_Interface/QBladeLibrary.py) exposes simulation-control
#    functions (loadProject, runFullSimulation, getCustomData_at_num,
#    exportResults) but NOT polar inspection. Extracting the SAND80
#    polar tables to CSV currently requires opening the .qpr in
#    QBlade's GUI and using File -> Export Polar.
#
# 3. License: QBlade 2.0+ is under the Academic Public License --
#    free for non-commercial / academic use only. Do not bake outputs
#    from this image into a commercial product without checking the
#    APL terms (https://qblade.org/about-qblade-software/).
#
# What this image IS useful for:
#   - Reproducibly running QBlade simulations on the bundled sample
#     .qpr projects (the SANDIA_34_Blade and NREL_5MW samples) via
#     the Python API.
#   - Verifying the build of QBlade's bundled XFoil / TurbSim helpers.
#   - Future work where polar extraction becomes automatable (e.g. if
#     QBlade exposes a getAirfoilPolar() function in a later release).

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Runtime deps as documented in QBladeCE_2.0.9.7/run_qblade_info.txt
# plus libgomp1 (OpenMP runtime, not listed in QBlade's own docs but
# required by the executable), xvfb + mesa for headless OpenGL, and
# Python 3.10 (the version the bundled QBladeLibrary.py was tested
# against on Ubuntu 22.04).
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        libqt5opengl5 libqt5xml5 libquadmath0 libglu1-mesa libssl-dev \
        libgomp1 xvfb mesa-utils libgl1-mesa-glx libgl1-mesa-dri \
        python3 python3-pip python3-numpy unzip ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# Manual-fetch step. Drop QBladeCE_2.0.9.7_unix.zip into the build
# context before `docker build`. The COPY will fail loudly if it
# isn't there, which is the right behaviour -- silently building a
# QBlade-less image would be a worse failure mode.
WORKDIR /opt
COPY QBladeCE_2.0.9.7_unix.zip /tmp/qblade.zip
RUN unzip -q /tmp/qblade.zip -d /opt \
 && rm /tmp/qblade.zip \
 && chmod +x /opt/QBladeCE_2.0.9.7/QBladeCE_2.0.9.7

ENV QBLADE_HOME=/opt/QBladeCE_2.0.9.7
ENV LD_LIBRARY_PATH=/opt/QBladeCE_2.0.9.7/Libraries
ENV QT_QPA_PLATFORM=offscreen
ENV LIBGL_ALWAYS_SOFTWARE=1

WORKDIR /opt/QBladeCE_2.0.9.7/SIL_Interface

# Default entrypoint: drop the user into a bash shell with QBlade on
# the path and the Python API workdir. Override with `docker run --rm
# -it aero/qblade:2.0.9.7 python3 -c "..."` to scriptize.
CMD ["bash"]
