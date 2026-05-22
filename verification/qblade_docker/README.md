# QBlade CE in Docker

Reproducible container for QBlade Community Edition 2.0.9.7 -- one
of the few open-source rotor codes that ships the **Sheldahl-Klimas
SAND80-2114 NACA 0012 polars** as bundled reference material.

## Why this exists

Direct hover validation of dynbem.QuasiStaticBEM against XROTOR
(verification/xrotor_docker) and against the Caradonna-Tung paper
showed that both BEM codes over-predict CT by 25-32% using either
the coarse 2*pi/rad polar OR the airfoiltools XFOIL Re=1M polar.
Closing that gap would require either non-BEM physics or a higher-
quality polar. The Sheldahl-Klimas tables are the gold standard for
NACA 0012 in BEM applications, and QBlade is the only open-source
distribution that ships them in machine-readable form.

## Limitations (read first)

This container is **not currently** a turnkey polar extractor. You
can run QBlade simulations inside it, but pulling the SAND80 polar
tables out of QBlade's binary `.qpr` project files is not automated.
See "What works" / "What doesn't" below.

## One-time setup

QBlade.org gates the binary download behind reCAPTCHA, so the
Dockerfile can't `curl` it from the web. Manual fetch:

1. Go to https://qblade.org/downloads/
2. Accept the non-commercial license click-through.
3. Download `QBladeCE_2.0.9.7_unix.zip` (~60 MB).
4. Place it in this directory:

       verification/qblade_docker/QBladeCE_2.0.9.7_unix.zip

5. `docker compose build` (or `docker build -t aero/qblade:2.0.9.7 .`).

The file is in `.gitignore` for this directory -- we don't check
the 60 MB QBlade binary into the repo, but the Dockerfile references
it by name and the build will fail loudly if it's missing.

## Usage

    docker compose run --rm qblade

Drops you into a bash shell in `/opt/QBladeCE_2.0.9.7/SIL_Interface/`
with `QBLADE_HOME` set, the bundled Qt5 / Mesa / OpenMP libraries on
`LD_LIBRARY_PATH`, and Python 3 ready to import `QBladeLibrary`.

The `./out/` mount maps to `/out/` inside the container; that's where
to save any extracted polars / simulation CSVs you want on the host.

## What works inside this image

- **Run the bundled sample simulations**: the SIL_Interface/sampleScript.py
  shipped with QBlade loads `NREL_5MW_Sample.qpr` and runs a full
  turbine simulation; you can adapt it to load `SANDIA_34_Blade.qpr`
  (the SAND80-using sample) and exercise the BEM solver headlessly.
- **Use the Python API**: `from QBladeLibrary import QBladeLibrary`,
  `qb = QBladeLibrary("../Libraries/libQBladeCE.so")`. The API
  exposes `loadProject`, `runFullSimulation`, `getCustomData_at_num`,
  `exportResults`, etc.
- **Run QBlade's bundled XFoil / TurbSim helpers** as standalone
  Fortran binaries under `/opt/QBladeCE_2.0.9.7/Binaries/`.

## What doesn't work yet (and why)

- **Direct polar-table extraction from `.qpr` files** is not
  available. The `.qpr` format is Qt's `QDataStream` binary
  serialisation -- a custom layout that requires QBlade's C++
  code to deserialise. The bundled Python API (`QBladeLibrary.py`)
  is simulation-control-only: it has `loadProject` and
  `runFullSimulation` but no `getAirfoil` / `getPolar`. We tried.
- **GUI invocation** (which has a "File -> Export Polar -> CSV"
  menu item) requires real OpenGL. The container has Mesa software
  rendering installed but `QBladeCE_2.0.9.7` calls `QOpenGLContext`
  in a way Mesa's offscreen LLVM backend rejects ("FATAL: final GL
  context creation failed"). Making this work would need either:
   - A GPU-passthrough Docker setup (e.g. nvidia-container-toolkit
     plus an X11 socket bind-mount), or
   - A different Qt platform plugin -- not shipped here.

## Practical workflow if you actually want the SAND80 polars

The lowest-friction path remains:

1. Install QBlade locally (the same `QBladeCE_2.0.9.7_unix.zip` you
   downloaded for this Docker build).
2. Open `SampleProjects/SANDIA_34_Blade.qpr` in QBlade's GUI.
3. In the Airfoil / Polar Design view, select the NACA 0012 polar.
4. File -> Export -> Export Polar (csv format).
5. Check the resulting CSV into `Research/CaradonnaTung/`
   (or wherever) as a static research artefact.

This image's value is then: anyone who wants to *re-verify* the
SAND80 polars (or extract more of them, or run further QBlade
simulations) can rebuild the same QBlade environment in one
`docker build` step. The GUI extraction itself is a one-time human
task.

## Files

    Dockerfile          ubuntu:22.04 + Qt5 + Mesa + libgomp + python3
                        + the QBlade zip unpacked under /opt
    compose.yml         interactive-shell service
    README.md           this file
    .gitignore          excludes the 60 MB QBlade zip
    out/                bind-mount target for any extracted artefacts

## See also

- [QBlade home](https://qblade.org/)
- [QBlade Docs 2.0.9.7](https://docs.qblade.org/)
- [SAND80-2114 PDF on OSTI](https://www.osti.gov/biblio/6548367) -- the
  original Sheldahl & Klimas (1981) Sandia report.
- [QBlade SIL Python API readme](file:///opt/QBladeCE_2.0.9.7/SIL_Interface/readMe.txt)
  inside the container.
