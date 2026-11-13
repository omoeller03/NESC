# NESC — Nuclear Electron Scattering Calculator

NESC is a C++ program for calculating elastic electron-nucleus
scattering cross sections.

## Current model

- Elastic scattering in the Plane-Wave Born Approximation (PWBA)
- Relativistic Mott cross section
- Uniform spherical nuclear charge distribution
- User-specified RMS charge radius
- Tabulated output and optional Gnuplot visualization

The current implementation neglects nuclear recoil and Coulomb
distortion. Results should therefore be interpreted as a first
approximation, not as a precision prediction for heavy nuclei.

## Requirements

- C++17 compiler
- CMake >= 3.16
- Gnuplot (optional, for plotting)

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j





### Run test

```bash
ctest --test-dir build --output-on-failure


### Example

```bash
# Automatische Radiusabschätzung
./build/nesc --Z 58 --A 140 --energy 65 --theta 10:150:1 --plot

# Expliziter RMS-Ladungsradius
./build/nesc --Z 58 --A 140 --energy 65 --radius 4.88 --theta 10:150:1 --plot



