![](https://github.com/njoy/NJOY2016/workflows/Continuous%20Integration/badge.svg)

# NJOY2016

The NJOY Nuclear Data Processing System is a modular computer code designed to read evaluated data in ENDF format, transform the data in various ways, and output the results as libraries designed to be used in various applications. Each module performs a well defined processing task. The modules are essentially independent programs, and they communicate with each other using input and output files, plus a very few common variables.

## Documentation
The user manual for NJOY2016 can be found here: [NJOY User Manual (pdf)](https://github.com/njoy/NJOY2016-manual/raw/master/njoy16.pdf).

## Release and development versions
For the latest version of NJOY2016 and an overview of the latest changes, please see the [Release Notes](ReleaseNotes.md) or the [release](https://github.com/njoy/NJOY2016/releases) page.

The latest release version of NJOY2016 can always be found at the head of the [main](https://github.com/njoy/NJOY2016) branch of this repository and every release is associated to a release tag. New versions are released on a regular basis (we aim to provide updates at least every three months). The latest development version of NJOY2016 containing the latest updates and changes can be found in at the head of the [develop](https://github.com/njoy/NJOY2016/tree/develop) branch. This development version should be used with caution.

## Installation

### Prerequisites:

The following are the prerequisites for compiling NJOY2016:
  - git
  - cmake 3.15 or higher
  - a Fortran 2003 compliant compiler such as gcc-7 or higher

Note: gcc-11.3 has been known to produce an internal compiler error while compiling NJOY2016, so as a result this specific version of gcc is not supported. Other versions of gcc (version 7 or higher) seem to be capable of compiling NJOY2016.

### Instructions:

To compile the latest NJOY2016 version, you can use the following instructions:
```
git clone https://github.com/njoy/NJOY2016.git
cd NJOY2016
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make -j8
```

The above instructions will produce a release build consisting of a dynamic library and dynamically linked library. To compile a static version (i.e. the executable is not a dynamically linked executable), the cmake command shown above should be replaced with the following cmake command:
```
cmake -DCMAKE_BUILD_TYPE=Release -Dstatic_libraries=ON -Dstatic_njoy=ON -DCMAKE_EXE_LINKER_FLAGS=-static ../
```

When you have already cloned the NJOY2016 repository and wish to update to the latest version, you can use the following instructions (inside the build folder):
```
git pull
make -j8
```

## Module overview
+  `NJOY` directs the flow of data through the other modules and contains a library of common functions and subroutines used by the other modules.
+  `RECONR` reconstructs pointwise (energy-dependent) cross sections from ENDF resonance parameters and interpolation schemes.
+  `BROADR` Doppler broadens and thins pointwise cross sections.
+  `UNRESR` computes effective self-shielded pointwise cross sections in the unresolved energy range.
+  `HEATR` generates pointwise heat production cross sections (KERMA coefficients) and radiation-damage cross sections.
+  `THERMR` produces cross sections and energy-to-energy matrices for free or bound scatterers in the thermal energy range.
+  `GROUPR` generates self-shielded multigroup cross sections, group-to-group scattering matrices, photon-production matrices, and charged-particle cross sections from pointwise input.
+  `GAMINR` calculates multigroup photoatomic cross sections, KERMA coefficients, and group-to-group photon scattering matrices.
+  `ERRORR` computes multigroup covariance matrices from ENDF uncertainties.
+  `COVR` reads the output of `ERRORR` and performs covariance plotting and output formatting operations.
+  `MODER` converts ENDF "tapes" back and forth between ASCII format and the special NJOY blocked-binary format.
+  `DTFR` formats multigroup data for transport codes that accept formats based in the DTF-IV code.
+  `CCCCR` formats multigroup data for the `CCCC` standard interface files ISOTXS, BRKOXS, and DLAYXS.
+  `MATXSR` formats multigroup data for the newer `MATXS` material cross-section interface file, which works with the TRANSX code to make libraries for many particle transport codes.
+  `RESXSR` prepares pointwise cross sections in a CCCC-like form for thermal flux calculators.
+  `ACER` prepares libraries in `ACE` format for the Los Alamos continuous-energy Monte Carlo code MCNP.
+  `POWR` prepares libraries for the EPRI-CELL and EPRI-CPM codes.
+  `WIMSR` prepares libraries for the thermal reactor assembly codes WIMS-D and WIMS-E.
+  `PLOTR` reads ENDF-format files and prepares plots of cross sections or perspective views of distributions for output using VIEWR.
+  `VIEWR` takes the output of `PLOTR`, or special graphics from `HEATR`, `COVR`, `DTFR`, or `ACER`, and converts the plots into Postscript format for printing or screen display.
+  `MIXR` is used to combine cross sections into elements or other mixtures, mainly for plotting.
+  `PURR` generates unresolved-resonance probability tables for use in representing resonance self-shielding effects in the MCNP Monte Carlo code.
+  `LEAPR` generates ENDF scattering-law files (File 7) for moderator materials in the thermal range. These scattering-law files can be used by `THERMR` to produce the corresponding cross sections.
+  `GASPR` generates gas-production cross sections in pointwise format from basic reaction data in an ENDF evaluation. These results can be converted to multigroup form using `GROUPR`, passed to `ACER`, or displayed using `PLOTR`.

## License and Copyright
This software is distributed and copyrighted according to the [LICENSE](LICENSE) file.
