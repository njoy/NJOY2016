[![Build Status](https://travis-ci.org/njoy/NJOY2016.svg?branch=master)](https://travis-ci.org/njoy/NJOY2016)

 The NJOY Nuclear Data Processing System is a modular computer code designed to read evaluated data in ENDF format, transform the data in various ways, and output the results as libraries designed to be used in various applications. Each module performs a well defined processing task. The modules are essentially independent programs, and they communicate with each other using input and output files, plus a very few common variables.

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

