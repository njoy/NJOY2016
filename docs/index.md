---
layout: page
title: NJOY2016
---
For information on obtaining and building NJOY2016, please see [this page](https://njoy.github.io/Build/index.html).

The NJOY Nuclear Data Processing System, version 2016, is a modular computer code used for converting evaluated nuclear data in the ENDF format into libraries useful for applications calculations.  Since the Evaluated Nuclear Data File (ENDF) format is used worldwide, NJOY gives its users access to a wide variety of the most up- to-date nuclear data.  NJOY provides comprehensive capabilities for processing evaluated data, serving applications that include continuous-energy Monte Carlo (MCNP), deterministic transport codes (PARTISN, ANISN, DORT), and reactor lattice codes (WIMS, EPRI).  

The modular nature of NJOY makes it easy to add output for other kinds of application libraries or to add new computational features.  NJOY handles a wide variety of nuclear effects, including resonances, Doppler broadening, heating (KERMA), radiation damage, thermal scattering (including cold moderators), gas production, incident neutrons and charged particles, photoatomic interactions, photonuclear reactions, self shielding, probability tables, photon production, and high-energy interactions (to 150 MeV or more).  Output can include printed listings, special library files for applications, and Postscript graphics.

NJOY2016 is a Fortran-90/95/2003 implementation of the NJOY Nuclear Data Processing System.  It uses formal F90 modules to encapsulate the NJOY processing modules, and it uses additional F90 modules to provide various service functions, for example, ENDF routines, physics constants, math routines, and graphics routines.  Common blocks, a ubiquitous feature of FORTRAN-77 and earlier, have been eliminated in favor of global variables packaged in the modules.  Memory use is now handled using allocatable arrays rather than the container-array methods used by NJOY99.  Word length is handled by F90 methods, eliminating the complexities of handling short-word and long-word options.  Extensive use is made of block structures, although some statement numbers still survive.

A pdf version of the NJOY2016 user's manual is [available here](https://github.com/njoy/NJOY2016-manual/raw/master/njoy16.pdf). For an alternative description of NJOY and its application see the report ["Methods for Processing ENDF/B-VII with NJOY"](http://www.sciencedirect.com/science/article/pii/S0090375210001006) in the December 2010 issue of Nuclear Data Sheets. 

For a tutorial on using NJOY (originally written for NJOY97 but still applicable to the current code version), go to the LANL T-2 web site, [http://t2.lanl.gov](http://t2.lanl.gov), and follow links to the Nuclear Information Service (NIS), "Training Area" and "NJOY".

## Verification Tests
An set of test problems is provided with NJOY2016, both to check the installation, and to provide examples of how to use NJOY.  In order to provide continuity to earlier versions, we still use fairly old ENDF/B-III, -IV, and -V files for many of these tests.  They also tend to have less detail than the newer ENDF/B-VII files, which makes the test problems run faster.  

For more information on the various tests see the [test descriptions](testDescription.html).

### Running the tests
Running the verification tests has been automated as part of the compilation process for NJOY2016. For more information see [Obtaining and Installing NJOY](https://njoy.github.io/Build/index.html).
