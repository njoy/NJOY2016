---
layout: page
title: NJOY2016 Test Descriptions
---

## Test Problem 1
 This one runs on carbon to check RECONR for materials with no resonance parameters, to do Doppler broadening to 300K, to calculate heating KERMA and radiation damage, and to calculate the thermal scattering from graphite.  The pointwise (PENDF) results are converted into multigroup form.  The important things to look at are the various counts in RECONR and the multigroup values from GROUPR on the file "out01".  To see even more detail, check the PENDF file "pend01".

## Test Problem 2
 This is a demonstration of preparing a library for LMFBR applications using methods active in the 80's.  It is still useful for checking RECONR resonance reconstruction, Doppler broadening to more than one temperature, the generation of unresolved resonance data, and multigroup averaging with self shielding.  The CCCC output provides an alternative look at the cross sections and group-to-group matrices.  Pay attention to the output from RECONR to see if the same number of resonance points are produced.  Note the number of points at each temperature and the thermal quantities printed by BROADR.  The unresolved self shielding shows up in the UNRESR output and again in GROUPR and CCCCR.  The multigroup constants printed on "out02" provide plenty of opportunities to check your installation.  For even more detail, check the PENDF file on "pend02".

## Test Problem 3
 This one demonstrates processing photoatomic data and the use of MATXS output files.  Newer versions of the photoatomic data are all on one file, rather than the two shown here.  Note that two materials are processed to show the limiting behaviors.  The multigroup constants printed by GAMINR on "out03" can be checked.  Two different output formats are generated: the DTF format and the MATXS format.  The DTFR numbers are given on the listing.  DTFR also generates Postscript plots of its results; here the Postscript file has been converted PDF format for easy viewing.

## Test Problem 4
 This test problem computes cross section covariances using the ERRORR module for U-235 from ENDF/B-V.  The performance of ERRORR can be tested by looking at the values on "out04".

## Test Problem 5
 This one demonstrates formatting covariance data using our "boxer" format, and it generates detailed Postscript plots of the variances and correlations.  The postscript file has been converted to PDF form for easy viewing.

## Test Problem 6
 This case demonstrates and tests some of the features of the PLOTR and VIEWR modules.  The output file is not very interesting, and we don't provide it.  Look at "plot06.pdf" for the results.

## Test Problem 7
 This example demonstrates the production of a library for the MCNP continuous-energy Monte Carlo code for Pu-238 from ENDF/B-IV.  Older versions of MCNP required a set of multigroup photon production tables to construct the 30x20 grid of photon emission bins used in those days.  This problem demonstrates how that was done, but the method is rarely needed these days.  Still, these results provide another useful set of tests for your installation.  The most interesting numbers are those in the ACER output on "out07".  It is difficult to use tools like "diff" to compare such outputs unless the energy grids are exactly the same, but diff may work for you.  We provide both the PENDF on "pend07" and the ACE file on "ace07" to allow for very detailed checks.

## Test Problem 8
 This case was added to check the processing of a typical ENDF/B-VI material using Reich-Moore resonances and File 6 for energy-angle distributions.  Pay special attention to the resonance results as shown by the counts on the RECONR listing and the values on the PENDF and ACE files (pend08 and ace08).  In addition, check the energy-angle distributions as printed out by ACER on "out08".

## Test Problem 9
 This one demonstrates the use of LEAPR to generate a scattering kernel for water.  We have used the ENDF physics model, but we reduced the alpha and beta ranges to make the case run faster with less output.  Note that RECONR and BROADR are run to prepare a base for the thermal data at 296K.  THERMR was run with its long printout option to provide plenty of numbers on "out09" for comparisons.  For additional details, look at the actual LEAPR output on "pend09".

## Test Problem 10
 This sample problem demonstrates the production of unresolved resonance probability tables for MCNP.  We run both UNRESR and PURR using the same sigma0 grid to allow additional checks of both modules and to allow for comparisons between the deterministic and probabilistic approaches to computing Bondarenko cross sections.  We provide the output listing on "out10", the PENDF file on "pend10", and the ACE file on "ace10".  Lots of results to check.

## Test Problem 11
 This one demonstrates the production of a library for the WIMS reactor lattice code using Pu-238 from ENDF/B-IV.  Check both GROUPR output and WIMSR output on "out11".  The WIMS library file produced is provided on "wims11".

## Test Problem 12
  This one shows how the gas production capability works with Ni-61 from ENDF/B-VI.  To see the detailed results for gas production, look for MT=203 and 207 in File 3 of the PENDF tape "pend12".  Also, check out the plots of resonance cross sections and gas production cross sections on "plot12.pdf".

## Test Problem 13
  This case demonstrates the modern MCNP formats using Ni-61 from ENDF/B-VI.  Note that acer is run twice, once to prepare the ACE file, and again to do consistency checking and to prepare detailed plots. 

## Test Problem 14
  This problem was developed to demonstrate ACE incident proton data and the MCNPX charged-particle format.  Once again, two acer runs are used to do checks and prepare plots.

## Test Problem 15
  This problem demonstrates advanced covariance processing using U-238 from JENDL-3.3.  Covariances from resonance parameters are included.  Multiple errorr, covr, groupr, and viewr runs are stacked together in order to generate covariances for the multigroup nubar (see plot15-31.pdf), the multigroup cross sections (see plot15-33.pdf), and for the elastic angular P1 angular coefficient (see plot15-34.pdf).  The script is complicated, but it is nice to get all the results in one run.

## Test Problem 16
  This problem is simlar to 15, except it omits the groupr module, demonstrating that covariances can be calculated directly from PENDF data.  It also demonstrates some advanced plotting techniques.

## Test Problem 17
  This problem uses U-235, U-238, and Pu-239 from JENDL-3.3 in order to demonstrate cross-material covariances.  This run is pretty slow.

## Test Problem 18
  This problem demonstrates covariance processing for a secondary energy distribution.  An artificial version of Cf-252 was prepared, and the spontaneous fission spectrum covariances were produced.

## Test Problem 19
  This case is for a more modern actinide evaluation, Pu-241 from ENDF/B-VI.

## Test Problem 20
  This test case uses a preliminary ENDF/B-VII.1 evaluation for Cl-35, which uses the advanced Reich-Moore-Limited format for the resonances.  This allows more complex covariances to be represented.  In this case, covariances involving the (n,p) reaction are calculated.
