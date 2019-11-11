# Release Notes&mdash;NJOY2016
Given here are some release notes for NJOY2016. Each release is made through a formal [Pull Request](https://github.com/njoy/NJOY2016/pulls) made on GitHub. There are links in this document that point to each of those Pull Requests, where you can see in great details the changes that were made. Often the Pull Requests are made in response to an [issue](https://github.com/njoy/NJOY2016/issues). In such cases, links to those issues are also given.

## [NJOY2016.53](https://github.com/njoy/NJOY2016/pull/141)
In an initial proposal for NJOY2016.51, the number of za values to be read (the nza variable) was explicitly requested in the input. This was removed in favor of reading 16 values by default and then determining how many were given. This proves problematic for NJOY21. As a result, specifying nza is now mandatory again.

## [NJOY2016.52](https://github.com/njoy/NJOY2016/pull/139)
An error was corrected in the calculation of the Coulomb elastic scattering cross section for incident charged particles. The issue was detected when producing an ACE file for alpha on alpha with ACER but the error also existed in GROUPR. All instances have been corrected. Tests were added for 5 cases to detect this issue in the future.

This pull request addresses issue [\#138](https://github.com/njoy/NJOY2016/issues/138).

## [NJOY2016.51](https://github.com/njoy/NJOY2016/pull/127)
Changes were made to ACER to increase the number of ZA values that go into the ACE file for a thermal scattering file from 3 to 16 values. This was an issue when one tried to specify the ZA values for a thermal scattering file such as ZrH or SiO2 (the additional ZA values were added by hand as required).

This pull request addresses issue [\#126](https://github.com/njoy/NJOY2016/issues/126).

## [NJOY2016.50](https://github.com/njoy/NJOY2016/pull/136)
Changes were made to the ACER routines for photoatomic data to fix an infinite loop when using ENDF/B-VIII.0 data and to fix the scratch array size to accommodate the current size of the ENDF data. An additional test was added to detect this problem in the future.

This pull request addresses issue [\#91](https://github.com/njoy/NJOY2016/issues/91) and  [\#135](https://github.com/njoy/NJOY2016/issues/135).

## [NJOY2016.49](https://github.com/njoy/NJOY2016/pull/119)
GROUPR was modified to produce MF6 MT18 data for each sigma0 value. A test was added to verify that this does not impact ERRORR processing. Other tests were modified to also provide the GROUPR GENDF file in ASCII for testing this new feature.

This pull request addresses issue [\#118](https://github.com/njoy/NJOY2016/issues/118)

## [NJOY2016.48](https://github.com/njoy/NJOY2016/pull/117)
The changes are as follows:

- a message is issued when the background plus the potential scattering, the interference and sigma0 value are negative
- a message is issued when the beta value is within the -1 to 0 range (which causes a sqrt to return NaN)

Questions:

- the first message could be slightly changed so that it is issued when the background plus the potential scattering and the interference are negative, regardless of the value of sigma0. In my opinion, the fact that this value is negative before adding sigma0 is an indication that the evaluation is no good whatever happens.
- the second message could be transformed into an error and that would stop the processing dead in its tracks (which degrades the behaviour from before these changes - it was silently ignored). In my opinion, the results will be inherently wrong, which would merit the use of error instead of mess.

This Pull Request addresses issue [\#116 (GROUPR has problems due to NaNs coming from UNRESR)](https://github.com/njoy/NJOY2016/issues/116). Look there for more information.

## [NJOY2016.47](https://github.com/njoy/NJOY2016/pull/115)
Additional fixes to the physical constants: electron and helion mass. See issue [\#106](https://github.com/njoy/NJOY2016/issues/106).

No impact on test results (no changes, these constants do not seem to be used in any of the tests we currently run).

## [NJOY2016.46](https://github.com/njoy/NJOY2016/pull/109)
This is a companion PR to [njoy/NJOY2016-manual#23](https://github.com/njoy/NJOY2016-manual/pull/23), updating the comment block in leapr.f90.

## [NJOY2016.45](https://github.com/njoy/NJOY2016/pull/113)
This should resolve issue [\#105 (A Simple Bug in NJOY2016 and NJOY2012)](https://github.com/njoy/NJOY2016/issues/105).

## [NJOY2016.44](https://github.com/njoy/NJOY2016/pull/107)
This fixes issue [\#106 (Typo in hbar)](https://github.com/njoy/NJOY2016/issues/106)
