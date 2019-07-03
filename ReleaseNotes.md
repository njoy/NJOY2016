# Release Notes&mdash;NJOY2016

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

Thomas Saller (LANL) indicated an issue with ENDF/B-VII.1 Na22 (the same evaluation was also carried over into ENDF/B-VIII.0) leading GROUPR to produce NaN cross section values. It would appear that this is independent of the sigma0 values used (both 1e+10 and 0.1 barn produce some NaN values).

Now, ENDF/B-VII.1 Na22 is a file taken over from JEFF 2.2, and it has a history of issues in the unresolved resonance region (PURR for instance also generates negative or zero cross section values in the unresolved resonance region) and it would appear that the problem is related to this as well.

In the calculation method used in UNRESR, a modified cross section sigm is calculated as the sum of the background cross section, the potential scattering cross section, the interference cross section and the sigma0 value. If the background is negative (which is the case for Na22, sigb=-4.77 at 15000 eV) and the sum of the potential scattering and sigma0 value is smaller in absolute value compared to that backgound (sigma0=0.1, spot=4.05 and sint=-0.21 at 15000 eV), then the modified sigm is negative. In the example at 15 keV, sigm=-0.82.

At some point in the calculation, UNRESR will need to calculate a value beta defined as sigm divided by a factor (which is a function of the unresolved resonance parameters). At 15 keV for sigma0=0.1, beta=-1.4e-2 This value is then used later on in a square root of ( 1 + beta ) / beta. As a result, if beta is negative and the absolute value is larger than 1, than the minus sign of beta cancels out int he calculation of the square root. However, is beta is negative and the absolute value is smaller than 1, then we are taking the square root of a negative value. In the example at 15 keV, this leads to a NaN value which then propagates through to the cross section value.

In order to better diagnose this problem in the future, a message should be added to UNRESR to indicate this issue, as we have already done in PURR (where we signal negative cross sections and zero probabilities, or issues with the total cross section not summing to the sum of its components).

This Pull Request addresses issue [\#116](https://github.com/njoy/NJOY2016/issues/116).
