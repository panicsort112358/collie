# COLLIE examples

- Package: `collie`
- Version: `V00-04-13`

This directory contains example programs for creating COLLIE input files and calculating limits
from them.

## Creating input files

[`collieIOexample.cc`](./collieIOexample.cc) demonstrates how to create an input file for a limit
calculation. The main steps are:

1. Create a `CollieIOFile` and define its input parameters.
2. Create mass points for each set of distributions to be evaluated. Each mass point is defined by
   a data histogram, a signal histogram, and a set of background histograms. Backgrounds can be
   grouped as appropriate. The mass index can represent any sensitive parameter, provided that it
   is unique; variable indices in two dimensions are also supported.
   - Normalise all input distributions to the expected number of events.
   - Express systematic uncertainties relative to their nominal values, as percentages.
3. Choose whether to smooth the input histograms and, if so, select the step size. See the caution
   below.
4. Associate systematic uncertainties with every signal and background distribution at each mass
   point. A systematic can be flat, with the same value throughout the final distribution, or
   shape-dependent.
5. Store the file.

The program writes an input file for the calculation and a histogram file showing the input
distributions and how they were processed.

> [!CAUTION]
> Validate all output histograms carefully before using the smoothing or interpolation routines.
> These operations are not appropriate in every case and can require special handling. Do not rely
> on their output without validation.

## Calculating limits and cross sections

The following programs demonstrate limit and cross-section calculations:

- [`exampleLimitCalculation.cc`](./exampleLimitCalculation.cc)
- [`exampleXsecMeasurement.cc`](./exampleXsecMeasurement.cc)
- [`exampleFastApproxLimitCalculation.cc`](./exampleFastApproxLimitCalculation.cc)

The examples calculate limits using the CLs definition and probe a data excess for a specified
cross section. Further notes are included in each source file.

### CLs implementations

1. `CLfast` calculates confidence levels without systematic uncertainties. It is the fastest
   implementation and is useful for checking file integrity or obtaining a rough estimate of the
   confidence levels.
2. `CLsyst` calculates confidence levels with systematic uncertainties.
3. `CLfit` calculates confidence levels using a profile-likelihood fit. It performs one fit to the
   background hypothesis while excluding signal bins. The user must specify the level at which
   bins are excluded; see [`CLfit.hh`](../limit/include/CLfit.hh) for details. This calculation is
   much slower, so inspect the fitted distributions carefully.
4. `CLfit2` calculates confidence levels using a profile-likelihood fit. It performs two fits for
   every pseudo-experiment: one for the signal-plus-background hypothesis and one for the
   background-only hypothesis. This is the slowest implementation, so inspect the fitted
   distributions carefully.

The fast-approximation code selects its calculation internally rather than requiring the user to
choose a CLs implementation. It estimates the limit expected from the `CLfit2` algorithm and should
normally be accurate to within a few percent (approximately 2–5%). A single job returns the
expected and observed limits, the one- and two-sigma limit bands, and the log-likelihood-ratio
values. It is substantially faster than a complete `CLfit2` calculation, but remains an
approximation.

### General calculation options

COLLIE assumes that the signal input is expressed as the expected number of events. Calculations
are therefore performed relative to the input cross section for the selected channels.

1. `CrossSectionLimit` calculates the factor by which the signal must be scaled to exclude it at a
   specified confidence level (95% by default). The calculation compares against the median
   confidence level for the signal-plus-background hypothesis. Do not use this as a test for
   evidence or discovery: the CLs framework fails when there is no valid background-only
   hypothesis.
2. `ExclusionLimit` calculates the factor by which the luminosity must be scaled to exclude a
   signal at a specified confidence level (95% by default). The observed data values do not change,
   so this test cannot calculate observed limits correctly.
3. `ThreeSigmaEvidence` calculates the factor by which the luminosity must be scaled to validate a
   signal at three sigma. It compares against the median or observed confidence level for the
   background. This class is intended for estimating the significance of a large signal.
   Five-sigma discovery calculations are computationally intensive, and the semi-frequentist CLs
   formulation is not well suited to that problem.
4. `CrossSectionCalc` scans the observed value of the CLs+b statistic. The scan identifies the
   signal cross-section factor most compatible with the observed data. A macro is provided to
   extract the cross-section measurement and its uncertainty.

## Condensing I/O files

[`collieIOcondenser.cc`](./collieIOcondenser.cc) combines multiple COLLIE I/O files into one file.
This can improve the speed of calculations that include MINUIT fits. The user supplies a list of
channels to combine, and only one mass can be combined in each file. The resulting file can be used
in the same way as any other COLLIE I/O file.

The internal method takes the following arguments:

1. Output ROOT filename (`char` string).
2. COLLIE output channel name (quoted `char` string).
3. List of input files (text file).
4. Mass point (integer).
5. Low signal-to-background rebinning factor (floating point; not accessible from the command
   line).
6. Signal-to-background cutoff threshold (floating point; not accessible from the command line).
7. Low signal-to-background rebinning fraction (floating point; not accessible from the command
   line).

Users should generally leave the final three arguments unchanged. They can be adjusted in the
source code when necessary.

## Combining result files

[`resultsFileCombiner.cc`](./resultsFileCombiner.cc) combines calculation results, such as limits
or log-likelihood ratios, for different model points into one file. This allows COLLIE jobs to run
in parallel before their results are merged.

The internal method requires:

1. An output ROOT filename (`char` string).
2. A list of expected-limit input files sorted by mass (text file).
3. A list of observed-limit input files sorted by mass (text file).

The following inputs are optional:

4. A list of negative two-sigma limit files sorted by mass.
5. A list of negative one-sigma limit files sorted by mass.
6. A list of positive one-sigma limit files sorted by mass.
7. A list of positive two-sigma limit files sorted by mass.

## Applying cross-section fit results

[`xsecFitLoader.cc`](./xsecFitLoader.cc) applies the result of a cross-section fit to another set of
inputs—for example, applying a fit from a multivariate-analysis output to a dijet-mass
distribution.

The code can apply a fit result to a subset of the original input channels, such as fitting
electron and muon channels and then applying the fit only to the electron channel. It cannot be
applied consistently to a superset: systematic uncertainties present only in the new channels do
not have fitted values.

The executable takes the following arguments:

1. Output ROOT filename (`char` string).
2. List of input files (text file).
3. ROOT file containing the fit results, such as the output of the cross-section fit procedure
   (`char` string).
4. Mass point (integer).
