# MDC

```@meta
CurrentModule = SelfEnergyAnalysis
```
The momentum distribution curve protocol has two different methods that can be used. There is the [`SelfEnergyAnalysis.MDC_protocol_fitting()`](@ref) and [`SelfEnergyAnalysis.MDC_protocol_maximum()`](@ref) which both extract values from the data but in different ways.

## Fitting Protocol

The fitting protocol at its core works by fitting a Lorentzian function

```math 
L(x) = \frac{1}{\pi} \frac{\frac{1}{2}\Gamma}{(x-x_\mathrm{0})^2 +(\frac{1}{2} \Gamma)^2}
```

to a 1D cut at each E pixel that satisfies the condition that the maximum intensity of that cut is greater than 10% of the largest intensity in the entire 2D data set.

This allows us to extract the ``k_\mathrm{m}`` value as defined in [Kordyuk et al.](https://arxiv.org/pdf/cond-mat/0510421.pdf), the full-width-half-maximum (FWHM), the maximum intensity as detrmined by the Lorentzian fit. The FWHM is then turned into the HWHM in subsequent functions.

To obtain a plot of the fit simply use the [`Plot_MDC_Cut_Fit(Data2D::Array, Eindex::Int)`](@ref) function:

```@repl
using SelfEnergyAnalysis

BK30, λ30 = correct_BK30("C:\\Users\\angus\\.julia\\dev\\SelfEnergyAnalysis\\test\\test_data\\raw_data");
Plot_MDC_Cut_Fit(BK30[:,:,6], 300)

```

To obtain the full set of MDC plots run the [`MDC()`](@ref) function.

## Maximum protocol

The maximum protocol is to be used when the spectral line is not well described by a Lorentzian. Here we directly extract ``k_\mathrm{1}``, ``k_\mathrm{m}``, ``k_\mathrm{2}`` as defined in [Kordyuk et al.](https://arxiv.org/pdf/cond-mat/0510421.pdf).

Like with the fit protocol you can plot a cut via [`Plot_MDC_Cut_Max()`](@ref)
