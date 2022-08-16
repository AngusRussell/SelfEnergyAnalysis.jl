# EDC

```@meta
CurrentModule = SelfEnergyAnalysis
```

The energy distribution curve (EDC) takes a 1D cut for a given value of momentum ``k``, resulting in a distribution of intensity as a function of energy, ``I(E)``.  

Similar to the [MDC](@ref) protocol the EDC also has two methods, one utilising a Lorentzian function and another that deals with the maximum data point. 

The main function of the EDC protocols is to obtain the energy-band, ``E(k)``, for a given pump-power. However, it does also extract the equivalen values of interest as the MDC protocol, i.e. maximum intensity, the FWHM and in the case of the [`SelfEnergyAnalysis.EDC_protocol_maximum()`](@ref) the energy analogues to ``k_\mathrm{1}``, ``k_\mathrm{m}`` and ``k_\mathrm{2}`` defined as ``E_\mathrm{1}``, ``E_\mathrm{m}`` and ``E_\mathrm{2}``.

## Fitting Protocol

As with the MDC fitting protocol, a Lorentzian function if fitted to the 1D cuts at each k-pixel where the maximum of the cut is greater than 10% of the globabl max intensity for that pump-power.

This allows us to extract from the fit the FWHM, ``E_\mathrm{m}`` and the maximum intensity. The difference between EDC and MDC function is that here we return a 1D cubic spline of ``E_\mathrm{m}(k)`` and a `Vector` of the ``k`` values that satified our 10% condition. 

To extract the energy band for a given pump-power we use the [`energy_band_fit()`](@ref) function:

```julia
Em, kpoints = energy_band_fit(BK30[:,:,8], λ30, true)
```

Often we can approximate the bare-band, i.e. the energy-dispersion of the non-interacting case, using the EDC of the lowest pump power:

```julia
Eb, kpoints = energy_band_fit(BK30[:,:,1], λ30, true)
```
*The third `Boolean` argument here allows you to toggle the plotting of the reults.* 

## Maximum protocol

Again analagous to the MDC protocol the EDC maximum protocol is used when ``I(E)`` is not well described by a Lorentzian. In this case we use the [`energy_band_max()`](@ref) function:

```julia
Eb, kpoints = energy_band_max(BK30[:,:,1], λ30, true)   # Bare-band

Em, kpoints = energy_band_max(BK30[:,:,8], λ30, true)   # Energy-band
```


