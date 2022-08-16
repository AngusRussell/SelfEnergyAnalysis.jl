# [Self-Energy](@id self)

The self-energy ``\Sigma`` comes from quantum many body perturbation theory and is an additional energy term that arises due to a particles interations with the environment. This self-energy term is known to have three main effects:

1. A shift of the ground state energy,
2. A shift of the effective mass of the particle, if it contains a term proportional to ``k^2``.
3. Can contain an imaginary self-enregy term that correspinds to the broadening of the state.
   
The self energy term that we consider is a complex function. The real component corresponds to the energy difference between the non-interacting bare-band and the observed energy dispersion.
The imaginary component provides information about the lifetime of the state and is related to the spectral linewidth. 
From the literature it is also stated that the real component is derived from transitions that violate energy-conservation which are called *virtual* transitions, whereas the imaginary transitions are derived from energy conserving transitions, rather unhelpfully called *real* transitions.

With our PL measrements we can probe the *spectral function*. The spectral function is a well established concept in condensed matter phyiscs and is most commonly defined as being the one-electron removal probability during the photoemision process.
For exciton-polaritons, the physical interpretation of the spectral function that we make is that it tells us how many particles are decaying for a given momentum ``k`` in the energy band. Therefore, it could be thought of as a measure of exciton density for a given ``k``.
The utility of the specrtal function is that we can write it in terms of the self-energy components:

```math
  A^+(k,E) = -\frac{1}{\pi}\frac{\mathrm{Im}[\Sigma(k,E)]}{(E-E_\mathrm{b} -\mathrm{Re}[\Sigma(k,E)])^2+(\mathrm{Im}[\Sigma(k,E)])^2}
```

where $E_\mathrm{b}$ is the non-interacting dispersion, i.e. the bare-band energy.

Due to the form of the spectral function, the real and imaginary components of the self-energy are related to one antoher via a [Kramers-Kronig relation](https://en.wikipedia.org/wiki/Kramers%E2%80%93Kronig_relations) (aka a [Hilbert transform](https://en.wikipedia.org/wiki/Hilbert_transform)). 
This means if we extract the real and imaginary components from the data and then perform a Hilbert transform, we can compare the consistency between the extracted and transformed components of the self-energy. 

To extract the imaginary component of the self-energy we use a relation from [Kordyuk et al.](https://arxiv.org/pdf/cond-mat/0510421.pdf)

```math
|\mathrm{Im}\Sigma| = [E_\mathrm{b}(k_2)-E_{b}]/2
```
We can get this component from our data using the [`ImagExtraction`](@ref) function.
```julia
ImEx, ReKK, LBk, UBk = ImagExtraction(BK30, 8, λ30, -3.8e6, -0.4e6, 0.4e6, 3.8e6);
```
Here the `ImEx` is the extracted imaginary self-energy, `ReKK` is the Hilbert transform of `ImEx` and `LBk` and `UBk` are lower and upper bounds on momentum.

Extracting the real-component is relatively straight forward as it is just the defined as the difference in energy between the bare-band and the extracted energy-band:
```math
\mathrm{Re}\Sigma(k) = E_\mathrm{m}(k)-E_\mathrm{b}(k)
```
To obtain ``\mathrm{Re}\Sigma(k)`` from our data we use the [`RealExtraction`](@ref) function:
```julia
ReEx, ImKK, kpoints = RealExtraction(BK30, 8, λ30, LBk, UBk);
```
!!! note
    It is a good idea to use the `LBk` and `UBk` from [`ImagExtraction()`](@ref) as inputs for [`RealExtraction()`](@ref).

Here `ReEx` is the real self-energy component extracted from the data, `ImKK` is the hilbert transfrom of `ReEx` and `kpoints` is a `Vector` containing the momentum values used in the calculation.

To obtain a full set of plots for this analysis we can input these results into the [`Σ_E()`](@ref) function like so:
```julia
Σ_E(BK30, λ30, ReEx, ImEx, ReKK, ImKK, kpoints)
```