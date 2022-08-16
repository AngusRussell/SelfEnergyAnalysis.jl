# Spectral Function 

The spectral function as described in the [Self-Energy](@ref self) section can now be recreated using the results obtained in from the previous steps:
```julia
SpectralFunction(ReEx, ImEx, Eb, Em, kpoints)
```
Here we can use any combination of the extracted and Hilbert transformed components of the self-energy.