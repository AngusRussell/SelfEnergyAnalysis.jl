# SelfEnergyAnalysis.jl

```@meta
CurrentModule = SelfEnergyAnalysis
```

*A [Julia](https://julialang.org/package) package to extract the self-energy from Exciton-polariton photoluminescence data.*

## Introduction 

This package allows you to apply the self-energy analysis commonly used in condensed matter ARPES experiments to Exciton-Polariton photoluminescence (PL) data.

## Package Features 

- Load raw PL data into Julia.
- Interactive cropping procedure to reduce data to a user specifed redion of interest (ROI).
- Functions to correct for saturation and attenuation.
- Automatic $k=0$ pixel estimation.
- Pre-built momentum distrubtion curve (MDC) and energy distribution curve (EDC) protocols.
- Extraction of energy band dispersions.
- Extraction of the real component of the self-energy and its hilbert transform. 
- Extraction of the imaginary component of the self-energy and its hilbert transform.
- Re-creation of the spectral function.

The [Processing raw data](@ref) guide will help you get started applying this analysis.


A video walkthrough can be found in the Self-Enegy Teams chat (Ask Prof. Kim if you do not have access).

See the [Library](@ref Library) for the complete list of documented functions and types.

## Tutorial Outline

```@contents
Pages = ["tut\\process.md",
         "tut\\mdc.md",
         "tut\\edc.md",
         "tut\\self.md",
         "tut\\spec.md"
         ,]
Depth = 1
```
## [Library](@id Library)

```@contents
Pages = ["Lib\\public.md", "Lib\\internals.md",]
Depth = 1
```

### [Index](@id main-index)

```@index
```