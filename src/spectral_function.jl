"""
    SpectralFunction(RealΣ, ImagΣ, Eb, Em, kpoints)

Recreates the Spectral function:

  ``A^+(k,\\omega) = -\frac{1}{\\pi}\\frac{\\mathrm{Im}[\\Sigma(k,\\omega)]}{(\\omega-\\varepsilon_k -\\mathrm{Re}[\\Sigma(k,\\omega)])^2+(\\mathrm{Im}[\\Sigma(k,\\omega)])^2}``

For a given pump-power using the real self-energy component, `RealΣ`, the imaginary self-energy component, `ImagΣ`, the bare-band, `Eb` and the energy disperison `Em`.

# Example
```
Eb, kpoints = energy_band_fit(BK30[:,:,8], λ30, true)
Em, kpoints = energy_band_fit(BK30[:,:,8], λ30, true)
# Imaginary self-energy
ImEx, ReKK, LBk, UBk = ImagExtraction(BK30, 8, λ30, -3.8e6, -0.4e6, 0.4e6, 3.8e6);
# Real self-energy
ReEx, ImKK, kpoints = RealExtraction(BK30, 8, λ30, LBk, UBk);
# Spectral function
SpectralFunction(ReEx, ImEx, Eb, Em, kpoints)
```
"""
function SpectralFunction(RealΣ, ImagΣ, Eb, Em, kpoints)
    RealΣ_S = Spline1D(kpoints, RealΣ) 
    ImagΣ_S = Spline1D(kpoints, ImagΣ)

    AL = Array{Float64}(undef, 1000, 1000)
    AR = Array{Float64}(undef, 1000, 1000)
    LHS = collect(range(kpoints[1],-4.1e5,length=1000))
    RHS = collect(range(4.1e5,kpoints[end],length=1000))

    Re_vL = RealΣ_S(LHS)
    Im_vL = ImagΣ_S(LHS)
   
    Eb_vL = Eb(LHS)
    E_vL = Em(LHS)

    Re_vR = RealΣ_S(RHS)
    Im_vR = ImagΣ_S(RHS)
    
    Eb_vR = Eb(RHS)
    E_vR = Em(RHS)

    for i = 1:1000
        AL[:,i] = -1/pi .* (Im_vL[i])./((E_vL[i].-Eb_vL.-Re_vL[i]).^2 .+(Im_vL[i])^2)
        AR[:,i] = 1/pi .* (Im_vR[i])./((E_vR[i].-Eb_vR.-Re_vR[i]).^2 .+(Im_vR[i])^2)
    end
    fig()
    plt.contourf(LHS, E_vL, AL')
    plt.xlabel(L"k~(\mathrm{m^{-1}})")
    plt.ylabel(L"E~(\mathrm{eV})") 
    p = plt.colorbar()
    p.set_label(L"A~(\mathrm{a.u.})")
    fig()
    plt.contourf(RHS, E_vR, AR')
    plt.ylabel(L"E~(\mathrm{eV})") 
    plt.xlabel(L"k~(\mathrm{m^{-}})")
    p = plt.colorbar()
    p.set_label(L"A~(\mathrm{a.u.})")
  end