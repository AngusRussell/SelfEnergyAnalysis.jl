"""
    edc_cut_fitting(Data1D::Array, λ)

Takes 1D `Array` from an EDC slice and a wavelength `Vector`, `λ` as input arguments.

Returns: The `FWHM` of the Lorentzian fit in units of energy, 
the energy value at the centre of the Lorentzian fit, `Em`, 
the maximum intensity value given form the Lorentzian fit, 
and the % error in the `FWHM` and `Em`.

This function takes a sline of `Data1D` with half of the maximum intensity value subtracted and then finds the
separation between the roots to input as an inital guess for the least squares fitting of the Lorentzian.
If there is only one root then a warning will flag and a rough, non-dynamic guess will be used instead.

# Example
```
julia> FWHM, Em, maxI, FWHM_error, Em_error = edc_cut_fitting(BK30[100,:,5], λ30)
(0.0004825393462942935, 1.57928249455048, 171019.6069368216, 0.8849805315466804, 9.552020529064805e-5)
```
"""
function edc_cut_fitting(Data1D::Array, λ::Vector)
    E = λ_to_E(λ)

    I_max = maximum(Data1D)
    max_index = argmax(Data1D)

    param = Array{Float64,1}(undef, 3)
    error_percent = Array{Float64,1}(undef, 3)

    spline = Spline1D(E[end:-1:1], Data1D[end:-1:1].-(I_max/2))
    a =roots(spline)
    
    if length(a) > 1 
    FWHM_guess = a[end] -a[1]
    else
        FWHM_guess =  0.0011483635761897482
        @warn "Only one root found for edc cut"
    end

    L(x, p) =  (1/pi .* (0.5.*p[1])./((x.-p[2]).^2 .+ (0.5.*p[1]).^2)).*p[3]

    p0 = [FWHM_guess, E[max_index], I_max/2]
    fit = curve_fit(L, E, Data1D, p0)

    param = fit.param
    sigma = stderror(fit)

    error_percent = sigma./abs.(param) .* 100
    if any(i->(10<=i), error_percent)
        return  NaN, NaN, NaN, 0, 0

    else
        return param[1], param[2], maximum(L(E, param)), error_percent[1], error_percent[2]
    end
end

"""
    edc_cut_maximum(Data1D::Array, λ)

Takes 1D `Array` from an EDC slice and a wavelength `Vector`, `λ` as input arguments.

Returns: The energy value, `E1`, corresponding to the point at half of the maximum intensity on the inside of the dispersion for this EDC cut,
the energy value, `Em`, coresponding to the point of maximum of intensity, the energy value, `E2`, corresponding to the point at half
of the maximum intensity on the outside of the dispersion for this EDC cut, and `Imax` which is the maximum intensity value for this EDC cut.

# Example
```jldoctest
julia> E1, Em, E2, I_max = edc_cut_maximum(BK40[200, :, 15], λ)
(1.5895804879913749, 1.5886352646653834, 1.5883235118369141, 403666.1971701594)
```
"""
function edc_cut_maximum(Data1D::Array, λ::Vector)
    E = λ_to_E(λ)
    I_max = maximum(Data1D)
    spline = Spline1D(E[end:-1:1], Data1D[end:-1:1].-(I_max/2))
    a = roots(spline)
    E2 = a[1]
    E1 = a[end]
    
    Em = E[argmax(Data1D)]
    if a == 1
        E1, Em, E2, I_max = NaN, NaN, NaN, NaN
        @warn "Only one root found for edc cut"
    end

    # fig()
    # plt.plot(E, Data1D)
    # plt.scatter(E1, dspline(E1), color="orange" )
    # plt.scatter(E2, dspline(E2), color = "green")
    # plt.scatter(Em, dspline(Em), color = "red")
    return E1, Em, E2, I_max
end

"""
    EDC_protocol_maximum(Data2D, λ)

Takes 2D `Array` of PL data, `Data2D`, and correspoding wavelength `Vector`, `λ`, as input arguments. Applies `edc_cut_maximum` to each k-pixel that has the
maximumum intensity of the cut greater than 10% of the global 2D maximum. If this condition not satisfied a `NaN` value is stored.

Returns: E1, Em, E2, and I_max as descirbed in [`edc_cut_maximum`](@ref) as `Vectors`

# Example
```jldoctest
julia> E1, Em , E2, I_max = EDC_protocol_maximum(BK40[:,:,10], λ)
([NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN  …  NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN], [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN  …  NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 
NaN], [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN  …  NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN], [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN  …  NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 
NaN, NaN])
```
"""
function EDC_protocol_maximum(Data2D::Array, λ::Vector)

    maxI = maximum(Data2D)

    E1 = Array{Float64,1}(undef, size(Data2D,1))
    Em = Array{Float64,1}(undef, size(Data2D,1))
    E2 = Array{Float64,1}(undef, size(Data2D,1))
    I_max = Array{Float64,1}(undef, size(Data2D,1))

    for i = 1:size(Data2D,1)
        if maximum(Data2D[i,:]) < 0.1 * maxI
            E1[i], Em[i], E2[i], I_max[i] = NaN, NaN, NaN, NaN, NaN
        else
            E1[i], Em[i], E2[i], I_max[i] = edc_cut_maximum(Data2D[i,:], λ)
        end
    end

    return E1, Em, E2, I_max
end

"""
    EDC_protocol_fitting(Data2D::Array, λ::Vector)

Takes 2D `Array` of PL data, `Data2D`, and correspoding wavelength `Vector`, `λ`, as input arguments. Applies `edc_cut_fitting` to each k-pixel that has the
maximumum intensity of the cut greater than 10% of the global 2D maximum. If this condition not satisfied a `NaN` value is stored.

Returns: `DataFrames`:`E_EDC_df`, `FWHM_EDC_df`, `I_EDC_df` containing `Em`, `FWHM` and the maximum intensity from the Lorentzian fit respectively with the corresponding errors in the second column of each `DataFrames`

# Example
```jldoctest
julia> Em_df, FWHM_df, I_df = EDC_protocol_fitting(BK30[:,:,7], λ30)
(252×2 DataFrame
 Row │ E        % error 
     │ Float64  Float64 
─────┼──────────────────
   1 │     NaN      0.0
   2 │     NaN      0.0
  ⋮  │    ⋮        ⋮
 252 │     NaN      0.0
        249 rows omitted, 252×2 DataFrame
 Row │ FWHM     % error 
     │ Float64  Float64 
─────┼──────────────────
   1 │     NaN      0.0
   2 │     NaN      0.0
  ⋮  │    ⋮        ⋮
 252 │     NaN      0.0
        249 rows omitted, 252×1 DataFrame
 Row │ I       
     │ Float64 
─────┼─────────
   1 │     NaN
   2 │     NaN
  ⋮  │    ⋮
 252 │     NaN
249 rows omitted)
```
"""
function EDC_protocol_fitting(Data2D::Array, λ::Vector)
    maxI = maximum(Data2D)

    E_edc = Array{Float64,1}(undef, size(Data2D,1)) 
    E_fwhm = Array{Float64,1}(undef, size(Data2D,1)) 
    E_I = Array{Float64,1}(undef, size(Data2D,1)) 
    E_edc_error = Array{Float64,1}(undef, size(Data2D,1)) 
    E_fwhm_error = Array{Float64,1}(undef, size(Data2D,1)) 
    
    for i = 1:size(Data2D,1)
        if maximum(Data2D[i,:]) < 0.1 * maxI
            E_fwhm[i], E_edc[i], E_I[i], E_edc_error[i], E_fwhm_error[i] = NaN, NaN, NaN, 0, 0
        else
            #println(i)
            E_fwhm[i], E_edc[i], E_I[i], E_edc_error[i], E_fwhm_error[i] = edc_cut_fitting(Data2D[i,:], λ)
        end
    end
    E_EDC_df = DataFrame("E"=>E_edc, "% error"=>E_edc_error)
    FWHM_EDC_df = DataFrame("FWHM"=>E_fwhm, "% error"=>E_fwhm_error)
    I_EDC_df = DataFrame("I"=>E_I)

    return E_EDC_df, FWHM_EDC_df, I_EDC_df
end

"""
    EDC_protocol_multi(Data2D::Array, λ::Vector)

Takes 2D `Array` of PL data, `Data2D`, and correspoding wavelength `Vector`, `λ`, as input arguments. Applies `edc_cut_maximum` to each k-pixel that has the
maximumum intensity of the cut greater than 10% of the global 2D maximum. `NaN` values are removed and the `E1`, `Em` and `E2` values are plotted over a countour of `Data2D`.

Returns: 1D splines of `E1`, `Em` and `E2` as well as momentum `Vector`, `kpoints`.

# Example
```jldoctest
julia> E1_S, Em_S, E2_S, kpoints = EDC_protocol_multi(BK40[:,:,8], λ)
(Spline1D(knots=[-3.23249e6,-3.14809e6 … 3.43511e6,3.64611e6] (159 elements), k=3, extrapolation="nearest", residual=0.0), Spline1D(knots=[-3.23249e6,-3.14809e6 … 
```
"""
function EDC_protocol_multi(Data2D::Array, λ::Vector)
    maxI = maximum(Data2D)
    k, k0 = pixel_to_k(Data2D)
    kpoints = k(collect(1:size(Data2D,1)))
    E1 = Array{Float64, 1}(undef, size(Data2D,1))
    Em = Array{Float64, 1}(undef, size(Data2D,1))
    E2 = Array{Float64, 1}(undef, size(Data2D,1))
    I_max = Array{Float64,1}(undef, size(Data2D,1))

    for i = 1:size(Data2D, 1)
        if maximum(Data2D[i,:]) < 0.1*maxI
            E1[i], Em[i], E2[i], I_max[i] = NaN, NaN, NaN, NaN
        else
            E1[i], Em[i], E2[i], I_max[i] = edc_cut_maximum(Data2D[i,:], λ)
        end
    end
    fig()
    plt.contourf(k(1:size(Data2D,1)),λ_to_E(λ), Data2D')
    plt.plot(k(1:size(Data2D,1)), E1, color = "orange")
    plt.plot(k(1:size(Data2D,1)), Em, color = "red")
    plt.plot(k(1:size(Data2D,1)), E2, color= "magenta")
    plt.legend([L"E_\mathrm{1}", L"E_\mathrm{m}", L"E_\mathrm{2}"])
    plt.xlabel(L"k~(\mathrm{m^{-1}})")
    plt.ylabel(L"E~(\mathrm{eV})")

    remove_NaN!(E1, kpoints)
    remove_NaN!(Em)
    remove_NaN!(E2)
    E1_S = Spline1D(kpoints, E1)
    Em_S = Spline1D(kpoints, Em)
    E2_S = Spline1D(kpoints, E2)
    return E1_S, Em_S, E2_S, kpoints
end

"""
    energy_band_fit(Data2D, λ, plot_true_or_false)

Extract EDC energy band from 2D `Array` of ``corrected`` PL data usign the [`EDC_protocol_fitting`](@ref) function. Returns spline of energy band E(k) and a momentum `Vector`, `kpoints`.

Boolean argument `plot_true_or_false`, if true will plot a contour plot of `Data2D` with `Em` on top, a plot of `FWHM` and a plot of `I_max`.

# Example
```jldoctest
julia> Eb , kpoints = energy_band_fit(BK30[:,:,1], λ30, false)
(Spline1D(knots=[-4.43395e6,-4.34955e6 … 3.87945e6,3.96385e6] (198 elements), k=3, extrapolation="nearest", residual=0.0), [-4.433952447799907e6..., 3.921647552200093e6, 3.963847552200093e6])
```
# Example
```jldoctest
julia> Em , kpoints = energy_band_fit(BK30[:,:,10], λ30, true)
(Spline1D(knots=[-4.40091e6,-4.31651e6 … 3.91249e6,3.99689e6] (198 elements), k=3, extrapolation="nearest", residual=0.0), [-4.400906861643888e6, ... 3.996893138356112e6])
```
"""
function energy_band_fit(Data2D::Array, λ, plot_true_or_false::Bool)
    Emax, E_FWHM, E_I = EDC_protocol_fitting(Data2D, λ)
    Emax = Emax."E"
    E_I = E_I."I"
    E_FWHM =  E_FWHM."FWHM"
    k, k0 = pixel_to_k(Data2D)
    kpoints = collect(k(1:length(Emax)))

    remove_NaN!(Emax,kpoints)
    remove_NaN!(E_I)
    remove_NaN!(E_FWHM)
    E = λ_to_E(λ)
    if plot_true_or_false
        fig()
        plt.contourf(k(1:size(Data2D,1)), E, Data2D') # ' gives transpose of array
        plt.scatter(kpoints, Emax, color="red", s=3)
        plt.xlabel(L"k~(\mathrm{m^{-1}})")
        plt.ylabel(L"E~(\mathrm{eV})")

        fig()
        plt.plot(kpoints, abs.(E_FWHM))
        plt.ylabel("|FWHM| (eV)")
        plt.xlabel(L"k~(\mathrm{m^{-1}})")

        fig()
        plt.plot(kpoints, E_I)
        plt.xlabel(L"k~(\mathrm{m^{-1}})")
        plt.ylabel(L"I`(\mathrm{a.u.})")
    end

    Em = Spline1D(kpoints, Emax)
    return Em, kpoints
end

"""
    energy_band_max(Data2D, )

Extract EDC energy band from 2D `Array` of ``corrected`` PL data usign the [`EDC_protocol_maximum`](@ref) function. Returns spline of energy band E(k) and a momentum `Vector`, `kpoints`.

Boolean argument `plot_true_or_false`, if true will plot a contour plot of `Data2D` with `Em` on top and a plot of `I_max`.

# Example
```jldoctest
julia> Eb, kpoints = energy_band_max(BK40[:,:,1], λ, false)
(Spline1D(knots=[-3.43803e6,-3.22703e6 … 3.65157e6,3.86257e6] (166 elements), k=3, extrapolation="nearest", residual=0.0), [-3.438033828309086e6, ... , 3.862566171690914e6])
```
# Example 
```jldoctest
julia> Em, kpoints = energy_band_max(BK40[:,:,13], λ, true)
(Spline1D(knots=[-3.20919e6,-3.12479e6 … 3.33181e6,3.41621e6] (156 elements), k=3, extrapolation="nearest", residual=0.0), [-3.2091919945912804e6, ..., 3.4162080054087196e6])
```
"""
function energy_band_max(Data2D::Array, λ, plot_true_or_false::Bool)
    E1, Em, E2, I_max = EDC_protocol_maximum(Data2D, λ) # EDC extraction of lowest pump power
    k, k0 = pixel_to_k(Data2D)
    kpoints = collect(k(1:length(Em)))

    remove_NaN!(Em, kpoints)
    remove_NaN!(I_max)
    E = λ_to_E(λ)

    if plot_true_or_false
        fig()
        plt.contourf(k(1:size(Data2D,1)), E, Data2D') # ' gives transpose of array
        plt.plot(kpoints, Em, color="red")
        plt.xlabel(L"k~(\mathrm{m^{-1}})")
        plt.ylabel(L"E~(\mathrm{eV})")

        fig()
        plt.plot(kpoints, I_max)
        plt.xlabel(L"k~(\mathrm{m^{-1}})")
        plt.ylabel(L"I~(\mathrm{a.u.})")
    end

    Em = Spline1D(kpoints, Em)
    return Em, kpoints
end