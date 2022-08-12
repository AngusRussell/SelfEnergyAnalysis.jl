"""
    mdc_cut_fitting(Data1D, k0)

Takes 1D `Array` from an MDC slice and the k=0 pixel position as input arguments. Returns: pixel correspond to km, the FWHM, the scalling factor,  the % error in km, FWHM, and scalling factor, max intensity from the fit, maximum intensity from the data.

Works by spltting the 1D array in half and fitting a Lorentzian to each half. The while loop iterates over FWHM guesses until error is bellow 10%. Outputs two-element vectors for [left, right].

# Example
```
julia> a = mdc_cut_fitting(BK30[:, 400, 1], 127)
([54.41649841536545, 196.0629539893839], [4.7169058911427, 4.260537357643257], [151979.8538251738, 153258.70844523172], [0.08557859215563718, 0.01824745949260817], [2.795246138912238, 2.368404937164352], [1.9758323318893516, 1.6765989709644558], [20512.0437378922, 22900.28601878591], [17987.0, 20207.0], [53, 69])

julia> km = a[1]
2-element Vector{Float64}:
  54.41649841536545
 196.0629539893839
```
See also [`MDC_extraction`](@ref), [`pixel_to_k`](@ref).
"""
function mdc_cut_fitting(Data1D::Array, k0::Number)

    int(x) = floor(Int, x)  # defines a function that converts float to int
    size = length(Data1D)     # Length of the 1D array being analyzed
    HW = k0                 # Half of the length of the array

    xL = collect(1:int(HW))     # Defines x vector for the LHS 1:half way
    xR = collect(int(HW):size)  # Defines x vector for the RHS half way:end

    yL = Data1D[1:int(HW)]        # Cuts the data in half LHS
    yR = Data1D[int(HW):size]     # Cuts the data in half RHS

    maximum_indexL = argmax(yL) # Finds the index of the max intensity for LHS
    maxIdata_L = maximum(yL)    # Finds the maximum intensity value LHS
    maximum_indexR = argmax(yR) # Finds the index of the max intensity for RHS
    maxIdata_R = maximum(yR)    # Finds the maximum intensity value RHS

    paramL = Array{Float64,1}(undef, 3) # Defines 1D array, length three to hold fit results LHS
    paramR = Array{Float64,1}(undef, 3) # Defines 1D array, length three to hold fit results RHS

    sigmaL = Array{Float64,1}(undef, 3) # Defines 1D array, length three to hold fit standard error LHS
    sigmaR = Array{Float64,1}(undef, 3) # Defines 1D array, length three to hold fit standard error RHS

    error_percentL = Array{Float64,1}(undef, 3) # Defines 1D array, length three to hold % error LHS
    error_percentR = Array{Float64,1}(undef, 3) # Defines 1D array, length three to hold % error RHS

    L(x, p) =  (1/pi .* (0.5.*p[1])./((x.-p[2]).^2 .+ (0.5.*p[1]).^2)).*p[3]  # Lorentcain function

    FWHM_guess = 0.5 #initial guess

    while true

        p0L = [FWHM_guess, maximum_indexL, maxIdata_L/2] # [FWHM, x0, scalling factor]
        p0R = [FWHM_guess, maximum_indexR, maxIdata_R/2] # [FWHM, x0, scalling factor]

        fitL = curve_fit(L, xL, yL, p0L) # fit L to x and y points, with p0 initial guesses
        fitR = curve_fit(L, xR, yR, p0R)

        paramL = fitL.param

        paramR = fitR.param

        sigmaL = stderror(fitL)
        sigmaR = stderror(fitR)

        error_percentL = sigmaL./abs.(paramL) .* 100
        error_percentR = sigmaR./abs.(paramR) .* 100

        if error_percentL[1] <10 && error_percentR[1]<10 && error_percentL[2]<10 && error_percentR[2]<10 || FWHM_guess >100
            break
        end
        FWHM_guess += 0.5
    end

    if FWHM_guess >100
        return [NaN,NaN], [NaN,NaN], [NaN,NaN], [0,0], [0,0], [0,0], [NaN,NaN], [NaN,NaN], [NaN,NaN]
    else
        x0L = paramL[2]

        x0L_error = error_percentL[2]
        x0R = paramR[2]
        x0R_error = error_percentR[2]

        maxI_L = L(paramL[2], paramL)
        maxI_R = L(paramR[2], paramR)

        FWHM_L = paramL[1]
        FWHM_L_error = error_percentL[1]
        FWHM_R = paramR[1]
        FWHM_R_error = error_percentR[1]
        # fig()
        # plt.scatter(xL,yL)
        # plt.plot(xL, L(xL, paramL))

        return [x0L, x0R], [FWHM_L, FWHM_R], [paramL[3], paramR[3]],  [x0L_error, x0R_error], [FWHM_L_error, FWHM_R_error], [error_percentL[3], error_percentR[3]],[maxI_L, maxI_R], [maxIdata_L,maxIdata_R], [maximum_indexL, maximum_indexR]
    end
end

"""
    mdc_cut_maximum(Data1D, k0)

Takes 1D `Array` from an MDC slice and the k=0 pixel position. 
Returns: pixel correspond to km, k1 and k2 as defined in [Kordyuk et al.](https://arxiv.org/pdf/cond-mat/0510421.pdf)

Works by spltting the 1D array in half and taking the maximum of `Data1D`.
A spline is then created of `Data1D` with half of its maximum value subtracted from each element.
The roots of this spline are then found to give k1 and k2. This function is best used for PL data sets where the spectral lines
are not well described by a Lorentzian function. If they are, one should intead use [`mdc_cut_fitting`](@ref).

# Example
```
julia> km, k1, k2 = mdc_cut_maximum(BK40[:,200, 5], 200)
([154, 250], [158.26734665230455, 245.55269763749718], [150.9186355501067, 251.4735060057741])
```
See also [`mdc_cut_fitting`](@ref), [`MDC_extraction`](@ref), [`pixel_to_k`](@ref).
"""
function mdc_cut_maximum(Data1D::Array, k0::Number)

    int(x) = floor(Int, x)  # defines a function that converts float to int
    size = length(Data1D)     # Length of the 1D array being analyzed
    HW = k0                 # Half of the length of the array

    xL = collect(1:int(HW))     # Defines x vector for the LHS 1:half way
    xR = collect(int(HW):size)  # Defines x vector for the RHS half way:end

    yL = Data1D[1:int(HW)]        # Cuts the data in half LHS
    yR = Data1D[int(HW):size]     # Cuts the data in half RHS

    maximum_indexL = argmax(yL) # Finds the index of the max intensity for LHS
    maxIdata_L = maximum(yL)    # Finds the maximum intensity value LHS
    maximum_indexR = argmax(yR) + k0 # Finds the index of the max intensity for RHS

    maxIdata_R = maximum(yR)    # Finds the maximum intensity value RHS


    poi_L = Spline1D(xL, yL.-(0.5*maxIdata_L))
    x0_L = roots(poi_L)
    x0_L = [x0_L[1], x0_L[end]]


    poi_R = Spline1D(xR, yR.-(0.5*maxIdata_R))
    x0_R = roots(poi_R)
    x0_R = [x0_R[1], x0_R[end]]


    return [maximum_indexL, maximum_indexR], [x0_L[2], x0_R[1]], [x0_L[1], x0_R[2]]
end


"""
    MDC_protocol_maximum(Data2D::Array)

Takes a 2D `Array` of PL data, sweeps through the energy pixels and applies the `mdc_cut_maximum` function to each 1D MDC. 

Returns `km`, `k1` and `k2` as defined in [Kordyuk et al.](https://arxiv.org/pdf/cond-mat/0510421.pdf) as `Arrays` with dimensions [length(E), 2] where column 1, i.e. [:,1], is the negative k-values and column 2, i.e. [:,2], is the positvie k-values.

# Exampple
```
julia> BK40, λ = correct_BK40();

julia> km, k1, k2 = MDC_protocol_maximum(BK40[:,:,5])
([NaN NaN; NaN 203.63971749241153; … ; 202.0 NaN; NaN NaN], [NaN NaN; NaN 216.6846069923789; … ; 188.0771116187305 NaN; NaN NaN], [NaN NaN; NaN 216.6846069923789; … ; 188.0771116187305 NaN; NaN NaN])
```
"""
function MDC_protocol_maximum(Data2D::Array)
    num_E_pixels = size(Data2D,2)
    maxI = maximum(Data2D)

    # Initialise Arrays to hold results
    km = Array{Float64, 2}(undef, num_E_pixels, 2)
    k1 =Array{Float64, 2}(undef, num_E_pixels, 2)
    k2 =Array{Float64, 2}(undef, num_E_pixels, 2)


    k, k0 = pixel_to_k(Data2D)

    # Loop through E values
    for i = 1:num_E_pixels
            MDC_i = Data2D[:,i]
        if maximum(MDC_i) < 0.1*maxI # id maximum of cut is less than 10% of max intensity, dont do mdc
            km[i,:], k1[i,:], k2[i,:] = [NaN,NaN], [NaN,NaN], [NaN,NaN]
        else
            km[i,:], k1[i,:], k2[i,:] = mdc_cut_maximum(Data2D[:,i], k0)

        end

    end

    
    km[:,2] = km[end:-1:1,2]
    k1[:,2] = k1[end:-1:1,2]
    k2[:,2] = k2[end:-1:1,2]
    return km, k1, k2

end

"""
    MDC_protocol_fitting(Data2D, cutoff)

Take 2D `Array` of ``corrected`` PL values and perform MDC analysis, extracting the quadruplete and corresponding errors. Returns DataFrames for for each parameter extracted from `mdc_cut_fitting`.

Optional `cutoff` parameter changes the thershold below which the MDC cuts are not applied. The default value is 10% of the global maximum.
# Example
```
julia> x0_df, FWHM_df, scale_factor_df, maximum_df= MDC_fitting_protocol(BK30[:,:,2])
(820×4 DataFrame
 Row │ Left x0  Right x0  Error Left  Error Right 
     │ Float64  Float64   Float64     Float64     
─────┼────────────────────────────────────────────
   1 │     NaN       NaN         0.0          0.0
   2 │     NaN       NaN         0.0          0.0
   3 │     NaN       NaN         0.0          0.0
   4 │     NaN       NaN         0.0          0.0
   5 │     NaN       NaN         0.0          0.0
   6 │     NaN       NaN         0.0          0.0
   7 │     NaN       NaN         0.0          0.0
  ⋮  │    ⋮        ⋮          ⋮            ⋮
 815 │     NaN       NaN         0.0          0.0
 816 │     NaN       NaN         0.0          0.0
 817 │     NaN       NaN         0.0          0.0
 818 │     NaN       NaN         0.0          0.0
 819 │     NaN       NaN         0.0          0.0
 820 │     NaN       NaN         0.0          0.0
                                  807 rows omitted, 820×4 DataFrame
 Row │ Left FWHM  Right FWHM  Error Left  Error Right 
     │ Float64    Float64     Float64     Float64     
─────┼────────────────────────────────────────────────
   1 │       NaN         NaN         0.0          0.0
   2 │       NaN         NaN         0.0          0.0
   3 │       NaN         NaN         0.0          0.0
   4 │       NaN         NaN         0.0          0.0
   5 │       NaN         NaN         0.0          0.0
   6 │       NaN         NaN         0.0          0.0
   7 │       NaN         NaN         0.0          0.0
  ⋮  │     ⋮          ⋮           ⋮            ⋮
 815 │       NaN         NaN         0.0          0.0
 816 │       NaN         NaN         0.0          0.0
 817 │       NaN         NaN         0.0          0.0
 818 │       NaN         NaN         0.0          0.0
 819 │       NaN         NaN         0.0          0.0
 820 │       NaN         NaN         0.0          0.0
                                      807 rows omitted, 820×4 DataFrame
 Row │ Left SF  Right SF  Error Left  Error Right 
     │ Float64  Float64   Float64     Float64     
─────┼────────────────────────────────────────────
   1 │     NaN       NaN         0.0          0.0
   2 │     NaN       NaN         0.0          0.0
   3 │     NaN       NaN         0.0          0.0
   4 │     NaN       NaN         0.0          0.0
   5 │     NaN       NaN         0.0          0.0
   6 │     NaN       NaN         0.0          0.0
   7 │     NaN       NaN         0.0          0.0
  ⋮  │    ⋮        ⋮          ⋮            ⋮
 815 │     NaN       NaN         0.0          0.0
 816 │     NaN       NaN         0.0          0.0
 817 │     NaN       NaN         0.0          0.0
 818 │     NaN       NaN         0.0          0.0
 819 │     NaN       NaN         0.0          0.0
 820 │     NaN       NaN         0.0          0.0
                                  807 rows omitted, 820×2 DataFrame
 Row │ Left Max I  Right Max I 
     │ Float64     Float64     
─────┼─────────────────────────
   1 │        NaN          NaN
   2 │        NaN          NaN
   3 │        NaN          NaN
   4 │        NaN          NaN
   5 │        NaN          NaN
   6 │        NaN          NaN
   7 │        NaN          NaN
  ⋮  │     ⋮            ⋮
 815 │        NaN          NaN
 816 │        NaN          NaN
 817 │        NaN          NaN
 818 │        NaN          NaN
 819 │        NaN          NaN
 820 │        NaN          NaN
               807 rows omitted)

```
See also [`mdc_cut_fitting`](@ref), [`correct_data`](@ref), [`pixel_to_k`](@ref).
"""
function MDC_protocol_fitting(Data2D, cutoff::Float64=0.1)
    num_E_pixels = size(Data2D,2)
    maxI = maximum(Data2D)


    # Initialise Arrays to hold results
    fit_x0 = Array{Float64, 2}(undef, num_E_pixels, 2)
    fit_FWHM =Array{Float64, 2}(undef, num_E_pixels, 2)
    fit_scale_factor = Array{Float64, 2}(undef, num_E_pixels, 2)
    fit_maximum = Array{Float64, 2}(undef, num_E_pixels, 2)
    data_maximum = Array{Float64, 2}(undef, num_E_pixels, 2)
    data_x0 = Array{Float64, 2}(undef, num_E_pixels, 2)
    error_FWHM = Array{Float64, 2}(undef, num_E_pixels, 2)
    error_x0 = Array{Float64, 2}(undef, num_E_pixels, 2)
    error_scale_factor = Array{Float64, 2}(undef, num_E_pixels, 2)

    k, k0 = pixel_to_k(Data2D)

    # Loop through E values
    for i = 1:num_E_pixels
            MDC_i = Data2D[:,i]
        if maximum(MDC_i) < cutoff*maxI # id maximum of cut is less than 10% of max intensity, dont do mdc
            fit_x0[i,:], fit_FWHM[i,:], fit_scale_factor[i,:], error_x0[i,:], error_FWHM[i,:], error_scale_factor[i,:], fit_maximum[i,:], data_maximum[i,:], data_x0[i,:] = [NaN,NaN], [NaN,NaN], [NaN,NaN], [0,0], [0,0], [0,0], [NaN,NaN], [NaN,NaN], [NaN,NaN]
        else
            fit_x0[i,:], fit_FWHM[i,:], fit_scale_factor[i,:], error_x0[i,:], error_FWHM[i,:], error_scale_factor[i,:], fit_maximum[i,:], data_maximum[i,:], data_x0[i,:] = mdc_cut_fitting(Data2D[:,i], k0)
        
        end

    end
    # Create dataframes for results, note RHS values are flipped so that htey are in ascending k
    x0_df = DataFrame("Left x0"=>fit_x0[:,1], "Right x0"=>fit_x0[end:-1:1,2], "Error Left" => error_x0[:,1], "Error Right" => error_x0[:,2])
    FWHM_df = DataFrame("Left FWHM"=>fit_FWHM[:,1], "Right FWHM" =>fit_FWHM[end:-1:1,2], "Error Left" => error_FWHM[:,1], "Error Right" => error_FWHM[:,2])
    scale_factor_df = DataFrame("Left SF"=>fit_scale_factor[:,1], "Right SF"=>fit_scale_factor[end:-1:1,2], "Error Left" => error_scale_factor[:,1], "Error Right" => error_scale_factor[end:-1:1,2])
    maximum_df = DataFrame("Left Max I"=>fit_maximum[:,1], "Right Max I"=>fit_maximum[:,2])
    return x0_df, FWHM_df, scale_factor_df, maximum_df
end 


"""
    MDC(Data2D, λ, plot_true_or_false)

Take 2D PL data `Array` and performs the mdc protocol with Lorentazian fitting. Returns `Vectors`: 
`kL`, `kR`, `EmL`, `EmR`, `HWHM_L`, `HWHM_R` which can then be used in the Imaginary self-energy protocols.

The Boolean argument `plot_true_or_false` will plot the results of `MDC_protocol_fitting` if `true`.

# Example
```
julia> kL, kR, EmL, EmR, HWHM_L, HWHM_R = MDC(BK30[:,:,8], λ30, true)
([NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN  …  NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN], [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN  …  NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN], [1.5926313998361783, 1.5926122982432607, 1.592593197717957, 1.5925740982602512, 1.592554999870126, 1.5925359025475647, 1.5925168062925505, 1.5924977111050667, 1.5924786169850962, 1.5924595239326225  …  1.5775073811557567, 1.5774891386045762, 1.5774708971065752, 1.5774526566617348, 1.577434417270037, 1.577416178931463, 1.5773979416459927, 1.577379705413609, 1.5773614702342915, 1.5773432361080215], [1.5773432361080215, 1.5773614702342915, 1.577379705413609, 1.5773979416459927, 1.577416178931463, 1.577434417270037, 1.5774526566617348, 1.5774708971065752, 1.5774891386045762, 1.5775073811557567  …  1.5924595239326225, 1.5924786169850962, 1.5924977111050667, 1.5925168062925505, 1.5925359025475647, 1.592554999870126, 1.5925740982602512, 1.592593197717957, 1.5926122982432607, 1.5926313998361783], [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN  …  NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN], [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN  …  NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN])
```
"""
function MDC(Data2D, λ, plot_true_or_false::Bool)
    x0_df, FWHM_df, scale_factor_df, maximum_df= MDC_protocol_fitting(Data2D) # MDC on Data2D
    x_pixels = collect(1:size(Data2D,1))
    E = λ_to_E(λ)
    k, x0 = pixel_to_k(Data2D)
    Left_x0 = x0_df."Left x0"
    Right_x0 = x0_df."Right x0"
    kL = k(Left_x0)
    kR = k(Right_x0)

    if plot_true_or_false
        fig()
        plt.contourf(k(x_pixels), E, Data2D', levels=300)
        plt.scatter(kL, E, s=0.5, color="magenta")
        plt.scatter(kR, E[end:-1:1], s=0.5, color="red")
        plt.xlabel(L"k~(\mathrm{m^{-1}})", labelpad=10, size=16, weight="bold", color="black")
        plt.ylabel(L"E~(\mathrm{eV})",labelpad=10, size=16, weight="bold", color="black")

        fig()
        plt.plot(x0_df."Error Left")
        plt.plot(x0_df."Error Right")
        plt.legend(["LHS", "RHS"])
        plt.xlabel("E Pixel", labelpad=10, size=16, weight="bold", color="black")
        plt.ylabel("Error (%)",labelpad=10, size=16, weight="bold", color="black")

        fig()
        plt.plot(abs.(FWHM_df."Left FWHM"))
        plt.plot(abs.(FWHM_df."Right FWHM"[end:-1:1]))
        plt.legend(["LHS", "RHS"])
        plt.xlabel("E Pixel", labelpad=10, size=16, weight="bold", color="black")
        plt.ylabel("k pixel length",labelpad=10, size=16, weight="bold", color="black")

        fig()
        plt.plot(FWHM_df."Error Left")
        plt.plot(FWHM_df."Error Right")
        plt.legend(["LHS", "RHS"])
        plt.xlabel("E Pixel", labelpad=10, size=16, weight="bold", color="black")
        plt.ylabel("Error (%)",labelpad=10, size=16, weight="bold", color="black")

        fig()
        plt.plot(maximum_df."Left Max I", E)
        plt.plot(maximum_df."Right Max I",E)
        plt.legend(["LHS", "RHS"])
        plt.ylabel(L"E~(\mathrm{eV})", labelpad=10, size=16, weight="bold", color="black")
        plt.xlabel(L"I~(\mathrm{a.u.})",labelpad=10, size=16, weight="bold", color="black")
    end

    EmL = λ_to_E(λ)
    EmR = EmL[end:-1:1]
    Right_FWHM = FWHM_df."Right FWHM"
    HWHM = abs.(Right_FWHM).*0.5
    HWHM_R = HWHM .* 4.22e4

    Left_FWHM = FWHM_df."Left FWHM"
    HWHM = abs.(Left_FWHM).*0.5
    HWHM_L = HWHM .* 4.22e4
    return kL, kR, EmL, EmR, HWHM_L, HWHM_R
end

"""
Plot_MDC_Cut(Data2D::Array, E_index::Int)

Takes 2D PL data `Array` and plots the MDC cut for a given energy pixel `E_index`.

# Example
```
julia> Plot_MDC_Cut(BK30[:,:,5], 200)
[210.67258734120222, 215.9525458388472]
PyObject Text(185.77777777777777, 0.5, 'I~(\\mathrm{a.u.})')
```
"""
function Plot_MDC_Cut(Data2D::Array, E_index::Int)
    k , k0 = pixel_to_k(Data2D)
    int(x) = floor(Int, x)  # defines a function that converts float to int

    Data1D = Data2D[: ,E_index]
    size = length(Data1D)     # Length of the 1D array being analyzed
    HW = k0                 # Half of the length of the array

    xL = collect(1:int(HW))     # Defines x vector for the LHS 1:half way
    xR = collect(int(HW):size)  # Defines x vector for the RHS half way:end

    yL = Data1D[1:int(HW)]        # Cuts the data in half LHS
    yR = Data1D[int(HW):size]     # Cuts the data in half RHS

    maximum_indexL = argmax(yL) # Finds the index of the max intensity for LHS
    maxIdata_L = maximum(yL)    # Finds the maximum intensity value LHS
    maximum_indexR = argmax(yR) + k0 # Finds the index of the max intensity for RHS

    maxIdata_R = maximum(yR)    # Finds the maximum intensity value RHS

    fig()
    plt.scatter(xL, yL, s=10, color = "blue")
    plt.scatter(maximum_indexL, maxIdata_L, color = "red", s=15)
    poi_L = Spline1D(xL, yL.-(0.5*maxIdata_L))
    x0_L = roots(poi_L)
    x0_L = [x0_L[1], x0_L[end]]
    y0_L =  Spline1D(xL, yL)(x0_L)
    plt.scatter(x0_L[1], y0_L[1], color= "green")
    plt.scatter(x0_L[2], y0_L[2], color= "orange")
    HWHM1 = collect(range(x0_L[1],maximum_indexL, length=50))
    HWHM2 = collect(range(maximum_indexL,x0_L[2], length=50))
    y_vL = fill(y0_L[1], 1:50)
    max_v = fill(maximum_indexL, 1:length(yL))
    plt.plot(HWHM1, y_vL, color="green", linestyle = "--")
    plt.plot(HWHM2, y_vL, color="orange", linestyle = "--")
    plt.plot(max_v, yL, color="red", linestyle = "-.")


    plt.scatter(xR, yR, s=10, color = "purple")
    plt.scatter(maximum_indexR, maxIdata_R, color = "red", s=15)
    poi_R = Spline1D(xR, yR.-(0.5*maxIdata_R))
    x0_R = roots(poi_R)
    println(x0_R)
    x0_R = [x0_R[1], x0_R[end]]
    y0_R =  Spline1D(xR, yR)(x0_R)
    plt.scatter(x0_R[1], y0_R[1], color= "orange")
    plt.scatter(x0_R[2], y0_R[2], color= "green")
    HWHM1_R = collect(range(x0_R[1],maximum_indexR, length=50))
    HWHM2_R = collect(range(maximum_indexR,x0_R[2], length=50))
    y_vR = fill(y0_R[1], 1:50)
    max_v_R = fill(maximum_indexR, 1:length(yR))
    plt.plot(HWHM1_R, y_vR, color="orange", linestyle = "--")
    plt.plot(HWHM2_R, y_vR, color="green", linestyle = "--")
    plt.plot(max_v_R, yR, color="red", linestyle = "-.")
    #plt.xlim([0,175])
    plt.ylim([-100, maximum(Data1D)*1.1])
    plt.xlabel(L"k~(\mathrm{pixels})")
    plt.ylabel(L"I~(\mathrm{a.u.})")
    #return [x0L, x0R], [FWHM_L, FWHM_R], [paramL[3], paramR[3]],  [0, 0], [0, 0], [0, 0],[maxI_L, maxI_R], [maxIdata_L,maxIdata_R], [maximum_indexL, maximum_indexR]

end