"""
    ImagExtraction(Data3D::Array, PumpPower, λ, [LBL, UBL, LBR, UBR])

Takes a 3D `Array` of PL data and extracts the _Imaginary_ component of the self-energy (Σ) for a pump power specified by `PumpPower` and
then calculates the Hilbert transform of the _Imaginary_ Σ to obtain the _Real_ compnent of Σ.

The optional arguments `LBL` and `UBL` set the lower and upper bounds of the negative MDC extracted momentum values respectively. 
The optional arguments `LBR` and `UBR` set the lower and upper bounds of the positive MDC extracted momentum values respectively. 

Returns `Vectors` `ImΣ_v8`, `ReΣ8`, `LBk`, `UBk`, which is the Imaginary self energy calculated from Eq.8 from [Kordyuk et al.](https://arxiv.org/pdf/cond-mat/0510421.pdf), 
the Hilbert transform of the extracted imaginary and the lower and upper momentum bounds.

# Example
```
julia> ImEx, ReKK, LBk, UBk =  Imag_extraction(BK30, 8, λ30)
Input lower bound for negative k:
-3.8e6
Input upper bound for negative k:
0.2e6
Happy with negative k crop? (y/n) or (skip) for mirrored
y
Input lower bound for positive k:
0.2e6
Input upper bound for positve k:
3.8e6
Happy with positive k crop? (y/n)
y
Happy with crop? (y/n)
y
([-0.00046139736483497, -0.0004594412177453874, -0.0004588813620354898, -0.0004594624369877758, -0.00046092908188474595, -0.0004630259360088984, -0.0004654976386427354, -0.0004680888290687556, -0.0004705441465694569, -0.0004726243636280423  …  0.00047649874693501113, 0.00047510432961822693, 0.0004737206502134636, 0.0004723380864263395, 0.0004709470159624756, 0.000469537816527493, 0.0004681008658270103, 0.00046662654156664845, 0.0004651052214520287, 0.00046352728318876935], [-0.0002152060070957816, -0.00021571726727061358, -0.00021854766917208025, -0.0002207725461813776, -0.00022328890431777722, -0.0002247923267341889, -0.0002261074133700307, -0.0002263255357303291, -0.00022630607107366771, -0.00022545373563568745  …  -0.00023934525575251508, -0.00023992731583510118, -0.00024014717269984319, -0.00024043636700954224, -0.0002404050210121233, -0.00024040936448398013, -0.0002401077097114925, -0.00023970922610572646, -0.0002389998061889069, -0.0002374854826986771], -3.8012317525328365e6, 3.80151318244036e6)
```
# Example 
```
julia> ImEx, ReKK, LBk, UBk = Imag_extraction(BK30, 8, λ30, -3.8e6, -0.2e6, 0.2e6, 3.8e6)
([-0.00046139736483497, -0.0004594412177453874, -0.0004588813620354898, -0.0004594624369877758, -0.00046092908188474595, -0.0004630259360088984, -0.0004654976386427354, -0.0004680888290687556, -0.0004705441465694569, -0.0004726243636280423  …  0.00047649874693501113, 0.00047510432961822693, 0.0004737206502134636, 0.0004723380864263395, 0.0004709470159624756, 0.000469537816527493, 0.0004681008658270103, 0.00046662654156664845, 0.0004651052214520287, 0.00046352728318876935], [-0.0002152060070957816, -0.00021571726727061358, -0.00021854766917208025, -0.0002207725461813776, -0.00022328890431777722, -0.0002247923267341889, -0.0002261074133700307, -0.0002263255357303291, -0.00022630607107366771, -0.00022545373563568745  …  -0.00023934525575251508, -0.00023992731583510118, -0.00024014717269984319, -0.00024043636700954224, -0.0002404050210121233, -0.00024040936448398013, -0.0002401077097114925, -0.00023970922610572646, -0.0002389998061889069, -0.0002374854826986771], -3.8012317525328365e6, 3.80151318244036e6)
```
"""
function ImagExtraction(Data3D::Array, PumpPower::Int, λ::Vector, LBL::Number=NaN, UBL::Number=NaN, LBR::Number=NaN, UBR::Number=NaN)
    
    Eb, kpoints = energy_band_fit(Data3D[:,:,1], λ, false)

    km_L, km_R, Em_L, Em_R, HWHM_L, HWHM_R = MDC(Data3D[:,:,PumpPower], λ, false)
    ## Remove NaN values

    remove_NaN!(km_L, Em_L)
    remove_NaN!(HWHM_L)
    remove_NaN!(km_R, Em_R)
    remove_NaN!(HWHM_R)

    # Crop data set
    if isnan(LBL) && isnan(UBL) && isnan(LBR) && isnan(UBR)
        LBL, UBL, LBR, UBR = kCropGUI([km_L;km_R], [Em_L;Em_R], Eb)
    end
    LBLi, UBLi = closest_index(km_L, LBL), closest_index(km_L, UBL)
    LBRi, UBRi = closest_index(km_R, LBR), closest_index(km_R, UBR)

    km_L, Em_L, HWHM_L = km_L[LBLi:UBLi], Em_L[LBLi:UBLi], HWHM_L[LBLi:UBLi]
    km_R, Em_R, HWHM_R = km_R[LBRi:UBRi], Em_R[LBRi:UBRi], HWHM_R[LBRi:UBRi]
    LBk = km_L[1]
    UBk = km_R[end]

    ###### Imaginary Calculation ########
    ImSL6 = Eb(km_L) .- Eb(km_L.+HWHM_L) # -3.8e6 m⁻¹ Lowerbound, upper bound removes bad fitting near k=0
    ImSR6 = Eb(km_R) .- Eb(km_R.-HWHM_R)  # lower bound removes bad fitting near k=0, upperbound 3.8e6 m⁻¹
    ImSL7 = -Eb(km_L) .+ Eb(km_L.-HWHM_L) # -3.8e6 m⁻¹ Lowerbound, upper bound removes bad fitting near k=0
    ImSR7 = -Eb(km_R) .+ Eb(km_R.+HWHM_R)  # lower bound removes bad fitting near k=0, upperbound 3.8e6 m⁻¹
    ImSL8 = (Eb(km_L.-HWHM_L) .- Eb(km_L.+HWHM_L))./2 # -3.8e6 m⁻¹ Lowerbound, upper bound removes bad fitting near k=0
    ImSR8 = (Eb(km_R.+HWHM_R) .- Eb(km_R.-HWHM_R))./2 # lower bound removes bad fitting near k=0, upperbound 3.8e6 m⁻¹

    ## Plot and save data ##
    fig() # Plot ImΣ(km)
    plt.plot(km_L, ImSL6.*1000,color="black")
    plt.plot(km_L, ImSL7.*1000,color="dimgrey")
    plt.plot(km_L, ImSL8.*1000,color="darkgray")
    plt.plot(km_R, ImSR6.*1000,color="black")
    plt.plot(km_R, ImSR7.*1000,color="dimgrey")
    plt.plot(km_R, ImSR8.*1000,color="darkgray")
    plt.legend([L"|\Sigma''|=\varepsilon(k_\mathrm{m}) -\varepsilon(k_\mathrm{1})", L"|\Sigma''|=\varepsilon(k_\mathrm{2}) -\varepsilon(k_\mathrm{m})", L"|\Sigma''|=[\varepsilon(k_\mathrm{2}) -\varepsilon(k_\mathrm{1})]/2"])
    plt.xlabel(L"k~(\mathrm{m^{-1}})") # the L prefix to the string enables latex notation
    plt.ylabel(L"E~(\mathrm{meV})")
    # concatenate LHS & RHS
    Im8 = vcat(-ImSL8,ImSR8)  
    Im6 = vcat(-ImSL6, ImSR6)
    Im7 = vcat(-ImSL7, ImSR7)
    km = vcat(km_L, km_R)
    #########################################
    ImSp8 = Spline1D(km,Im8) # 1D spline to connect LHS and RHS
    ImSp6 = Spline1D(km,Im6)
    ImSp7 = Spline1D(km,Im7)
    ### Plot of Spline
    fig()
    plt.plot(range(LBk,UBk,length=10000), ImSp6(range(LBk,UBk,length=10000)).*1000, color="black")
    plt.plot(range(LBk,UBk,length=10000), ImSp7(range(LBk,UBk,length=10000)).*1000, color="dimgrey")
    plt.plot(range(LBk,UBk,length=10000), ImSp8(range(LBk,UBk,length=10000)).*1000, color="darkgray")
    plt.legend([L"|\Sigma''|=\varepsilon(k_\mathrm{m}) -\varepsilon(k_\mathrm{1})", L"|\Sigma''|=\varepsilon(k_\mathrm{2}) -\varepsilon(k_\mathrm{m})", L"|\Sigma''|=[\varepsilon(k_\mathrm{2}) -\varepsilon(k_\mathrm{1})]/2"])
    plt.xlabel(L"k~(\mathrm{m^{-1}})") # the L prefix to the string enables latex notation
    plt.ylabel(L"E~(\mathrm{meV})")

    # Vector of values from spline 
    ImΣ_v6 =  ImSp6(range(LBk,UBk,length=10000))
    ImΣ_v7 =  ImSp7(range(LBk,UBk,length=10000))
    ImΣ_v8 =  ImSp8(range(LBk,UBk,length=10000))

    add_tails!(ImΣ_v6, 2000) # Adds 2000 values to either end of vector
    add_tails!(ImΣ_v7, 2000) # Adds 2000 values to either end of vector
    add_tails!(ImΣ_v8, 2000) # Adds 2000 values to either end of vector
    ImΣ6  = Array{Float64, 2}(undef,1, length(ImΣ_v6)) # Initialises row vector that is neccesarry form for the Hilbert.jl package
    ImΣ7  = Array{Float64, 2}(undef,1, length(ImΣ_v7)) # Initialises row vector that is neccesarry form for the Hilbert.jl package
    ImΣ8  = Array{Float64, 2}(undef,1, length(ImΣ_v8)) # Initialises row vector that is neccesarry form for the Hilbert.jl package

    ImΣ6[1,:] = ImΣ_v6 # Row vector takes values from column vector
    ImΣ7[1,:] = ImΣ_v7
    ImΣ8[1,:] = ImΣ_v8 # Row vector takes values from column vector

    ReΣ6 = -imag.(hilbert(ImΣ6))[1,:] # Hilbert transform of ImΣ to get ReΣ
    remove_tails!(ReΣ6, 2000)
    remove_tails!(ImΣ_v6,2000)

    ReΣ7 = -imag.(hilbert(ImΣ7))[1,:] # Hilbert transform of ImΣ to get ReΣ
    remove_tails!(ReΣ7, 2000)
    remove_tails!(ImΣ_v7,2000)

    ReΣ8 = -imag.(hilbert(ImΣ8))[1,:] # Hilbert transform of ImΣ to get ReΣ
    remove_tails!(ReΣ8, 2000)
    remove_tails!(ImΣ_v8,2000)


    kpoints = collect(range(LBk,UBk,length=10000))
    fig()
    plt.plot(kpoints, ReΣ6.*1000, color="darkgreen") # Plot ReΣ wihtout the tails.
    plt.plot(kpoints, ReΣ7.*1000, color="mediumseagreen") # Plot ReΣ wihtout the tails.
    plt.plot(kpoints, ReΣ8.*1000, color="springgreen") # Plot ReΣ wihtout the tails.
    plt.legend([L"KK[\mathrm{Eq}(6)]", L"KK[\mathrm{Eq}(7)]", L"KK[\mathrm{Eq}(8)]"])
    plt.xlabel(L"k~(\mathrm{m^{-1}})") # the L prefix to the string enables latex notation
    plt.ylabel(L"E~(\mathrm{meV})")
 
    return ImΣ_v8, ReΣ8, LBk, UBk
end


"""
    RealExtraction(Data3D ,PumpPower, λ, [LBk, UBk])

Takes a 3D `Array` of PL data and extracts the _Real_ component of the self-energy (Σ) for a pump power specified by `PumpPower` and
then calculates the Hilbert transform of the _Real_ Σ to obtain the _Imaginary_ compnent of Σ.

The optional arguments `LBk` and `UBk` will set the momentum boundaries and can be used to align the momentum axis with ImΣ and ReΣ 
obtained from [`ImagExtraction`](@ref). If no arguments passed the boundaries obtained from the EDC protocol in calculating the energy-band `Em` will be used.

# Example
```
julia> ReEx, ImKK, kpoints = RealExtraction(BK30, 7, λ30, -3.8e6, 3.8e6)
([2.9664571152387964e-5, 2.9790523568884453e-5, 2.991958608378198e-5, 3.0051750320447823e-5, 3.018700790069495e-5, 3.0325350447668598e-5, 3.0466769583847864e-5, 3.061125693304412e-5, 3.075880411729237e-5, 
3.090940276018195e-5  …  7.428860256530889e-5, 7.421338523982968e-5, 7.41368548877741e-5, 7.405845043484405e-5, 7.397766063887978e-5, 7.389402208768381e-5, 7.380707329729397e-5, 7.371635278574651e-5, 7.36213990681911e-5, 7.352175066244193e-5], [-9.468758379541248e-5, -9.5039629607532e-5, -9.531386639313429e-5, -9.559182989949589e-5, -9.584337325568394e-5, -9.609789825802171e-5, -9.633595080357406e-5, -9.65765037015618e-5, -9.68046854523619e-5, -9.703499732407897e-5  …  8.607905586474244e-5, 8.587189292395486e-5, 8.566794428284767e-5, 8.546092738835575e-5, 8.525634350529895e-5, 8.504403153449403e-5, 8.483309280865156e-5, 8.46044317211318e-5, 8.43755075515656e-5, 8.408264358960797e-5], [-3.8e6, -3.799239923992399e6, -3.7984798479847983e6, -3.797719771977198e6, -3.796959695969597e6, -3.796199619961996e6, -3.795439543954395e6, -3.794679467946795e6, -3.793919391939194e6, -3.793159315931593e6  …  3.793159315931593e6, 3.793919391939194e6, 3.794679467946795e6, 3.795439543954395e6, 3.796199619961996e6, 3.796959695969597e6, 3.797719771977198e6, 3.7984798479847983e6, 3.799239923992399e6, 3.8e6])
```
# Example
```
julia> ReEx, ImKK, kpoints = RealExtraction(BK30, 7, λ30)
([0.00010422385218955554, 6.371453207476563e-6, -4.65585557229975e-5, -5.456219781829752e-5, -5.236149671028478e-5, -4.8452370591034466e-5, -4.722816447189082e-5, -4.537699598095024e-5, -4.018479436718181e-5, -3.8168634289625913e-5  …  8.74556669829829e-5, 7.55151892055661e-5, 6.48101694484815e-5, 2.5642976922535254e-5, 2.8391609154576614e-5, 7.254151970736977e-5, 0.00031220609081006323, 0.0003674518442831065, 0.00045268712396273614, 0.0005444240512653131], [-0.0003310594241268489, -0.00030251690891883827, -0.0003644797499028658, -0.00040517010003138275, -0.000438188501598904, -0.0004555005937439928, -0.00047358316854034864, -0.0004896552316938953, -0.0005051446221573488, -0.0005187357350299954  …  -0.0008875739245731761, -0.0009280475980466711, -0.0009403797114826156, -0.0010139781913010586, -0.0010963812015614923, -0.0012696450494998042, -0.0013086108212050165, -0.0012125797282928044, -0.0012345959100698075, -0.001131267272622972], [-4.488655670649381e6, -4.446455670649381e6, -4.404255670649381e6, -4.362055670649381e6, -4.319855670649381e6, -4.277655670649381e6, -4.235455670649381e6, -4.1932556706493814e6, -4.1510556706493814e6, -4.1088556706493814e6  …  3.7403443293506186e6, 3.7825443293506186e6, 3.8247443293506186e6, 3.8669443293506186e6, 3.9091443293506186e6, 3.9513443293506186e6, 3.9935443293506186e6, 4.0357443293506186e6, 4.0779443293506186e6, 4.1201443293506186e6])
```
"""
function RealExtraction(Data3D::Array ,PumpPower::Int, λ, LBk::Number=NaN, UBk::Number=NaN)
    
    Eb, kpoints = energy_band_fit(Data3D[:,:,1], λ, false)
    Em, kpoints = energy_band_fit(Data3D[:,:, PumpPower], λ, false)

    if isnan(LBk) == false && isnan(UBk) == false
        kpoints = collect(range(LBk,UBk,length=10000))
    end

    Real_v = Em(kpoints) - Eb(kpoints) # Calculate ReΣ
    add_tails!(Real_v, 10000)

    real = Array{Float64,2}(undef, 1,length(Real_v))
    real[1,:] = Real_v
    
    Imag_KK = imag.(hilbert(real))[1,:] # hilbert/KKR of ReΣ
    remove_tails!(Imag_KK, 10000)
    remove_tails!(Real_v, 10000)

    fig()
    plt.title(L"\mathrm{Re}\Sigma_\mathrm{Ex}", weight = "bold")
    plt.plot(kpoints, Real_v.*1000, color= "black")
    plt.xlabel(L"k~(\mathrm{m^{-1}})")
    plt.ylabel(L"E~(\mathrm{meV})")

    fig()
    plt.title(L"\mathrm{Im}\Sigma_\mathrm{KK}", weight = "bold")
    plt.plot(kpoints, Imag_KK.*1000, color = "green")
    plt.xlabel(L"k~(\mathrm{m^{-1}})")
    plt.ylabel(L"E~(\mathrm{meV})")


    return Real_v, Imag_KK, kpoints
end

"""
    Σ_E(Data3D, λ, ReEx, ImEx, ReKK, ImKK, kpoints)

Takes the results of (`ImagExtraction`)[@ref] and (`RealExtraction`)[@ref] and plots all the different Σ components as a function of energy and momentum.

# Example
```
julia> ImEx, ReKK, LBk, UBk = Imag_extraction(BK30, 8, λ30, -3.8e6, -0.2e6, 0.2e6, 3.8e6);

julia> ReEx, ImKK, kpoints = RealExtraction(BK30, 8, λ30, LBk, UBk);

julia> Σ_E(BK30, λ30, ReEx, ImEx, ReKK, ImKK, kpoints)
PyObject <matplotlib.legend.Legend object at 0x000000006D660310>
```
"""
function Σ_E(Data3D::Array, λ::Vector, ReEx::Vector, ImEx::Vector, ReKK::Vector, ImKK::Vector, kpoints::Vector)
    Em, kbb = energy_band_fit(Data3D[:,:,1], λ, false)
    Re_Ex_S = Spline1D(kpoints, ReEx)
    Im_Ex_S = Spline1D(kpoints, ImEx)
    Re_KK_S = Spline1D(kpoints, ReKK)
    Im_KK_S = Spline1D(kpoints, ImKK)

    LHS = collect(range(kpoints[1],0,length=10000))
    RHS = collect(range(0,kpoints[end],length=10000))

    fig()
   # plt.title("Extracted Σ")
    plt.plot(kpoints, Re_Ex_S(kpoints).*1000, color = "blue")
    plt.plot(kpoints, Im_Ex_S(kpoints).*1000, color = "red")
    plt.xlabel(L"k~(\mathrm{m^{-1}})") # the L prefix to the string enables latex notation
    plt.ylabel(L"E~(\mathrm{meV})")
    plt.legend([L"ReΣ_\mathrm{Ex}", L"ImΣ_\mathrm{Ex}"])

    fig()
    #plt.title("KK Σ")
    plt.plot(kpoints, Re_KK_S(kpoints).*1000, color = "cyan")
    plt.plot(kpoints, Im_KK_S(kpoints).*1000, color = "magenta")
    plt.xlabel(L"k~(\mathrm{m^{-1}})") # the L prefix to the string enables latex notation
    plt.ylabel(L"E~(\mathrm{meV})")
    plt.legend([L"ReΣ_\mathrm{KK}", L"ImΣ_\mathrm{KK}"])


    fig()
    # LHS Real
    plt.plot(Re_Ex_S(LHS).*1000, Em(LHS), color = "blue")
    plt.plot(Re_KK_S(LHS).*1000, Em(LHS), color = "cyan")
    plt.xlabel(L"\Sigma~(\mathrm{meV})") # the L prefix to the string enables latex notation
    plt.ylabel(L"E~(\mathrm{eV})")
    plt.legend([L"ReΣ_\mathrm{EX}", L"ReΣ_\mathrm{KK}"])

    fig()
    # RHS Real
    plt.plot(Re_Ex_S(RHS).*1000, Em(RHS), color = "blue")
    plt.plot(Re_KK_S(RHS).*1000, Em(RHS), color = "cyan")
    plt.xlabel(L"\Sigma~(\mathrm{meV})") # the L prefix to the string enables latex notation
    plt.ylabel(L"E~(\mathrm{eV})")
    plt.legend([L"ReΣ_\mathrm{EX}", L"ReΣ_\mathrm{KK}"])

    fig()
    # LHS Imag
    plt.plot(Im_Ex_S(LHS).*1000, Em(LHS), color = "red")
    plt.plot(Im_KK_S(LHS).*1000, Em(LHS), color = "magenta")
    plt.xlabel(L"\Sigma~(\mathrm{meV})") # the L prefix to the string enables latex notation
    plt.ylabel(L"E~(\mathrm{eV})")
    plt.legend([L"ImΣ_\mathrm{EX}", L"ImΣ_\mathrm{KK}"])

    fig()
    # RHS Imag
    plt.plot(Im_Ex_S(RHS).*1000, Em(RHS), color = "red")
    plt.plot(Im_KK_S(RHS).*1000, Em(RHS), color = "magenta")
    plt.xlabel(L"\Sigma~(\mathrm{meV})") # the L prefix to the string enables latex notation
    plt.ylabel(L"E~(\mathrm{eV})")
    plt.legend([L"ImΣ_\mathrm{EX}", L"ImΣ_\mathrm{KK}"])

    fig()
    plt.plot(Re_Ex_S(LHS).*1000, Em(LHS), color = "blue")
    plt.plot(Im_KK_S(LHS).*1000, Em(LHS), color = "magenta")
    plt.xlabel(L"\Sigma~(\mathrm{meV})") # the L prefix to the string enables latex notation
    plt.ylabel(L"E~(\mathrm{eV})")
    plt.legend([L"ReΣ_\mathrm{EX}", L"ImΣ_\mathrm{KK}"])

    fig()
    # RHS Imag
    plt.plot(Re_Ex_S(RHS).*1000, Em(RHS), color = "blue")
    plt.plot(Im_KK_S(RHS).*1000, Em(RHS), color = "magenta")
    plt.xlabel(L"\Sigma~(\mathrm{meV})") # the L prefix to the string enables latex notation
    plt.ylabel(L"E~(\mathrm{eV})")
    plt.legend([L"ReΣ_\mathrm{EX}", L"ImΣ_\mathrm{KK}"])


    fig()
    # LHS Imag
    plt.plot(Im_Ex_S(LHS).*1000, Em(LHS), color = "red")
    plt.plot(Re_KK_S(LHS).*1000, Em(LHS), color = "cyan")
    plt.xlabel(L"\Sigma~(\mathrm{meV})") # the L prefix to the string enables latex notation
    plt.ylabel(L"E~(\mathrm{eV})")
    plt.legend([L"ImΣ_\mathrm{EX}", L"ReΣ_\mathrm{KK}"])

    fig()
    # RHS Imag
    plt.plot(Im_Ex_S(RHS).*1000, Em(RHS), color = "red")
    plt.plot(Re_KK_S(RHS).*1000, Em(RHS), color = "cyan")
    plt.xlabel(L"\Sigma~(\mathrm{meV})") # the L prefix to the string enables latex notation
    plt.ylabel(L"E~(\mathrm{eV})")
    plt.legend([L"ImΣ_\mathrm{EX}", L"ReΣ_\mathrm{KK}"])

end