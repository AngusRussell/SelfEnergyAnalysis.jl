function LowerPolaritonDispersion(k, p)
    # Define constants
# ħ = 6.5821e-16 #units of eV
ħ = 1.054571817e-34 # J⋅sd
m_e = 9.1093837e-31 # kg
e = 1.602176634e-19 # C
m_x = 0.2 * m_e
m_c =  3.75e-5 * m_e
E_x0 = 1.5935 # eV
#E_c0 = 1.58 # eV
Ω = 5.75e-3 # eV

E_x(k) = p[2] .+ (ħ^2 .* k.^2)./(2 * m_x *e)
E_LP(k, p) = 0.5 .* (E_x(k) .+ (p[1] .+(ħ^2 .* k.^2)./(2 * m_c *e))) .- (((E_x(k) .- (p[1] .+(ħ^2 .* k.^2)./(2 * m_c *e))).^2 .+(2 * Ω)^2).^0.5) .* 0.5

return E_LP(k, p)
end

function LowerPolaritonDispersion_fit()

    Eb = p1_bare_band()
    x = collect(range(-3.45e6,3.6e6, length=200))
    y = Eb(x)
    p0 = [1.5913, 1.5953]#, 3.75e-5]

    fit = curve_fit(LowerPolaritonDispersion, x, y, p0)

    param = fit.param


    sigma = stderror(fit)
    error = sigma ./ param .*100
    # fig()
    # plt.plot(x, LowerPolaritonDispersion(x,param))
    fig()
    plt.scatter(x,y,s=10, marker="o", facecolors="none", edgecolors="blue")
    plt.plot(x, LowerPolaritonDispersion(x, param), color="red")
    plt.xlabel(L"k_∥~(\mathrm{m^{-1}})", weight = "bold")
    plt.ylabel(L"E~(\mathrm{eV})",  weight = "bold")
    plt.legend([L"E_\mathrm{LP}~\mathrm{Fit}",L"P=0.5~\mathrm{mW}"], loc = "upper center")

    return param, error
end


function LPD(k)
    # Define constants

    ħ = 1.054571817e-34 # J⋅sd
    m_e = 9.1093837e-31 # kg
    e = 1.602176634e-19 # C
    m_x = 0.2 * m_e
    m_c =  3.75e-5 * m_e # from fit
    E_x0 = 1.5959 # eV
    E_c0 = 1.5926 # eV : from fit
    Ω = 5.75e-3 # eV
    p = [E_c0, m_c]

    E_x(k) = E_x0 .+ (ħ^2 .* k.^2)./(2 * m_x *e)
    E_LP(k, p) = 0.5 .* (E_x(k) .+ (p[1] .+(ħ^2 .* k.^2)./(2 * p[2]  *e))) .- (((E_x(k) .- (p[1] .+(ħ^2 .* k.^2)./(2 * p[2]  *e))).^2 .+(2 * Ω)^2).^0.5) .* 0.5

return E_LP(k, p)
end