function SpectralFunction(Re_S, Im_S, Eb, EDC_spline, kpoints, Figname)
    
    AL = Array{Float64}(undef, 1000, 1000)
  #  AR = Array{Float64}(undef, 10000, 10000)
    LHS = collect(range(kpoints[1],-4.1e5,length=1000))
   # RHS = collect(range(10,kpoints[end],length=10000))

    Re_vL = Re_S(LHS)
    Im_vL = Im_S(LHS)
   
    Eb_vL = Eb(LHS)
    E_vL = EDC_spline(LHS)
    # Re_vR = Re_S(RHS)
    # Im_vR = Im_S(RHS)
    # Eb_vR = Eb(RHS)
    # E_vR = EDC_spline(RHS)

    for i = 1:1000
        AL[:,i] = -1/pi .* (Im_vL[i])./((E_vL[i].-Eb_vL.-Re_vL[i]).^2 .+(Im_vL[i])^2)
        #AR[:,i] = 1/pi .* (Im_vR[i])./((E_vR[i].-Eb_vR.-Re_vR[i]).^2 .+(SpIm_vR[i])^2)
    end
    fig()
    plt.contourf(LHS, E_vL, AL')
    plt.xlabel(L"k~(\mathrm{m^{-1}})")
    plt.ylabel(L"E~(\mathrm{eV})") 
    p = plt.colorbar()
    p.set_label(L"A~(\mathrm{a.u.})")
    plt.savefig(Figname)
    plt.close()
end