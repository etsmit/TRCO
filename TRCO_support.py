#support functions for TRCO_play.ipynb
#eventual python replacement for tcal_calc.pro



#function getFluxCalib, src, freq, coeffs=coeffs, specindex=specindex














def getApEff(elev, freq, coeffs=None):

    if coeffs == None:
        eff_long = 0.71

    eff_long=0.71
    rms = 230

    #theres some coefficient math here that is skipped in tcal_calc

    #where does 4.19... come from?
    arg = -(4.19169e-8 * rms * freq)**2

    return eff_long*np.exp( arg )

















#function getTau, freqs, coeffs=coeffs








#function ElevFromAirMass, A







#pro Ta2Flux, tau=tau, ap_eff=ap_eff






#function quickTatm, freqs, TempK






#function samplerToIDX, str











#function cvrtFlux2Ta, flux=flux, tau=tau, ap_eff=ap_eff










#function cvrtTa2Flux, Ta, tau=tau, ap_eff=ap_eff












#function blankMask, arr
