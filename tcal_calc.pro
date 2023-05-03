;calculate TCal and uhh
;ETS 12/8/22
;ETS added calibrate section from TCal_Calibrate.pro



;==================================================
;
;functions here
;for procedures scroll down to line X
;
;==================================================


;+
; Returns a vector of fluxes for a given vector of frequencies for
; standard calibrators.
;
; <p>Can also be used to generate fluxes for non-standard calibrators.
;
; <p><b>Note: </b>If coeffs is given, then src and specindex are ignored.  If
; specindex is given, then src is ignored. src is only used if both
; specindex and coeffs is not used.
;
; <p>Recognized source names are: 3C48, 3C123, 3C147, 3C161, 3C218, 3C227,
; 3C249.1, VIRA, 3C286, 3C295, 3C309.1, 3C348, 3C353, NGC7027.
;
; <p>Calibrator coefficients are from the Ott et all or Peng list of
; calibrators.
;
; @param src {in}{required}{type=string} source name.  Must be one in
; the enclosed catalog.  This is ignored if coeffs and specindex are
; provided.
; @param freq {in}{required} list of frequencies for which fluxes are
; required.
; @keyword coeffs {in}{optional} For a non standard calibrator,
; polynomial coefficients to use to determine the flux of a source.
; log(S) = coeff[0] + coeff[1]*log(freq) + coeff[2]*log(freq)^2
; @keyword specindex {in}{optional} For a non standard calibrator,
; spectral index coefficients to use to determine the flux of a
; source.  S = specindex[0] * (freq/specindex[1])^(specindex[2])
; @returns a vector of fluxes at the given frequencies.
;
; @file_comments scalUtils is a collection of routines that return
; various quantities needed for calibration.  Users will need to look
; over and maybe modify getTau and getFluxCalib before using any of
; the other scal or getVctr routines.  Contact the contributor for
; additional details.
;
; <p><B>See also</B> the scal User's Guide found in the
; documentation for <a href="scal.html">scal.pro</a>
;
; <p><B>Contributed By: Ron Maddalena, NRAO-GB</B>
;
; @version $Id: scalUtils.pro,v 1.4 2009/12/02 18:57:35 bgarwood Exp $
;-
function getFluxCalib, src, freq, coeffs=coeffs, specindex=specindex

    srcNames = ['J0133-3629', '3C48',    'FornaxA',  '3C123',   'J0444-2809', '3C138',   'PictorA',  'TaurusA', '3C147',   '3C161',   '3C196',   '3C218',    '3C227',  '3C249.1',  '3C274',  '3C286',  '3C295',  '3C309.1',  '3C348',  '3C353',  '3C380',  'CygnusA',  '3C444',  'CassA',  'NGC7027']
    a =        [1.044,        1.3253,    2.2175,     1.8017,    0.971,        1.0088,    1.938,      2.9516,    1.4516,    1.25,      1.2872,    1.7795,     6.757,    2.537,      2.4466,   1.2481,   1.4701,   2.617,      1.8298,   1.8627,    1.232,   3.3498,     1.1064,   3.3584,    1.322   ]
    b =        [-0.6619,     -0.7553,   -0.6606,    -0.7884,   -0.8938,      -0.4981,   -0.747,     -0.2173,   -0.6961,    0.726,    -0.853,    -0.9176,    -2.801,   -0.565,     -0.8116,  -0.4507,  -0.7658,  -0.437,     -1.0247,  -0.6938,   -0.7909, -1.0022,    -1.0052,  -0.7518,   -0.134   ]
    c =        [-0.2252,     -0.1914,    0,         -0.1035,   -0.1176,      -0.1552,   -0.0739,    -0.0473,   -0.2007,   -0.2286,   -0.1534,   -0.0843,     0.2969,  -0.0404,    -0.0483,  -0.1798,  -0.278,   -0.0373,    -0.0951,  -0.0998,    0.0947, -0.2246,    -0.075,   -0.0347,    0       ]
    d =        [0,            0.0498,    0,         -0.0248,    0,           -0.0102,    0,         -0.0674,    0.064,     0,        -0.02,     -0.0139,     0,        0,          0,        0.0357,  -0.0347,   0,          0,       -0.0732,    0.0976,  0.0227,    -0.0767,  -0.0705,    0       ]
    e =        [0,            0,         0,          0.009,     0,            0.0223,    0,          0,        -0.0464,    0,         0.0201,    0.0295,     0,        0,          0,        0,        0.0399,   0,          0,        0,        -0.1794,  0.0425,     0,        0,         0       ]
    f =        [0,            0,         0,          0,         0,            0,         0,          0,         0.0289,    0,         0,         0,          0,        0,          0,        0,        0,        0,          0,        0,        -0.1566,  0,          0,        0,         0       ]
    MhzGhz =   [1000,         1000,      1000,       1000,      1000,         1000,      1000,       1000,      1000,      1,         1000,      1000,       1,        1,          1000,     1000,     1000,     1,          1000,     1000,      1000,    1000,       1000,     1000,      1       ]
    f1 =       [0.2,          0.05,      0.2,        0.05,      0.2,          0.2,       0.2,        0.05,      0.05,      1.4,       0.05,      0.05,       1.4,    1.4,          0.05,     0.05,     0.05,     1.4,        0.2,      0.2,       0.05,    0.05,       0.2,      0.2,      10.5     ]
    f2 =       [4,           50,         0.5,       50,         2,           50,         4,          4,        50,        10.6,      50,        12,          4.8,    4.8,         12,       50,       50,       32,         12,        0.4,      50,      12,         12,        4,        43.2     ]

    ; f1, f2 = Frequency limits in GHz
    ; MHzGhz = 1000 if the b,c,d values are for freqs in GHz, = 1 if in MHz
    ; a,b,c = coefficients for logS = a + b*log(freq) + c*log(freq)^2 + ... f*log(freq)^5

    ; srcNames = ['3C48',  '3C123', '3C147', '3C161', '3C218', '3C227', '3C249.1',  'VIRGOA', '3C286', '3C295', '3C309.1', '3C348', '3C353', 'NGC7027']
    ; a =        [2.465,     2.525,   2.806,   1.250,   4.729,   6.757,     2.537,   4.484,   0.956,   1.490,     2.617,   3.852,   3.148,     -9.625592]
    ; b =        [-0.004,    0.246,  -0.140,   0.726,  -1.025,  -2.801,    -0.565,  -0.603,   0.584,   0.756,    -0.437,  -0.361,  -0.157,    5.002555]
    ; c =        [-0.1251, -0.1638, -0.1031, -0.2286,  0.0130,  0.2969,   -0.0404, -0.0280, -0.1644, -0.2545,   -0.0373, -0.1053, -0.0911,     -0.6042999  ]
    ; a =        [2.69116,     2.525,   2.806,   1.250,   4.729,   6.757,     2.537,   4.484,   0.956,   1.490,     2.617,   3.852,   3.148,     -9.625592]
    ; b =        [-0.124817,    0.246,  -0.140,   0.726,  -1.025,  -2.801,    -0.565,  -0.603,   0.584,   0.756,    -0.437,  -0.361,  -0.157,    5.002555]
    ; c =        [-0.109415, -0.1638, -0.1031, -0.2286,  0.0130,  0.2969,   -0.0404, -0.0280, -0.1644, -0.2545,   -0.0373, -0.1053, -0.0911,     -0.6042999  ]

    if (n_elements(coeffs) eq 0 and n_elements(specindex) eq 0) then begin
        ; Use the table coeffs
        n = n_elements(freq)
        i = where(srcNames eq strupcase(strtrim(src,2)))
        if i ge 0 then begin
            if min(freq)/1000. LT f1(i) OR max(freq)/1000. GT f2(i) then begin
                print, 'Frequency range is beyond that for ', src, '...  Calculating flux anyway'
            endif

            cvrt=replicate(MHzGhz(i), n)
            aa = replicate(a(i), n)
            bb = replicate(b(i), n)
            cc = replicate(c(i), n)
            dd = replicate(d(i), n)
            ee = replicate(e(i), n)
            ff = replicate(f(i), n)

            logfreq = alog10(freq/cvrt)
            return, 10^(aa + bb*logfreq + cc*logfreq*logfreq + dd*logfreq^3 + ee*logfreq^4 + ff*logfreq^5)
        endif
        if (n_elements(coeffs) eq 0) then begin
            print, src, " is not in the calibration table ... Using S = 1.0"
            return, replicate(1.0, n)
        endif
    endif

    if (n_elements(coeffs) ne 0) then begin
        ; Use the user-supplied coeffs
        logfreq = alog10(freq)
        flux = 0
        for i=0, n_elements(coeffs)-1 do begin
            flux = flux + coeffs[i]*logfreq^i
        end
        return, 10^flux
    endif

    if (n_elements(specindex) ne 0) then begin
        ; Use the user-supplied spectral index, flux, and frequency
        return, specindex[0] * (freq/specindex[1])^(specindex[2])
    endif

end

;+
; Returns a vector of aperture efficiencies.
;
; <p><B>Note: </B>You cannot supply a surface rms without also
; supplying a long wavelength efficiency.
;
; @param freq {in}{required}{type=vector} list of frequencies in MHz for which an opacity is needed
; @param elev {in}{required}{type=float} elevation in degrees of the observation
; @keyword coeffs {in}{optional}{type=vector} coeffs[0] = the long wavelength efficiency (Default : 0.72)
; coeffs[1] = Surface rms in microns (Default : 220 microns)
; @returns vector of aperture efficiences at freq
;
; @examples
; <pre>
; a = getApEff(45.0, freqs) ; returns the PTCS  model for ap_eff
; a = getApEff(45.0, freqs, coeffs=[0.69]) ; uses 184 microns but a long-wavelength eff of 69%
; a = getApEff(45.0, freqs, coeffs=[0.73, 250]) ; uses 250 microns and a long-wavelength eff of 73%
; </pre>
;
;-
function getApEff, elev, freq, coeffs=coeffs

    ; Default is 71% efficiency at long wavelengths and the current best model for
    ; the surface, which has a minimum surface RMS of 220 microns near the rigging angle
    ; and then deteriates away from that angle.

    if (n_elements(coeffs) eq 0) then begin
        eff_long=0.71
    end

    eff_long=0.71
    ; rms = 415.364254 - 7.10557065*elev + 0.0656154394*elev*elev
    rms = 230.

    if (n_elements(coeffs) ge 1) then eff_long = coeffs[0]
    if (n_elements(coeffs) ge 2) then begin
        rms = 0
        for i=1, n_elements(coeffs)-1 do begin
            rms = rms + coeffs[i]*elev^(i-1)
        end
    endif

    ; print, 'D', eff_long, rms, freq, eff_long*exp(-(4.19169e-8*rms*freq)^2)
    return, eff_long*exp(-(4.19169e-8*rms*freq)^2)
end

;+
; Returns a vector of atmospheric zenith opacities.
; @param freqs {in}{required} list of frequencies in MHz for which an opacity is needed
; @keyword coeffs {in}{optional} polynomial coefficients tau = coeff[0] + coeff[1]*freq + coeff[2]*freq^2 + ....
; @returns vector of atmospheric zenith opacities at freqs
;
; @examples
; <pre>
; a = getTau(freqs, coeffs=[0.01]) ; an opacity that is constant with freq
; a = getTau(freqs, coeffs=[0.0234, 0.4567, 0.0045])
; </pre>
;
;-
function getTau, freqs, coeffs=coeffs

    n = n_elements(freqs)
    if (n_elements(coeffs) eq 0) then return, replicate(0.0, n)

    tau = replicate(0.0, n)
    for i=0, n_elements(coeffs)-1 do begin
        tau = tau + coeffs[i]*(freqs/1000.)^i
    end
    ; print, 'C', tau
    return, tau
end

;+
; Estimate the airmass as a function of elevation in degrees.
;
; @param elev{in}{required}{float} elevation in degrees.
; @returns airmass
;-
function AirMass, elev
    if (elev LT 28) then begin
        ; print, 'A', -0.023437  + 1.0140 / sin( (!pi/180.)*(elev + 5.1774 / (elev + 3.3543) ) )
        return, -0.023437  + 1.0140 / sin( (!pi/180.)*(elev + 5.1774 / (elev + 3.3543) ) )
    endif else begin
        ; print, 'B', 1./sin(!pi*elev/180.)
        return, 1./sin(!pi*elev/180.)
    endelse
end

function ElevFromAirMass, A
    ; Reverse calculation
    if (A GT 2.12488) then begin
       B=-(3.3543 + (180/!pi)*asin(1.0140/(A + 0.023437)))
       return, ((-B+SQRT(B*B-4*5.1774))/2-3.3543)
    endif else begin
        ; print, 'B', 1./sin(!pi*elev/180.)
        return, (180/!pi)*asin(1./A)
    endelse
end

;+
; Converts data in buffer 0 from Ta to Jy.
;
; @keyword tau {in}{required} atmospheric zenith opacity encoded as a
; vector.  See the documentation for getTau for the format of the
; vector.
; @keyword ap_eff {in}{required} aperture efficiency encoded as a
; vector.  See the documentation for getApEff for the format of the
; vector.
;-
pro Ta2Flux, tau=tau, ap_eff=ap_eff

    elev=!g.s[0].elevation
    num_chan = n_elements(getdata(0))
    freqs = chantofreq(!g.s[0],seq(0,num_chan-1))/1.e6

    if n_elements(tau) eq 0 then begin
        tau = getForecastedTau(dateToMJD(!g.s[0].timestamp),!g.s[0].center_frequency/1.e6)
    endif else begin
        if n_elements(tau) eq 1 then begin
            if tau eq -1 then tau = getForecastedTau(dateToMJD(!g.s[0].timestamp),!g.s[0].center_frequency/1.e6)
        endif
    endelse
    tauVctr = getTau(freqs, coeffs=tau)

    effVctr = getApEff(elev, freqs, coeffs=ap_eff)
    setdata, getdata(0) * exp(tauVctr*AirMass(elev))/(2.8 * effVctr )
    !g.s[0].units = "Jy"

end

;+
; Calculates an estimate to Tatm from ground air temperature and
; frequencies.
;
; <p>Only appropriate for freqs < 50 GHz.
;
; <p>The results of Maddalena & Johnson (2005, Bulletin of the American Astronomical
; Society, Vol. 37, p.1438).   The rms uncertainty in my model is 3.5 K
;
; @param freqs {in}{required}{type=float} list of frequencies in MHz
; for which an opacity is needed
; @param TempK {in}{requried}{type=float} ground temperature in K
;
;-
function quickTatm, freqs, TempK
    f = freqs/1000.
    A = 259.691860 - 1.66599001*f + 0.226962192*f^2 - 0.0100909636*f^3 + 0.00018402955*f^4 - 0.00000119516*f^5
    B = 0.42557717 + 0.03393248*f + 0.000257983*f^2 - 0.0000653903*f^3 + 0.00000157104*f^4 - 0.00000001182*f^5
    return, A + B*(TempK-273.15)
end

function dateToMjd, dateString
    ; Converts a date string to an MJD
    year  = strmid(dateString,0,4)
    month = strmid(dateString,5,2)
    day   = strmid(dateString,8,2)
    hour  = strmid(dateString,11,2)
    minute= strmid(dateString,14,2)
    second= strmid(dateString,17)

    ;now convert from julian day to mjd
    return, julday(month, day, year, hour, minute, second) - 2400000.5d0
end

function getForecastedTau, mjd, freqMHz
    ; Gets the forecasted opeacity from the weather data base for the list of MJD's and frequqncies freqMHz
    if freqMHz LE 2000 then return, 0.009
    spawn, '/users/rmaddale/bin/getForecastValues -timeList ' + string(mjd) + ' -freqList ' + string( freqMHz/1000., /print) + ' | grep = | cut -f3 -d" "', tau
    return, tau
end

function samplerToIDX, str

    ;
    aa=[   "A1_0","A1_1","A1_2","A1_3","A1_4","A1_5","A1_6","A1_7","A2_0","A2_1","A2_2","A2_3","A2_4","A2_5","A2_6","A2_7"]
    aa=[aa,"B1_0","B1_1","B1_2","B1_3","B1_4","B1_5","B1_6","B1_7","B2_0","B2_1","B2_2","B2_3","B2_4","B2_5","B2_6","B2_7"]
    aa=[aa,"C1_0","C1_1","C1_2","C1_3","C1_4","C1_5","C1_6","C1_7","C2_0","C2_1","C2_2","C2_3","C2_4","C2_5","C2_6","C2_7"]
    aa=[aa,"D1_0","D1_1","D1_2","D1_3","D1_4","D1_5","D1_6","D1_7","D2_0","D2_1","D2_2","D2_3","D2_4","D2_5","D2_6","D2_7"]
    aa=[aa,"E1_0","E1_1","E1_2","E1_3","E1_4","E1_5","E1_6","E1_7","E2_0","E2_1","E2_2","E2_3","E2_4","E2_5","E2_6","E2_7"]
    aa=[aa,"F1_0","F1_1","F1_2","F1_3","F1_4","F1_5","F1_6","F1_7","F2_0","F2_1","F2_2","F2_3","F2_4","F2_5","F2_6","F2_7"]
    aa=[aa,"G1_0","G1_1","G1_2","G1_3","G1_4","G1_5","G1_6","G1_7","G2_0","G2_1","G2_2","G2_3","G2_4","G2_5","G2_6","G2_7"]
    aa=[aa,"H1_0","H1_1","H1_2","H1_3","H1_4","H1_5","H1_6","H1_7","H2_0","H2_1","H2_2","H2_3","H2_4","H2_5","H2_6","H2_7"]

    idx = where(aa eq str)

    if idx LT 0 then begin
        aa=[   "A1","A2","A3","A4","A5","A5","A7","A8","A9","A10","A11","A12","A13","A14","A15","A16"]
        aa=[aa,"B1","B2","B3","B4","B5","B5","B7","B8","B9","B10","B11","B12","B13","B14","B15","B16"]
        aa=[aa,"C1","C2","C3","C4","C5","C5","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16"]
        aa=[aa,"D1","D2","D3","D4","D5","D5","D7","D8","D9","D10","D11","D12","D13","D14","D15","D16"]

        idx = where(aa eq str)
    endif

    return, idx
end

function cvrtFlux2Ta, flux=flux, tau=tau, ap_eff=ap_eff

    ; Converts a vector of fluxes to a vector of Ta, correcting for atmosphere Tau
    ; and aperture effeciency
    ;
    ; Uses the data in DC0 to get the frequency vector.  If flux not provided,
    ; will try to use the source name and catalog of fluxes.  If Tau not given
    ; uses weather database.  If ap_eff not provided, will use the above model

    elev=!g.s[0].elevation
    num_chan = n_elements(getdata(0))
    freqs = chantofreq(!g.s[0],seq(0,num_chan-1))/1.e6

    if n_elements(flux) eq 0 then flux = getfluxcalib(!g.s[0].source, freqs)

    if n_elements(tau) eq 0 then begin
        tau = getForecastedTau(dateToMJD(!g.s[0].timestamp),!g.s[0].center_frequency/1.e6)
    endif else begin
        if n_elements(tau) eq 1 then begin
            if tau eq -1 then tau = getForecastedTau(dateToMJD(!g.s[0].timestamp),!g.s[0].center_frequency/1.e6)
        endif
    endelse
    tauVctr = getTau(freqs, coeffs=tau)

    effVctr = getApEff(elev, freqs, coeffs=ap_eff)

    print, "Tau, ApEff, Center Freq = ", mean(tauVctr), mean(effVctr), mean(freqs)
    return, flux*exp(-tauVctr*AirMass(elev))*(2.8 * effVctr)

end

function cvrtTa2Flux, Ta, tau=tau, ap_eff=ap_eff

    ; Converts a vector of Ta to a vector of fluxes, correcting for atmosphere Tau
    ; and aperture effeciency
    ;
    ; Uses the data in DC0 to get the frequency vector.  If Tau not given
    ; uses weather database.  If ap_eff not provided, will use the above model

    elev=!g.s[0].elevation
    num_chan = n_elements(getdata(0))
    freqs = chantofreq(!g.s[0],seq(0,num_chan-1))/1.e6

    if n_elements(tau) eq 0 then begin
        tau = getForecastedTau(dateToMJD(!g.s[0].timestamp),!g.s[0].center_frequency/1.e6)
    endif else begin
        if n_elements(tau) eq 1 then begin
            if tau eq -1 then tau = getForecastedTau(dateToMJD(!g.s[0].timestamp),!g.s[0].center_frequency/1.e6)
        endif
    endelse

    tauVctr = getTau(freqs, coeffs=tau)

    effVctr = getApEff(elev, freqs, coeffs=ap_eff)
    return, Ta*exp(tauVctr*AirMass(elev))/(2.8 * effVctr)

end




;=======================================================================
;
;tcal procedures here
;
;=======================================================================




pro tcal_calc,onscan,offscan,ifnum=ifnum,plnum=plnum,fdnum=fdnum,fileout=fileout


write=1
if (n_elements(fdnum) eq 0) then fdnum = 0
if (n_elements(fileout) eq 0) then write = 0
;set frequencies i guess
gettp,onscan,ifnum=ifnum,plnum=plnum,fdnum=fdnum

;get signal scan, cal on
onsource_calon_chunk=getchunk(scan=onscan,ifnum=ifnum,plnum=plnum,fdnum=fdnum,sig='T',cal='T',count=nons_con_chunk)
for ons_con_int=0,nons_con_chunk-1 do accum,dc=onsource_calon_chunk[ons_con_int]
data_free,onsource_calon_chunk
ave,/quiet
onsource_calon_data=getdata()

;get signal scan, cal off
onsource_caloff_chunk=getchunk(scan=onscan,ifnum=ifnum,plnum=plnum,fdnum=fdnum,sig='T',cal='F',count=nons_coff_chunk)
for ons_coff_int=0,nons_coff_chunk-1 do accum,dc=onsource_caloff_chunk[ons_coff_int]
data_free,onsource_caloff_chunk
ave,/quiet
onsource_caloff_data=getdata()

;get ref scan, cal on
offsource_calon_chunk=getchunk(scan=offscan,ifnum=ifnum,plnum=plnum,fdnum=fdnum,sig='T',cal='T',count=noffs_con_chunk)
for offs_con_int=0,noffs_con_chunk-1 do accum,dc=offsource_calon_chunk[offs_con_int]
data_free,offsource_calon_chunk
ave,/quiet
offsource_calon_data=getdata()

;get ref scan, cal off
offsource_caloff_chunk=getchunk(scan=offscan,ifnum=ifnum,plnum=plnum,fdnum=fdnum,sig='T',cal='F',count=noffs_coff_chunk)
for offs_coff_int=0,noffs_coff_chunk-1 do accum,dc=offsource_caloff_chunk[offs_coff_int]
data_free,offsource_caloff_chunk
ave,/quiet
offsource_caloff_data=getdata()


num_chan = n_elements(offsource_caloff_data)
freqs = chantofreq(!g.s[0],seq(0,num_chan-1))/1.e6
freqs = freqs[sort(freqs)]
;print,freqs
fluxS_Vctr = getFluxCalib(!g.s[0].source,freqs)
ApEff = getApEff(!g.s[0].elevation,mean(freqs))

; To Convert flux to temperature scale, (TA', i.e. TA * e(tau*Airmass))
; Run the procedure below, note that it is best to do this after stitching together the IF chunks, etc
; From Memo #302 ap_eff=0.352 TA'/Sv
;S_to_T = ApEff/0.352
;fluxT_Vctr = fluxS_Vctr*S_to_T
;multiply Temperature measurements by fluxT_Vctr

calcounts_onsource=onsource_calon_data-onsource_caloff_data
calcounts_offsource=offsource_calon_data-offsource_caloff_data
sourcecounts_calon=onsource_calon_data-offsource_calon_data
sourcecounts_caloff=onsource_caloff_data-offsource_caloff_data

calonsource=calcounts_onsource/sourcecounts_calon
caloffsource=calcounts_offsource/sourcecounts_caloff
meancalonsource=mean(calonsource[0.1*num_chan:0.9*num_chan],/nan)
meancaloffsource=mean(caloffsource[0.1*num_chan:0.9*num_chan],/nan)
sdcalonsource=stddev(calonsource[0.1*num_chan:0.9*num_chan],/nan)
sdcaloffsource=stddev(caloffsource[0.1*num_chan:0.9*num_chan],/nan)
diffmean=abs(((meancalonsource-meancaloffsource)/meancaloffsource)*100.0)
diffsd=abs(((sdcalonsource-sdcaloffsource)/sdcaloffsource)*100.0)
Tcal=(calonsource+caloffsource)/2.0




if (diffmean gt 1.0) or (diffsd gt 5.0) then begin
    print,'difference in means of on and off scans - ',diffmean
    print,'difference in standard deviations of on and off scans - ',diffsd
    print,'differences between on-source and off-source cal measurements are large, this is probably not a good dataset for determining Tcals'
endif

Tsys_calon=offsource_calon_data/sourcecounts_calon
Tsys_caloff=offsource_caloff_data/sourcecounts_caloff

flist=strjoin(string(freqs/1000.0,format='(f7.3)'),' ')


;spawn,'~rmaddale/bin/getForecastValues -type AtmTsys -freqList '+flist+' -elev '+string(!g.s[0].elevation)+' -timeList '+string(!g.s[0].mjd),result
;Tsky=strmid(result,5,10,/reverse_offset)
;Tbg=2.73
;Tspill=3.0
;Trx=Tsys_caloff-Tsky-Tbg-Tspill

;uncalibrated vector Tcal and Tsys
Tcal=reverse(Tcal)
;Trx=reverse(Trx)
Tsys_caloff=reverse(Tsys_caloff)

;setdata, Tcal & !g.s[0].units = 'Tcal (K)'
;copy,0,1
;setdata, Tsys_caloff & !g.s[0].units = 'Tsys (K)'
;copy,0,2
;print,'Tcal in buffer 1;  Tsys in buffer 2'
;copy, 1,0
;show


;TCal_Calibrate section

scannos= [onscan,offscan]

ApEffArr=[-1]
AveElevArr=[-1]
gettp,scannos[0]
fluxS_Vctr = getFluxCalib(!g.s[0].source,freqs);*1000.0)
for ijk=0,n_elements(scannos)-1 do begin
gettp,scannos[ijk],/quiet
ApEff = getApEff(!g.s[0].elevation,mean(freqs))
ApEffArr=[ApEffArr,ApEff]
AveElevArr=[AveElevArr,!g.s[0].elevation]
endfor
unfreeze
ApEff = mean(ApEffArr[1:*])
AveEl=mean(AveElevArr[1:*])
print,'ApEff:',ApEff
print,'AveEl:',AveEl


; To Convert flux to temperature scale, (TA', i.e. TA * e(tau*Airmass))
; Run the procedure below, note that it is best to do this after stitching together the IF chunks, etc
; From Memo #302 ap_eff=0.352 TA'/Sv
S_to_T = ApEff/0.352
fluxT_Vctr = fluxS_Vctr*S_to_T

; Need to adjust this calibration factor for tau and airmass
; Tau is based on the last scanno provided
; Airmass is based on an average of the scan elevations
; Atmospheric Tsky is also based on an average of the scan elevations

TCal_Cal=TCal*fluxT_Vctr
TSys_Cal=TSys_caloff*fluxT_Vctr
flist=strjoin(string(freqs/1000.0,format='(f8.3)'),' ')


flist=strjoin(string((findgen(round(max(freqs)-min(freqs)))+min(freqs)+1)/1000.0,format='(f8.3)'))

;print,flist
print,string(!g.s[0].mjd)
spawn,'~rmaddale/bin/getForecastValues -type Opacity -freqList '+flist+' -timeList '+string(!g.s[0].mjd),result
Tau=strmid(result,5,10,/reverse_offset)
;print,'Tau:',Tau
AM=Airmass(AveEl)
print,'AM:',AM


;check 2GHz and below taus are replaced by a flat value evaluated at 2GHz
if freqs[0] lt 2000 then begin
	print,'replacing some taus'
	;if partial freq range under 2
	tau[0:(where(tau gt -9,count))[0]] = tau[(where(tau gt -9,count))[0]]
endif
if freqs[n_elements(freqs)-1] lt 2000 then begin
	;if all freq range under 2
	print,'replacing all taus'
	spawn,'~rmaddale/bin/getForecastValues -type Opacity -freqList 2 -timeList '+string(!g.s[0].mjd),result
	tau[0:*] = strmid(result,5,10,/reverse_offset)
endif

;print,Tau

print, freqs[0]
print, freqs[n_elements(freqs)-1]
print, freqs[1] - freqs[0]

TCal_Cal=TCal_Cal/exp(Tau*AM)
TSys_Cal=TSys_Cal/exp(Tau*AM)

spawn,'~rmaddale/bin/getForecastValues -type AtmTsys -freqList '+flist+' -elev '+string(AveEl)+' -timeList '+string(!g.s[0].mjd),result
Tsky=strmid(result,5,10,/reverse_offset)
Tbg=2.73
Tspill=3.0
Trx=TSys_Cal-Tsky-Tbg-Tspill

maxnu=max(freqs)
minnu=min(freqs)
dnu=abs(!g.s[0].frequency_interval/1d9)
nnu=ceil((maxnu-minnu)/dnu)+1
Freqout=(findgen(nnu)*dnu)+minnu
Trxint=interpol(Trx,Freqs,Freqout)
TCalint=interpol(TCal_Cal,Freqs,Freqout)
TSysint=interpol(TSys_Cal,Freqs,Freqout)


setdata, TCal_Cal & !g.s[0].units = 'Tcal (K)'
print,!g.s[0].reference_frequency
print,!g.s[0].reference_channel
print,(freqs[n_elements(freqs) - 1 - !g.s[0].reference_channel])*1e6

!g.s[0].reference_frequency = (freqs[n_elements(freqs) - 1 - !g.s[0].reference_channel])*1e6
show
copy,0,1
setdata, TSys_Cal & !g.s[0].units = 'Tsys (K)'
copy,0,2
setdata, TRx & !g.s[0].units = 'Trx (K)'
copy,0,3
print,'Tcal in buffer 1;  Tsys in buffer 2; Trx in buffer 3'
copy, 1,0
show,0



if write eq 0 then begin
    print,'Freq (GHz)','TCal (K)','TSys (K)',format="(A14,A14,A14)"
    print,'----------','--------','--------',format="(A14,A14,A14)"
    for ijk=0,n_elements(freqs)-1 do begin
	print,freqs[ijk]/1000.0,Tcal[ijk],Tsys_caloff[ijk],format="(f14.4,f14.3,f14.3)"
    endfor
endif else begin
    if write eq 1 then begin
	openw,lun, fileout, /get_lun
	printf,lun,'Freq (GHz)','TCal (K)','TSys (K)',format="(A14,A14,A14)"
	printf,lun,'----------','--------','--------',format="(A14,A14,A14)"
	for ijk=0,n_elements(freqs)-1 do begin
	    printf,lun,freqs[ijk]/1000.0,Tcal[ijk],Tsys_caloff[ijk],format="(f14.4,f14.3,f14.3)"
	endfor
	close,lun
	free_lun, lun
    endif
endelse
end




;load database vector tcal from the rcvr fits file in /home/gbtdata/
pro load_db_tcal, infile, colname, buf, ext, tc_freq = tc_freq, tc_temp = tc_temp
;infile: TCAL fits file /home/gbtdata/[project]/[RcvrX_Y]/Z.fits
;colname: column name in fits, e.g. 'LOW_CAL_TEMP'
;buf: buffer to save to
;ext: extension number in fits file. corresponds to diff pols. use fv to pick the right one

tcal_data = readfits(infile,hdr,exten_no=ext)

tc_freq = tbget(hdr,tcal_data,'FREQUENCY')

tc_temp = tbget(hdr,tcal_data,colname)

setdata, tc_temp
;!g.s[0].reference_frequency = tc_freq[!g.s[0].reference_channel]
!g.s[0].reference_frequency = tc_freq[n_elements(tc_freq) - 1 - !g.s[0].reference_channel]
;!g.s[0].frequency_resolution = tc_freq[1] - tc_freq[0]
show
copy,0,buf


end
