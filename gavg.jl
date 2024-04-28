# function EXPSATO computes (2.10) from Sato and Fehler
#    Inputs:  eps  =  RMS velocity perturbation
#             a    =  correlation distance
#             m    =  wavenumber
#    Returns: P(m) =  PSDF (Power Spectral Density Function)
#
function EXPSATO(eps::Float64, a::Float64, m::Float64)
    kappa::Float64 = 0.5
    gamma_k1::Float64 = gamma(kappa)
    gamma_k2::Float64 = gamma(kappa + 0.5)
    expon::Float64 = kappa + 1.5
    expsato::Float64 = (8. * pi ^ 1.5 * eps ^ 2 * a ^ 3 * gamma_k2) / (
        gamma_k1 * (1. + a ^ 2 * m ^ 2) ^ expon)
    return expsato
end

# XSATO computes (4.50) from Sato and Fehler
#   Inputs:  psi  =  spherical coor. angle (radians)
#            zeta =  sph. coor angle from x3 (radians)
#            nu   =  density vs. velocity pert. scaling (see 4.48)
#            gam0 =  Pvel/Svel (often assumed = sqrt(3) )
#   Returns: xpp, xps, xsp, xss_psi, xss_zeta  =  from eqn. (4.50)
#
function XSATO(psi::Float64, zeta::Float64, nu::Float64, gam0::Float64)
    gam2::Float64 = gam0 ^ 2
    cpsi::Float64 = cos(psi)
    c2psi::Float64 = cos(2 * psi)
    spsi::Float64 = sin(psi)
    czeta::Float64 = cos(zeta)
    szeta::Float64 = sin(zeta)

    xpp::Float64 = (1. / gam2) * (nu * (-1. + cpsi + (2. / gam2) * spsi ^ 2)
        - 2. + (4. / gam2) * spsi ^ 2)
    xps::Float64 = -spsi * (nu * (1. - (2. / gam0) * cpsi) - (4. / gam0) * cpsi)
    xsp::Float64 = (1. / gam2) * spsi * czeta * (nu * (1. -  (2. / gam0) * cpsi)
        - (4. / gam0) * cpsi)
    
    xss_psi::Float64 = czeta * (nu * (cpsi - c2psi) - 2. * c2psi)
    xss_zeta::Float64 = szeta * (nu * (cpsi - 1.) + 2. * cpsi)
    return xpp, xps, xsp, xss_psi, xss_zeta
end

# GSATO computes (4.52) from Sato and Fehler (exponential autocor)
#   Inputs:  psi  =  spherical coor. angle (radians)
#            zeta =  sph. coor angle from x3 (radians)
#            el   =  S-wave wavenumber (=om/beta0)
#            nu   =  density vs. velocity pert. scaling (see 4.48)
#            gam0 =  Pvel/Svel (often assumed = sqrt(3) )
#            eps  =  RMS velocity perturbation
#            a    =  correlation distance
#   Returns: gpp,gps,gsp,gss  =  from eqn. (4.52)
#            spol =  S-to-S scattered S polarization (radians)
#                 =  0 for pure psi direction
function GSATO(psi::Float64, zeta::Float64, el::Float64, nu::Float64, gam0::Float64, eps::Float64, a::Float64)
    pi4::Float64 = 4.0 * pi
    el4::Float64 = el ^ 4
    gam2::Float64 = gam0 ^ 2

    xpp, xps, xsp, xss_psi, xss_zeta = XSATO(psi, zeta, nu, gam0)

    arg::Float64 = (2. * el / gam0) * sin(psi / 2.)
    gpp::Float64 = (el4 / pi4) * xpp ^ 2 * EXPSATO(eps, a, arg)
    if gpp < 1.0e-30 gpp = 0.0 end

    arg = (el / gam0) * sqrt(1. + gam2 - 2. * gam0 * cos(psi))
    gps::Float64 = (1. / gam0) * (el4 / pi4) * xps ^ 2 * EXPSATO(eps, a, arg)
    if gps < 1.0e-30 gps = 0.0 end
    gsp::Float64 = gam0 * (el4 / pi4) * xsp ^ 2 * EXPSATO(eps, a, arg)
    if gsp < 1.0e-30 gsp = 0.0 end

    arg = 2. * el * sin(psi / 2.)
    gss::Float64 = (el4 / pi4) * (xss_psi ^ 2 + xss_zeta ^ 2) * EXPSATO(eps, a, arg)
    if gss < 1.0e-30 gss = 0.0 end
    
    spol = atan(xss_zeta, xss_psi)
    return gpp, gps, gsp, gss, spol
end

# GAVGSATO2 computes average over solid angle of g functions
# in (4.52) of Sato and Fehler (exponential autocor)
#
# This version also does momentum scattering coef.
#
# Inputs:    el   =  S-wave wavenumber (=om/beta0)
#            nu   =  density vs. velocity pert. scaling (see 4.48)
#            gam0 =  Pvel/Svel (often assumed = sqrt(3) )
#            eps  =  RMS velocity perturbation
#            a    =  correlation distance
#   Returns: gpp0,gps0,gsp0,gss0  =  averaged over solid angle
#            gppm,gpsm,gspm,gssm  =  momentum scattering coef.
#
function GAVGSATO2(el::Float64, nu::Float64, gam0::Float64, eps::Float64, a::Float64)
    n::Int64 = 60
    pi4::Float64 = 4.0 * pi
    el4::Float64 = el ^ 4
    gam2::Float64 = gam0 ^ 2

    sum1::Float64 = 0.0
    sum2::Float64 = 0.0
    sumgpp::Float64 = 0.0
    sumgps::Float64 = 0.0
    sumgsp::Float64 = 0.0
    sumgss::Float64 = 0.0
    sumgpp2::Float64 = 0.0
    sumgps2::Float64 = 0.0
    sumgsp2::Float64 = 0.0
    sumgss2::Float64 = 0.0

    for i in 1:n
        psi::Float64 = float(i) * pi / float(n)
        for j in -n:n
            zeta::Float64 = float(j) * pi / float(n)
            gpp, gps, gsp, gss, spol = GSATO(psi, zeta, el, nu, gam0, eps, a)
            area::Float64 = sin(psi) * pi ^ 2 / float(n) ^ 2
            sum1 = sum1 + area
            sumgpp = sumgpp + gpp * area
            sumgps = sumgps + gps * area
            sumgsp = sumgsp + gsp * area
            sumgss = sumgss + gss * area
            fact::Float64 = 1.0 - cos(psi)
            sum2 = sum2 + fact * area
            sumgpp2 = sumgpp2 + gpp * fact * area
            sumgps2 = sumgps2 + gps * fact * area
            sumgsp2 = sumgsp2 + gsp * fact * area
            sumgss2 = sumgss2 + gss * fact * area
        end
    end
    println("test one = ", sum1/pi4)
    println("test one (2) = ", sum2/pi4)
    gpp0::Float64 = sumgpp / sum1
    gps0::Float64 = sumgps / sum1
    gsp0::Float64 = sumgsp / sum1
    gss0::Float64 = sumgss / sum1
    gppm::Float64 = sumgpp2 / sum2
    gpsm::Float64 = sumgps2 / sum2
    gspm::Float64 = sumgsp2 / sum2
    gssm::Float64 = sumgss2 / sum2
    return gpp0, gps0, gsp0, gss0, gppm, gpsm, gspm, gssm
end