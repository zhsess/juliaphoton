# Program to model seismic phases as energy packets
# Originally written 2003-2004 by Peter Shearer (pshearer@ucsd.edu)
# Converted to f90 in 2012 by Nick Mancinelli (n.j.mancinelli@gmail.com)
# Converted to Julia in 2024 by Hao Zhang (zhsess@gmail.com)
# Last modified: 2024-04-24

using SpecialFunctions

# read in parameters
vmodel::String = "./iasp91.smo"
zsource::Float64 = 0.0 # assumed source depth (km)
freq::Float64 = 1.0 # assumed freq. (Hz), affects attenuation and scattering
om::Float64 = 2.0 * pi * freq

nqlay::Int64 = 3 # number of Q layers (max=6)
zminq = Float64[0, 220, 5153.9] # min depth of Q layer
zmaxq = Float64[220, 2889, 9999.0] # max depth of Q layer
qalpha = Float64[227, 1383, 360] # Q alpha

nscatvol::Int64 = 4 # number of scattering volumes (to follow) (max=6)
zminscat = Float64[1, 200, 600, 1700] # min scat depth (km)
zmaxscat = Float64[200, 600, 1700, 2880] # max scat depth (km)
xmaxscat = Float64[999999, 999999, 999999, 999999] # max scat range from source (km)
pvelref = Float64[8.08, 9.13, 11.734, 13.131] # reference P & S velocity for layer
svelref = Float64[4.47, 4.93, 6.504, 7.055]
el = om ./ svelref
gam0 = pvelref ./ svelref
nu = Float64[0.8, 0.8, 0.8, 0.8] # relative size of density perturbation (0.8 often assumed)
eps = Float64[0.02, 0.02, 0.005, 0.005] # rms perturbation
alen = Float64[4.0, 4.0, 8.0, 8.0] # scale length (km)

# define variables
nface0::Int64 = 6
nlay0::Int64 = 655
erad::Float64 = 6371.0

z_s = zeros(Float64, nlay0)
alpha_s = zeros(Float64, nlay0)
beta_s = zeros(Float64, nlay0)
rho = zeros(Float64, nlay0)

z = zeros(Float64, nlay0)
alpha = zeros(Float64, nlay0)
beta = zeros(Float64, nlay0)

q = zeros(Float64, nlay0, 2)
slow = zeros(Float64, nlay0, 2)

iflag = zeros(Int64, nlay0)
iflagvol = zeros(Int64, nlay0)

dep::Float64 = 0.0
izsource::Int64 = 0
npts::Int64 = 0

scatprob = zeros(Float64, nface0, 2, 3)

# calculate flat earth transformation
function FLATTEN(z_s::Float64, vel_s::Float64)
    erad:: Float64 = 6371.0
    r::Float64 = erad - z_s
    z_f::Float64 = - erad * log(r / erad)
    vel_f::Float64 = vel_s * (erad / r)
    return z_f, vel_f
end

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
    # println("test one = ", sum1/pi4)
    # println("test one (2) = ", sum2/pi4)
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

# read velocity model and transform to flat earth
# value at center of earth is removed to avoid singularity

open(vmodel, "r") do f
    for i in 1:nlay0
        # Read the line and split it into variables
        line_ = split(readline(f))
        z_s[i] = parse(Float64, line_[1])
        dep = parse(Float64, line_[2])
        alpha_s[i] = parse(Float64, line_[3])
        beta_s[i] = parse(Float64, line_[4])
        rho[i] = parse(Float64, line_[5])

        if z_s[i] == zsource global izsource = i end
        if z_s[i] == erad
            global npts = i
            break
        end
        if beta_s[i] == 0.0 beta_s[i] = 0.001 end

        z[i], alpha[i] = FLATTEN(z_s[i], alpha_s[i])
        if beta_s[i] != 0.0
            z[i], beta[i] = FLATTEN(z_s[i], beta_s[i])
        else
            beta[i] = 0.0
        end
        iflag[i] = 0
        iflagvol[i] = 0

        q[i, 1] = 999999.
        for k in 1:nqlay
            if (z_s[i] >= zminq[k]) & (z_s[i] <= zmaxq[k])
                q[i, 1] = qalpha[k]
            end
        end
        q[i, 2] = (4/9) * q[i, 1]
        global npts = i
    end
end

vpmin = alpha[1]
vsmin = beta[1]
vpsource = alpha[izsource]
vssource = beta[izsource]

z_s[npts] = z_s[npts - 1]
alpha_s[npts] = alpha_s[npts - 1]
beta_s[npts] = beta_s[npts - 1]
q[npts, 1] = q[npts - 1, 1]
q[npts, 2] = q[npts - 1, 2]
z[npts], alpha[npts] = FLATTEN(z_s[npts], alpha_s[npts])
z[npts], beta[npts] = FLATTEN(z_s[npts], beta_s[npts])

println("finished reading model")
println("Depth points in model = $npts \n")

for i in 1:npts
    slow[i, 1] = 1.0 / alpha[i]
    if beta[i] != 0.0
        slow[i, 2] = 1.0 / beta[i]
    else
        slow[i, 2] = 1.0 / alpha[i]
    end
end

# set up scattering volume flags and other stuff
# Note that if scattering volumes overlap, iflagvol is set to last of scattering volumes
for iscat in 1:nscatvol
    gpp0, gps0, gsp0, gss0, gppm, gpsm, gspm, gssm = GAVGSATO2(el[iscat], nu[iscat], gam0[iscat], eps[iscat], alen[iscat])

    # println("Scattering Strengths for layer $iscat:")
    # println("gpp0 = $gpp0, gps0 = $gps0, gsp0 = $gsp0, gss0 = $gss0")
    # println("afp = ", 1. / (gpp0 + gps0), ", ", 1. / (gsp0 + gss0))
    # println("gppm = $gppm, gpsm = $gpsm, gspm = $gspm, gssm = $gssm")
    # println("moment afp =", 1. / (gppm + gpsm), ", ", 1. / (gspm + gssm), "\n")

    scatprob[iscat, 1, 1] = gpp0 + gps0
    scatprob[iscat, 1, 2] = gpp0
    scatprob[iscat, 1, 3] = gps0
    scatprob[iscat, 2, 1] = gsp0 + gss0
    scatprob[iscat, 2, 2] = gsp0
    scatprob[iscat, 2, 3] = gss0

    for i in 2:npts
        if (z_s[i-1] >= zminscat[iscat]) & (z_s[i] <= zmaxscat[iscat])
            iflagvol[i] = iscat
        end
    end

end

# line 565