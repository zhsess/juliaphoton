include("gavg.jl")

# calculate flat earth transformation
function FLATTEN(z_s::Float64, vel_s::Float64)
    erad:: Float64 = 6371.0
    r::Float64 = erad - z_s
    z_f::Float64 = - erad * log(r / erad)
    vel_f::Float64 = vel_s * (erad / r)
    return z_f, vel_f
end

# SPH_LOC finds location of second point on sphere, given range 
# and azimuth at first point.
#
# Inputs:  flat1  =  latitude of first point (degrees) 
#          flon2  =  longitude of first point (degrees)
#          del    =  angular separation between points (degrees)
#          azi    =  azimuth at 1st point to 2nd point, from N (deg.)
# Returns: flat2  =  latitude of second point (degrees)
#          flon2  =  longitude of second point (degrees)
#
function SPH_LOC(flat1::Float64, flon1::Float64, del::Float64, azi::Float64)
    if del == 0.0
        flat2::Float64 = flat1
        flon2::Float64 = flon1
        return flat2, flon2
    end
    reddeg::Float64 = pi / 180.0
    delr::Float64 = del * reddeg
    azr::Float64 = azi * reddeg
    theta1::Float64 = (90.0 - flat1) * reddeg
    phi1::Float64 = flon1 * reddeg
    ctheta2::Float64 = sin(delr) * sin(theta1) * cos(azr) + cos(theta1) * cos(delr)
    theta2 = acos(ctheta2)
    if theta1 == 0.0
        phi2 = azr
    elseif theta2 == 0.0
        phi2 = 0.0
    else
        sphi2::Float64 = sin(delr) * sin(azr) / sin(theta2)
        cphi2::Float64 = (cos(delr) - cos(theta1) * ctheta2) / (sin(theta1) * sin(theta2))
        phi2 = phi1 + atan(sphi2, cphi2)
    end
    flat2 = 90.0 - theta2 / reddeg
    flon2 = phi2 / reddeg
    if flon2 > 360.0 flon2 = flon2 - 360.0 end
    if flon2 < 0.0 flon2 = flon2 + 360.0 end
    return flat2, flon2
end

# SPH_AZI computes distance and azimuth between two points on sphere
#
# Inputs:  flat1  =  latitude of first point (degrees) 
#          flon2  =  longitude of first point (degrees)
#          flat2  =  latitude of second point (degrees)
#          flon2  =  longitude of second point (degrees)
# Returns: del    =  angular separation between points (degrees)
#          azi    =  azimuth at 1st point to 2nd point, from N (deg.)
#
# Note:  This routine is inaccurate for del less than about 0.5 degrees. 
#        For greater accuracy, use SPH_AZIDP or perform a separate
#        calculation for close ranges using Cartesian geometry.
#
function SPH_AZI(flat1::Float64, flon1::Float64, flat2::Float64, flon2::Float64)
    if (flat1 == flat2 && flon1 == flon2) ||(flat1 == 90.0 && flat2 == 90.0) || (flat1 == -90.0 && flat2 == -90.0)
        del::Float64 = 0.0
        azi::Float64 = 0.0
        return del, azi
    end
    raddeg::Float64 = pi / 180.0
    theta1::Float64 = (90.0 - flat1) * raddeg
    theta2::Float64 = (90.0 - flat2) * raddeg
    phi1::Float64 = flon1 * raddeg
    phi2::Float64 = flon2 * raddeg
    stheta1::Float64 = sin(theta1)
    ctheta1::Float64 = cos(theta1)
    stheta2::Float64 = sin(theta2)
    ctheta2::Float64 = cos(theta2)
    cang::Float64 = stheta1 * stheta2 * cos(phi2 - phi1) + ctheta1 * ctheta2
    ang::Float64 = acos(cang)
    del = ang / raddeg
    sang::Float64 = sqrt(1.0 - cang ^ 2)
    if (sang == 0.0) || (stheta1 == 0.0)
        del = 0.0
        azi = 0.0
        return del, azi
    end
    caz::Float64 = (ctheta2 - cang * ctheta1) / (sang * stheta1)
    saz::Float64 = - stheta2 * sin(phi1 - phi2) / sang
    az = atan(saz, caz)
    azi = az / raddeg
    if azi < 0.0 azi = azi + 360.0 end
    return del, azi
end

# RANDANG2 computes random scattering angle using Sato and Fehler
# equations.  This version works for multiple scattering layers
#   Input: iscat  = scattering layer number (max=6)
#          itype  = 1 for gpp
#                 = 2 for gps
#                 = 3 for gsp
#                 = 4 for gss
#            el   =  S-wave wavenumber (=om/beta0)
#            nu   =  density vs. velocity pert. scaling (see 4.48)
#            gam0 =  Pvel/Svel (often assumed = sqrt(3) )
#            eps  =  RMS velocity perturbation
#            a    =  correlation distance
#   Returns: psi  =  spherical coor. angle (radians)
#            zeta =  sph. coor angle from x3 axis (radians)
#            spol =  scattered S polarization (radians)
#                 =  0 for pure psi direction
#                    (only meaningful for gps and gss)
#
# Note:  Program computes g arrays upon first call.  Subsequent
#        calls use this array and will not be correct if 
#        el, nu, etc., are changed
mutable struct SavedValues
    firstcall::BitArray
    psi2::Array{Float64, 2}
    zeta2::Array{Float64, 2}
    sumg::Array{Float64, 3}
    nk::Int64
end
function RANDANG2(iscat::Int64, itype::Int64, el::Float64, nu::Float64, gam0::Float64, eps::Float64, a::Float64, sd::SavedValues)
    n::Int64 = 100
    npts0::Int64 = 66000
    nscat0::Int64 = 6
    g = zeros(Float64, npts0, 4, nscat0)
    spol2 = zeros(Float64, npts0, nscat0)

    idum::Int64 = 0
    if sd.firstcall[iscat]
        println("first call RANDANG2: setting up iscat arrays: $iscat")
        sd.firstcall[iscat] = false
        k = 0
        sumg1::Float64 = 0.0
        sumg2::Float64 = 0.0
        sumg3::Float64 = 0.0
        sumg4::Float64 = 0.0
        for i in 0:n
            for j in 0:n
                k = k + 1
                sd.psi2[k, iscat] = float(i) * pi / float(n)
                sd.zeta2[k, iscat] = float(j) * pi / float(n)
                g[k, 1, iscat], g[k, 2, iscat], g[k, 3, iscat], g[k, 4, iscat], spol2[k, iscat] = GSATO(
                    sd.psi2[k, iscat], sd.zeta2[k, iscat], el, nu, gam0, eps, a)
                area::Float64 = sin(sd.psi2[k, iscat]) * pi ^ 2 / float(n) ^ 2
                sumg1 = sumg1 + area * g[k, 1, iscat]
                sumg2 = sumg2 + area * g[k, 2, iscat]
                sumg3 = sumg3 + area * g[k, 3, iscat]
                sumg4 = sumg4 + area * g[k, 4, iscat]
                sd.sumg[k, 1, iscat] = sumg1
                sd.sumg[k, 2, iscat] = sumg2
                sd.sumg[k, 3, iscat] = sumg3
                sd.sumg[k, 4, iscat] = sumg4
            end
        end
        sd.nk = k
        for k in 1:sd.nk
            sd.sumg[k, 1, iscat] = sd.sumg[k, 1, iscat] / sumg1
            sd.sumg[k, 2, iscat] = sd.sumg[k, 2, iscat] / sumg2
            sd.sumg[k, 3, iscat] = sd.sumg[k, 3, iscat] / sumg3
            sd.sumg[k, 4, iscat] = sd.sumg[k, 4, iscat] / sumg4
        end
        println("Finish setting up arrays")
    end

    fran::Float64 = rand()
    k1::Int64 = 1
    k2::Int64 = sd.nk
    for it in 1:16
        k::Int64 = floor((k1 + k2) / 2)
        if fran < sd.sumg[k, itype, iscat]
            k2 = k
        elseif fran > sd.sumg[k, itype, iscat]
            k1 = k
        else
            k2 = k
        end
    end

    zeta = sd.zeta2[k, iscat]
    psi = sd.psi2[k, iscat]
    spol::Float64 = 0
    if itype == 4 spol = spol2[k2, iscat] end

    return psi, zeta, spol
end

# UPDATE_P gets new ray index number for P/S or S/P conversion
#  Inputs:    p  =  table of ray parameters for both P and S
#            np  =  number of p values in p array
#            iw1 =  incident wave type (1=P, 2=S)
#            iw2 =  scattered wave type (1=P, 2=S)
#            ip1  =  ray index for incident wave
#  Returns:  ip2  =  ray index for scattered wave
#
#  Note:  Selects closest ray parameter SMALLER than target p to avoid
#         possibility of non-existent ray (for p close to material slowness)
#
function UPDATE_P(p::Array{Float64, 2}, np::Int64, iw1::Int64, iw2::Int64, ip1::Int64)
    nray0::Int64 = 50000
    ptarg::Float64 = p[ip1, iw1]
    frac::Float64 = (ptarg - p[1, iw2]) / (p[np, iw2] - p[1, iw2])
    ip2::Int64 = floor(frac * (np - 1)) + 1
    if p[ip2, iw2] < ptarg ip2 = ip2 - 1 end
    if ip2 < 1 ip2 = ip2 + 1 end
    if ip2 > np ip2 = ip2 - 1 end
    if (ip2 < 1) || (ip2 > np) || (p[ip2, iw2] > ptarg)
        println("***ERROR in UPDATE_P")
        ip2 = 1
        return ip2
    end
    return ip2
end